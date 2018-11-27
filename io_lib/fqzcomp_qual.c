/*
 * Format:
 *
 * (General flags & meta-data)
 * byte    CRAM major version (3 or 4)
 * byte    Bit flags (4=do_rle, 2=do_dedup, 1=stored_qmap)
 * byte    Max quality sym count (eg 4, 8, 40)
 * (If flags&1; stored_qmap)
 * byte       Nsym: number of quality symbols
 * byte[nsym] Quality alphabet; eg to replace char !,-,3,E with int 0,1,2,3.
 * (Qual meta-data)
 * byte    q_qctxbits<<4  +  q_qctxshift
 * byte    q_2ctxbits<<4  +  q_mctxbits
 * byte[]  Array for q_dmap of size q_mctxbits
 * (Run length meta-data; iff do_rle set)
 * byte    r_qctxbits<<4  +  r_qctxshift
 * byte    r_pctxbits<<4  +  r_pctxshift
 * byte    r_lctxbits<<4  +  r_mctxbits
 * byte[]  Array for r_dmap of size r_mctxbits
 *
 * Models contexts must be no more than 16-bit in size, but smaller is better
 * for speed / memory..
 *
 * Arrays are 256 elements wide and must consists of a series of values between
 * 0 and 1<<mctxbits (eg 1<<q_mctxbits) in numerically ascending order.
 * They are stored as a series of 1<<mctxbits values indicating the run
 * length for that element.  Any remaining values (up to 256) are copied
 * from the last value.
 *
 * Quality context for q_{n} is computed using:
 * q_{n-1} = previous qual
 * q_{n-2} = qual 2 ago
 * q_{n-3} = qual 3 ago
 *
 * qlast = (qlast<<q_qctxshift) + q_{n-1}
 * ctx = qlast & ((1<<q_qctxbits)-1)
 * if (q_2ctxbits)
 *     ctx = (ctx<<1) | (q_{n-1} == q_{n-3});
 * ctx <<= q_mctxbits
 * ctx += q_dmap[MIN(255, delta)]
 * delta += (q_{n-1} != q_{n-2});
 *
 * Fixed params.
 *   MAXR = 8 (nsym for run-lengths, keep adding MAXR to run until we read 0)
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "cram_block_compression.h"
#include "fqzcomp_qual.h"

// Define this for version 3.1 of CRAM only.  In the version 4.0
// prototype we switch the orientation of quality values to their
// original direction.  The gain is sufficient that even with v3 it's
// best to do this and store additional data in this codec to record
// this fact (ie a duplicate of the BAM_FREVERSE flag).

static const char *name(void) {
    return "fqzcomp-qual";
}

#define QMAX 128
#define QBITS 12
#define QSIZE (1<<QBITS)

#define NSYM 256
#include "c_simple_model.h"

#undef NSYM
#define NSYM QMAX
//#include "c_escape_model.h"
#include "c_simple_model.h"

#undef NSYM
#define NSYM 2
#include "c_simple_model.h"

// Fast tuning for small slice sizes, but not giving much up on large ones.
#define MAXR 8
#undef NSYM
#define NSYM MAXR
#include "c_simple_model.h"

static int store_array(unsigned char *out, int *array, int size, int bits) {
    int i, j, k;
    for (i = j = k = 0; i < size && j < (1<<bits); j++) {
	int run_len = i;
	while (i < size && array[i] == j)
	    i++;
	run_len = i-run_len;
	out[k++] = run_len;
    }
    while (j < (1<<bits))
	out[k++] = 0, j++;

    return k;
}

static int read_array(unsigned char *in, int *array, int bits) {
    int i, j, k;
    for (i = j = k = 0; i < (1<<bits); i++) {
	int run_len = in[k++];
	while (run_len && j < 256)
	    run_len--, array[j++] = i;
    }

    // Copy last element to end
    i = array[j-1];
    while (j < 256)
	array[j++] = i;
    return k;
}

unsigned char *compress_block_fqz2f(int vers,
				    int level,
				    cram_slice *s,
				    unsigned char *in,
				    size_t in_size,
				    size_t *out_size) {
    //approx sqrt(delta), must be sequential
    static int dsqr[] = {
	0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
	6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
    };

    //approx sqrt(delta)/2
    static int dsqr2[] = {
	0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3,
	4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
	7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
//	0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
//	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3,
//	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
//	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    };

    unsigned char *comp = (unsigned char *)malloc(in_size*1.1+1000);
    unsigned char *comp2 = (unsigned char *)malloc(in_size*1.1+1000);
    if (!comp || !comp2) {
	free(comp);
	free(comp2);
	return NULL;
    }

    int comp_idx = 0;
    size_t i, j;
    ssize_t rec = 0;
    unsigned int last = 0, qlast = 0, rlast = 0;
    RangeCoder rc, rc2;
    unsigned char q1 = 1, q2 = 0;
    int run_len = 0;

    uint32_t qhist[256] = {0}, nsym, max_sym, qrle = 0;
    int do_rle = 0; // RLE not helping much at the moment, so avoid complexity and memory.
    int do_dedup = 1;
    for (last = i = 0; i < in_size; i++) {
	qhist[in[i]]++;
	if (last == in[i])
	    qrle++;
	last = in[i];
    }
    for (i = max_sym = nsym = 0; i < 256; i++)
	if (qhist[i])
	    max_sym = i, nsym++;

    int store_qmap = (nsym <= 8 && nsym*2 < max_sym);

    if (100.*qrle/in_size < 75)
	do_rle = 0;

    // If nsym is significant lower than max_sym, store a lookup table and
    // transform prior to compression.  Eg 4 or 8 binned data.
    comp[comp_idx++] = vers;
    comp[comp_idx++] = (do_rle<<2)|(do_dedup<<1)|(store_qmap);
    comp[comp_idx++] = max_sym;

    // optimal in Q40 100k/slice
    int q_qctxbits =10;
    int q_qctxshift=5; // qmax 64, although we can store up to 128 if needed
    int q_2ctxbits =0; // 0 or 1 only
    int q_pctxbits =4;
    int q_pctxshift=2;
    int q_mctxbits =2;
    int q_mctxshift=1;

    if (nsym <= 4) {
	// NovaSeq
	q_qctxbits =10;
	q_qctxshift=2; // qmax 64, although we can store up to 128 if needed
	q_pctxbits =2;
	q_pctxshift=5;
	if (in_size < 5000000) {
	    q_pctxbits =1;
	    q_pctxshift=6;
	    q_mctxbits =2;
	    q_mctxshift=1;
	}
    } else if (nsym <= 8) {
	// HiSeqX
	q_qctxbits =9;
	q_qctxshift=3;
	q_pctxbits =3;
	q_pctxshift=5;
	q_mctxbits =2;
	q_mctxshift=2;
	//q_2ctxbits =1;
	if (in_size < 5000000) {
	    q_qctxbits =6;
	    q_pctxbits =2;
	    q_pctxshift=5;
	}
    }

    if (in_size < 300000) {
	q_qctxbits=q_qctxshift;
	q_2ctxbits=0;
	q_mctxbits=2;
    }

    for (i = 0; i < sizeof(dsqr)/sizeof(*dsqr); i++)
	if (dsqr[i] > (1<<q_mctxbits)-1)
	    dsqr[i] = (1<<q_mctxbits)-1;

    if (store_qmap) {
	comp[comp_idx++] = nsym;
	int comp_idx_start = comp_idx;
	for (i = 0; i < 256; i++)
	    if (qhist[i])
		comp[comp_idx] = i, qhist[i] = comp_idx++ -comp_idx_start;
	max_sym = nsym;
    } else {
	nsym = 255;
	for (i = 0; i < 256; i++)
	    qhist[i] = i;
    }

    comp[comp_idx++] = (q_qctxbits<<4)|q_qctxshift;
    comp[comp_idx++] = (q_pctxbits<<4)|q_pctxshift;
    comp[comp_idx++] = (q_2ctxbits<<7)|(q_mctxbits<<3)|q_mctxshift;

    comp_idx += store_array(comp+comp_idx, dsqr, sizeof(dsqr)/sizeof(*dsqr), q_mctxbits);

    int r_qctxbits =8;
    int r_qctxshift=6;
    int r_pctxbits =3;
    int r_pctxshift=1;
    int r_lctxbits =2;
    int r_mctxbits =3;

    for (i = 0; i < sizeof(dsqr2)/sizeof(*dsqr2); i++)
	if (dsqr2[i] > (1<<r_mctxbits)-1)
	    dsqr2[i] = (1<<r_mctxbits)-1;

    if (do_rle) {
	if (nsym <= 4) {
	    if (in_size < 5000000) {
		r_qctxbits=2;
		r_lctxbits=2;
	    } else {
		r_qctxbits=5;
		r_lctxbits=3;
	    }
	    r_qctxshift=2;
	} else if (nsym <= 8) {
	    r_qctxbits =6;
	    r_qctxshift=3;
	    r_lctxbits =2;
	}

	comp[comp_idx++] = (r_qctxbits<<4)|r_qctxshift;
	comp[comp_idx++] = (r_pctxbits<<4)|r_pctxshift;
	comp[comp_idx++] = (r_lctxbits<<4)|r_mctxbits;

	comp_idx += store_array(comp+comp_idx, dsqr2, sizeof(dsqr2)/sizeof(*dsqr2), r_mctxbits);
    }

    SIMPLE_MODEL(QMAX,_) *model_qual;

    model_qual = malloc(sizeof(*model_qual) * (1<<16));
    if (!model_qual)
	return NULL;

    for (i = 0; i < (1<<16); i++)
	SIMPLE_MODEL(QMAX,_init)(&model_qual[i], max_sym+1);

    SIMPLE_MODEL(256,_) model_len[4];
    for (i = 0; i < 4; i++)
	SIMPLE_MODEL(256,_init)(&model_len[i],256);

    SIMPLE_MODEL(2,_) model_strand;
    SIMPLE_MODEL(2,_init)(&model_strand,2);

    SIMPLE_MODEL(MAXR,_) *model_run = NULL;
    if (do_rle) {
	model_run = malloc(sizeof(*model_run)*(1<<16));
	if (!model_run)
	    return NULL;
	for (i = 0; i < (1<<16); i++)
	    SIMPLE_MODEL(MAXR,_init)(&model_run[i],MAXR);
    }

    int delta = 0, delta2 = 0, j2 = 0;
    int last_len = 0;
    SIMPLE_MODEL(2,_) model_dup;
    SIMPLE_MODEL(2,_init)(&model_dup,2);

    RC_SetOutput(&rc, (char *)comp+comp_idx+4);
    RC_StartEncode(&rc);

    RC_SetOutput(&rc2, (char *)comp2);
    RC_StartEncode(&rc2);

    int ndup0 = 0, ndup1 = 0;

    if (vers == 3) {
	// Pass 1, reverse all seqs if necessary.
	i = rec = j = 0;
	while (i < in_size) {
	    int len = rec < s->hdr->num_records-1
		? s->crecs[rec].len
		: in_size - i;

	    if (s->crecs[rec].flags & BAM_FREVERSE) {
		// Reverse complement sequence - note: modifies buffer
		int I,J;
		unsigned char *cp = in+i;
		for (I = 0, J = len-1; I < J; I++, J--) {
		    unsigned char c;
		    c = cp[I];
		    cp[I] = cp[J];
		    cp[J] = c;
		}
	    }

	    i += len;
	    rec++;
	}
	rec = 0;
    }

    // We encode in original orientation for V3.1 (automatically done for V4.0).
    // Therefore during decode we'll need to do two passes, first to decode
    // and perform delta, and second to reverse back again.
    // These two can be merged with a bit of cleverness, but we do it simply for now.
    int nrun = 0;
    for (run_len = i = j = 0; i < in_size; i++, j--) {
	if (j == 0) {
	    // Quality buffer maybe longer than sum of reads if we've
	    // inserted a specific base + quality pair.
	    int len = rec < s->hdr->num_records-1
		? s->crecs[rec].len
		: in_size - i;
	    SIMPLE_MODEL(256,_encodeSymbol)(&model_len[0], &rc, (len>> 0) & 0xff);
	    SIMPLE_MODEL(256,_encodeSymbol)(&model_len[1], &rc, (len>> 8) & 0xff);
	    SIMPLE_MODEL(256,_encodeSymbol)(&model_len[2], &rc, (len>>16) & 0xff);
	    SIMPLE_MODEL(256,_encodeSymbol)(&model_len[3], &rc, (len>>24) & 0xff);

	    if (vers == 3) {
		// no need to reverse complement for V4.0 as the core format
		// already has this feature.
		if (s->crecs[rec].flags & BAM_FREVERSE)
		    SIMPLE_MODEL(2,_encodeSymbol)(&model_strand, &rc, 1);
		else
		    SIMPLE_MODEL(2,_encodeSymbol)(&model_strand, &rc, 0);
	    }

	    rec++;
	    j = len;
	    delta = 0;
	    qlast = last = 0;

	    if (do_dedup) {
		// Possible dup of previous read?
		if (i && len == last_len && !memcmp(in+i-last_len, in+i, len)) {
		    SIMPLE_MODEL(2,_encodeSymbol)(&model_dup, &rc, 1);
		    i += len-1;
		    j = 1;
		    ndup1++;
		    continue;
		}
		SIMPLE_MODEL(2,_encodeSymbol)(&model_dup, &rc, 0);
		ndup0++;

		last_len = len;
	    }
	}

	unsigned char q = in[i];
        //unsigned char q = in[i] & (QMAX-1);
	//assert(in[i] < QMAX && in[i] >= 0);
	if ((do_rle && q != q1) && i > 0) {
	    // Every symbol is sym+rep_count, even the A+0 case.
	    // Rep count is based on symbol and history itself
	    int looped = 0;
	    rlast = (rlast<<r_qctxshift) + qhist[q1];
	    do {
		int r = run_len>MAXR-1?MAXR-1:run_len;
		nrun++;

		int ctx = rlast & ((1<<r_qctxbits)-1);
		ctx <<= r_pctxbits; ctx |= MIN((1<<r_pctxbits)-1, j2>>r_pctxshift);
		ctx <<= r_lctxbits; ctx |= MIN((1<<r_lctxbits)-1,looped);
		ctx <<= r_mctxbits; ctx |= dsqr2[MIN(sizeof(dsqr2)/sizeof(*dsqr2)-1,delta2)];

		SIMPLE_MODEL(MAXR,_encodeSymbol)(&model_run[ctx], &rc2, r);
		run_len -= MAXR-1;
		looped++;
	    } while (run_len >= 0);
	    run_len = 0;
	} else if (i>0) {
	    run_len++;
	}

	if (q != q1 || !do_rle) {
	    SIMPLE_MODEL(QMAX,_encodeSymbol)(&model_qual[last], &rc, qhist[q]);

	    delta2 = delta;
	    j2 = j;

	    qlast = (qlast<<q_qctxshift) + qhist[q];
	    last = qlast & ((1<<q_qctxbits)-1);
	    if (q_2ctxbits) {
		last <<= 1; last |= (q == q2);
	    }
	    if (q_pctxbits)
		last = (last<<q_pctxbits) | MIN((1<<q_pctxbits)-1, j>>q_pctxshift);
	    if (q_mctxbits) {
		last <<= q_mctxbits;
		last += dsqr[MIN(sizeof(dsqr)/sizeof(*dsqr)-1,delta>>q_mctxshift)];
	    }
	    q2 = q1;
	}

	delta += (q1 != q);
	q1 = q;
    }

    // Run length for last symbol
    if (do_rle) {
	int looped=0;
	rlast = (rlast<<r_qctxshift) + qhist[q1];
	do {
	    int r = run_len>MAXR-1?MAXR-1:run_len;

	    int ctx = rlast & ((1<<r_qctxbits)-1);
	    ctx <<= r_pctxbits; ctx |= MIN((1<<r_pctxbits)-1, j2>>r_pctxshift);
	    ctx <<= r_lctxbits; ctx |= MIN((1<<r_lctxbits)-1,looped);
	    ctx <<= r_mctxbits; ctx |= dsqr2[MIN(sizeof(dsqr2)/sizeof(*dsqr2)-1,delta2)];

	    SIMPLE_MODEL(MAXR,_encodeSymbol)(&model_run[ctx], &rc2, r);
	    run_len -= MAXR-1;
	    looped++;
	} while (run_len >= 0);
    }

    //SIMPLE_MODEL(QMAX,_encodeSymbol)(&model_qual[last], &rc, QMAX-1);
    RC_FinishEncode(&rc);
    RC_FinishEncode(&rc2);

    if (vers == 3) {
	// Pass 3, un-reverse all seqs if necessary.
	i = rec = j = 0;
	while (i < in_size) {
	    int len = rec < s->hdr->num_records-1
		? s->crecs[rec].len
		: in_size - i;

	    if (s->crecs[rec].flags & BAM_FREVERSE) {
		// Reverse complement sequence - note: modifies buffer
		int I,J;
		unsigned char *cp = in+i;
		for (I = 0, J = len-1; I < J; I++, J--) {
		    unsigned char c;
		    c = cp[I];
		    cp[I] = cp[J];
		    cp[J] = c;
		}
	    }

	    i += len;
	    rec++;
	}
    }

    comp[comp_idx++] = (RC_OutSize(&rc) >> 0) & 0xff;
    comp[comp_idx++] = (RC_OutSize(&rc) >> 8) & 0xff;
    comp[comp_idx++] = (RC_OutSize(&rc) >>16) & 0xff;
    comp[comp_idx++] = (RC_OutSize(&rc) >>24) & 0xff;

    memcpy(comp + comp_idx + RC_OutSize(&rc), comp2, RC_OutSize(&rc2));
    free(comp2);

    *out_size = comp_idx + RC_OutSize(&rc) + RC_OutSize(&rc2);

    free(model_qual);
    if (do_rle)
	free(model_run);

    return comp;
}

unsigned char *uncompress_block_fqz2f(cram_slice *s,
				      unsigned char *in,
				      size_t in_size,
				      size_t *out_size) {
    int dsqr[256] = {0};
    int dsqr2[256] = {0};

    unsigned char *uncomp = NULL;
    RangeCoder rc, rc2;
    size_t i, j, rec = 0, len = *out_size, in_idx = 0;
    unsigned char q1 = 1, q2 = 0;
    unsigned int last = 0, qlast = 0, rlast = 0;
    unsigned int run_len = 0;

    int vers       = in[in_idx++];
    int flags      = in[in_idx++];
    int do_rle     = flags&4;
    int do_dedup   = flags&2;
    int store_qmap = flags&1;
    int max_sym    = in[in_idx++];

    int nsym = store_qmap ? in[in_idx++] : 255;
    int qmap[256];

    if (nsym <= 8) {
	for (i = 0; i < nsym; i++)
	    qmap[i] = in[in_idx++];
	max_sym = nsym;
    } else {
	for (i = 0; i < 256; i++)
	    qmap[i] = i;
    }

    int q_qctxbits = in[in_idx]>>4;
    int q_qctxshift= in[in_idx++]&15;
    int q_pctxbits = in[in_idx]>>4;
    int q_pctxshift= in[in_idx++]&15;
    int q_2ctxbits = in[in_idx]>>7;
    int q_mctxbits = (in[in_idx]>>3)&15;
    int q_mctxshift= in[in_idx++]&7;

    if (q_qctxbits+q_pctxbits+q_2ctxbits+q_mctxbits > 16)
	return NULL;

    in_idx += read_array(in+in_idx, dsqr, q_mctxbits);

    int r_qctxbits =0;
    int r_qctxshift=0;
    int r_pctxbits =0;
    int r_pctxshift=0;
    int r_lctxbits =0;
    int r_mctxbits =0;
    if (do_rle) {
	r_qctxbits = in[in_idx]>>4;
	r_qctxshift= in[in_idx++]&15;
	r_pctxbits = in[in_idx]>>4;
	r_pctxshift= in[in_idx++]&15;
	r_lctxbits = in[in_idx]>>4;
	r_mctxbits = in[in_idx++]&15;

	in_idx += read_array(in+in_idx, dsqr2, r_mctxbits);

	if (r_qctxbits+r_pctxbits+r_lctxbits+r_mctxbits > 16)
	    return NULL;
    }

    SIMPLE_MODEL(QMAX,_) *model_qual;
    model_qual = malloc(sizeof(*model_qual) * (1<<16));
    if (!model_qual)
	return NULL;
    for (i = 0; i < (1<<16); i++)
	SIMPLE_MODEL(QMAX,_init)(&model_qual[i],max_sym+1);

    SIMPLE_MODEL(256,_) model_len[4];
    for (i = 0; i < 4; i++)
	SIMPLE_MODEL(256,_init)(&model_len[i],256);

    SIMPLE_MODEL(MAXR,_) *model_run = NULL;
    if (do_rle) {
	model_run = malloc(sizeof(*model_run)*(1<<16));
	if (!model_run)
	    return NULL;
	for (i = 0; i < (1<<16); i++)
	    SIMPLE_MODEL(MAXR,_init)(&model_run[i],MAXR);
    }

    SIMPLE_MODEL(2,_) model_strand;
    SIMPLE_MODEL(2,_init)(&model_strand,2);

    int delta = 0, delta2 = 0, j2 = 0, rev = 0;

    SIMPLE_MODEL(2,_) model_dup;
    SIMPLE_MODEL(2,_init)(&model_dup,2);

    uncomp = (unsigned char *)malloc(*out_size);
    if (!uncomp)
	return NULL;

    uint32_t cq_len = 
	  (in[in_idx+0]<< 0)
	| (in[in_idx+1]<< 8)
	| (in[in_idx+2]<<16)
	| (in[in_idx+3]<<24);
    in_idx+=4;

    RC_SetInput(&rc, (char *)in+in_idx);
    RC_StartDecode(&rc);

    RC_SetInput(&rc2, (char *)in+in_idx+cq_len);
    RC_StartDecode(&rc2);

    int nrec = 1000;
    char *rev_a = malloc(nrec);
    int *len_a = malloc(nrec * sizeof(int));
    if (!rev_a || !len_a)
	return NULL;

    for (rec = i = j = 0; i < len; i++, j--) {
	if (rec >= nrec) {
	    nrec *= 2;
	    rev_a = realloc(rev_a, nrec);
	    len_a = realloc(len_a, nrec*sizeof(int));
	    if (!rev_a || !len_a)
		return NULL;
	}

	if (j == 0) {
	    int len;
	    len  = SIMPLE_MODEL(256,_decodeSymbol)(&model_len[0], &rc);
	    len |= SIMPLE_MODEL(256,_decodeSymbol)(&model_len[1], &rc)<<8;
	    len |= SIMPLE_MODEL(256,_decodeSymbol)(&model_len[2], &rc)<<16;
	    len |= SIMPLE_MODEL(256,_decodeSymbol)(&model_len[3], &rc)<<24;
	    if (vers == 3) {
		rev = SIMPLE_MODEL(2,_decodeSymbol)(&model_strand, &rc);
		rev_a[rec] = rev;
		len_a[rec] = len;
	    }

	    if (do_dedup) {
		if (SIMPLE_MODEL(2,_decodeSymbol)(&model_dup, &rc)) {
		    // Dup of last line
		    memcpy(uncomp+i, uncomp+i-len, len);
		    i += len-1;
		    j = 1;
		    rec++;
		    continue;
		}
	    }

	    rec++;
	    j = len;
	    delta = 0;
	    qlast = last = 0;
	}

	unsigned char q, Q;
	
	if (run_len) {
	    q = q1;
	    run_len--;
	} else {
	    Q = SIMPLE_MODEL(QMAX,_decodeSymbol)(&model_qual[last], &rc);
	    q = qmap[Q];

	    qlast = (qlast<<q_qctxshift) + Q;
	    last = qlast & ((1<<q_qctxbits)-1);
	    if (q_2ctxbits) {
		last <<= 1; last |= (q == q2);
	    }
	    if (q_pctxbits)
		last = (last<<q_pctxbits) | MIN((1<<q_pctxbits)-1, j>>q_pctxshift);
	    if (q_mctxbits) {
		last <<= q_mctxbits;
		last += dsqr[MIN(255, delta>>q_mctxshift)];
	    }

	    q2 = q1;

	    delta2 = delta;
	    j2 = j;

	    run_len = 0;
	    if (do_rle) {
		int r;
		int looped = 0;
		rlast = (rlast<<r_qctxshift) + Q;
		do {
		    int ctx = rlast & ((1<<r_qctxbits)-1);
		    ctx <<= r_pctxbits; ctx |= MIN((1<<r_pctxbits)-1, j2>>r_pctxshift);
		    ctx <<= r_lctxbits; ctx |= MIN((1<<r_lctxbits)-1,looped);
		    ctx <<= r_mctxbits; ctx |= dsqr2[MIN(255,delta2)];

		    r = SIMPLE_MODEL(MAXR,_decodeSymbol)(&model_run[ctx], &rc2);
		    run_len += r;
		    looped++;
		} while (r == MAXR-1);
	    }
	}

	delta += (q1 != q);
	q1 = q;
        uncomp[i] = q;
    }
    rev_a[rec] = rev;
    len_a[rec] = len;

    if (vers == 3) {
	for (i = rec = 0; i < len; i += len_a[rec++]) {
	    if (!rev_a[rec])
		continue;

	    int I, J;
	    unsigned char *cp = uncomp+i;
	    for (I = 0, J = len_a[rec]-1; I < J; I++, J--) {
		unsigned char c;
		c = cp[I];
		cp[I] = cp[J];
		cp[J] = c;
	    }
	}
    }

    RC_FinishDecode(&rc);
    RC_FinishDecode(&rc2);
    free(model_qual);
    if (do_rle)
	free(model_run);
    free(rev_a);
    free(len_a);

    return uncomp;
}


static cram_compressor c = {
    'q', //FOUR_CC("FQZq"),
    1<<DS_QS, // quality only
    1.0,
    name,
    //compress_block_fqz2_BIN,

    //compress_block_fqz2,
    //uncompress_block_fqz2,

    compress_block_fqz2f,
    uncompress_block_fqz2f,

    //compress_block,
    //uncompress_block,
};

cram_compressor *cram_compressor_init(void) {
    return &c;
}

char *fqz_compress(int vers, cram_slice *s, char *in, size_t uncomp_size,
		   size_t *comp_size, int level) {
    return (char *)compress_block_fqz2f(vers, level, s, (unsigned char *)in,
					uncomp_size, comp_size);
}

char *fqz_decompress(char *in, size_t comp_size, size_t *uncomp_size) {
    return (char *)uncompress_block_fqz2f(NULL, (unsigned char *)in,
					  comp_size, uncomp_size);
}

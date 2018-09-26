#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "cram_block_compression.h"
#include "fqzcomp_qual.h"

// Comment out if you wish to not de-dup the last quality string.
// Generally it makes little difference, but some data sets (eg tophat
// alignments) benefit greatly from duplication between adjacent rows.
#define DEDUP

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

#define MAXR 9
#undef NSYM
#define NSYM MAXR
#include "c_simple_model.h"

// Fqzcomp -q2 equiv
unsigned char *compress_block_fqz2f(int vers,
				    int level,
				    cram_slice *s,
				    unsigned char *in,
				    size_t in_size,
				    size_t *out_size) {
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
    unsigned int last = 0;
    RangeCoder rc, rc2;
    unsigned char q1 = 1, q2 = 0;
    int run_len = 0;

    uint32_t qhist[256] = {0}, nsym, max_sym;
    for (i = 0; i < in_size; i++)
	qhist[in[i]]++;
    for (i = max_sym = nsym = 0; i < 256; i++)
	if (qhist[i])
	    max_sym = i, nsym++;

    // If nsym is significant lower than max_sym, store a lookup table and
    // transform prior to compression.  Eg 4 or 8 binned data.
    comp[comp_idx++] = vers;
    comp[comp_idx++] = max_sym;
    if (nsym <= 8 && nsym*2 < max_sym) {
	comp[comp_idx++] = nsym;
	for (i = 0; i < 256; i++)
	    if (qhist[i])
		comp[comp_idx] = i, qhist[i] = comp_idx++ -3;
	max_sym = nsym;
    } else {
	comp[comp_idx++] = 0;
    }

    SIMPLE_MODEL(QMAX,_) *model_qual;

    model_qual = malloc(sizeof(*model_qual) * QSIZE*16);
    if (!model_qual)
	return NULL;
    for (i = 0; i < QSIZE*16; i++)
	SIMPLE_MODEL(QMAX,_init)(&model_qual[i], max_sym+1);

    SIMPLE_MODEL(256,_) model_len[4];
    for (i = 0; i < 4; i++)
	SIMPLE_MODEL(256,_init)(&model_len[i],256);

    SIMPLE_MODEL(2,_) model_strand;
    SIMPLE_MODEL(2,_init)(&model_strand,2);

    SIMPLE_MODEL(MAXR,_) *model_run;
    model_run = malloc(sizeof(*model_run)*(QMAX<<11));
    if (!model_run)
	return NULL;
    for (i = 0; i < QMAX<<11; i++)
	SIMPLE_MODEL(MAXR,_init)(&model_run[i],MAXR);

    int delta = 5, delta2 = 5, j2 = 0;
#ifdef DEDUP
    int last_len = 0;
    SIMPLE_MODEL(2,_) model_dup;
    SIMPLE_MODEL(2,_init)(&model_dup,2);
#endif

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
    int nswitch = 0;
    for (i = j = 0; i < in_size; i++, j--) {
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
	    delta = 5;
	    last = 0;

#ifdef DEDUP
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
#endif
	}

	unsigned char q = in[i];
        //unsigned char q = in[i] & (QMAX-1);
	//assert(in[i] < QMAX && in[i] >= 0);
	if (q != q1 && i > 0) {
	    // Every symbol is sym+rep_count, even the A+0 case.
	    // Rep count is based on symbol and history itself
	    int looped = 0;
	    do {
		int r = run_len>MAXR-1?MAXR-1:run_len;
		nrun++;

		int ctx = q1 & (QMAX-1);
		ctx <<= 4; ctx |= ((j2/16)&15);
		ctx <<= 2; ctx |= looped;
		ctx <<= 2; ctx |= ((delta2/16) & 0x7);

		SIMPLE_MODEL(MAXR,_encodeSymbol)(&model_run[ctx], &rc2, r);
		run_len -= MAXR-1;
		looped++;
		if (looped>3) looped=3;
	    } while (run_len >= 0);
	    run_len = 0;
	} else if (i>0) {
	    run_len++;
	}

	if (q != q1) {
	    nswitch++;
	    if (nsym <= 8)
		SIMPLE_MODEL(QMAX,_encodeSymbol)(&model_qual[last], &rc, qhist[q]);
	    else
		SIMPLE_MODEL(QMAX,_encodeSymbol)(&model_qual[last], &rc, q);

	    delta2 = delta;
	    j2 = j;

	    last = ((q1<<6) | q) & (QSIZE-1);
	    last |= (q == q2)<<QBITS;
	    last |= (MIN(7*8, delta*4)&0xf8) << (QBITS-2);
	    q2 = q1;
	}

	delta += (q1 != q);
	q1 = q;
    }

    // Run length for last symbol
    {
	int looped=0;
	do {
	    int r = run_len>MAXR-1?MAXR-1:run_len;

	    int ctx = q1 & (QMAX-1);
	    ctx <<= 4; ctx |= ((j2/16)&15);
	    ctx <<= 2; ctx |= looped;
	    ctx <<= 2; ctx |= ((delta2/16) & 0x7);

	    SIMPLE_MODEL(MAXR,_encodeSymbol)(&model_run[ctx], &rc2, r);
	    run_len -= MAXR-1;
	    looped++;
	    if (looped>3) looped=3;
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

//    fprintf(stderr, "%d switches, %d runs\n", nswitch, nrun);
//    fprintf(stderr, "comp_idx %d\n", (int)comp_idx+4);
//    fprintf(stderr, "rc  size %d\n", (int)RC_OutSize(&rc));
//    fprintf(stderr, "rc2 size %d\n", (int)RC_OutSize(&rc2));
//    fprintf(stderr, "dup = %d + %d\n", ndup0, ndup1);

    comp[comp_idx++] = (RC_OutSize(&rc) >> 0) & 0xff;
    comp[comp_idx++] = (RC_OutSize(&rc) >> 8) & 0xff;
    comp[comp_idx++] = (RC_OutSize(&rc) >>16) & 0xff;
    comp[comp_idx++] = (RC_OutSize(&rc) >>24) & 0xff;

    memcpy(comp + comp_idx + RC_OutSize(&rc), comp2, RC_OutSize(&rc2));
    free(comp2);

    *out_size = comp_idx + RC_OutSize(&rc) + RC_OutSize(&rc2);

    free(model_qual);
    free(model_run);

    return comp;
}

unsigned char *uncompress_block_fqz2f(cram_slice *s,
				      unsigned char *in,
				      size_t in_size,
				      size_t *out_size) {
    unsigned char *uncomp = NULL;
    RangeCoder rc, rc2;
    size_t i, j, rec = 0, len = *out_size, in_idx = 0;
    unsigned char q1 = 1, q2 = 0;
    unsigned int last = 0;
    unsigned int run_len = 0;

    int vers = in[in_idx++];
    int max_sym = in[in_idx++];
    int nsym = in[in_idx++];

    int qmap[256];
    if (nsym && nsym <= 8) {
	for (i = 0; i < nsym; i++)
	    qmap[i] = in[in_idx++];
	max_sym = nsym;
    } else {
	nsym = QMAX;
    }

    SIMPLE_MODEL(QMAX,_) *model_qual;
    model_qual = malloc(sizeof(*model_qual) * QSIZE*16);
    if (!model_qual)
	return NULL;
    for (i = 0; i < QSIZE*16; i++)
	SIMPLE_MODEL(QMAX,_init)(&model_qual[i],max_sym+1);

    SIMPLE_MODEL(256,_) model_len[4];
    for (i = 0; i < 4; i++)
	SIMPLE_MODEL(256,_init)(&model_len[i],256);

    SIMPLE_MODEL(MAXR,_) *model_run;
    model_run = malloc(sizeof(*model_run)*(QMAX<<11));
    if (!model_run)
	return NULL;
    for (i = 0; i < QMAX<<11; i++)
	SIMPLE_MODEL(MAXR,_init)(&model_run[i],MAXR);

    SIMPLE_MODEL(2,_) model_strand;
    SIMPLE_MODEL(2,_init)(&model_strand,2);

    int delta = 5, delta2 = 5, j2 = 0, rev = 0;
#ifdef DEDUP
    SIMPLE_MODEL(2,_) model_dup;
    SIMPLE_MODEL(2,_init)(&model_dup,2);
#endif

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

#ifdef DEDUP
	    if (SIMPLE_MODEL(2,_decodeSymbol)(&model_dup, &rc)) {
		// Dup of last line
		memcpy(uncomp+i, uncomp+i-len, len);
		i += len-1;
		j = 1;
		rec++;
		continue;
	    }
#endif

	    rec++;
	    j = len;
	    delta = 5;
	    last = 0;
	}

	unsigned char q;
	
	if (run_len) {
	    q = q1;
	    run_len--;
	} else {
	    q = SIMPLE_MODEL(QMAX,_decodeSymbol)(&model_qual[last], &rc);
	    if (nsym <= 8) q = qmap[q]; // remove conditional here by always filling qmap.

	    last = ((q1<<6) | q) & (QSIZE-1);
	    last += (q == q2)<<QBITS;
	    last += (MIN(7*8, delta*4)&0xf8) << (QBITS-2);

	    q2 = q1;

	    delta2 = delta;
	    j2 = j;

	    int r;
	    run_len = 0;
	    int looped = 0;
	    do {
		int ctx = q;
		ctx <<= 4; ctx |= ((j2/16)&15);
		ctx <<= 2; ctx |= looped;
		ctx <<= 2; ctx |= ((delta2/16) & 0x7);

		r = SIMPLE_MODEL(MAXR,_decodeSymbol)(&model_run[ctx], &rc2);
		run_len += r;
		looped++;
		if (looped>3) looped=3;
	    } while (r == MAXR-1);
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

// We use generic maps to turn 0-M into 0-N where N <= M
// before adding these into the context.  These are used
// for positions, running-diffs and quality values.
//
// This can be used as a simple divisor, eg pos/24 to get
// 2 bits of positional data for each quarter along a 100bp
// read, or it can be tailored for specific such as noting
// the first 5 cycles are poor, then we have stability and
// a gradual drop off in the last 20 or so.  Perhaps we then
// map pos 0-4=0, 5-79=1, 80-89=2, 90-99=3.
//
// We don't need to specify how many bits of data we are
// using (2 in the above example), as that is just implicit
// in the values in the map.  Specify not to use a map simply
// disables that context type (our map is essentially 0-M -> 0).

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#ifdef TEST_MAIN
#include <fcntl.h>
/*
 * Standalone hackery.
 * Fixed sizes, just for testing and benchmarking
 */
typedef struct {
    int num_records;
} cram_hdr;

typedef struct {
    int len;
    int qual;
    int flags;
} cram_crec;

typedef struct {
    cram_hdr *hdr;
    cram_crec crecs[1000000];
} cram_slice;

static cram_hdr fixed_hdr;
static cram_slice fixed_slice = {0};
#define BAM_FREVERSE 0
#define BAM_FREAD2 128

cram_slice *fake_slice(size_t buf_len, int len) {
    fixed_slice.hdr = &fixed_hdr;
    fixed_hdr.num_records = (buf_len+len-1) / len;
    assert(fixed_hdr.num_records <= 1000000);
    int i, tlen = 0;
    for (i = 0; i < fixed_hdr.num_records; i++) {
	fixed_slice.crecs[i].len = len;
	fixed_slice.crecs[i].qual = tlen;
	fixed_slice.crecs[i].flags = 0;
	tlen += len;
    }

    return &fixed_slice;
}

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#else
#include "cram_block_compression.h"
#include "fqzcomp_qual.h"
#endif

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

#define NSYM 2
#include "c_simple_model.h"

#undef NSYM
#define NSYM 256
#include "c_simple_model.h"

#undef NSYM
#define NSYM QMAX
//#include "c_escape_model.h"
#include "c_simple_model.h"
//#include "c_cdf_model.h"
//#include "c_cdf16_model.h"

// An array of 0,0,0, 1,1,1,1, 3, 5,5
// is turned into a run-length of 3x0, 4x1, 0x2, 1x4, 0x4, 2x5,
// which then becomes 3 4 0 1 0 2.
//
// NB: size > 255 therefore means we need to repeatedly read to find
// the actual run length.
// Alternatively we could bit-encode instead of byte encode, eg BETA.
static int store_array(unsigned char *out, unsigned int *array, int size) {
    int i, j, k;
    for (i = j = k = 0; i < size; j++) {
	int run_len = i;
	while (i < size && array[i] == j)
	    i++;
	run_len = i-run_len;

	int r;
	do {
	    r = MIN(255, run_len);
	    out[k++] = r;
	    run_len -= r;
	} while (r == 255);
    }
    while (i < size)
	out[k++] = 0, j++;

    // RLE on out.
    //    1 2 3 3 3 3 3 4 4    5
    // => 1 2 3 3 +3... 4 4 +0 5
    int last = -1;
    for (i = j = 0; j < k; i++) {
	out[i] = out[j++];
	if (out[i] == last) {
	    int n = j;
	    while (j < k && out[j] == last)
		j++;
	    out[++i] = j-n;
	} else {
	    last = out[i];
	}
    }
    k = i;

//    fprintf(stderr, "Store_array %d => %d {", size, k);
//    for (i = 0; i < k; i++)
//	fprintf(stderr, "%d,", out[i]);
//    fprintf(stderr, "}\n");
    return k;
}

static int read_array(unsigned char *in, unsigned int *array, int size) {
#if 0
    unsigned char A[1024];
    int i, j, k, last = -1, nb = 0;

    size = MIN(1024, size);

    // Remove level one of run-len encoding
    for (i = j = k = 0; k < size; i++) {
	k += in[i];
	A[j++] = in[i];
	int run = in[i];
	if (in[i] == last) {
	    int r = in[++i];
	    while (r-- && k < size) {
		A[j++] = run;
		k += run;
	    }
	}
	last = run;
    }
    nb = i;

    // Now expand inner level of run-length encoding
    in = A;
    for (i = j = k = 0; j < size; i++) {
	int run_len = 0;
	int run_part;
	do {
	    run_part = in[k++];
	    run_len += run_part;
	} while (run_part == 255);

	while (run_len && j < size)
	    run_len--, array[j++] = i;
    }

    return nb;
#else
    int i, j, k, last = -1, r2 = 0;

    for (i = j = k = 0; j < size; i++) {
	int run_len;
	if (r2) {
	    run_len = last;
	} else {
	    run_len = 0;
	    int r, loop = 0;
	    do {
		r = in[k++];
		//fprintf(stderr, "%d,", r);
		if (++loop == 3)
		    run_len += r*255, r=255;
		else
		    run_len += r;
	    } while (r == 255);
	}
	if (r2 == 0 &&  run_len == last) {
	    r2 = in[k++];
	    //fprintf(stderr, "%d,", r2);
	} else {
	    if (r2) r2--;
	    last = run_len;
	}

	while (run_len && j < size)
	    run_len--, array[j++] = i;
    }
    //fprintf(stderr, "}\n");

    return k;
#endif
}

// FIXME: how to auto-tune these rather than trial and error?
static int strat_opts[][10] = {
    {10, 5, 4,-1, 2, 1, 0, 9, 10, 14},  // basic options (level < 7)
    //{10, 5, 7, 0, 2, 1, 0, 15, 9, 14},  // e.g. HiSeq 2000
    {9, 5, 7, 0, 2, 0, 7, 15, 0, 14},  // e.g. HiSeq 2000, ~2.1% smaller than above
    {12, 6, 2, 0, 2, 3, 0, 9, 12, 14},  // e.g. MiSeq
    {12, 6, 0, 0, 0, 0, 0, 12, 0, 0},   // e.g. IonTorrent; adaptive O1
    {0, 0, 0, 0, 0,  0, 0, 0, 0, 0},    // custom
};

#ifdef __SSE__
#   include <xmmintrin.h>
#else
#   define _mm_prefetch(a,b)
#endif

unsigned char *compress_block_fqz2f(int vers,
				    int level,
				    cram_slice *s,
				    unsigned char *in,
				    size_t in_size,
				    size_t *out_size) {
    //approx sqrt(delta), must be sequential
    int dsqr[] = {
	0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
	6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
    };
    unsigned int qtab[256]  = {0};
    unsigned int ptab[1024] = {0};
    unsigned int dtab[256]  = {0};
    unsigned int stab[256]  = {0};

    int comp_idx = 0;
    size_t i, j;
    ssize_t rec = 0;
    unsigned int last = 0, qlast = 0;
    RangeCoder rc;
    unsigned char q1 = 1;

    int strat = vers>>8;
    if (strat > 4) strat = 4;
    vers &= 0xff;

    unsigned char *comp = (unsigned char *)malloc(in_size*1.1+1000);
    if (!comp)
	return NULL;

    // Compute quality length per sequence.
    // This isn't just s->crecs[i].len as some times we emit extra QS
    // records, eg for feature code B.  Instead look at the .qual field.
    unsigned int *q_len = (unsigned int *)malloc(s->hdr->num_records * sizeof(*q_len));
    if (!q_len) {
	free(comp);
	return NULL;
    }
    for (i = 0; i < s->hdr->num_records; i++) {
	q_len[i] = (i < s->hdr->num_records-1)
	    ? s->crecs[i+1].qual - s->crecs[i].qual
	    : in_size - s->crecs[i].qual;
    }

#define NP 128
    uint32_t qhist[256] = {0}, nsym, max_sym;
    uint32_t qhist1[256][NP] = {{0}}, qhist2[256][NP] = {{0}};
    int do_dedup = 0;
    int do_rev = (vers==3); // V3.0 doesn't store qual in original orientation
    for (i = 0; i < in_size; i++)
	qhist[in[i]]++;

    // Compute possible strand dependence and whether dedup helps.
    int dir = 0;
    int last_len = 0;
    uint64_t t1[NP] = {0}, t2[NP] = {0};
    for (rec = i = j = 0; i < in_size; i++, j--) {
	if (j == 0) {
	    if (rec < s->hdr->num_records) {
		j = q_len[rec];
		dir = s->crecs[rec].flags & BAM_FREAD2 ? 1 : 0;
		if (i > 0 && j == last_len
		    && !memcmp(in+i-last_len, in+i, j))
		    do_dedup++; // cache which records are dup?
		last_len = j;
	    } else {
		j = in_size - i;
		dir = 0;
	    }
	    rec++;
	}
	if (dir)
	    qhist2[in[i]][j&(NP-1)]++, t2[j&(NP-1)]++;
	else
	    qhist1[in[i]][j&(NP-1)]++, t1[j&(NP-1)]++;
    }
    do_dedup = ((rec+1)/(do_dedup+1) < 500);
    rec = 0;
    last_len = 0;

    double e2 = 0;
    for (j = 0; j < NP; j++) {
	if (!t1[j] || !t2[j]) continue;
	for (i = 0; i < 256; i++) {
	    if (!qhist[i]) continue;
	    e2 += qhist1[i][j] > qhist2[i][j]
		? qhist1[i][j]/(t1[j]+1.) - qhist2[i][j]/(t2[j]+1.)
		: qhist2[i][j]/(t2[j]+1.) - qhist1[i][j]/(t1[j]+1.);
	}
    }
    int do_strand = (e2 > 77);

    for (i = max_sym = nsym = 0; i < 256; i++)
	if (qhist[i])
	    max_sym = i, nsym++;

    int store_qmap = (nsym <= 8 && nsym*2 < max_sym);

    // Check for fixed length.
    int first_len = q_len[0];
    for (i = 1; i < s->hdr->num_records; i++) {
	if (q_len[i] != first_len)
	    break;
    }
    int fixed_len = (i == s->hdr->num_records);

    int store_qtab = 0; // unused by current encoder
    int q_qctxbits = strat_opts[strat][0];
    int q_qctxshift= strat_opts[strat][1];
    int q_pctxbits = strat_opts[strat][2];
    int q_pctxshift= strat_opts[strat][3];
    int q_dctxbits = strat_opts[strat][4];
    int q_dctxshift= strat_opts[strat][5];
    int q_qloc     = strat_opts[strat][6];
    int q_sloc     = strat_opts[strat][7];
    int q_ploc     = strat_opts[strat][8];
    int q_dloc     = strat_opts[strat][9];

    if (strat >= 4)
	goto manually_set; // used in TEST_MAIN for debugging

    if (q_pctxshift < 0)
	q_pctxshift = MAX(0, log((double)s->crecs[0].len/(1<<q_pctxbits))/log(2)+.5);

    if (nsym <= 4) {
	// NovaSeq
	q_qctxshift=2; // qmax 64, although we can store up to 128 if needed
	if (in_size < 5000000) {
	    q_pctxbits =2;
	    q_pctxshift=5;
	}
    } else if (nsym <= 8) {
	// HiSeqX
	q_qctxbits =MIN(q_qctxbits,9);
	q_qctxshift=3;
	if (in_size < 5000000)
	    q_qctxbits =6;
    }

    if (in_size < 300000) {
	q_qctxbits=q_qctxshift;
	q_dctxbits=2;
    }

 manually_set:
    // If nsym is significant lower than max_sym, store a lookup table and
    // transform prior to compression.  Eg 4 or 8 binned data.
    comp[comp_idx++] = 5; // FMT 5: this code
    comp[comp_idx++] =
	(store_qtab<<7)|
	((q_dctxbits>0)<<6)|
	((q_pctxbits>0)<<5)|
	(do_rev<<4)|
	(do_strand<<3)|
	(fixed_len<<2)|
	(do_dedup<<1)|
	(store_qmap);
    comp[comp_idx++] = max_sym;

//    fprintf(stderr, "strat %d 0x%x%x%x%x%x%x%x%x%x%x\n", strat,
//	    q_qctxbits, q_qctxshift,
//	    q_pctxbits, q_pctxshift,
//	    q_dctxbits, q_dctxshift,
//	    q_qloc    , q_sloc,
//	    q_ploc    , q_dloc);

    for (i = 0; i < sizeof(dsqr)/sizeof(*dsqr); i++)
	if (dsqr[i] > (1<<q_dctxbits)-1)
	    dsqr[i] = (1<<q_dctxbits)-1;

    comp[comp_idx++] = (q_qctxbits<<4)|q_qctxshift;
    comp[comp_idx++] = (q_qloc<<4)|q_sloc;
    comp[comp_idx++] = (q_ploc<<4)|q_dloc;

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

    // Produce ptab from pctxshift.
    if (q_qctxbits) {
	for (i = 0; i < 256; i++) {
	    qtab[i] = i;
	    // Alternative mappings:
	    //qtab[i] = i<max_sym ? i : max_sym;
	    //qtab[i] = i > 30 ? MIN(max_sym,i)-15 : i/2;  // eg for 9827 BAM
	    //qtab[i] = pow(i,1.5)*0.12;
	}

	// Eg map Q8 to Q4 and then compress with 2 bit instead of 3.
//	int j = 0;
//	qtab[0] = j;
//	qtab[1] = j++;
//	qtab[2] = j;
//	qtab[3] = j++;
//	qtab[4] = j;
//	qtab[5] = j++;
//	qtab[6] = j;
//	q_qctxshift = 2;

	// custom qtab
	if (store_qtab)
	    comp_idx += store_array(comp+comp_idx, qtab, 256);
    }

    if (q_pctxbits) {
	for (i = 0; i < 1024; i++)
	    ptab[i] = MIN((1<<q_pctxbits)-1, i>>q_pctxshift);

//	// Eg a particular novaseq run, 2.8% less, -x 0xa2302108ae
//	i = j = 0;
//	while (i <= 10)  {ptab[i++] = j;}  j++;
//	while (i <= 12)  {ptab[i++] = j;}  j++;
//	while (i <= 22)  {ptab[i++] = j;}  j++;
//	while (i <= 24)  {ptab[i++] = j;}  j++;
//	while (i <= 27)  {ptab[i++] = j;}  j++;
//	while (i <= 37)  {ptab[i++] = j;}  j++;
//	while (i <= 83)  {ptab[i++] = j;}  j++;
//	while (i <= 95)  {ptab[i++] = j;}  j++;
//	while (i <= 127) {ptab[i++] = j;}  j++;
//	while (i <= 143) {ptab[i++] = j;}  j++;
//	while (i <= 1023){ptab[i++] = j;}  j++;

	// custom ptab
	comp_idx += store_array(comp+comp_idx, ptab, 1024);

	for (i = 0; i < 1024; i++)
	    ptab[i] <<= q_ploc;
    }

    if (q_dctxbits) {
	for (i = 0; i < 256; i++)
	    dtab[i] = dsqr[MIN(sizeof(dsqr)/sizeof(*dsqr)-1, i>>q_dctxshift)];

	// custom dtab
	comp_idx += store_array(comp+comp_idx, dtab, 256);

	for (i = 0; i < 256; i++)
	    dtab[i] <<= q_dloc;
    }

    if (do_strand)
	stab[1] = 1<<q_sloc;



    SIMPLE_MODEL(QMAX,_) *model_qual;

    model_qual = malloc(sizeof(*model_qual) * (1<<16));
    if (!model_qual)
	return NULL;

    for (i = 0; i < (1<<16); i++)
	SIMPLE_MODEL(QMAX,_init)(&model_qual[i], max_sym+1);

    SIMPLE_MODEL(256,_) model_len[4];
    for (i = 0; i < 4; i++)
	SIMPLE_MODEL(256,_init)(&model_len[i],256);

    SIMPLE_MODEL(2,_) model_revcomp;
    SIMPLE_MODEL(2,_init)(&model_revcomp,2);

    SIMPLE_MODEL(2,_) model_strand;
    SIMPLE_MODEL(2,_init)(&model_strand,2);

    int delta = 0;
    SIMPLE_MODEL(2,_) model_dup;
    SIMPLE_MODEL(2,_init)(&model_dup,2);

    RC_SetOutput(&rc, (char *)comp+comp_idx);
    RC_StartEncode(&rc);

    int ndup0 = 0, ndup1 = 0;

    if (do_rev) {
	// Pass 1, reverse all seqs if necessary.
	i = rec = j = 0;
	while (i < in_size) {
	    int len = rec < s->hdr->num_records-1
		? q_len[rec] : in_size - i;

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
    int read2 = 0;
    for (i = j = 0; i < in_size; i++, j--) {
	if (j == 0) {
	    // Quality buffer maybe longer than sum of reads if we've
	    // inserted a specific base + quality pair.
	    int len = rec < s->hdr->num_records
		? q_len[rec] : in_size - s->crecs[s->hdr->num_records-1].qual;
	    if (!fixed_len || rec == 0) {
		SIMPLE_MODEL(256,_encodeSymbol)(&model_len[0], &rc, (len>> 0) & 0xff);
		SIMPLE_MODEL(256,_encodeSymbol)(&model_len[1], &rc, (len>> 8) & 0xff);
		SIMPLE_MODEL(256,_encodeSymbol)(&model_len[2], &rc, (len>>16) & 0xff);
		SIMPLE_MODEL(256,_encodeSymbol)(&model_len[3], &rc, (len>>24) & 0xff);
	    }

	    if (do_rev) {
		// no need to reverse complement for V4.0 as the core format
		// already has this feature.
		if (s->crecs[rec].flags & BAM_FREVERSE)
		    SIMPLE_MODEL(2,_encodeSymbol)(&model_revcomp, &rc, 1);
		else
		    SIMPLE_MODEL(2,_encodeSymbol)(&model_revcomp, &rc, 0);
	    }

	    if (do_strand) {
		read2 = (s->crecs[rec].flags & BAM_FREAD2) ? 0 : 1;
		SIMPLE_MODEL(2,_encodeSymbol)(&model_strand, &rc, read2);
	    }

	    rec++;
	    j = len;
	    delta = last = qlast = q1 = 0;

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
	SIMPLE_MODEL(QMAX,_encodeSymbol)(&model_qual[last], &rc, qhist[q]);

	qlast = (qlast<<q_qctxshift) + qtab[qhist[q]];
	last = (qlast & ((1<<q_qctxbits)-1)) << q_qloc;
	last += ptab[MIN(j,1024)]; //limits max pos
	last += stab[read2];
	last += dtab[delta];

	last &= 0xffff;
        //_mm_prefetch((const char *)&model_qual[last], _MM_HINT_T0);

	delta += (q1 != q);
	q1 = q;
    }

    RC_FinishEncode(&rc);

    if (do_rev) {
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

    *out_size = comp_idx + RC_OutSize(&rc);

//    fprintf(stderr, "%d / %d %d %d %d %d / %d %d %d %d %d %d %d %d %d %d = %d to %d\n",
//	    nsym, do_rev, do_strand, fixed_len, do_dedup, store_qmap,
//	    q_qctxbits,
//	    q_qctxshift,
//	    q_pctxbits,
//	    q_pctxshift,
//	    q_dctxbits,
//	    q_dctxshift,
//	    q_qloc,
//	    q_sloc,
//	    q_ploc,
//	    q_dloc,
//	    (int)in_size, (int)*out_size);

    free(model_qual);
    free(q_len);

    return comp;
}

unsigned char *uncompress_block_fqz2f(cram_slice *s,
				      unsigned char *in,
				      size_t in_size,
				      size_t *out_size) {
    unsigned int qtab[256]  = {0};
    unsigned int ptab[1024] = {0};
    unsigned int dtab[256]  = {0};
    unsigned int stab[256]  = {0};

    unsigned char *uncomp = NULL;
    RangeCoder rc;
    size_t i, j, rec = 0, len = *out_size, in_idx = 0;
    unsigned char q1 = 1;
    unsigned int last = 0, qlast = 0;

    int vers       = in[in_idx++];
    if (vers != 5) {
	fprintf(stderr, "This version of fqzcomp only supports format 5\n");
	return NULL;
    }

    int flags      = in[in_idx++];
    int have_qtab  = flags&128;
    int have_dtab  = flags&64;
    int have_ptab  = flags&32;
    int do_rev     = flags&16;
    int do_strand  = flags&8;
    int fixed_len  = flags&4;
    int do_dedup   = flags&2;
    int store_qmap = flags&1;
    int max_sym    = in[in_idx++];

    int q_qctxbits = in[in_idx]>>4;
    int q_qctxshift= in[in_idx++]&15;
    int q_qloc     = in[in_idx]>>4;
    int q_sloc     = in[in_idx++]&15;
    int q_ploc     = in[in_idx]>>4;
    int q_dloc     = in[in_idx++]&15;

    int qmap[256];
    if (store_qmap) {
	int nsym = in[in_idx++];
	for (i = 0; i < nsym; i++)
	    qmap[i] = in[in_idx++];
	max_sym = nsym;
    } else {
	for (i = 0; i < 256; i++)
	    qmap[i] = i;
    }

    // Produce ptab from pctxshift.
    if (q_qctxbits) {
	if (have_qtab)
	    in_idx += read_array(in+in_idx, qtab, 256);
	else
	    for (i = 0; i < 256; i++)
		qtab[i] = i;
    }

    if (have_ptab)
	in_idx += read_array(in+in_idx, ptab, 1024);
    for (i = 0; i < 1024; i++)
	ptab[i] <<= q_ploc;

    if (have_dtab)
        in_idx += read_array(in+in_idx, dtab, 256);
    for (i = 0; i < 256; i++)
	dtab[i] <<= q_dloc;

    if (do_strand)
	stab[1] = 1<<q_sloc;

    SIMPLE_MODEL(QMAX,_) *model_qual;
    model_qual = malloc(sizeof(*model_qual) * (1<<16));
    if (!model_qual)
	return NULL;
    for (i = 0; i < (1<<16); i++)
	SIMPLE_MODEL(QMAX,_init)(&model_qual[i],max_sym+1);

    SIMPLE_MODEL(256,_) model_len[4];
    for (i = 0; i < 4; i++)
	SIMPLE_MODEL(256,_init)(&model_len[i],256);

    SIMPLE_MODEL(2,_) model_revcomp;
    SIMPLE_MODEL(2,_init)(&model_revcomp,2);

    SIMPLE_MODEL(2,_) model_strand;
    SIMPLE_MODEL(2,_init)(&model_strand,2);

    int delta = 0, rev = 0;

    SIMPLE_MODEL(2,_) model_dup;
    SIMPLE_MODEL(2,_init)(&model_dup,2);

    uncomp = (unsigned char *)malloc(*out_size);
    if (!uncomp)
	return NULL;

    RC_SetInput(&rc, (char *)in+in_idx);
    RC_StartDecode(&rc);

    int nrec = 1000;
    char *rev_a = malloc(nrec);
    int *len_a = malloc(nrec * sizeof(int));
    if (!rev_a || !len_a)
	return NULL;

    int last_len = 0, read2 = 0;
    for (rec = i = j = 0; i < len; i++, j--) {
	if (rec >= nrec) {
	    nrec *= 2;
	    rev_a = realloc(rev_a, nrec);
	    len_a = realloc(len_a, nrec*sizeof(int));
	    if (!rev_a || !len_a)
		return NULL;
	}

	if (j == 0) {
	    int len = last_len;
	    if (!fixed_len || rec == 0) {
		len  = SIMPLE_MODEL(256,_decodeSymbol)(&model_len[0], &rc);
		len |= SIMPLE_MODEL(256,_decodeSymbol)(&model_len[1], &rc)<<8;
		len |= SIMPLE_MODEL(256,_decodeSymbol)(&model_len[2], &rc)<<16;
		len |= SIMPLE_MODEL(256,_decodeSymbol)(&model_len[3], &rc)<<24;
		last_len = len;
	    }
	    if (do_rev) {
		rev = SIMPLE_MODEL(2,_decodeSymbol)(&model_revcomp, &rc);
		rev_a[rec] = rev;
		len_a[rec] = len;
	    }
	    if (do_strand)
		read2 = SIMPLE_MODEL(2,_decodeSymbol)(&model_strand, &rc);

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
	    delta = last = qlast = q1 = 0;
	}

	unsigned char q, Q;
	
	Q = SIMPLE_MODEL(QMAX,_decodeSymbol)(&model_qual[last], &rc);
	q = qmap[Q];

	qlast = (qlast<<q_qctxshift) + qtab[Q];
	last = (qlast & ((1<<q_qctxbits)-1)) << q_qloc;
	last += ptab[MIN(j,1024)]; //limits max pos
	last += stab[read2];
	last += dtab[delta];

	last &= 0xffff;

	delta += (q1 != q) * (delta<255);
	q1 = q;
        uncomp[i] = q;
    }
    rev_a[rec] = rev;
    len_a[rec] = len;

    if (do_rev) {
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
    free(model_qual);
    free(rev_a);
    free(len_a);

    return uncomp;
}


#ifndef TEST_MAIN
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

#else // TEST_MAIN
#include <unistd.h>

#define BS 1024*1024
static unsigned char *load(char *fn, size_t *lenp) {
    unsigned char *data = NULL;
    uint64_t dsize = 0;
    uint64_t dcurr = 0;
    signed int len;

    build_rcp_freq();

    int fd = open(fn, O_RDONLY);
    if (!fd) {
	perror(fn);
	return NULL;
    }

    do {
	if (dsize - dcurr < BS) {
	    dsize = dsize ? dsize * 2 : BS;
	    data = realloc(data, dsize);
	}

	len = read(fd, data + dcurr, BS);
	if (len > 0)
	    dcurr += len;
    } while (len > 0);

    if (len == -1) {
	perror("read");
    }
    close(fd);

    *lenp = dcurr;
    return data;
}

#define BLK_SIZE 100*1000000

int main(int argc, char **argv) {
    unsigned char *in, *out;
    size_t in_len, out_len;
    int decomp = 0, vers = 4;

    if (argc > 1 && strcmp(argv[1], "-d") == 0) {
	decomp = 1;
	argv++;
	argc--;
    }
    while (argc > 1 && strcmp(argv[1], "-b") == 0) {
	vers += 256;
	argv++;
	argc--;
    }
    while (argc > 2 && strcmp(argv[1], "-x") == 0) {
	uint64_t x = strtol(argv[2], NULL, 0);
	argv+=2;
	argc-=2;
	int z;
	for (z=0; z<10; z++) {
	    strat_opts[4][9-z] = x&15;
	    x>>=4;
	}
	vers = 4*256 + 4;
    }
    in = load(argc > 1 ? argv[1] : "/dev/stdin", &in_len);
    if (!in)
	exit(1);

    int rec_len = 100;
    if (argc > 2)
	rec_len = atoi(argv[2]);
    int blk_size = BLK_SIZE; // MAX
    if (argc > 3)
	blk_size = atoi(argv[3]);
    if (blk_size > BLK_SIZE)
	blk_size = BLK_SIZE;

    if (decomp) {
	rec_len = *(int *)in;
	unsigned char *in2 = in+4;
	in_len -= 4;
	while (in_len > 0) {
	    size_t out_len = *(size_t *)in2;
	    size_t in2_len = *(size_t *)(in2+8);
	    in2 += 16;
	    out = uncompress_block_fqz2f(fake_slice(out_len, rec_len),
					 in2, in_len-16, &out_len);
	    if (write(1, out, out_len) < 0) return 1;
	    free(out);
	    in2 += in2_len;
	    in_len -= in2_len+16;
	}
    } else {
	unsigned char *in2 = in;
	long t_out = 0;
	if (write(1, &rec_len, 4) < 0) return 1;
	out = NULL;
	while (in_len > 0) {
	    size_t in2_len = in_len <= blk_size ? in_len : blk_size;
	    out = compress_block_fqz2f(vers, 0, fake_slice(in2_len, rec_len),
				       in2, in2_len, &out_len);
	    //fprintf(stderr, "%d to %d\n", (int)in2_len, (int)out_len);
	    if (write(1, &in2_len, 8)  < 0) return 1;
	    if (write(1, &out_len, 8)  < 0) return 1;
	    if (write(1, out, out_len) < 0) return 1;
	    in_len -= in2_len;
	    in2 += in2_len;
	    t_out += out_len+16;
	}
	free(out);
	fprintf(stderr, "Total output = %ld\n", t_out);
    }

    free(in);

    return 0;
}
#endif

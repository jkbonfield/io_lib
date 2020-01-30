/*
 * Copyright (c) 2013 Genome Research Ltd.
 * Author(s): James Bonfield
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 * 
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *    Institute nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific
 *    prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <ctype.h>

#include "io_lib/cram.h"
#include "io_lib/os.h"

cram_stats *cram_stats_create(void) {
    return calloc(1, sizeof(cram_stats));
}

void cram_stats_add(cram_stats *st, int64_t val) {
    st->nsamp++;

    if (val < MAX_STAT_VAL && val >= 0) {
	st->freqs[val]++;
    } else {
	HashItem *hi;

	if (!st->h) {
	    st->h = HashTableCreate(2048, HASH_DYNAMIC_SIZE|HASH_NONVOLATILE_KEYS|HASH_INT_KEYS);
	}

	if ((hi = HashTableSearch(st->h, (char *)(size_t)val, 8))) {
	    hi->data.i++;
	} else {
	    HashData hd;
	    hd.i = 1;
	    HashTableAdd(st->h, (char *)(size_t)val, 8, hd, NULL);
	}
    }
}

void cram_stats_del(cram_stats *st, int64_t val) {
    st->nsamp--;

    if (val < MAX_STAT_VAL && val >= 0) {
	st->freqs[val]--;
	assert(st->freqs[val] >= 0);
    } else if (st->h) {
	HashItem *hi;

	if ((hi = HashTableSearch(st->h, (char *)(size_t)val, 8))) {
	    if (--hi->data.i == 0)
		HashTableDel(st->h, hi, 0);
	} else {
	    fprintf(stderr, "Failed to remove val %"PRId64" from cram_stats\n", val);
	    st->nsamp++;
	}
    } else {
	fprintf(stderr, "Failed to remove val %"PRId64" from cram_stats\n", val);
	st->nsamp++;
    }
}

void cram_stats_dump(cram_stats *st) {
    int i;
    fprintf(stderr, "cram_stats:\n");
    for (i = 0; i < MAX_STAT_VAL; i++) {
	if (!st->freqs[i])
	    continue;
	fprintf(stderr, "\t%d\t%d\n", i, st->freqs[i]);
    }
    if (st->h) {
	HashIter *iter=  HashTableIterCreate();
	HashItem *hi;

	while ((hi = HashTableIterNext(st->h, iter))) {
	    fprintf(stderr, "\t%d\t%d\n", (int)(size_t)hi->key,
		    (int)hi->data.i);
	}
	HashTableIterDestroy(iter);
    }
}

/*
 * Computes entropy from integer frequencies for various encoding methods and
 * picks the best encoding.
 *
 * FIXME: we could reuse some of the code here for the actual encoding
 * parameters too. Eg the best 'k' for SUBEXP or the code lengths for huffman.
 *
 * Returns the best codec to use.
 */
enum cram_encoding cram_stats_encoding(cram_fd *fd, cram_stats *st) {
    int nvals, i, ntot = 0, max_val = 0, min_val = INT_MAX;
    int *vals = NULL, *freqs = NULL, vals_alloc = 0;

    //cram_stats_dump(st);

    /* Count number of unique symbols */
    for (nvals = i = 0; i < MAX_STAT_VAL; i++) {
	if (!st->freqs[i])
	    continue;
	if (nvals >= vals_alloc) {
	    vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
	    vals  = realloc(vals,  vals_alloc * sizeof(int));
	    freqs = realloc(freqs, vals_alloc * sizeof(int));
	    if (!vals || !freqs) {
		if (vals)  free(vals);
		if (freqs) free(freqs);
		return E_HUFFMAN; // Cannot do much else atm
	    }
	}
	vals[nvals] = i;
	freqs[nvals] = st->freqs[i];
	ntot += freqs[nvals];
	if (max_val < i) max_val = i;
	if (min_val > i) min_val = i;
	nvals++;
    }
    if (st->h) {
	HashIter *iter=  HashTableIterCreate();
	HashItem *hi;
	int i;

	while ((hi = HashTableIterNext(st->h, iter))) {
	    if (nvals >= vals_alloc) {
		vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
		vals  = realloc(vals,  vals_alloc * sizeof(int));
		freqs = realloc(freqs, vals_alloc * sizeof(int));
		if (!vals || !freqs)
		    return E_HUFFMAN; // Cannot do much else atm
	    }
	    i = (size_t)hi->key;
	    vals[nvals]=i;
	    freqs[nvals] = hi->data.i;
	    ntot += freqs[nvals];
	    if (max_val < i) max_val = i;
	    if (min_val > i) min_val = i;
	    nvals++;
	}
	HashTableIterDestroy(iter);
    }

    st->nvals = nvals;
    assert(ntot == st->nsamp);
    free(vals);
    free(freqs);

    // Crude and simple alternative.
    return nvals > 1 ? E_EXTERNAL : E_HUFFMAN;


#ifdef RANDOMISER
    // RANDOMISER
    switch(random()%10) {
    case 0:  return E_HUFFMAN;
    case 1:  return E_HUFFMAN;
    //case 1:  return E_BETA; // Java doesn't support E_BETA for BYTE vals
    default: return E_EXTERNAL;
    }
#endif

    // Single value items are best served in HUFFMAN as this takes
    // zero bits to store (it's only needs the compression header).
    if (nvals <= 1) {
	if (fd->verbose > 1)
	    fprintf(stderr, "0 values => 0 bits\n");

	return E_HUFFMAN;
    }

    if (fd->verbose > 1)
	fprintf(stderr, "Range = %d..%d, nvals=%d, ntot=%d\n",
		min_val, max_val, nvals, ntot);

    // Massive cull of old entropy analysis.
    // If you're interested in code to compare huffman vs beta
    // bit streams, see svn revision r4051
    return E_EXTERNAL;
}

void cram_stats_free(cram_stats *st) {
    if (st->h)
	HashTableDestroy(st->h, 0);
    free(st);
}

/*
 * Analyse c->slice[?]->qual_blk to identify appropriate
 * quality encoding method.
 * Does it benefit from bit packing?
 * Does it benefit from RLE pre-packing?
 * Does it benefit from RLE post-packing?
 */
static void hist(unsigned char *data, int dlen, int hist[256]) {
    // FIXME: speed up; see hist4 in rans
    int i;
    for (i = 0; i < dlen; i++)
	hist[data[i]]++;
}

void cram_stats_qual(cram_container *c,
		     int *nval, int val[256],
		     int *nrle, int rle[256]) {
    int i, last = -1;

    // Number of distinct symbols
    // We end up with val[256] being a lookup table to go from A,B,C to 0,1,2
    memset(val, 0, 256*sizeof(*val));
    for (i = 0; i < c->curr_slice; i++)
	hist(c->slice[i].qual_blk->data, c->slice[i].qual_blk->byte, val);

    *nval = 0;
    unsigned char map[256];
    for (i = 0; i < 256; i++)
	if (val[i])
	    map[i] = (*nval)++;

    // Their RLE costs (packed or unpacked).
    // We end up with rle[256] being a lookup table of symbols to RLE
    for (i = 0; i < 256; i++) rle[i] = -99;
    if (*nval <= 1) {
	// NOP: it's a constant.
    } else if (*nval <= 2) {
	// 8 per byte
	unsigned char *d = c->slice[0].qual_blk->data;
	int len = c->slice[0].qual_blk->byte;
	for (i = 0; i < len-7; i += 8) {
	    unsigned char p =
		(map[d[i+0]]<<7) +
		(map[d[i+1]]<<6) +
		(map[d[i+2]]<<5) +
		(map[d[i+3]]<<4) +
		(map[d[i+4]]<<3) +
		(map[d[i+5]]<<2) +
		(map[d[i+6]]<<1) +
		(map[d[i+7]]);
	    rle[p] += (p == last) ? 1 : -2;
	    last = p;
	}
    } else if (*nval <= 4) {
	// 4 per byte
	unsigned char *d = c->slice[0].qual_blk->data;
	int len = c->slice[0].qual_blk->byte;
	for (i = 0; i < len-3; i += 4) {
	    unsigned char p =
		(map[d[i+0]]<<6) +
		(map[d[i+1]]<<4) +
		(map[d[i+2]]<<2) +
		(map[d[i+3]]);
	    rle[p] += (p == last) ? 1 : -2;
	    last = p;
	}
    } else if (*nval <= 16) {
	// 2 per byte
	unsigned char *d = c->slice[0].qual_blk->data;
	int len = c->slice[0].qual_blk->byte;
	for (i = 0; i < len-1; i += 2) {
	    unsigned char p = (map[d[i+0]]<<4) + map[d[i+1]];
	    rle[p] += (p == last) ? 1 : -2;
	    last = p;
	}
    } else {
	// unpacked
	unsigned char *d = c->slice[0].qual_blk->data;
	int len = c->slice[0].qual_blk->byte;
	for (i = 0; i < len; i++) {
	    rle[d[i]] += (d[i] == last) ? 1 : -2;
	    last = d[i];
	}
    }
//    for (i = 0; i < 256; i++)
//	if (rle[i] > 0)
//	    fprintf(stderr, "rle[%d]=%d\n", i, rle[i]);

    *nrle = 0;
    for (i = 0; i < 256; i++) {
	if (rle[i] > 0)
	    (*nrle)++;
    }
}

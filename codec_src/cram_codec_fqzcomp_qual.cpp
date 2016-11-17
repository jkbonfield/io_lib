#define __STDC_LIMIT_MACROS 1

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "io_lib/cram_block_compression.h"
#include "clr.cdr"
#include "simple_model.h"

// Comment out if you wish to not de-dup the last quality string.
#define DEDUP

static const char *name(void) {
    return "fqzcomp-qual";
}

#ifndef QMAX
#  define QMAX 64
#endif
#define QMAX_B 9
#define QBITS 12
#define QSIZE (1<<QBITS)

// unsigned char *compress_block(int level,
// 			      cram_slice *s,
//                               unsigned char *in,
//                               size_t in_size,
//                               size_t *out_size) {
//     unsigned char *comp = (unsigned char *)malloc(in_size*1.1+1000);
//     size_t i;
//     unsigned int last = 0;
//     RangeCoder rc;
//     unsigned char q1 = 0, q2 = 0, q3 = 0;
//     SIMPLE_MODEL<QMAX> *model_qual;
//     model_qual = new SIMPLE_MODEL<QMAX>[QSIZE];
//     int rle = 0;
// 
//     rc.output((char *)comp);
//     rc.StartEncode();
// 
// #if 0
//     // Reverse complement
//     int rec = 0;
//     i = 0;
//     while (i < in_size) {
// 	int len = s->crecs[rec].len;
// 	if (s->crecs[rec].flags & BAM_FREVERSE) {
// 	    int I, J;
// 	    unsigned char *cp = in+i;
// 	    for (I = 0, J = len-1; I < J; I++, J--) {
// 		unsigned char c;
// 		c = cp[I];
// 		cp[I] = cp[J];
// 		cp[J] = c;
// 	    }
// 	}
// 	rec++;
// 	i+=len;
//     }
// #endif
// 
//     //fprintf(stderr, "%.1000s\n", in);
// 
//     for (i = 0; i < in_size; i++) {
// 	unsigned char q = in[i] < QMAX ? in[i] : QMAX-1;
//         //unsigned char q = in[i] & (QMAX-1);
// 	//assert(in[i] < QMAX && in[i] >= 0);
//         model_qual[last].encodeSymbol(&rc, q);
//         //last = ((MAX(q1, q2)<<6) + q) & (QSIZE-1);
//         last = (((q1>>2)<<6) + q) & (QSIZE-1); // quicker learning
// 
//         //last += (q2==q3) << QBITS;
//         //q3 = q2;
// 
//         //_mm_prefetch((const char *)&model_qual[last], _MM_HINT_T0);
//         q2 = q1; q1 = q;
//     }
// 
//     model_qual[last].encodeSymbol(&rc, QMAX-1); /* terminator */
//     rc.FinishEncode();
// 
//     *out_size = rc.size_out();
// 
//     //fprintf(stdout, "%d -> %d\n", (int)in_size, (int)*out_size);
// 
//     delete[] model_qual;
// 
//     return comp;
// }
// 
// unsigned char *uncompress_block(cram_slice *s,
// 				unsigned char *in,
//                                 size_t in_size,
//                                 size_t *out_size) {
//     int block_size, data_size;
//     unsigned char *uncomp = NULL;
//     RangeCoder rc;
//     size_t i, len = *out_size;
//     unsigned char q1 = 0, q2 = 0;
//     unsigned int last = 0;
//     SIMPLE_MODEL<QMAX> *model_qual;
//     model_qual = new SIMPLE_MODEL<QMAX>[QSIZE];
// 
//     uncomp = (unsigned char *)malloc(*out_size);
// 
//     rc.input((char *)in);
//     rc.StartDecode();
// 
//     for (i = 0; i < len; i++) {
//         unsigned char q = model_qual[last].decodeSymbol(&rc);
//         uncomp[i] = q;
// 
//         //last = ((MAX(q1, q2)<<6) + q) & (QSIZE-1);
//         last = (((q1>>2)<<6) + q) & (QSIZE-1); // quicker learning
// 
//         q2 = q1; q1 = q;
//     }
// 
//     rc.FinishDecode();
//     delete[] model_qual;
// 
//     return uncomp;
// }
// 
// // Fqzcomp -q2 equiv
// unsigned char *compress_block_fqz2(int level,
// 				   cram_slice *s,
// 				   unsigned char *in,
// 				   size_t in_size,
// 				   size_t *out_size) {
//     unsigned char *comp = (unsigned char *)malloc(in_size*1.1+1000);
//     size_t i, j, rec = 0;
//     unsigned int last = 0;
//     RangeCoder rc;
//     unsigned char q1 = 0, q2 = 0, q3 = 0;
//     SIMPLE_MODEL<QMAX> *model_qual;
//     model_qual = new SIMPLE_MODEL<QMAX>[QSIZE*16];
//     SIMPLE_MODEL<256>  model_len[4];
//     int rle = 0, delta = 5;
// 
//     rc.output((char *)comp);
//     rc.StartEncode();
// 
//     for (i = j = 0; i < in_size; i++, j--) {
// 	if (j == 0) {
// 	    int len = s->crecs[rec].len;
// 	    model_len[0].encodeSymbol(&rc, (len>> 0) & 0xff);
// 	    model_len[1].encodeSymbol(&rc, (len>> 8) & 0xff);
// 
// 	    rec++;
// 	    j = len;
// 	    delta = 5;
// 	}
// 
//         unsigned char q = in[i] & (QMAX-1);
// 	assert(in[i] < QMAX && in[i] >= 0);
//         model_qual[last & (QSIZE*16-1)].encodeSymbol(&rc, q);
//         last = ((MAX(q1, q2)<<6) + q) & (QSIZE-1);
// 	last  += (q1==q2) << QBITS;
// 
// 	delta += (q1>q)*(q1-q);
// 	last  += (MIN(7*8, delta)&0xf8) << (QBITS-2);
// 	//if (i%16 == 0) delta/=2;
// 
//         q2 = q1; q1 = q;
//     }
// 
//     model_qual[last].encodeSymbol(&rc, QMAX-1); /* terminator */
//     rc.FinishEncode();
// 
//     *out_size = rc.size_out();
// 
//     //fprintf(stdout, "%d -> %d\n", (int)in_size, (int)*out_size);
// 
//     delete[] model_qual;
// 
//     return comp;
// }
// 
// unsigned char *uncompress_block_fqz2(cram_slice *s,
// 				     unsigned char *in,
// 				     size_t in_size,
// 				     size_t *out_size) {
//     int block_size, data_size;
//     unsigned char *uncomp = NULL;
//     RangeCoder rc;
//     size_t i, j, rec = 0, len = *out_size;
//     unsigned char q1 = 0, q2 = 0;
//     unsigned int last = 0;
//     SIMPLE_MODEL<QMAX> *model_qual;
//     model_qual = new SIMPLE_MODEL<QMAX>[QSIZE*16];
//     SIMPLE_MODEL<256>  model_len[4];
//     int delta = 5;
// 
//     uncomp = (unsigned char *)malloc(*out_size);
// 
//     rc.input((char *)in);
//     rc.StartDecode();
// 
//     for (i = j = 0; i < len; i++, j--) {
// 	if (j == 0) {
// 	    int len;
// 	    len  = model_len[0].decodeSymbol(&rc);
// 	    len |= model_len[1].decodeSymbol(&rc)<<8;
// 	    rec++;
// 	    j = len;
// 	    delta = 5;
// 	}
// 
//         unsigned char q = model_qual[last].decodeSymbol(&rc);
//         uncomp[i] = q;
// 
//         last = ((MAX(q1, q2)<<6) + q) & (QSIZE-1);
// 	last  += (q1==q2) << QBITS;
// 
// 	delta += (q1>q)*(q1-q);
// 	last  += (MIN(7*8, delta)&0xf8) << (QBITS-2);
// 
// 	//if (i%16 == 0) delta/=2;
//         q2 = q1; q1 = q;
//     }
// 
//     rc.FinishDecode();
//     delete[] model_qual;
// 
//     return uncomp;
// }


// Fqzcomp -q2 equiv
unsigned char *compress_block_fqz2f(int level,
				    cram_slice *s,
				    unsigned char *in,
				    size_t in_size,
				    size_t *out_size) {
    unsigned char *comp = (unsigned char *)malloc(in_size*1.1+1000);
    size_t i, j;
    ssize_t rec = 0;
    unsigned int last = 0;
    RangeCoder rc;
    unsigned char q1 = 0, q2 = 0;
    SIMPLE_MODEL<QMAX> *model_qual;
    model_qual = new SIMPLE_MODEL<QMAX>[QSIZE*16];
    SIMPLE_MODEL<256>  model_len[4];
    SIMPLE_MODEL<2>  model_strand;
    int delta = 5;
#ifdef DEDUP
    int last_len = 0;
    SIMPLE_MODEL<2>  model_dup;
#endif

    rc.output((char *)comp);
    rc.StartEncode();


    for (i = j = 0; i < in_size; i++, j--) {
	if (j == 0) {
	    // Quality buffer maybe longer than sum of reads if we've
	    // inserted a specific base + quality pair.
	    // FIXME: how to handle this?
	    int len = rec < s->hdr->num_records-1
		? s->crecs[rec].len
		: in_size - i;
	    model_len[0].encodeSymbol(&rc, (len>> 0) & 0xff);
	    model_len[1].encodeSymbol(&rc, (len>> 8) & 0xff);

	    if (s->crecs[rec].flags & BAM_FREVERSE) {
		model_strand.encodeSymbol(&rc, 1);
		// Reverse complement sequence - note: modifies buffer
		int I,J;
		unsigned char *cp = in+i;
		for (I = 0, J = len-1; I < J; I++, J--) {
		    unsigned char c;
		    c = cp[I];
		    cp[I] = cp[J];
		    cp[J] = c;
		}
	    } else {
		model_strand.encodeSymbol(&rc, 0);
	    }

	    rec++;
	    j = len;
	    delta = 5;
	    last = 0; // reset last too?

#ifdef DEDUP
	    // Possible dup of previous read?
	    if (i && len == last_len && !memcmp(in+i-last_len, in+i, len)) {
		model_dup.encodeSymbol(&rc, 1);
		i += len-1;
		j = 1;
		continue;
	    }
	    model_dup.encodeSymbol(&rc, 0);

	    last_len = len;
#endif
	}

	unsigned char q = in[i] < QMAX ? in[i] : QMAX-1;
        //unsigned char q = in[i] & (QMAX-1);
	//assert(in[i] < QMAX && in[i] >= 0);
        model_qual[last].encodeSymbol(&rc, q);
        last = ((MAX(q1, q2)<<6) + q) & (QSIZE-1);
	last  += (q1==q2) << QBITS;

	delta += (q1>q)*(q1-q);
	last  += (MIN(7*8, delta)&0xf8) << (QBITS-2);
	//if (i%16 == 0) delta/=2;

        q2 = q1; q1 = q;
    }

    model_qual[last].encodeSymbol(&rc, QMAX-1); /* terminator */
    rc.FinishEncode();

    *out_size = rc.size_out();

    //fprintf(stdout, "%d -> %d\n", (int)in_size, (int)*out_size);

    delete[] model_qual;

    return comp;
}

unsigned char *uncompress_block_fqz2f(cram_slice *s,
				      unsigned char *in,
				      size_t in_size,
				      size_t *out_size) {
    unsigned char *uncomp = NULL;
    RangeCoder rc;
    size_t i, j, rec = 0, len = *out_size;
    unsigned char q1 = 0, q2 = 0;
    unsigned int last = 0;
    SIMPLE_MODEL<QMAX> *model_qual;
    model_qual = new SIMPLE_MODEL<QMAX>[QSIZE*16];
    SIMPLE_MODEL<256>  model_len[4];
    SIMPLE_MODEL<2> model_strand;
    int delta = 5, rev = 0;
    unsigned char *seq_st = NULL, *seq_en = NULL;
#ifdef DEDUP
    SIMPLE_MODEL<2> model_dup;
#endif

    uncomp = (unsigned char *)malloc(*out_size);

    rc.input((char *)in);
    rc.StartDecode();

    for (i = j = 0; i < len; i++, j--) {
	if (j == 0) {
	    // rev complement last read if necessary
	    if (rev) {
		for (; seq_st < seq_en; seq_st++, seq_en--) {
		    unsigned char c = *seq_st;
		    *seq_st = *seq_en;
		    *seq_en = c;
		}
	    }

	    int len;
	    len  = model_len[0].decodeSymbol(&rc);
	    len |= model_len[1].decodeSymbol(&rc)<<8;
	    rev = model_strand.decodeSymbol(&rc);

#ifdef DEDUP
	    if (model_dup.decodeSymbol(&rc)) {
		// Dup of last line
		memcpy(uncomp+i, uncomp+i-len, len);
		i += len-1;
		j = 1;
		continue;
	    }
#endif

	    seq_st = uncomp+i;
	    seq_en = seq_st + len-1;

	    rec++;
	    j = len;
	    delta = 5;
	    last = 0;
	}

        unsigned char q = model_qual[last].decodeSymbol(&rc);
        uncomp[i] = q;

        last = ((MAX(q1, q2)<<6) + q) & (QSIZE-1);
	last  += (q1==q2) << QBITS;

	delta += (q1>q)*(q1-q);
	last  += (MIN(7*8, delta)&0xf8) << (QBITS-2);

	//if (i%16 == 0) delta/=2;
        q2 = q1; q1 = q;
    }

    if (rev) {
	for (; seq_st < seq_en; seq_st++, seq_en--) {
	    unsigned char c = *seq_st;
	    *seq_st = *seq_en;
	    *seq_en = c;
	}
    }

    rc.FinishDecode();
    delete[] model_qual;

    return uncomp;
}


// // Maps to values to 0-7
// static void squash_qual(unsigned char *data,
// 			size_t l) {
//     size_t i, j;
// 
//     unsigned int map[256] = {0};
// 
//     for (i = 0; i < l; i++)
// 	map[data[i]]++;
// 
//     for (i = j = 0; i < 256; i++)
// 	if (map[i])
// 	    map[i] = j++;
// 
//     //fprintf(stderr, "%d symbols\n", j);
//     assert(j <= 8); // binned only. Maybe return this as a value
// 
//     for (i = 0; i < l; i++)
// 	data[i] = map[data[i]];
// }
// 
// 
// // Fqzcomp -q2 equiv
// unsigned char *compress_block_fqz2_BIN(int level,
// 				       cram_slice *s,
// 				       unsigned char *in,
// 				       size_t in_size,
// 				       size_t *out_size) {
//     unsigned char *comp = (unsigned char *)malloc(in_size*1.1+1000);
//     size_t i;
//     unsigned int last = 0, last_q = 0;
//     RangeCoder rc;
//     unsigned char q1 = 0, q2 = 0;
//     SIMPLE_MODEL<QMAX_B> *model_qual;
//     model_qual = new SIMPLE_MODEL<QMAX_B>[QSIZE*16];
//     int rle = 0;
//     int delta = 5;
//     int rec = 0, j;
// 
//     rc.output((char *)comp);
//     rc.StartEncode();
// 
// #if 1
//     // Reverse complement
//     rec = 0;
//     i = 0;
//     while (i < in_size) {
// 	int len = s->crecs[rec].len;
// 	if (s->crecs[rec].flags & BAM_FREVERSE) {
// 	    int I, J;
// 	    unsigned char *cp = in+i;
// 	    for (I = 0, J = len-1; I < J; I++, J--) {
// 		unsigned char c;
// 		c = cp[I];
// 		cp[I] = cp[J];
// 		cp[J] = c;
// 	    }
// 	}
// 	rec++;
// 	i+=len;
//     }
// #endif
// 
//     squash_qual(in, in_size); // binned only
// 
//     rec = 0;
//     for (i = j = 0; i < in_size; i++, j--) {
// 	if (j == 0) {
// 	    j = s->crecs[rec++].len;
// 	    delta = 5;
// 	}
// 
//         unsigned char q = in[i] & (QMAX-1);
// 	assert(in[i] < QMAX_B && in[i] >= 0); //binned only
//         model_qual[last].encodeSymbol(&rc, q);
// 
// 	last = last_q & (QSIZE-1);
// 	last_q = (last_q<<4) | q;  //binned only
// 
// 	// delta saves 3-4%, but adds 14% cpu
// 	delta += (q1>q)*(q1-q);
// 	last  += (MIN(7*8, delta)&0xf8) << (QBITS-2);
// 
// 	//_mm_prefetch((const char *)&model_qual[last], _MM_HINT_T0);
// 
// //	last = ((MAX(q1, q2)<<6) + q) & ((1<<QBITS)-1);
// //	last  += (q1==q2) << QBITS;
// //
// //	// delta saves 3-4%, but adds 14% cpu
// //	delta += (q1>q)*(q1-q);
// //	last  += (MIN(7*8, delta)&0xf8) << (QBITS-2);
// //
// //	//_mm_prefetch((const char *)&model_qual[last], _MM_HINT_T0);
// //	q2 = q1; q1 = q;
// //	//if (i%8 == 0) delta/=2;
//     }
// 
//     model_qual[last].encodeSymbol(&rc, QMAX_B-1); /* terminator */
//     rc.FinishEncode();
// 
//     *out_size = rc.size_out();
// 
//     //fprintf(stdout, "%d -> %d\n", (int)in_size, (int)*out_size);
// 
//     delete[] model_qual;
// 
//     return comp;
// }



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

extern "C" {
cram_compressor *cram_compressor_init(void) {
    return &c;
}
}

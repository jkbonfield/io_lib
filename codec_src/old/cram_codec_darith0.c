//gcc -g -fPIC -shared cram_codec_squash.c -o libcram_codec_squash.so -I.. -I. -I/nfs/sam_scratch/jkb/opt/squash/include/squash-0.7 -L/nfs/sam_scratch/jkb/opt/squash/lib/ -lsquash0.7 -Wl,-rpath,/nfs/sam_scratch/jkb/opt/squash/lib/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include "io_lib/cram_block_compression.h"

#ifndef ARITH_ORDER
#    define ARITH_ORDER 0
#endif

static const char *name(void) {
  return "Dynamic arith-0";
}

/*-----------------------------------------------------------------------------
 * Arithmetic coder. clr.cdr from fqz_comp. By Eugene Shelwein
 */

/*
 * Note it is up to the calling code to ensure that no overruns on input and
 * output buffers occur.
 *
 * Call the input() and output() functions to set and query the current
 * buffer locations.
 */

#define  TOP	   (1<<24)

typedef unsigned char uc;

typedef struct {
    uint64_t low;
    uint32_t range, code;

    uc *in_buf;
    uc *out_buf;
} RngCoder;

static inline void RC_input(RngCoder *rc, char *in) {
    rc->out_buf = rc->in_buf = (uc *)in;
}

static inline void RC_output(RngCoder *rc, char *out) {
    rc->in_buf = rc->out_buf = (uc *)out;
}

static inline int RC_size_out(RngCoder *rc) {
    return rc->out_buf - rc->in_buf;
}

static inline int RC_size_in(RngCoder *rc) {
    return rc->in_buf - rc->out_buf;
}

static inline void RC_StartEncode(RngCoder *rc) { 
    rc->low=0;
    rc->range=(uint32_t)-1; 
}

static inline void RC_StartDecode(RngCoder *rc) { 
    int i;

    rc->low=0;  
    rc->range=(uint32_t)-1;
    for (i = 0; i < 8; i++)
	rc->code = (rc->code<<8) | *rc->in_buf++;
}

static inline void RC_FinishEncode(RngCoder *rc) { 
    int i;
    for (i = 0; i < 8; i++) {
	*rc->out_buf++ = rc->low>>56;
	rc->low<<=8;
    }
}

static inline void RC_FinishDecode(RngCoder *rc) {}

static inline void RC_Encode(RngCoder *rc, uint32_t cumFreq, uint32_t freq,
	       uint32_t totFreq) {
    rc->low  += cumFreq * (rc->range/= totFreq);
    rc->range*= freq;

    if (cumFreq + freq > totFreq)
	abort();

    while (rc->range<TOP) {
	if ( (uc)((rc->low^(rc->low+rc->range))>>56) ) 
	    rc->range = (((uint32_t)rc->low | (TOP-1))-(uint32_t)rc->low);
	*rc->out_buf++ = rc->low>>56, rc->range<<=8, rc->low<<=8;
    }
}

static inline uint32_t RC_GetFreq(RngCoder *rc, uint32_t totFreq) {
    return rc->code / (rc->range /= totFreq);
}

static inline void RC_Decode(RngCoder *rc, uint32_t cumFreq, uint32_t freq,
	       uint32_t totFreq) {
    uint32_t temp = cumFreq*rc->range;
    rc->low  += temp;
    rc->code -= temp;
    rc->range*= freq;
 
    while (rc->range<TOP) {
	if ( (uc)((rc->low^(rc->low+rc->range))>>56) ) 
	    rc->range = (((uint32_t)rc->low | (TOP-1))-(uint32_t)rc->low);
	rc->code = (rc->code<<8) | *rc->in_buf++, rc->range<<=8, rc->low<<=8;
    }
}



/*
 *--------------------------------------------------------------------------
 * A simple frequency model.
 *
 * This keeps a list of symbols and their frequencies, approximately
 * sorted by symbol frequency. We allow for a single symbol to periodically
 * move up the list when emitted, effectively doing a single step of
 * bubble sort periodically. This means it's largely the same complexity
 * irrespective of alphabet size.
 * It's more efficient on strongly biased distributions than random data.
 *
 * There is no escape symbol, so the model is tailored to relatively
 * stationary samples (although we do have occasional normalisation to
 * avoid frequency counters getting too high).
 *--------------------------------------------------------------------------
 */

// Shrinking this to 1<<10 gives 2-3% smaller qualities, but 50% longer
#ifndef MAX_FREQ
#    define MAX_FREQ (1<<16)-16
#endif

#ifndef MODEL_STEP
#    define MODEL_STEP 4
#endif

#define MODEL_NSYM 256

typedef struct {
    uint16_t Symbol;
    uint16_t Freq;
} SymFreqs;

typedef struct {
    uint32_t TotFreq;  // Total frequency

    // Array of Symbols approximately sorted by Freq. 
    SymFreqs sentinel, F[MODEL_NSYM+1];
} SimpleModel;

void SM_Init(SimpleModel *sm) {
    int i;

    for (i=0; i<MODEL_NSYM; i++ ) {
	sm->F[i].Symbol = i;
	sm->F[i].Freq   = 1;
    }

    sm->TotFreq         = MODEL_NSYM;
    sm->sentinel.Symbol = 0;
    sm->sentinel.Freq   = MAX_FREQ; // Always first; simplifies sorting.

    sm->F[MODEL_NSYM].Freq = 0; // terminates normalize() loop. See below.
}

static inline void SM_normalize(SimpleModel *sm) {
    SymFreqs *s;

    /* Faster than F[i].Freq for 0 <= i < MODEL_NSYM */
    sm->TotFreq=0;
    for (s = sm->F; s->Freq; s++) {
	s->Freq     -= s->Freq>>1;
	sm->TotFreq += s->Freq;
    }
}

static inline void SM_encodeSymbol(SimpleModel *sm, RngCoder *rc, uint16_t sym) {
    SymFreqs *s = sm->F;
    uint32_t AccFreq  = 0;

    while (s->Symbol != sym)
	AccFreq += s++->Freq;

    RC_Encode(rc, AccFreq, s->Freq, sm->TotFreq);
    s->Freq     += MODEL_STEP;
    sm->TotFreq += MODEL_STEP;

    if (sm->TotFreq > MAX_FREQ)
	SM_normalize(sm);

    /* Keep approx sorted */

    if (s[0].Freq > s[-1].Freq) {
	SymFreqs t = s[0];
	s[0] = s[-1];
	s[-1] = t;
    }
}

static inline uint16_t SM_decodeSymbol(SimpleModel *sm, RngCoder *rc) {
    SymFreqs* s = sm->F;
    uint32_t freq = RC_GetFreq(rc, sm->TotFreq);
    uint32_t AccFreq;

    for (AccFreq = 0; (AccFreq += s->Freq) <= freq; s++);
    AccFreq -= s->Freq;

    RC_Decode(rc, AccFreq, s->Freq, sm->TotFreq);
    s->Freq     += MODEL_STEP;
    sm->TotFreq += MODEL_STEP;

    if (sm->TotFreq > MAX_FREQ)
	SM_normalize(sm);

    /* Keep approx sorted */
    if (s[0].Freq > s[-1].Freq) {
	SymFreqs t = s[0];
	s[0] = s[-1];
	s[-1] = t;
	return t.Symbol;
    }

    return s->Symbol;
}

/*-----------------------------------------------------------------------------
 * Memory to memory compression functions.
 */
unsigned char *arith_compress(unsigned char *in, size_t in_size,
			      size_t *out_size, int order) {
    unsigned char *out_buf = malloc(2*in_size+100);
    unsigned char *cp = out_buf;
    SimpleModel *ctx;
    int _out_size;
    RngCoder rc;
    int order_mask = (1<<(8*order))-1;
    unsigned int last, i;

    if (!out_buf)
	return NULL;

    ctx = (SimpleModel *)malloc((1<<(order*8)) * sizeof(*ctx));
    if (!ctx) {
	free(out_buf);
	return NULL;
    }

    {
	unsigned int *cp_f = (unsigned int *)&ctx[0];
	unsigned int *cp_t = (unsigned int *)&ctx[1];
	unsigned int *end = (unsigned int *)&ctx[1<<(order*8)];
	assert(sizeof(*ctx)%4 == 0);
	SM_Init(&ctx[0]);
	while (cp_t != end)
	    *cp_t++ = *cp_f++;
    }
    //for (i = 0; i < (1<<(order*8)); i++) SM_Init(&ctx[i]);

    RC_output(&rc, out_buf+5);
    RC_StartEncode(&rc);
    for (last = i = 0; i < in_size; i++) {
	unsigned char c = in[i];
	SM_encodeSymbol(&ctx[last], &rc, c);
	last = ((last<<8) | c) & order_mask;
    }
    RC_FinishEncode(&rc);
    *out_size = RC_size_out(&rc) + 5;

    *cp++ = order;
    *cp++ = (in_size>> 0) & 0xff;
    *cp++ = (in_size>> 8) & 0xff;
    *cp++ = (in_size>>16) & 0xff;
    *cp++ = (in_size>>24) & 0xff;

    free(ctx);

    return out_buf;
}

unsigned char *arith_uncompress(unsigned char *in, size_t in_size,
				size_t *out_size) {
    unsigned char *out_buf = NULL, *cp;
    SimpleModel *ctx;
    int _out_size;
    RngCoder rc;
    int order_mask, order;
    unsigned int last, i, in_sz, out_sz;

    order = in[0];
    order_mask = (1<<(8*order))-1;
    out_sz = ((in[1])<<0) | ((in[2])<<8) | ((in[3])<<16) | ((in[4])<<24);
    out_buf = malloc(out_sz);
    if (!out_buf)
	return NULL;

    ctx = (SimpleModel *)malloc((1<<(order*8)) * sizeof(*ctx));
    if (!ctx) {
	free(out_buf);
	return NULL;
    }

    {
	unsigned int *cp_f = (unsigned int *)&ctx[0];
	unsigned int *cp_t = (unsigned int *)&ctx[1];
	unsigned int *end = (unsigned int *)&ctx[1<<(order*8)];
	assert(sizeof(*ctx)%4 == 0);
	SM_Init(&ctx[0]);
	while (cp_t != end)
	    *cp_t++ = *cp_f++;
    }
    //for (i = 0; i < (1<<(order*8)); i++) SM_Init(&ctx[i]);

    RC_input(&rc,in+5);
    RC_StartDecode(&rc);
    for (last = i = 0; i < out_sz; i++) {
	unsigned char c = SM_decodeSymbol(&ctx[last], &rc);
	out_buf[i] = c;
	last = ((last<<8) | c) & order_mask;
    }
    RC_FinishDecode(&rc);

    *out_size = out_sz;

    free(ctx);

    return out_buf;
}

unsigned char *compress_block(int level,
			      cram_slice *s,
			      unsigned char *in,
			      size_t in_size,
			      size_t *out_size) {
    return arith_compress(in, in_size, out_size, ARITH_ORDER);
}

unsigned char *uncompress_block(unsigned char *in,
				// cram_slice *s,
				size_t in_size,
				size_t *out_size) {
    return arith_uncompress(in, in_size, out_size);
}

static cram_compressor c = {
    '0'+ARITH_ORDER,
    0, // all data series
    1.0,
    name,
    compress_block,
    uncompress_block,
};

cram_compressor *cram_compressor_init(void) {
    return &c;
}



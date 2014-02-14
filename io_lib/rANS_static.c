/*-------------------------------------------------------------------------- */
/* rans_byte.h from https://github.com/rygorous/ryg_rans */

// Simple byte-aligned rANS encoder/decoder - public domain - Fabian 'ryg' Giesen 2014
//
// Not intended to be "industrial strength"; just meant to illustrate the general
// idea.

#ifndef RANS_BYTE_HEADER
#define RANS_BYTE_HEADER

#include <stdint.h>

// READ ME FIRST:
//
// This is designed like a typical arithmetic coder API, but there's three
// twists you absolutely should be aware of before you start hacking:
//
// 1. You need to encode data in *reverse* - last symbol first. rANS works
//    like a stack: last in, first out.
// 2. Likewise, the encoder outputs bytes *in reverse* - that is, you give
//    it a pointer to the *end* of your buffer (exclusive), and it will
//    slowly move towards the beginning as more bytes are emitted.
// 3. Unlike basically any other entropy coder implementation you might
//    have used, you can interleave data from multiple independent rANS
//    encoders into the same bytestream without any extra signaling;
//    you can also just write some bytes by yourself in the middle if
//    you want to. This is in addition to the usual arithmetic encoder
//    property of being able to switch models on the fly. Writing raw
//    bytes can be useful when you have some data that you know is
//    incompressible, and is cheaper than going through the rANS encode
//    function. Using multiple rANS coders on the same byte stream wastes
//    a few bytes compared to using just one, but execution of two
//    independent encoders can happen in parallel on superscalar and
//    Out-of-Order CPUs, so this can be *much* faster in tight decoding
//    loops.
//
//    This is why all the rANS functions take the write pointer as an
//    argument instead of just storing it in some context struct.

// --------------------------------------------------------------------------

// L ('l' in the paper) is the lower bound of our normalization interval.
// Between this and our byte-aligned emission, we use 31 (not 32!) bits.
// This is done intentionally because exact reciprocals for 31-bit uints
// fit in 32-bit uints: this permits some optimizations during encoding.
#define RANS_BYTE_L (1u << 23)  // lower bound of our normalization interval

// State for a rANS encoder. Yep, that's all there is to it.
typedef uint32_t RansState;

// Initialize a rANS encoder.
static inline void RansEncInit(RansState* r)
{
    *r = RANS_BYTE_L;
}

// Encodes a single symbol with range start "start" and frequency "freq".
// All frequencies are assumed to sum to "1 << scale_bits", and the
// resulting bytes get written to ptr (which is updated).
//
// NOTE: With rANS, you need to encode symbols in *reverse order*, i.e. from
// beginning to end! Likewise, the output bytestream is written *backwards*:
// ptr starts pointing at the end of the output buffer and keeps decrementing.
static inline void RansEncPut(RansState* r, uint8_t** pptr, uint32_t start, uint32_t freq, uint32_t scale_bits)
{
    // renormalize
    uint32_t x = *r;
    uint32_t x_max = ((RANS_BYTE_L >> scale_bits) << 8) * freq; // this turns into a shift.
    if (x >= x_max) {
        uint8_t* ptr = *pptr;
        do {
            *--ptr = (uint8_t) (x & 0xff);
            x >>= 8;
        } while (x >= x_max);
        *pptr = ptr;
    }

    // x = C(s,x)
    *r = ((x / freq) << scale_bits) + (x % freq) + start;
}

// Flushes the rANS encoder.
static inline void RansEncFlush(RansState* r, uint8_t** pptr)
{
    uint32_t x = *r;
    uint8_t* ptr = *pptr;

    ptr -= 4;
    ptr[0] = (uint8_t) (x >> 0);
    ptr[1] = (uint8_t) (x >> 8);
    ptr[2] = (uint8_t) (x >> 16);
    ptr[3] = (uint8_t) (x >> 24);

    *pptr = ptr;
}

// Initializes a rANS decoder.
// Unlike the encoder, the decoder works forwards as you'd expect.
static inline void RansDecInit(RansState* r, uint8_t** pptr)
{
    uint32_t x;
    uint8_t* ptr = *pptr;

    x  = ptr[0] << 0;
    x |= ptr[1] << 8;
    x |= ptr[2] << 16;
    x |= ptr[3] << 24;
    ptr += 4;

    *pptr = ptr;
    *r = x;
}

// Returns the current cumulative frequency (map it to a symbol yourself!)
static inline uint32_t RansDecGet(RansState* r, uint32_t scale_bits)
{
    return *r & ((1u << scale_bits) - 1);
}

// Advances in the bit stream by "popping" a single symbol with range start
// "start" and frequency "freq". All frequencies are assumed to sum to "1 << scale_bits",
// and the resulting bytes get written to ptr (which is updated).
static inline void RansDecAdvance(RansState* r, uint8_t** pptr, uint32_t start, uint32_t freq, uint32_t scale_bits)
{
    uint32_t mask = (1u << scale_bits) - 1;

    // s, x = D(x)
    uint32_t x = *r;
    x = freq * (x >> scale_bits) + (x & mask) - start;

    // renormalize
    if (x < RANS_BYTE_L) {
        uint8_t* ptr = *pptr;
        do x = (x << 8) | *ptr++; while (x < RANS_BYTE_L);
        *pptr = ptr;
    }

    *r = x;
}

// --------------------------------------------------------------------------

// That's all you need for a full encoder; below here are some utility
// functions with extra convenience or optimizations.

// Description for a symbol
typedef struct {
    uint32_t start;     // Start of range
    uint32_t freq;      // Frequency of symbol (=size of range)
    uint32_t rcp_freq;  // Reciprocal frequency
    uint32_t rcp_shift; // Reciprocal shift
} RansSymbol;

// Initializes a symbol to start "start" and frequency "freq"
static inline void RansSymbolInit(RansSymbol* s, uint32_t start, uint32_t freq)
{
    s->start = start;
    s->freq = freq;
    if (freq < 2) // 0 is unsupported (div by zero!) and 1 requires special care
        s->rcp_freq = s->rcp_shift = 0;
    else {
        // Alverson, "Integer Division using reciprocals"
        uint32_t shift = 0;
        //while ((freq >> shift) > 1)
	while (freq > (1u << shift))
            shift++;

        s->rcp_freq = (uint32_t) (((1ull << (shift + 31)) + freq-1) / freq);
        s->rcp_shift = shift - 1;
    }
}

// Encodes a given symbol. This is faster than straight RansEnc since we can do
// multiplications instead of a divide.
static inline void RansEncPutSymbol(RansState* r, uint8_t** pptr, RansSymbol const* sym, uint32_t scale_bits)
{
    // renormalize
    uint32_t x = *r;
    uint32_t x_max = ((RANS_BYTE_L >> scale_bits) << 8) * sym->freq; // this turns into a shift.
    if (x >= x_max) {
        uint8_t* ptr = *pptr;
        do {
            *--ptr = (uint8_t) (x & 0xff);
            x >>= 8;
        } while (x >= x_max);
        *pptr = ptr;
    }

    // x = C(s,x)
    if (sym->freq == 1)
        x = (x << scale_bits) + sym->start;
    else {
        // written strangely so the compiler generates a multiply high, but only
        // uses 32-bit shifts.
        uint32_t q = (uint32_t) (((uint64_t)x * sym->rcp_freq) >> 32) >> sym->rcp_shift;
	//assert(q == x/sym->freq);
        uint32_t p = x - q * sym->freq;
        x = (q << scale_bits) + sym->start + p;
    }

    *r = x;
}

// Equivalent to RansDecAdvance that takes a symbol.
static inline void RansDecAdvanceSymbol(RansState* r, uint8_t** pptr, RansSymbol const* sym, uint32_t scale_bits)
{
    RansDecAdvance(r, pptr, sym->start, sym->freq, scale_bits);
}

// Advances in the bit stream by "popping" a single symbol with range start
// "start" and frequency "freq". All frequencies are assumed to sum to "1 << scale_bits".
// No renormalization or output happens.
static inline void RansDecAdvanceStep(RansState* r, uint32_t start, uint32_t freq, uint32_t scale_bits)
{
    uint32_t mask = (1u << scale_bits) - 1;

    // s, x = D(x)
    uint32_t x = *r;
    *r = freq * (x >> scale_bits) + (x & mask) - start;
}

// Equivalent to RansDecAdvanceStep that takes a symbol.
static inline void RansDecAdvanceSymbolStep(RansState* r, RansSymbol const* sym, uint32_t scale_bits)
{
    RansDecAdvanceStep(r, sym->start, sym->freq, scale_bits);
}


// Renormalize.
static inline void RansDecRenorm(RansState* r, uint8_t** pptr)
{
    // renormalize
    uint32_t x = *r;
    if (x < RANS_BYTE_L) {
        uint8_t* ptr = *pptr;
        do x = (x << 8) | *ptr++; while (x < RANS_BYTE_L);
        *pptr = ptr;
    }

    *r = x;
}

#endif // RANS_BYTE_HEADER

/*-------------------------------------------------------------------------- */
/*
 * Example wrapper to use the rans_byte.h functions included above.
 *
 * This demonstrates how to use, and unroll, an order-0 and order-1 frequency
 * model.
 */

/*
 * Copyright (c) 2014 Genome Research Ltd.
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
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2014
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>

#define TF_SHIFT 14
#define TOTFREQ (1<<TF_SHIFT)

#define ABS(a) ((a)>0?(a):-(a))
#define BLK_SIZE 2000000

// Room to allow for expanded BLK_SIZE on worst case compression.
#define BLK_SIZE2 ((int)(1.05*BLK_SIZE))

/*-----------------------------------------------------------------------------
 * Memory to memory compression functions.
 *
 * These are original versions without any manual loop unrolling. They
 * are easier to understand, but can be up to 2x slower.
 */

unsigned char *rans_compress_O0(unsigned char *in, unsigned int in_size,
				unsigned int *out_size) {
    unsigned char *out_buf = malloc(1.05*in_size + 257*257*3 + 4);
    unsigned char *cp, *out_end;
    RansSymbol syms[256];
    RansState rans0, rans1, rans2, rans3;
    uint8_t* ptr;
    int F[256], C[256], T = 0, i, j, n, tab_size, rle;
    unsigned char c;

    if (!out_buf)
	return NULL;

    ptr = out_end = out_buf + (int)(1.05*in_size) + 257*257*3 + 4;

    // Compute statistics
    memset(F, 0, 256*sizeof(int));
    memset(C, 0, 256*sizeof(int));
    for (i = 0; i < in_size; i++) {
	F[c = in[i]]++;
	T++;
    }

    // Normalise so T[i] == 65536
    for (n = j = 0; j < 256; j++)
	if (F[j])
	    n++;

    for (j = 0; j < 256; j++) {
	if (!F[j])
	    continue;
	if ((F[j] *= ((double)TOTFREQ-n)/T) == 0)
	    F[j] = 1;
    }

    // Encode statistics.
    cp = out_buf+4;
    for (T = rle = j = 0; j < 256; j++) {
	C[j] = T;
	T += F[j];
	if (F[j]) {
	    if (rle) {
		rle--;
	    } else {
		*cp++ = j;
		if (!rle && j && F[j-1]) {
		    for(rle=j+1; rle<256 && F[rle]; rle++)
			;
		    rle -= j+1;
		    *cp++ = rle;
		}
	    }
	    if (F[j]<128) {
		*cp++ = F[j];
	    } else {
		*cp++ = 128 | (F[j]>>8);
		*cp++ = F[j]&0xff;
	    }
	    RansSymbolInit(&syms[j], C[j], F[j]);
	}
    }
    *cp++ = 0;
    tab_size = cp-out_buf;

    RansEncInit(&rans0);
    RansEncInit(&rans1);
    RansEncInit(&rans2);
    RansEncInit(&rans3);

    switch (i=(in_size&3)) {
    case 3: RansEncPutSymbol(&rans2, &ptr, &syms[in[in_size-(i-2)]], TF_SHIFT);
    case 2: RansEncPutSymbol(&rans1, &ptr, &syms[in[in_size-(i-1)]], TF_SHIFT);
    case 1: RansEncPutSymbol(&rans0, &ptr, &syms[in[in_size-(i-0)]], TF_SHIFT);
    case 0:
	break;
    }
    for (i=(in_size &~3); i>0; i-=4) {
	unsigned char c3 = in[i-1];
	unsigned char c2 = in[i-2];
	unsigned char c1 = in[i-3];
	unsigned char c0 = in[i-4];

	RansEncPutSymbol(&rans3, &ptr, &syms[c3], TF_SHIFT);
	RansEncPutSymbol(&rans2, &ptr, &syms[c2], TF_SHIFT);
	RansEncPutSymbol(&rans1, &ptr, &syms[c1], TF_SHIFT);
	RansEncPutSymbol(&rans0, &ptr, &syms[c0], TF_SHIFT);
    }

    RansEncFlush(&rans3, &ptr);
    RansEncFlush(&rans2, &ptr);
    RansEncFlush(&rans1, &ptr);
    RansEncFlush(&rans0, &ptr);

    // Finalise block size and return it
    *out_size = (out_end - ptr) + tab_size;

    cp = out_buf;
    *cp++ = (in_size>> 0) & 0xff;
    *cp++ = (in_size>> 8) & 0xff;
    *cp++ = (in_size>>16) & 0xff;
    *cp++ = (in_size>>24) & 0xff;

    memmove(out_buf + tab_size, ptr, out_end-ptr);

    return out_buf;
}

typedef struct {
    struct {
	int F;
	int C;
    } fc[256];
    unsigned char *R;
} ari_decoder;

unsigned char *rans_uncompress_O0(unsigned char *in, unsigned int in_size,
				  unsigned int *out_size) {
    /* Load in the static tables */
    unsigned char *cp = in + 4;
    int i, j, x, out_sz, rle;
    char *out_buf;
    ari_decoder D;
    RansSymbol syms[256];

    memset(&D, 0, sizeof(D));

    out_sz = ((in[0])<<0) | ((in[1])<<8) | ((in[2])<<16) | ((in[3])<<24);
    out_buf = malloc(out_sz);
    if (!out_buf)
	return NULL;

    //fprintf(stderr, "out_sz=%d\n", out_sz);

    // Precompute reverse lookup of frequency.
    j = *cp++;
    x = 0;
    rle = 0;
    do {
	if ((D.fc[j].F = *cp++) >= 128) {
	    D.fc[j].F &= ~128;
	    D.fc[j].F = ((D.fc[j].F & 127) << 8) | *cp++;
	}
	D.fc[j].C = x;

	RansSymbolInit(&syms[j], D.fc[j].C, D.fc[j].F);

	/* Build reverse lookup table */
	if (!D.R) D.R = (unsigned char *)malloc(TOTFREQ);
	memset(&D.R[x], j, D.fc[j].F);

	x += D.fc[j].F;

	if (!rle && *cp == j+1) {
	    j = *cp++;
	    rle = *cp++;
	} else if (rle) {
	    rle--;
	    j++;
	} else {
	    j = *cp++;
	}
    } while(j);

    assert(x < TOTFREQ);

    RansState rans0, rans1, rans2, rans3;
    uint8_t *ptr = cp;
    RansDecInit(&rans0, &ptr);
    RansDecInit(&rans1, &ptr);
    RansDecInit(&rans2, &ptr);
    RansDecInit(&rans3, &ptr);

    int out_end = (out_sz&~3);
    for (i=0; i < out_end; i+=4) {
	unsigned char c0 = D.R[RansDecGet(&rans0, TF_SHIFT)];
	unsigned char c1 = D.R[RansDecGet(&rans1, TF_SHIFT)];
	unsigned char c2 = D.R[RansDecGet(&rans2, TF_SHIFT)];
	unsigned char c3 = D.R[RansDecGet(&rans3, TF_SHIFT)];

	out_buf[i+0] = c0;
	out_buf[i+1] = c1;
	out_buf[i+2] = c2;
	out_buf[i+3] = c3;

	RansDecAdvanceSymbolStep(&rans0, &syms[c0], TF_SHIFT);
	RansDecAdvanceSymbolStep(&rans1, &syms[c1], TF_SHIFT);
	RansDecAdvanceSymbolStep(&rans2, &syms[c2], TF_SHIFT);
	RansDecAdvanceSymbolStep(&rans3, &syms[c3], TF_SHIFT);

	RansDecRenorm(&rans0, &ptr);
	RansDecRenorm(&rans1, &ptr);
	RansDecRenorm(&rans2, &ptr);
	RansDecRenorm(&rans3, &ptr);
    }

    switch(out_sz&3) {
	unsigned char c;
    case 0:
	break;
    case 1:
	c = D.R[RansDecGet(&rans0, TF_SHIFT)];
	RansDecAdvanceSymbol(&rans0, &ptr, &syms[c], TF_SHIFT);
	out_buf[out_end] = c;
	break;

    case 2:
	c = D.R[RansDecGet(&rans0, TF_SHIFT)];
	RansDecAdvanceSymbol(&rans0, &ptr, &syms[c], TF_SHIFT);
	out_buf[out_end] = c;

	c = D.R[RansDecGet(&rans1, TF_SHIFT)];
	RansDecAdvanceSymbol(&rans1, &ptr, &syms[c], TF_SHIFT);
	out_buf[out_end+1] = c;
	break;

    case 3:
	c = D.R[RansDecGet(&rans0, TF_SHIFT)];
	RansDecAdvanceSymbol(&rans0, &ptr, &syms[c], TF_SHIFT);
	out_buf[out_end] = c;

	c = D.R[RansDecGet(&rans1, TF_SHIFT)];
	RansDecAdvanceSymbol(&rans1, &ptr, &syms[c], TF_SHIFT);
	out_buf[out_end+1] = c;

	c = D.R[RansDecGet(&rans2, TF_SHIFT)];
	RansDecAdvanceSymbol(&rans2, &ptr, &syms[c], TF_SHIFT);
	out_buf[out_end+2] = c;
	break;
    }
    
    *out_size = out_sz;

    if (D.R) free(D.R);

    return (unsigned char *)out_buf;
}

unsigned char *rans_compress_O1(unsigned char *in, unsigned int in_size,
				unsigned int *out_size) {
    unsigned char *out_buf = malloc(1.05*in_size + 257*257*3 + 4);
    unsigned char *cp = out_buf, *out_end;
    unsigned int last, tab_size;
    RansSymbol syms[256][256];

    if (!out_buf)
	return NULL;

    out_end = out_buf + (int)(1.05*in_size) + 257*257*3 + 4;
    cp = out_buf+4;

    int F[256][256], C[256][256], T[256], i, j;
    unsigned char c;

    memset(F, 0, 256*256*sizeof(int));
    memset(C, 0, 256*256*sizeof(int));
    memset(T, 0, 256*sizeof(int));
    //for (last = 0, i=in_size-1; i>=0; i--) {
    //	F[last][c = in[i]]++;
    //	T[last]++;
    //	last = c;
    //}

    for (last=i=0; i<in_size; i++) {
	F[last][c = in[i]]++;
	T[last]++;
	last = c;
    }
    F[0][in[1*(in_size>>3)]]++;
    F[0][in[2*(in_size>>3)]]++;
    F[0][in[3*(in_size>>3)]]++;
    F[0][in[4*(in_size>>3)]]++;
    F[0][in[5*(in_size>>3)]]++;
    F[0][in[6*(in_size>>3)]]++;
    F[0][in[7*(in_size>>3)]]++;
    T[0]+=7;

    // Normalise so T[i] == 65536
    for (i = 0; i < 256; i++) {
	int t = T[i], t2, n;

	if (t == 0)
	    continue;

	for (n = j = 0; j < 256; j++)
	    if (F[i][j])
		n++;

	for (t2 = j = 0; j < 256; j++) {
	    if (!F[i][j])
		continue;
	    if ((F[i][j] *= ((double)TOTFREQ-n)/t) == 0)
		F[i][j] = 1;
	    t2 += F[i][j];
	}

	assert(t2 <= TOTFREQ);
    }

    for (i = 0; i < 256; i++) {
	unsigned int x = 0;
	int rle = 0;

	if (!T[i])
	    continue;

	*cp++ = i;
	for (j = 0; j < 256; j++) {
	    C[i][j] = x;
	    if (F[i][j]) {
		//fprintf(stderr, "F[%d][%d]=%d, x=%d\n", i, j, F[i][j], x);
		x += F[i][j];
		if (rle) {
		    rle--;
		} else {
		    *cp++ = j;
		    if (!rle && j && F[i][j-1]) {
			for(rle=j+1; rle<256 && F[i][rle]; rle++)
			    ;
			rle -= j+1;
			*cp++ = rle;
		    }
		}
		if (F[i][j]<128) {
		    *cp++ = F[i][j];
		} else {
		    *cp++ = 128 | (F[i][j]>>8);
		    *cp++ = F[i][j]&0xff;
		}

		RansSymbolInit(&syms[i][j], C[i][j], F[i][j]);
	    }
	}
	*cp++ = 0;
	T[i] = x;
    }
    *cp++ = 0;
    tab_size = cp - out_buf;
    assert(tab_size < 257*257*3);

    RansState rans0, rans1, rans2, rans3, rans4, rans5, rans6, rans7;
    RansEncInit(&rans0);
    RansEncInit(&rans1);
    RansEncInit(&rans2);
    RansEncInit(&rans3);
    RansEncInit(&rans4);
    RansEncInit(&rans5);
    RansEncInit(&rans6);
    RansEncInit(&rans7);

    uint8_t* ptr = out_end;

    int isz4 = in_size>>3;
    int i0 = 1*isz4-2;
    int i1 = 2*isz4-2;
    int i2 = 3*isz4-2;
    int i3 = 4*isz4-2;
    int i4 = 5*isz4-2;
    int i5 = 6*isz4-2;
    int i6 = 7*isz4-2;
    int i7 = 8*isz4-2;

    int l0 = in[i0+1];
    int l1 = in[i1+1];
    int l2 = in[i2+1];
    int l3 = in[i3+1];
    int l4 = in[i4+1];
    int l5 = in[i5+1];
    int l6 = in[i6+1];
    int l7 = in[i7+1];

    // Deal with the remainder
    l7 = in[in_size-1];
    for (i7 = in_size-2; i7 > 8*isz4-2; i7--) {
	unsigned char c7 = in[i7];
	RansEncPutSymbol(&rans7, &ptr, &syms[c7][l7], TF_SHIFT);
	l7 = c7;
    }

    for (; i0 >= 0; i0--, i1--, i2--, i3--, i4--, i5--, i6--, i7--) {
	unsigned char c0 = in[i0];
	unsigned char c1 = in[i1];
	unsigned char c2 = in[i2];
	unsigned char c3 = in[i3];
	unsigned char c4 = in[i4];
	unsigned char c5 = in[i5];
	unsigned char c6 = in[i6];
	unsigned char c7 = in[i7];

	RansEncPutSymbol(&rans7, &ptr, &syms[c7][l7], TF_SHIFT);
	RansEncPutSymbol(&rans6, &ptr, &syms[c6][l6], TF_SHIFT);
	RansEncPutSymbol(&rans5, &ptr, &syms[c5][l5], TF_SHIFT);
	RansEncPutSymbol(&rans4, &ptr, &syms[c4][l4], TF_SHIFT);
	RansEncPutSymbol(&rans3, &ptr, &syms[c3][l3], TF_SHIFT);
	RansEncPutSymbol(&rans2, &ptr, &syms[c2][l2], TF_SHIFT);
	RansEncPutSymbol(&rans1, &ptr, &syms[c1][l1], TF_SHIFT);
	RansEncPutSymbol(&rans0, &ptr, &syms[c0][l0], TF_SHIFT);

	l0 = c0;
	l1 = c1;
	l2 = c2;
	l3 = c3;
	l4 = c4;
	l5 = c5;
	l6 = c6;
	l7 = c7;
    }

    RansEncPutSymbol(&rans7, &ptr, &syms[0][l7], TF_SHIFT);
    RansEncPutSymbol(&rans6, &ptr, &syms[0][l6], TF_SHIFT);
    RansEncPutSymbol(&rans5, &ptr, &syms[0][l5], TF_SHIFT);
    RansEncPutSymbol(&rans4, &ptr, &syms[0][l4], TF_SHIFT);
    RansEncPutSymbol(&rans3, &ptr, &syms[0][l3], TF_SHIFT);
    RansEncPutSymbol(&rans2, &ptr, &syms[0][l2], TF_SHIFT);
    RansEncPutSymbol(&rans1, &ptr, &syms[0][l1], TF_SHIFT);
    RansEncPutSymbol(&rans0, &ptr, &syms[0][l0], TF_SHIFT);

    RansEncFlush(&rans7, &ptr);
    RansEncFlush(&rans6, &ptr);
    RansEncFlush(&rans5, &ptr);
    RansEncFlush(&rans4, &ptr);
    RansEncFlush(&rans3, &ptr);
    RansEncFlush(&rans2, &ptr);
    RansEncFlush(&rans1, &ptr);
    RansEncFlush(&rans0, &ptr);

    *out_size = (out_end - ptr) + tab_size;

    cp = out_buf;
    *cp++ = (in_size>> 0) & 0xff;
    *cp++ = (in_size>> 8) & 0xff;
    *cp++ = (in_size>>16) & 0xff;
    *cp++ = (in_size>>24) & 0xff;

    memmove(out_buf + tab_size, ptr, out_end-ptr);

    return out_buf;
}

unsigned char *rans_uncompress_O1(unsigned char *in, unsigned int in_size,
				  unsigned int *out_size) {
    /* Load in the static tables */
    unsigned char *cp = in + 4;
    int i, j = -999, x, out_sz;
    char *out_buf;
    ari_decoder D[256];
    RansSymbol syms[256][256];
    
    memset(D, 0, 256*sizeof(*D));

    out_sz = ((in[0])<<0) | ((in[1])<<8) | ((in[2])<<16) | ((in[3])<<24);
    out_buf = malloc(out_sz);
    if (!out_buf)
	return NULL;

    //fprintf(stderr, "out_sz=%d\n", out_sz);

    i = *cp++;
    do {
	int rle = 0;
	j = *cp++;
	x = 0;
	do {
	    if ((D[i].fc[j].F = *cp++) >= 128) {
		D[i].fc[j].F &= ~128;
		D[i].fc[j].F = ((D[i].fc[j].F & 127) << 8) | *cp++;
	    }
	    D[i].fc[j].C = x;

	    //fprintf(stderr, "i=%d j=%d F=%d C=%d\n", i, j, D[i].fc[j].F, D[i].fc[j].C);

	    if (!D[i].fc[j].F)
		D[i].fc[j].F = TOTFREQ;

	    RansSymbolInit(&syms[i][j], D[i].fc[j].C, D[i].fc[j].F);

	    /* Build reverse lookup table */
	    if (!D[i].R) D[i].R = (unsigned char *)malloc(TOTFREQ);
	    memset(&D[i].R[x], j, D[i].fc[j].F);

	    x += D[i].fc[j].F;
	    if (!rle && *cp == j+1) {
		j = *cp++;
		rle = *cp++;
	    } else if (rle) {
		rle--;
		j++;
	    } else {
		j = *cp++;
	    }
	} while(j);

	i = *cp++;
    } while (i);

    // Precompute reverse lookup of frequency.

    RansState rans0, rans1, rans2, rans3, rans4, rans5, rans6, rans7;
    uint8_t *ptr = cp;
    RansDecInit(&rans0, &ptr);
    RansDecInit(&rans1, &ptr);
    RansDecInit(&rans2, &ptr);
    RansDecInit(&rans3, &ptr);
    RansDecInit(&rans4, &ptr);
    RansDecInit(&rans5, &ptr);
    RansDecInit(&rans6, &ptr);
    RansDecInit(&rans7, &ptr);

    int isz4 = out_sz>>3;
    int i0 = 0*isz4;
    int i1 = 1*isz4;
    int i2 = 2*isz4;
    int i3 = 3*isz4;
    int i4 = 4*isz4;
    int i5 = 5*isz4;
    int i6 = 6*isz4;
    int i7 = 7*isz4;
    int l0 = 0;
    int l1 = 0;
    int l2 = 0;
    int l3 = 0;
    int l4 = 0;
    int l5 = 0;
    int l6 = 0;
    int l7 = 0;
    for (; i0 < isz4; i0++, i1++, i2++, i3++, i4++, i5++, i6++, i7++) {
	unsigned char c0 = D[l0].R[RansDecGet(&rans0, TF_SHIFT)];
	unsigned char c1 = D[l1].R[RansDecGet(&rans1, TF_SHIFT)];
	unsigned char c2 = D[l2].R[RansDecGet(&rans2, TF_SHIFT)];
	unsigned char c3 = D[l3].R[RansDecGet(&rans3, TF_SHIFT)];
	unsigned char c4 = D[l4].R[RansDecGet(&rans4, TF_SHIFT)];
	unsigned char c5 = D[l5].R[RansDecGet(&rans5, TF_SHIFT)];
	unsigned char c6 = D[l6].R[RansDecGet(&rans6, TF_SHIFT)];
	unsigned char c7 = D[l7].R[RansDecGet(&rans7, TF_SHIFT)];

	out_buf[i0] = c0;
	out_buf[i1] = c1;
	out_buf[i2] = c2;
	out_buf[i3] = c3;
	out_buf[i4] = c4;
	out_buf[i5] = c5;
	out_buf[i6] = c6;
	out_buf[i7] = c7;

	RansDecAdvanceSymbolStep(&rans0, &syms[l0][c0], TF_SHIFT);
	RansDecAdvanceSymbolStep(&rans1, &syms[l1][c1], TF_SHIFT);
	RansDecAdvanceSymbolStep(&rans2, &syms[l2][c2], TF_SHIFT);
	RansDecAdvanceSymbolStep(&rans3, &syms[l3][c3], TF_SHIFT);
	RansDecAdvanceSymbolStep(&rans4, &syms[l4][c4], TF_SHIFT);
	RansDecAdvanceSymbolStep(&rans5, &syms[l5][c5], TF_SHIFT);
	RansDecAdvanceSymbolStep(&rans6, &syms[l6][c6], TF_SHIFT);
	RansDecAdvanceSymbolStep(&rans7, &syms[l7][c7], TF_SHIFT);

	RansDecRenorm(&rans0, &ptr);
	RansDecRenorm(&rans1, &ptr);
	RansDecRenorm(&rans2, &ptr);
	RansDecRenorm(&rans3, &ptr);
	RansDecRenorm(&rans4, &ptr);
	RansDecRenorm(&rans5, &ptr);
	RansDecRenorm(&rans6, &ptr);
	RansDecRenorm(&rans7, &ptr);

	l0 = c0;
	l1 = c1;
	l2 = c2;
	l3 = c3;
	l4 = c4;
	l5 = c5;
	l6 = c6;
	l7 = c7;
    }

    // Remainder
    for (; i7 < out_sz; i7++) {
	unsigned char c7 = D[l7].R[RansDecGet(&rans7, TF_SHIFT)];
	out_buf[i7] = c7;
	RansDecAdvanceSymbol(&rans7, &ptr, &syms[l7][c7], TF_SHIFT);
	l7 = c7;
    }
    
    *out_size = out_sz;

    for (i = 0; i < 256; i++)
	if (D[i].R) free(D[i].R);

    return (unsigned char *)out_buf;
}

/*-----------------------------------------------------------------------------
 * Simple interface to the order-0 vs order-1 encoders and decoders.
 */
unsigned char *rans_compress(unsigned char *in, unsigned int in_size,
			     unsigned int *out_size, int order) {
    return order
	? rans_compress_O1(in, in_size, out_size)
	: rans_compress_O0(in, in_size, out_size);
}

unsigned char *rans_uncompress(unsigned char *in, unsigned int in_size,
			       unsigned int *out_size, int order) {
    return order
	? rans_uncompress_O1(in, in_size, out_size)
	: rans_uncompress_O0(in, in_size, out_size);
}

#ifdef TEST_MAIN
/*-----------------------------------------------------------------------------
 * Main
 */
int main(int argc, char **argv) {
    int opt, order = 1;
    unsigned char in_buf[BLK_SIZE2+257*257*3];
    int decode = 0;
    FILE *infp = stdin, *outfp = stdout;
    struct timeval tv1, tv2;
    size_t bytes = 0;

    extern char *optarg;
    extern int optind, opterr, optopt;

    while ((opt = getopt(argc, argv, "o:d")) != -1) {
	switch (opt) {
	case 'o':
	    order = atoi(optarg);
	    break;

	case 'd':
	    decode = 1;
	    break;
	}
    }

    order = order ? 1 : 0; // Only support O(0) and O(1)

    if (optind < argc) {
	if (!(infp = fopen(argv[optind], "rb"))) {
	    perror(argv[optind]);
	    return 1;
	}
	optind++;
    }

    if (optind < argc) {
	if (!(outfp = fopen(argv[optind], "wb"))) {
	    perror(argv[optind]);
	    return 1;
	}
	optind++;
    }

    gettimeofday(&tv1, NULL);

    if (decode) {
	// Only used in some test implementations of RC_GetFreq()
	//RC_init();
	//RC_init2();

	for (;;) {
	    uint32_t in_size, out_size;
	    unsigned char *out;

	    order = fgetc(infp);
	    if (4 != fread(&in_size, 1, 4, infp))
		break;
	    if (in_size != fread(in_buf, 1, in_size, infp)) {
		fprintf(stderr, "Truncated input\n");
		exit(1);
	    }
	    out = rans_uncompress(in_buf, in_size, &out_size, order);
	    if (!out)
		abort();

	    fwrite(out, 1, out_size, outfp);
	    free(out);

	    bytes += out_size;
	}
    } else {
	for (;;) {
	    uint32_t in_size, out_size;
	    unsigned char *out;

	    in_size = fread(in_buf, 1, BLK_SIZE, infp);
	    if (in_size <= 0)
		break;

	    out = rans_compress(in_buf, in_size, &out_size, order);

	    fputc(order, outfp);
	    fwrite(&out_size, 1, 4, outfp);
	    fwrite(out, 1, out_size, outfp);
	    free(out);

	    bytes += in_size;
	}
    }

    gettimeofday(&tv2, NULL);

    fprintf(stderr, "Took %ld microseconds, %5.1f MB/s\n",
	    (long)(tv2.tv_sec - tv1.tv_sec)*1000000 +
	    tv2.tv_usec - tv1.tv_usec,
	    (double)bytes / ((long)(tv2.tv_sec - tv1.tv_sec)*1000000 +
			     tv2.tv_usec - tv1.tv_usec));
    return 0;
}
#endif

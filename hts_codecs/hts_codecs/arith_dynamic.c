/*
 * Copyright (c) 2019 Genome Research Ltd.
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
 *       Institute nor the names of its contributors may be used to endorse
 *       or promote products derived from this software without specific
 *       prior written permission.
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

// As per standard rANS_static but using optional RLE or bit-packing
// techniques prior to entropy encoding.  This is a significant
// reduction in some data sets.

// top bits in order byte
#define X_PACK  0x80    // Pack 2,4,8 or infinite symbols into a byte.
#define X_RLE   0x40    // Run length encoding with runs & lits encoded separately
#define X_CAT   0x20    // Nop; for tiny segments where rANS overhead is too big
#define X_NOSZ  0x10    // Don't store the original size; used by X4 mode
#define X_4     0x08    // For 4-byte integer data; rotate & encode 4 streams.
#define X_EXT   0x04    // External compression codec via magic num (gz, xz, bz2)
#define X_ORDER 0x03    // Mask to obtain order

#include <bzlib.h>

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <limits.h>
#ifndef NO_THREADS
#include <pthread.h>
#endif

#include "arith_dynamic.h"
#include "varint.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

/*-----------------------------------------------------------------------------
 * Memory to memory compression functions.
 *
 * These are original versions without any manual loop unrolling. They
 * are easier to understand, but can be up to 2x slower.
 */
#define MAGIC 8

unsigned int arith_compress_bound(unsigned int size, int order) {
    return (order == 0
	? 1.05*size + 257*3 + 4
	: 1.05*size + 257*257*3 + 4 + 257*3+4) +
	((order & X_PACK) ? 1 : 0) +
	((order & X_RLE) ? 1 + 257*3+4: 0) + 5;
}

#ifndef MODEL_256 // see fqzcomp_qual_fuzz.c
#define NSYM 256
#include "c_simple_model.h"
#endif

// Compresses in_size bytes from 'in' to *out_size bytes in 'out'.
//
// NB: The output buffer does not hold the original size, so it is up to
// the caller to store this.
unsigned char *arith_compress_O0(unsigned char *in, unsigned int in_size,
				 unsigned char *out, unsigned int *out_size) {
    int i, bound = arith_compress_bound(in_size,0)-5; // -5 for order/size

    if (!out) {
	*out_size = bound;
	out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
	return NULL;

    unsigned int m = 0;
    for (i = 0; i < in_size; i++)
	if (m < in[i])
	    m = in[i];
    m++;
    *out = m;

    SIMPLE_MODEL(256,_) byte_model;
    SIMPLE_MODEL(256,_init)(&byte_model, m);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out+1);
    RC_StartEncode(&rc);

    for (i = 0; i < in_size; i++)
	SIMPLE_MODEL(256, _encodeSymbol)(&byte_model, &rc, in[i]);

    RC_FinishEncode(&rc);

    // Finalise block size and return it
    *out_size = RC_OutSize(&rc)+1;

    return out;
}

unsigned char *arith_uncompress_O0(unsigned char *in, unsigned int in_size,
				   unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;
    int i;
    unsigned int m = in[0] ? in[0] : 256;

    SIMPLE_MODEL(256,_) byte_model;
    SIMPLE_MODEL(256,_init)(&byte_model, m);

    if (!out)
	out = malloc(out_sz);
    if (!out)
	return NULL;

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    for (i = 0; i < out_sz; i++)
	out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model, &rc);

    RC_FinishDecode(&rc);
    
    return out;
}


//-----------------------------------------------------------------------------

// Bit packing symbols to take 0, 1(8 sym), 2(4 sym), or 4(16) bits.
static uint8_t *a_pack(uint8_t *data, int64_t len,
		       uint8_t *out_meta, int *out_meta_len, uint64_t *out_len) {
    int p[256] = {0}, n;
    uint64_t i, j;
    uint8_t *out = malloc(len+1);
    if (!out)
	return NULL;

    // count syms
    for (i = 0; i < len; i++)
	p[data[i]]=1;
    
    for (i = n = 0; i < 256; i++) {
	if (p[i]) {
	    p[i] = n++; // p[i] is now the code number
	    out_meta[n] = i;
	}
    }
    j = n+1;

    // 1 value per byte
    if (n > 16 || len < j + len/2) {
	out_meta[0] = 255; // anything >= 16 will do
	*out_meta_len = 1;
	// FIXME shortcut this by returning data and checking later.
	memcpy(out, data, len);
	*out_len = len;
	return out;
    }

    // Work out how many values per byte to encode.
    int val_per_byte;
    if (n > 4)
	val_per_byte = 2;
    else if (n > 2)
	val_per_byte = 4;
    else if (n > 1)
	val_per_byte = 8;
    else
	val_per_byte = 0; // infinite

    // We store the actual number of symbols, and intuit from that
    // during decode the number of symbols per byte.
    out_meta[0] = n;

    *out_meta_len = j;
    j = 0;

    switch (val_per_byte) {
    case 2:
	for (i = 0; i < (len & ~1); i+=2)
	    out[j++] = (p[data[i]]<<0) | (p[data[i+1]]<<4);
	switch (len-i) {
	case 1: out[j++] = p[data[i]];
	}
	*out_len = j;
	return out;

    case 4: {
	for (i = 0; i < (len & ~3); i+=4)
	    out[j++] = (p[data[i]]<<0) | (p[data[i+1]]<<2) | (p[data[i+2]]<<4) | (p[data[i+3]]<<6);
	out[j] = 0;
	int s = len-i, x = 0;
	switch (s) {
	case 3: out[j] |= p[data[i++]] << x; x+=2;
	case 2: out[j] |= p[data[i++]] << x; x+=2;
	case 1: out[j] |= p[data[i++]] << x; x+=2;
	    j++;
	}
	*out_len = j;
	return out;
    }

    case 8: {
	for (i = 0; i < (len & ~7); i+=8)
	    out[j++] = (p[data[i+0]]<<0) | (p[data[i+1]]<<1) | (p[data[i+2]]<<2) | (p[data[i+3]]<<3)
		     | (p[data[i+4]]<<4) | (p[data[i+5]]<<5) | (p[data[i+6]]<<6) | (p[data[i+7]]<<7);
	out[j] = 0;
	int s = len-i, x = 0;
	switch (s) {
	case 7: out[j] |= p[data[i++]] << x++;
	case 6: out[j] |= p[data[i++]] << x++;
	case 5: out[j] |= p[data[i++]] << x++;
	case 4: out[j] |= p[data[i++]] << x++;
	case 3: out[j] |= p[data[i++]] << x++;
	case 2: out[j] |= p[data[i++]] << x++;
	case 1: out[j] |= p[data[i++]] << x++;
	    j++;
	}
	*out_len = j;
	return out;
    }

    case 0:
	*out_len = j;
	return out;
    }

    return NULL;
}

// expands the unpack meta data and returns the number of bytes read.
// nsym is number of symbols per byte, with the symbol values
// themselve sreturned in map.
//
// Returns number of bytes of data[] consumed on success,
//         zero on failure.
static uint8_t a_unpack_meta(uint8_t *data, uint32_t data_len,
			     uint64_t udata_len, uint8_t *map, int *nsym) {
    if (data_len <= 1)
	return 0;

    // Decode symbol map
    int j = 1, c = 0;
    int d = data[0];
    if (d <= 1)
	*nsym = 0;
    else if (d <= 2)
	*nsym = 8;
    else if (d <= 4)
	*nsym = 4;
    else if (d <= 16)
	*nsym = 2;
    else
	*nsym = 1;

    if (*nsym == 1)
	return 1;

    do {
	map[c++] = data[j++];
    } while (j-1 < d && c < 16 && j < data_len);

    return j-1 < d ? 0 : j;
}

static uint8_t *a_unpack(uint8_t *data, int64_t len, uint8_t *out, uint64_t out_len, int nsym, uint8_t *p) {
    //uint8_t *out;
    uint8_t c = 0;
    int64_t i, j = 0, olen;

    if (nsym == 1) {
	// raw data; FIXME: shortcut the need for malloc & memcpy here
	memcpy(out, data, len);
	return out;
    }

    // //SIMPLE decode, for clarity.
    // switch(nsym) {
    // case 8: {
    // 	unsigned char v = 0;
    // 	for (i = 0; i < out_len; i++) {
    // 	    if (i%8 == 0)
    // 		v = data[j++];
    // 	    out[i] = p[v & 1];
    // 	    v >>= 1;
    // 	}
    // 	break;
    // }
    // case 4: {
    // 	unsigned char v = 0;
    // 	for (i = 0; i < out_len; i++) {
    // 	    if (i%4 == 0)
    // 		v = data[j++];
    // 	    out[i] = p[v & 3];
    // 	    v >>= 2;
    // 	}
    // 	break;
    // }
    // case 2: {
    // 	unsigned char v = 0;
    // 	for (i = 0; i < out_len; i++) {
    // 	    if (i%2 == 0)
    // 		v = data[j++];
    // 	    out[i] = p[v & 15];
    // 	    v >>= 4;
    // 	}
    // 	break;
    // }
    // case 0:
    // 	memset(out, p[0], out_len);
    // 	break;
    // 
    // default:
    // 	return NULL;
    // }

    switch(nsym) {
    case 8: {
#ifdef ALLOW_UAC
	uint64_t map[256], x0, x1, x2, x3, x4, x5, x6, x7;
	int x;
	for (x = 0; x < 256; x++) {
	    map[x]=
		(((uint64_t)p[x>>0&1])<<0)+
		(((uint64_t)p[x>>1&1])<<8)+
		(((uint64_t)p[x>>2&1])<<16)+
		(((uint64_t)p[x>>3&1])<<24)+
		(((uint64_t)p[x>>4&1])<<32)+
		(((uint64_t)p[x>>5&1])<<40)+
		(((uint64_t)p[x>>6&1])<<48)+
		(((uint64_t)p[x>>7&1])<<56);
	}
#endif
	if ((out_len+7)/8 > len)
	    return NULL;
	olen = out_len & ~7;

#ifdef ALLOW_UAC
	for (i = 0; i < olen; i+=8) {
	    uint64_t w = map[data[j++]];
	    *(uint64_t *)&out[i] = w;
	}
#else
	for (i = j = 0; i < olen; i+=8) {
	    c = data[j++];
	    out[i+0] = p[(c>>0)&1];
	    out[i+1] = p[(c>>1)&1];
	    out[i+2] = p[(c>>2)&1];
	    out[i+3] = p[(c>>3)&1];
	    out[i+4] = p[(c>>4)&1];
	    out[i+5] = p[(c>>5)&1];
	    out[i+6] = p[(c>>6)&1];
	    out[i+7] = p[(c>>7)&1];
	}
#endif
	if (out_len != olen) {
	    c = data[j++];
	    while (i < out_len) {
		out[i++] = p[c & 1];
		c >>= 1;
	    }
	}
	break;
    }

    case 4: {
#ifdef ALLOW_UAC
	uint32_t map[256], x, y, z, _, P=0;
	for (x = 0; x < 4; x++)
	    for (y = 0; y < 4; y++)
		for (z = 0; z < 4; z++)
		    for (_ = 0; _ < 4; _++, P++)
			map[P] = (p[x]<<24)+(p[y]<<16)+(p[z]<<8)+(p[_]<<0);
#endif

//	// Yay a real world usage of Duff's Device!
//	// Sadly it's slower due to needing a reverse loop in this scenario.
//	i = out_len;
//	j += out_len/4;
//	c = data[j--];
//	switch(out_len % 4) {
//	    while (i) {
//	    case 0: out[--i] = p[(c>>0)&3];
//	    case 3: out[--i] = p[(c>>2)&3];
//	    case 2: out[--i] = p[(c>>4)&3];
//	    case 1: out[--i] = p[(c>>6)&3];
//		c = data[j--];
//	    }
//	}

	if ((out_len+3)/4 > len)
	    return NULL;
	olen = out_len & ~3;

#ifdef ALLOW_UAC
	for (; i < olen-12; i+=16) {
	    uint32_t w1 = map[data[j++]];
	    uint32_t w2 = map[data[j++]];
	    uint32_t w3 = map[data[j++]];
	    uint32_t w4 = map[data[j++]];
	    *(uint32_t *)&out[i   ] = w1;
	    *(uint32_t *)&out[i+ 4] = w2;
	    *(uint32_t *)&out[i+ 8] = w3;
	    *(uint32_t *)&out[i+12] = w4;
	}

	for (; i < olen; i+=4) {
	    uint32_t w = map[data[j++]];
	    *(uint32_t *)&out[i] = w;
	}
#else
	for (i = j = 0; i < olen; i+=4) {
	    c = data[j++];
	    out[i+0] = p[(c>>0)&3];
	    out[i+1] = p[(c>>2)&3];
	    out[i+2] = p[(c>>4)&3];
	    out[i+3] = p[(c>>6)&3];
	}
#endif
	if (out_len != olen) {
	    c = data[j++];
	    while (i < out_len) {
		out[i++] = p[c & 3];
		c >>= 2;
	    }
	}
	break;
    }

    case 2: {
#ifdef ALLOW_UAC
	uint16_t map[256], x, y;
	for (x = 0; x < 16; x++)
	    for (y = 0; y < 16; y++)
		map[x*16+y] = p[x]*256+p[y];
#endif

	if ((out_len+1)/2 > len)
	    return NULL;
	olen = out_len & ~1;
#ifdef ALLOW_UAC
	for (i = j = 0; i+2 < olen; i+=4) {
	    uint16_t w1 = map[data[j++]];
	    uint16_t w2 = map[data[j++]];
	    *(uint16_t *)&out[i  ] = w1;
	    *(uint16_t *)&out[i+2] = w2;
	}
	for (; i < olen; i+=2) {
	    uint16_t w1 = map[data[j++]];
	    *(uint16_t *)&out[i] = w1;
	}
#else
	for (i = j = 0; i < olen; i+=2) {
	    c = data[j++];
	    out[i+0] = p[(c>>0)&15];
	    out[i+1] = p[(c>>4)&15];
	}
#endif
	if (out_len != olen) {
	    c = data[j++];
	    out[i+0] = p[c&15];
	}
	break;
    }

    case 0:
	memset(out, p[0], out_len);
	break;

    default:
	return NULL;
    }

    return out;
}

//-----------------------------------------------------------------------------
unsigned char *arith_compress_O1(unsigned char *in, unsigned int in_size,
				 unsigned char *out, unsigned int *out_size) {
    int i, bound = arith_compress_bound(in_size,0)-5; // -5 for order/size

    if (!out) {
	*out_size = bound;
	out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
	return NULL;

    SIMPLE_MODEL(256,_) byte_model[256];
    unsigned int m = 0;
    if (1 || in_size > 1000) {
	for (i = 0; i < in_size; i++)
	    if (m < in[i])
		m = in[i];
	//fprintf(stderr, "%d max %d\n", in_size, m);
	m++;
    }
    *out = m;
    for (i = 0; i < 256; i++)
	SIMPLE_MODEL(256,_init)(&byte_model[i], m);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out+1);
    RC_StartEncode(&rc);

    uint8_t last = 0;
    for (i = 0; i < in_size; i++) {
	SIMPLE_MODEL(256, _encodeSymbol)(&byte_model[last], &rc, in[i]);
	last = in[i];
    }

    RC_FinishEncode(&rc);

    // Finalise block size and return it
    *out_size = RC_OutSize(&rc)+1;

    return out;
}

unsigned char *arith_uncompress_O1(unsigned char *in, unsigned int in_size,
				   unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;

    SIMPLE_MODEL(256,_) byte_model[256];
    unsigned int m = in[0] ? in[0] : 256, i;
    for (i = 0; i < 256; i++)
	SIMPLE_MODEL(256,_init)(&byte_model[i], m);
    
    if (!out)
	out = malloc(out_sz);
    if (!out)
	return NULL;

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    unsigned char last = 0;
    for (i = 0; i < out_sz; i++) {
	out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model[last], &rc);
	last = out[i];
    }

    RC_FinishDecode(&rc);
    
    return out;
}

//-----------------------------------------------------------------------------

// Disable O2 for now
#if 0

#if 0
unsigned char *arith_compress_O2(unsigned char *in, unsigned int in_size,
				 unsigned char *out, unsigned int *out_size) {
    fprintf(stderr, "WARNING: using undocumented O2 arith\n");

    int i, j;
    int bound = arith_compress_bound(in_size,0)-5; // -5 for order/size

    if (!out) {
	*out_size = bound;
	out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
	return NULL;

    unsigned int m = 0;
    if (1 || in_size > 1000) {
	for (i = 0; i < in_size; i++)
	    if (m < in[i])
		m = in[i];
	//fprintf(stderr, "%d max %d\n", in_size, m);
	m++;
    }
    *out = m;

    SIMPLE_MODEL(256,_) *byte_model;
    byte_model = malloc(256*256*sizeof(*byte_model));
    for (i = 0; i < 256; i++)
	for (j = 0; j < 256; j++)
	    SIMPLE_MODEL(256,_init)(&byte_model[i*256+j], m);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out+1);
    RC_StartEncode(&rc);

    unsigned char last1 = 0, last2 = 0;
    for (i = 0; i < in_size; i++) {
	SIMPLE_MODEL(256, _encodeSymbol)(&byte_model[last1*256 + last2], &rc, in[i]);
	last2 = last1;
	last1 = in[i];
    }

    free(byte_model);
    RC_FinishEncode(&rc);

    // Finalise block size and return it
    *out_size = RC_OutSize(&rc)+1;

    return out;
}
#else
unsigned char *arith_compress_O2(unsigned char *in, unsigned int in_size,
				 unsigned char *out, unsigned int *out_size) {
    fprintf(stderr, "WARNING: using undocumented O2 arith\n");

    int i, j;
    int bound = arith_compress_bound(in_size,0)-5; // -5 for order/size

    if (!out) {
	*out_size = bound;
	out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
	return NULL;

    unsigned int m = 0;
    if (1 || in_size > 1000) {
	for (i = 0; i < in_size; i++)
	    if (m < in[i])
		m = in[i];
	//fprintf(stderr, "%d max %d\n", in_size, m);
	m++;
    }
    *out = m;

    SIMPLE_MODEL(256,_) *byte_model;
    byte_model = malloc(256*256*sizeof(*byte_model));
    for (i = 0; i < 256; i++)
	for (j = 0; j < 256; j++)
	    SIMPLE_MODEL(256,_init)(&byte_model[i*256+j], m);
    SIMPLE_MODEL(256,_) byte_model1[256];
    for (i = 0; i < 256; i++)
	SIMPLE_MODEL(256,_init)(&byte_model1[i], m);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out+1);
    RC_StartEncode(&rc);

    unsigned char last1 = 0, last2 = 0;
    for (i = 0; i < in_size; i++) {
	// Use Order-1 is order-2 isn't sufficiently advanced yet (75+ symbols)
	if (byte_model[last1*256+last2].TotFreq <= m+75*16) {
	    SIMPLE_MODEL(256, _encodeSymbol)(&byte_model1[last1], &rc, in[i]);
	    SIMPLE_MODEL(256, _updateSymbol)(&byte_model[last1*256 + last2], &rc, in[i]);
	} else {
	    SIMPLE_MODEL(256, _encodeSymbol)(&byte_model[last1*256 + last2], &rc, in[i]);
	    //SIMPLE_MODEL(256, _updateSymbol)(&byte_model1[last1], &rc, in[i]);
	}
	last2 = last1;
	last1 = in[i];
    }

    free(byte_model);
    RC_FinishEncode(&rc);

    // Finalise block size and return it
    *out_size = RC_OutSize(&rc)+1;

    return out;
}
#endif

unsigned char *arith_uncompress_O2(unsigned char *in, unsigned int in_size,
				   unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;

    SIMPLE_MODEL(256,_) *byte_model;
    byte_model = malloc(256*256*sizeof(*byte_model));
    unsigned int m = in[0] ? in[0] : 256, i, j;
    for (i = 0; i < 256; i++)
	for (j = 0; j < 256; j++)
	    SIMPLE_MODEL(256,_init)(&byte_model[i*256+j], m);
    
    if (!out)
	out = malloc(out_sz);
    if (!out)
	return NULL;

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    unsigned char last1 = 0, last2 = 0;
    for (i = 0; i < out_sz; i++) {
	out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model[last1*256 + last2], &rc);
	last2 = last1;
	last1 = out[i];
    }

    free(byte_model);
    RC_FinishDecode(&rc);
    
    return out;
}

#endif // Disable O2
/*-----------------------------------------------------------------------------
 */

#undef NSYM
#define NSYM 258
#include "c_simple_model.h"
#define MAX_RUN 4

unsigned char *arith_compress_O0_RLE(unsigned char *in, unsigned int in_size,
				     unsigned char *out, unsigned int *out_size) {
    int i, bound = arith_compress_bound(in_size,0)-5; // -5 for order/size

    if (!out) {
	*out_size = bound;
	out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
	return NULL;

    unsigned int m = 0;
    for (i = 0; i < in_size; i++)
	if (m < in[i])
	    m = in[i];
    m++;
    *out = m;

    SIMPLE_MODEL(256,_) byte_model;
    SIMPLE_MODEL(256,_init)(&byte_model, m);

    SIMPLE_MODEL(NSYM,_) run_model[NSYM];
    for (i = 0; i < NSYM; i++)
	SIMPLE_MODEL(NSYM,_init)(&run_model[i], MAX_RUN);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out+1);
    RC_StartEncode(&rc);

    unsigned char last = 0;
    for (i = 0; i < in_size;) {
	//SIMPLE_MODEL(256, _encodeSymbol)(&byte_model, &rc, in[i]);
	SIMPLE_MODEL(256, _encodeSymbol)(&byte_model, &rc, in[i]);
	//fprintf(stderr, "lit %c (ctx %c)\n", in[i], last);
	int run = 0;
	last = in[i++];
	while (i < in_size && in[i] == last/* && run < MAX_RUN-1*/)
	    run++, i++;
	int rctx = last;
	do {
	    int c = run < MAX_RUN ? run : MAX_RUN-1;
	    SIMPLE_MODEL(NSYM, _encodeSymbol)(&run_model[rctx], &rc, c);
	    run -= c;

	    if (rctx == last)
		rctx = 256;
	    else
		rctx += (rctx < NSYM-1);
	    if (c == MAX_RUN-1 && run == 0)
		SIMPLE_MODEL(NSYM, _encodeSymbol)(&run_model[rctx], &rc, 0);
	} while (run);
    }

    RC_FinishEncode(&rc);

    // Finalise block size and return it
    *out_size = RC_OutSize(&rc)+1;

    //fprintf(stderr, "RLE %d to %d\n", in_size, *out_size);

    return out;
}

unsigned char *arith_uncompress_O0_RLE(unsigned char *in, unsigned int in_size,
				       unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;
    int i;
    unsigned int m = in[0] ? in[0] : 256;

    SIMPLE_MODEL(256,_) byte_model;
    SIMPLE_MODEL(256,_init)(&byte_model, m);

    SIMPLE_MODEL(NSYM,_) run_model[NSYM];
    for (i = 0; i < NSYM; i++)
	SIMPLE_MODEL(NSYM,_init)(&run_model[i], MAX_RUN);

    if (!out)
	out = malloc(out_sz);
    if (!out)
	return NULL;

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    for (i = 0; i < out_sz; i++) {
	unsigned char last;
	last = out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model, &rc);
	//fprintf(stderr, "lit %c\n", last);
	int run = 0, r = 0, rctx = out[i];
	do {
	    r = SIMPLE_MODEL(NSYM, _decodeSymbol)(&run_model[rctx], &rc);
	    if (rctx == last)
		rctx = 256;
	    else
		rctx += (rctx < NSYM-1);
	    //fprintf(stderr, "run %d (ctx %d, %d)\n", r, last, l);
	    run += r;
	} while (r == MAX_RUN-1 && run < out_sz);
	while (run-- && i+1 < out_sz)
	    out[++i] = last;
    }

    RC_FinishDecode(&rc);

    return out;
}

unsigned char *arith_compress_O1_RLE(unsigned char *in, unsigned int in_size,
				     unsigned char *out, unsigned int *out_size) {
    int i, bound = arith_compress_bound(in_size,0)-5; // -5 for order/size

    if (!out) {
	*out_size = bound;
	out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
	return NULL;

    unsigned int m = 0;
    for (i = 0; i < in_size; i++)
	if (m < in[i])
	    m = in[i];
    m++;
    *out = m;

    //SIMPLE_MODEL(256,_) byte_model;
    //SIMPLE_MODEL(256,_init)(&byte_model, m);

    SIMPLE_MODEL(256,_) byte_model[256];
    for (i = 0; i < 256; i++)
	SIMPLE_MODEL(256,_init)(&byte_model[i], m);

    SIMPLE_MODEL(NSYM,_) run_model[NSYM];
    for (i = 0; i < NSYM; i++)
	SIMPLE_MODEL(NSYM,_init)(&run_model[i], MAX_RUN);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out+1);
    RC_StartEncode(&rc);

    unsigned char last = 0;
    for (i = 0; i < in_size;) {
	//SIMPLE_MODEL(256, _encodeSymbol)(&byte_model, &rc, in[i]);
	SIMPLE_MODEL(256, _encodeSymbol)(&byte_model[last], &rc, in[i]);
	//fprintf(stderr, "lit %c (ctx %c)\n", in[i], last);
	int run = 0;
	last = in[i++];
	while (i < in_size && in[i] == last/* && run < MAX_RUN-1*/)
	    run++, i++;
	int rctx = last;
	do {
	    int c = run < MAX_RUN ? run : MAX_RUN-1;
	    SIMPLE_MODEL(NSYM, _encodeSymbol)(&run_model[rctx], &rc, c);
	    run -= c;

	    if (rctx == last)
		rctx = 256;
	    else
		rctx += (rctx < NSYM-1);
	    if (c == MAX_RUN-1 && run == 0)
		SIMPLE_MODEL(NSYM, _encodeSymbol)(&run_model[rctx], &rc, 0);
	} while (run);
    }

    RC_FinishEncode(&rc);

    // Finalise block size and return it
    *out_size = RC_OutSize(&rc)+1;

    //fprintf(stderr, "RLE %d to %d\n", in_size, *out_size);

    return out;
}

unsigned char *arith_uncompress_O1_RLE(unsigned char *in, unsigned int in_size,
				       unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;
    int i;
    unsigned int m = in[0] ? in[0] : 256;

    SIMPLE_MODEL(256,_) byte_model[256];
    for (i = 0; i < 256; i++)
	SIMPLE_MODEL(256,_init)(&byte_model[i], m);

    SIMPLE_MODEL(NSYM,_) run_model[NSYM];
    for (i = 0; i < NSYM; i++)
	SIMPLE_MODEL(NSYM,_init)(&run_model[i], MAX_RUN);

    if (!out)
	out = malloc(out_sz);
    if (!out)
	return NULL;

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    unsigned char last = 0;
    for (i = 0; i < out_sz; i++) {
	out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model[last], &rc);
	//fprintf(stderr, "lit %c (ctx %c)\n", out[i], last);
	last = out[i];
	int run = 0, r = 0, rctx = last;

	do {
	    r = SIMPLE_MODEL(NSYM, _decodeSymbol)(&run_model[rctx], &rc);
	    if (rctx == last)
		rctx = 256;
	    else
		rctx += (rctx < NSYM-1);
	    run += r;
	} while (r == MAX_RUN-1 && run < out_sz);
	while (run-- && i+1 < out_sz)
	    out[++i] = last;
    }

    RC_FinishDecode(&rc);

    return out;
}

/*-----------------------------------------------------------------------------
 * Simple interface to the order-0 vs order-1 encoders and decoders.
 *
 * Smallest is method, <in_size> <input>, so worst case 2 bytes longer.
 */
unsigned char *arith_compress_to(unsigned char *in,  unsigned int in_size,
				 unsigned char *out, unsigned int *out_size,
				 int order) {
    unsigned int c_meta_len;
    uint8_t *rle = NULL, *packed = NULL;

    if (!out) {
	*out_size = arith_compress_bound(in_size, order);
	if (!(out = malloc(*out_size)))
	    return NULL;
    }

    if (in_size%4 != 0 || in_size <= 20)
	order &= ~X_4;

    if (order & X_CAT) {
	out[0] = X_CAT;
	c_meta_len = 1 + u32tou7(&out[1], in_size);
	memcpy(out+c_meta_len, in, in_size);
	*out_size = in_size+c_meta_len;
    }

    if (order & X_4) {
	unsigned char *in4 = malloc(in_size);
	if (!in4)
	    return NULL;
	unsigned int len4 = in_size/4, i4[4];
	int i;

	for (i = 0; i < 4; i++)
	    i4[i] = i*len4;
	for (i = 0; i4[0] < len4; i4[0]++, i4[1]++, i4[2]++, i4[3]++, i+=4) {
	    in4[i4[0]] = in[i+0];
	    in4[i4[1]] = in[i+1];
	    in4[i4[2]] = in[i+2];
	    in4[i4[3]] = in[i+3];
	}

	unsigned int olen2;
	unsigned char *out2;
	c_meta_len = 1;
	*out = order;
	c_meta_len += u32tou7(out+c_meta_len, in_size);
	out2 = out+26;
	for (i = 0; i < 4; i++) {
	    // Brute force try all methods.
	    // FIXME: optimise this bit.  Maybe learn over time?
	    int j, best_j = 0, best_sz = INT_MAX;

	    // Works OK with read names. The first byte is the most important,
	    // as it has most variability (little-endian).  After that it's
	    // often quite predictable.
	    //
	    // Do we gain in any other context in CRAM? Aux tags maybe?
	    int m[][4] = {{3, 1,64,0},
			  {2, 1,0},
			  {2, 1,128},
			  {2, 1,128}};

//	    int m[][6] = {{4, 1,64,2,0},  //test of adding in an order-2 codec
//			  {3, 1,2,0},
//			  {3, 1,2,128},
//			  {3, 1,2,128}};

// Other possibilities for methods to try.
//	    int m[][10] = {{8, 1,128,129,64,65,192,193,4,0},
//			   {8, 1,128,129,64,65,192,193,4,0},
//			   {8, 1,128,129,64,65,192,193,4,0},
//			   {8, 1,128,129,64,65,192,193,4,0}};

//	    int m[][9] = {{5, 1,128,64,65,0},
//			  {5, 1,128,64,65,0},
//			  {5, 1,128,64,65,0},
//			  {5, 1,128,64,65,0}};

//	    int m[][6] = {{4, 0,1,128,64},
//			  {5, 0,1,128,65,193},
//			  {3, 0,1,128},
//			  {3, 0,1,128}};

//	    int m[][6] = {{4, 1,128,64,0},
//			  {4, 1,128,65,0},
//			  {2, 128,0},
//			  {2, 128,0}};

//	    int m[][6] = {{2, 64,0},
//			  {1, 0},
//			  {1, 128},
//			  {1, 128}};

//	    int m[][6] = {{1, 0},
//			  {2, 128,0},
//			  {1, 128},
//			  {1, 128}};

	    for (j = 1; j <= m[i][0]; j++) {
		olen2 = *out_size - (out2 - out);
		//fprintf(stderr, "order=%d m=%d\n", order&3, m[i][j]);
		if ((order&3) == 0 && (m[i][j]&1))
		    continue;

		arith_compress_to(in4+i*len4, len4, out2, &olen2, m[i][j] | X_NOSZ);
		if (best_sz > olen2) {
		    best_sz = olen2;
		    best_j = j;
		}
	    }
//	    if (best_j == 0) // none desireable
//		return NULL;
	    if (best_j != j-1) {
		olen2 = *out_size - (out2 - out);
		if (!arith_compress_to(in4+i*len4, len4, out2, &olen2, m[i][best_j] | X_NOSZ))
		    return NULL;
	    }
	    out2 += olen2;
	    c_meta_len += u32tou7(out+c_meta_len, olen2);
	}
	memmove(out+c_meta_len, out+26, out2-(out+26));
	free(in4);
	*out_size = c_meta_len + out2-(out+26);
	return out;
    }

    int do_pack = order & X_PACK;
    int do_rle  = order & X_RLE;
    int no_size = order & X_NOSZ;
    int do_ext  = order & X_EXT;

    out[0] = order;
    c_meta_len = 1;

    if (!no_size)
	c_meta_len += u32tou7(&out[1], in_size);

    order &= 0x3;

    // Format is compressed meta-data, compressed data.
    // Meta-data can be empty, pack, rle lengths, or pack + rle lengths.
    // Data is either the original data, bit-packed packed, rle literals or
    // packed + rle literals.

    if (do_pack && in_size) {
	// PACK 2, 4 or 8 symbols into one byte.
	int pmeta_len;
	uint64_t packed_len;
	packed = a_pack(in, in_size, out+c_meta_len, &pmeta_len, &packed_len);
	if (!packed || (pmeta_len == 1 && out[c_meta_len] == 1)) {
	    out[0] &= ~X_PACK;
	    do_pack = 0;
	    free(packed);
	    packed = NULL;
	} else {
	    in = packed;
	    in_size = packed_len;
	    c_meta_len += pmeta_len;

	    // Could derive this rather than storing verbatim.
	    // Orig size * 8/nbits (+1 if not multiple of 8/n)
	    int sz = u32tou7(out+c_meta_len, in_size);
	    c_meta_len += sz;
	    *out_size -= sz;
	}
    } else if (do_pack) {
	out[0] &= ~X_PACK;
    }

    if (do_rle && !in_size) {
	out[0] &= ~X_RLE;
    }

    *out_size -= c_meta_len;
    if (order && in_size < 8) {
	out[0] &= ~3;
	order  &= ~3;
    }

    if (do_ext) {
	// Use an external compression library instead.
	// For now, bzip2
	if (BZ_OK != BZ2_bzBuffToBuffCompress((char *)out+c_meta_len, out_size,
					      (char *)in, in_size, 9, 0, 30))
	    *out_size = in_size; // Didn't fit with bz2; force X_CAT below instead

//      // lzma doesn't help generally, at least not for the name tokeniser
//	size_t lzma_size = 0;
//	lzma_easy_buffer_encode(9, LZMA_CHECK_CRC32, NULL,
//				in, in_size, out+c_meta_len, &lzma_size,
//				*out_size);
//	*out_size = lzma_size;

    } else {
	if (do_rle) {
	    if (order == 0)
		arith_compress_O0_RLE(in, in_size, out+c_meta_len, out_size);
	    else
		arith_compress_O1_RLE(in, in_size, out+c_meta_len, out_size);
	} else {
	    //if (order == 2)
	    //	arith_compress_O2(in, in_size, out+c_meta_len, out_size);
	    //else
	    if (order == 1)
		arith_compress_O1(in, in_size, out+c_meta_len, out_size);
	    else
		arith_compress_O0(in, in_size, out+c_meta_len, out_size);
	}
    }

    if (*out_size >= in_size) {
	out[0] &= ~(3|X_EXT); // no entropy encoding, but keep e.g. PACK
	out[0] |= X_CAT | no_size;
	memcpy(out+c_meta_len, in, in_size);
	*out_size = in_size;
    }

    free(rle);
    free(packed);

    *out_size += c_meta_len;

    return out;
}

unsigned char *arith_compress(unsigned char *in, unsigned int in_size,
			      unsigned int *out_size, int order) {
    return arith_compress_to(in, in_size, NULL, out_size, order);
}

unsigned char *arith_uncompress_to(unsigned char *in,  unsigned int in_size,
				   unsigned char *out, unsigned int *out_size) {
    unsigned char *in_end = in + in_size;
    unsigned char *out_free = NULL;
    unsigned char *tmp_free = NULL;

    if (in_size == 0)
	return NULL;

    if (*in & X_4) {
	unsigned int ulen, olen, clen4[4];
	int c_meta_len = 1, i, j;

	// Decode lengths
	c_meta_len += u7tou32(in+c_meta_len, in_end, &ulen);
	if (!out) {
	    if (!(out_free = out = malloc(ulen)))
		return NULL;
	    *out_size = ulen;
	}
	if (ulen != *out_size || (ulen%4 != 0)) {
	    free(out_free);
	    return NULL;
	}

	for (i = 0; i < 4; i++)
	    c_meta_len += u7tou32(in+c_meta_len, in_end, &clen4[i]);

	//fprintf(stderr, "    x4 meta %d\n", c_meta_len); //c-size

	// Uncompress the 4 streams
	unsigned char *out4 = malloc(ulen);
	if (!out4) {
	    free(out_free);
	    return NULL;
	}
	for (i = 0; i < 4; i++) {
	    olen = ulen/4;
	    if (in_size < c_meta_len) {
		free(out_free);
		free(out4);
		return NULL;
	    }
	    if (!arith_uncompress_to(in+c_meta_len, in_size-c_meta_len, out4 + i*(ulen/4), &olen)
		|| olen != ulen/4) {
		free(out_free);
		free(out4);
		return NULL;
	    }
	    c_meta_len += clen4[i];
	}

	unsigned int i4[4] = {0*(ulen/4), 1*(ulen/4), 2*(ulen/4), 3*(ulen/4)};
	j = 0;
	while (j < ulen) {
	    out[j++] = out4[i4[0]++];
	    out[j++] = out4[i4[1]++];
	    out[j++] = out4[i4[2]++];
	    out[j++] = out4[i4[3]++];
	}
	free(out4);
	*out_size = ulen;
	return out;
    }

    int order = *in++;  in_size--;
    int do_pack = order & X_PACK;
    int do_rle  = order & X_RLE;
    int do_cat  = order & X_CAT;
    int no_size = order & X_NOSZ;
    int do_ext  = order & X_EXT;
    order &= 3;

    int sz = 0;
    unsigned int osz;
    if (!no_size)
	sz = u7tou32(in, in_end, &osz);
    else
	sz = 0, osz = *out_size;
    in += sz;
    in_size -= sz;

    if (osz >= INT_MAX)
	return NULL;

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
	// Limit maximum size to get fast turnaround on fuzzing test cases
	if (osz > 100000)
	    goto err;
#endif

    if (no_size && !out)
	return NULL; // Need one or the other

    if (!out) {
	*out_size = osz;
	if (!(out_free = out = malloc(*out_size)))
	    return NULL;
    } else {
	if (*out_size < osz)
	    return NULL;
	*out_size = osz;
    }

    uint32_t c_meta_size = 0;
    unsigned int tmp1_size = *out_size;
    unsigned int tmp2_size = *out_size;
    unsigned char *tmp1 = NULL, *tmp2 = NULL, *tmp = NULL;

    // Need In, Out and Tmp buffers with temporary buffer of the same size
    // as output.  Our entropy decode is either arithmetic (with/without RLE)
    // or external (bz2, gzip, lzma) but with an optional unPACK transform
    // at the end.
    //
    // To avoid pointless memcpy when unpacking we switch around which
    // buffers we're writing to accordingly.

    // Format is pack meta data if present, followed by compressed data.
    if (do_pack) {
	if (!(tmp_free = tmp = malloc(*out_size)))
	    goto err;
	tmp1 = tmp;  // uncompress
	tmp2 = out;  // unpack
    } else {
	// no pack
	tmp  = NULL;
	tmp1 = out;  // uncompress
	tmp2 = out;  // NOP
    }

    
    // Decode the bit-packing map.
    uint8_t map[16] = {0};
    int npacked_sym = 0;
    uint64_t unpacked_sz = 0; // FIXME: rename to packed_per_byte
    if (do_pack) {
	c_meta_size = a_unpack_meta(in, in_size, *out_size, map, &npacked_sym);
	if (c_meta_size == 0)
	    goto err;;

	unpacked_sz = osz;
	in      += c_meta_size;
	in_size -= c_meta_size;

	// New unpacked size.  We could derive this bit from *out_size
	// and npacked_sym.
	unsigned int osz;
	sz = u7tou32(in, in_end, &osz);
	in += sz;
	in_size -= sz;
	if (osz > tmp1_size)
	    goto err;;
	tmp1_size = osz;
    }

    //fprintf(stderr, "    meta_size %d bytes\n", (int)(in - orig_in)); //c-size

    // uncompress RLE data.  in -> tmp1
    if (in_size) {
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
	// Limit maximum size to get fast turnaround on fuzzing test cases
	if (tmp1_size > 100000)
	    goto err;
#endif
	if (do_cat) {
	    //fprintf(stderr, "    CAT %d\n", tmp1_size); //c-size
	    if (tmp1_size > in_size)
		goto err;;
	    if (tmp1_size > *out_size)
		goto err;;
	    memcpy(tmp1, in, tmp1_size);
	} else if (do_ext) {
	    if (BZ_OK != BZ2_bzBuffToBuffDecompress((char *)tmp1, &tmp1_size,
						    (char *)in, in_size, 0, 0))
		goto err;;
	} else {
	    // in -> tmp1
	    if (do_rle) {
		tmp1 = order == 1
		    ? arith_uncompress_O1_RLE(in, in_size, tmp1, tmp1_size)
		    : arith_uncompress_O0_RLE(in, in_size, tmp1, tmp1_size);
	    } else {
		//if (order == 2)
		//    tmp1 = arith_uncompress_O2(in, in_size, tmp1, tmp1_size)
		//else
		tmp1 = order == 1
		    ? arith_uncompress_O1(in, in_size, tmp1, tmp1_size)
		    : arith_uncompress_O0(in, in_size, tmp1, tmp1_size);
	    }
	    if (!tmp1)
		goto err;;
	}
    } else {
	tmp1 = NULL;
	tmp1_size = 0;
    }

    if (do_pack) {
	// Unpack bits via pack-map.  tmp1 -> tmp2
	if (npacked_sym == 1)
	    unpacked_sz = tmp1_size;
	//uint8_t *porig = unpack(tmp2, tmp2_size, unpacked_sz, npacked_sym, map);
	//memcpy(tmp3, porig, unpacked_sz);
	if (!a_unpack(tmp1, tmp1_size, tmp2, unpacked_sz, npacked_sym, map))
	    goto err;;
	tmp2_size = unpacked_sz;
    } else {
	tmp2_size = tmp1_size;
    }

    if (tmp)
	free(tmp);

    *out_size = tmp2_size;
    return tmp2;

 err:
    free(tmp_free);
    free(out_free);
    return NULL;
}

unsigned char *arith_uncompress(unsigned char *in, unsigned int in_size,
				unsigned int *out_size) {
    return arith_uncompress_to(in, in_size, NULL, out_size);
}


/*-----------------------------------------------------------------------------
 * Main
 */
// A simple test harness for compression and decompression via the command line.
//
// This also permits us to fuzz test the decoder on random (invalid) inputs.
#ifdef TEST_MAIN

#ifndef BLK_SIZE
#  define BLK_SIZE 1013*1047
#endif

// Room to allow for expanded BLK_SIZE on worst case compression.
#define BLK_SIZE2 ((105LL*BLK_SIZE)/100)

static unsigned char in_buf[BLK_SIZE2+257*257*3];

int main(int argc, char **argv) {
    int opt, order = 0;
    int decode = 0, test = 0;
    FILE *infp = stdin, *outfp = stdout;
    struct timeval tv1, tv2, tv3, tv4;
    size_t bytes = 0;

    extern char *optarg;
    extern int optind;

    while ((opt = getopt(argc, argv, "o:dt")) != -1) {
	switch (opt) {
	case 'o':
	    order = atoi(optarg);
	    break;

	case 'd':
	    decode = 1;
	    break;
	    
	case 't':
	    test = 1;
	    break;
	}
    }

    //order = order ? 1 : 0; // Only support O(0) and O(1)

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

    if (test) {
	size_t len, in_sz = 0, out_sz = 0;
	typedef struct {
	    unsigned char *blk;
	    uint32_t sz;
	} blocks;
	blocks *b = NULL, *bc = NULL, *bu = NULL;
	int nb = 0, i;
	
	while ((len = fread(in_buf, 1, BLK_SIZE, infp)) != 0) {
	    // inefficient, but it'll do for testing
	    b = realloc(b, (nb+1)*sizeof(*b));
	    bu = realloc(bu, (nb+1)*sizeof(*bu));
	    bc = realloc(bc, (nb+1)*sizeof(*bc));
	    b[nb].blk = malloc(len);
	    b[nb].sz = len;
	    memcpy(b[nb].blk, in_buf, len);
	    bc[nb].sz = arith_compress_bound(BLK_SIZE, order);
	    bc[nb].blk = malloc(bc[nb].sz);
	    bu[nb].sz = len;
	    bu[nb].blk = malloc(BLK_SIZE);
	    nb++;
	    in_sz += len;
	}
	fprintf(stderr, "Testing %d blocks\n", nb);

#ifndef NTRIALS
#define NTRIALS 10
#endif
	int trials = NTRIALS;
	while (trials--) {
	    // Warmup
	    for (i = 0; i < nb; i++) memset(bc[i].blk, 0, bc[i].sz);

	    gettimeofday(&tv1, NULL);

	    out_sz = 0;
	    for (i = 0; i < nb; i++) {
		unsigned int csz = bc[i].sz;
		bc[i].blk = arith_compress_to(b[i].blk, b[i].sz, bc[i].blk, &csz, order);
		assert(csz <= bc[i].sz);
		out_sz += 5 + csz;
	    }

	    gettimeofday(&tv2, NULL);
	    
	    // Warmup
	    for (i = 0; i < nb; i++) memset(bu[i].blk, 0, BLK_SIZE);

	    gettimeofday(&tv3, NULL);

	    for (i = 0; i < nb; i++)
		bu[i].blk = arith_uncompress_to(bc[i].blk, bc[i].sz, bu[i].blk, &bu[i].sz);

	    gettimeofday(&tv4, NULL);

	    for (i = 0; i < nb; i++) {
		if (b[i].sz != bu[i].sz || memcmp(b[i].blk, bu[i].blk, b[i].sz))
		    fprintf(stderr, "Mismatch in block %d, sz %d/%d\n", i, b[i].sz, bu[i].sz);
		//free(bc[i].blk);
		//free(bu[i].blk);
	    }

	    fprintf(stderr, "%5.1f MB/s enc, %5.1f MB/s dec\t %lld bytes -> %lld bytes\n",
		    (double)in_sz / ((long)(tv2.tv_sec - tv1.tv_sec)*1000000 +
				     tv2.tv_usec - tv1.tv_usec),
		    (double)in_sz / ((long)(tv4.tv_sec - tv3.tv_sec)*1000000 +
				     tv4.tv_usec - tv3.tv_usec),
		    (long long)in_sz, (long long)out_sz);
	}

	exit(0);
	
    }

    if (decode) {
	// Only used in some test implementations of RC_GetFreq()
	//RC_init();
	//RC_init2();

	for (;;) {
	    uint32_t in_size, out_size;
	    unsigned char *out;

	    if (4 != fread(&in_size, 1, 4, infp))
		break;
	    if (in_size > BLK_SIZE)
		exit(1);

	    if (in_size != fread(in_buf, 1, in_size, infp)) {
		fprintf(stderr, "Truncated input\n");
		exit(1);
	    }
	    out = arith_uncompress(in_buf, in_size, &out_size);
	    if (!out)
		exit(1);

	    fwrite(out, 1, out_size, outfp);
	    fflush(outfp);
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

	    if (in_size < 4)
		order &= ~1;

	    out = arith_compress(in_buf, in_size, &out_size, order);

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


// Encode, decode, compare and abort is different.
//
// This is designed for use within AFL to check all inputs can be round-tripped.
#ifdef TEST_MAIN2

#ifndef BLK_SIZE
#  define BLK_SIZE 1013*1047
#endif

// Room to allow for expanded BLK_SIZE on worst case compression.
#define BLK_SIZE2 ((105LL*BLK_SIZE)/100)

/*-----------------------------------------------------------------------------
 * Main
 */
static unsigned char in_buf[BLK_SIZE2+257*257*3];

int main(int argc, char **argv) {
    int opt, order = 0;
    int decode = 0, test = 0;
    FILE *infp = stdin, *outfp = stdout;
    struct timeval tv1, tv2, tv3, tv4;
    size_t bytes = 0;

    extern char *optarg;
    extern int optind;

    //order = order ? 1 : 0; // Only support O(0) and O(1)

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

    if (1) {
	size_t len, in_sz = 0, out_sz = 0;
	typedef struct {
	    unsigned char *blk;
	    uint32_t sz;
	} blocks;
	blocks *b = NULL, *bc = NULL, *bu = NULL;
	int nb = 0, i;

	while ((len = fread(in_buf, 1, BLK_SIZE, infp)) != 0) {
	    // inefficient, but it'll do for testing
	    b = realloc(b, (nb+1)*sizeof(*b));
	    bu = realloc(bu, (nb+1)*sizeof(*bu));
	    bc = realloc(bc, (nb+1)*sizeof(*bc));
	    b[nb].blk = malloc(len);
	    b[nb].sz = len;
	    memcpy(b[nb].blk, in_buf, len);
	    bc[nb].sz = arith_compress_bound(BLK_SIZE, order);
	    bc[nb].blk = malloc(bc[nb].sz);
	    bu[nb].sz = BLK_SIZE;
	    bu[nb].blk = malloc(BLK_SIZE);
	    nb++;
	    in_sz += len;
	}
	fprintf(stderr, "Testing %d blocks\n", nb);

	int O, order_map[10] = {0,1,64,65,128,129,192,193,4,128+4};
	for (O=0; O<40; O++) {
	    order = order_map[O%10];

	    order |= ((O/10)&1) ? X_NOSZ : 0;
	    order |= ((O/10)&2) ? X_4    : 0;

	    out_sz = 0;
	    for (i = 0; i < nb; i++) {
		unsigned int csz = bc[i].sz;
		bc[i].blk = arith_compress_to(b[i].blk, b[i].sz, bc[i].blk, &csz, order);
		assert(csz <= bc[i].sz);
		out_sz += 5 + csz;
		fprintf(stderr, "%08x C %d -> %d\n", order, b[i].sz, csz);
	    }
	    for (i = 0; i < nb; i++) {
		bu[i].blk = arith_uncompress_to(bc[i].blk, bc[i].sz,
						bu[i].blk, &bu[i].sz);
		fprintf(stderr, "%08x D %d -> %d\n", order, bc[i].sz, bu[i].sz);
	    }

	    for (i = 0; i < nb; i++) {
		if (b[i].sz != bu[i].sz || memcmp(b[i].blk, bu[i].blk, b[i].sz)) {
		    fprintf(stderr, "Mismatch in block %d, sz %d/%d\n", i, b[i].sz, bu[i].sz);
		    abort();
		}
	    }
	}
    }

    return 0;
}
#endif

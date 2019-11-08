/*
 * Copyright (c) 2013, 2014, 2015 Genome Research Ltd.
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

/*
 * FIXME: add checking of cram_external_type to return NULL on unsupported
 * {codec,type} tuples.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <stddef.h>
#include <htscodecs/varint.h>

#include "io_lib/cram.h"

/*
 * ---------------------------------------------------------------------------
 * Block bit-level I/O functions.
 * All defined static here to promote easy inlining by the compiler.
 */

#if 0
/* Get a single bit, MSB first */
static signed int get_bit_MSB(cram_block *block) {
    unsigned int val;

    if (block->byte > block->alloc)
	return -1;

    val = block->data[block->byte] >> block->bit;
    if (--block->bit == -1) {
	block->bit = 7;
	block->byte++;
	//printf("(%02X)", block->data[block->byte]);
    }

    //printf("-B%d-", val&1);

    return val & 1;
}
#endif

/*
 * Count number of successive 0 and 1 bits
 */
static int get_one_bits_MSB(cram_block *block) {
    int n = 0, b;
    if (block->byte >= block->uncomp_size)
        return -1;
    do {
	b = block->data[block->byte] >> block->bit;
	if (--block->bit == -1) {
	    block->bit = 7;
	    block->byte++;
	    if (block->byte == block->uncomp_size && (b&1))
	        return -1;
	}
	n++;
    } while (b&1);

    return n-1;
}

static int get_zero_bits_MSB(cram_block *block) {
    int n = 0, b;
    if (block->byte >= block->uncomp_size)
        return -1;
    do {
	b = block->data[block->byte] >> block->bit;
	if (--block->bit == -1) {
	    block->bit = 7;
	    block->byte++;
	    if (block->byte == block->uncomp_size && !(b&1))
	        return -1;
	}
	n++;
    } while (!(b&1));

    return n-1;
}

#if 0
/* Stores a single bit */
static void store_bit_MSB(cram_block *block, unsigned int bit) {
    if (block->byte >= block->alloc) {
	block->alloc = block->alloc ? block->alloc*2 : 1024;
	block->data = realloc(block->data, block->alloc);
    }

    if (bit)
	block->data[block->byte] |= (1 << block->bit);

    if (--block->bit == -1) {
	block->bit = 7;
	block->byte++;
	block->data[block->byte] = 0;
    }
}
#endif

#if 0
/* Rounds to the next whole byte boundary first */
static void store_bytes_MSB(cram_block *block, char *bytes, int len) {
    if (block->bit != 7) {
	block->bit = 7;
	block->byte++;
    }

    while (block->byte + len >= block->alloc) {
	block->alloc = block->alloc ? block->alloc*2 : 1024;
	block->data = realloc(block->data, block->alloc);
    }

    memcpy(&block->data[block->byte], bytes, len);
    block->byte += len;
}
#endif

/* Local optimised copy for inlining */
static inline int64_t get_bits_MSB(cram_block *block, int nbits) {
    uint64_t val = 0;
    int i;

#if 0
    // Fits within the current byte */
    if (nbits <= block->bit+1) {
	val = (block->data[block->byte]>>(block->bit-(nbits-1))) & ((1<<nbits)-1);
	if ((block->bit -= nbits) == -1) {
	    block->bit = 7;
	    block->byte++;
	}
	return val;
    }

    // partial first byte
    val = block->data[block->byte] & ((1<<(block->bit+1))-1);
    nbits -= block->bit+1;
    block->bit = 7;
    block->byte++;

    // whole middle bytes
    while (nbits >= 8) {
	val = (val << 8) | block->data[block->byte++];
	nbits -= 8;
    }

    val <<= nbits;
    val |= (block->data[block->byte]>>(block->bit-(nbits-1))) & ((1<<nbits)-1);
    block->bit -= nbits;
    return val;
#endif

#if 0
    /* Inefficient implementation! */
    //printf("{");
    for (i = 0; i < nbits; i++)
	//val = (val << 1) | get_bit_MSB(block);
	GET_BIT_MSB(block, val);
#endif

#if 1
    /* Combination of 1st two methods */
    if (nbits <= block->bit+1) {
	val = (block->data[block->byte]>>(block->bit-(nbits-1))) & ((1<<nbits)-1);
	if ((block->bit -= nbits) == -1) {
	    block->bit = 7;
	    block->byte++;
	}
	return val;
    }

//    /* Consume as many as possible from current byte */
//    val = block->data[block->byte] & ((1<<(block->bit+1))-1);
//    nbits -= block->bit+1;
//    block->bit = 7;
//    block->byte++;

    switch(nbits) {
//    case 15: GET_BIT_MSB(block, val);
//    case 14: GET_BIT_MSB(block, val);
//    case 13: GET_BIT_MSB(block, val);
//    case 12: GET_BIT_MSB(block, val);
//    case 11: GET_BIT_MSB(block, val);
//    case 10: GET_BIT_MSB(block, val);
//    case  9: GET_BIT_MSB(block, val);
    case  8: GET_BIT_MSB(block, val);
    case  7: GET_BIT_MSB(block, val);
    case  6: GET_BIT_MSB(block, val);
    case  5: GET_BIT_MSB(block, val);
    case  4: GET_BIT_MSB(block, val);
    case  3: GET_BIT_MSB(block, val);
    case  2: GET_BIT_MSB(block, val);
    case  1: GET_BIT_MSB(block, val);
	break;

    default:
	for (i = 0; i < nbits; i++)
	    //val = (val << 1) | get_bit_MSB(block);
	    GET_BIT_MSB(block, val);
    }
#endif

    //printf("=0x%x}", val);

    return val;
}

/*
 * Can store up to 24-bits worth of data encoded in an integer value
 * Possibly we'd want to have a less optimal store_bits function when dealing
 * with nbits > 24, but for now we assume the codes generated are never
 * that big. (Given this is only possible with 121392 or more
 * characters with exactly the correct frequency distribution we check
 * for it elsewhere.)
 */
static int store_bits_MSB(cram_block *block, uint64_t val, int nbits) {
    /* fprintf(stderr, " store_bits: %02x %d\n", val, nbits); */

    /*
     * Use slow mode until we tweak the huffman generator to never generate
     * codes longer than 24-bits.
     */
    unsigned int mask;

    if (block->byte+8 >= block->alloc) {
	if (block->byte) {
	    block->alloc *= 2;
	    block->data = realloc(block->data, block->alloc + 8);
	    if (!block->data)
		return -1;
	} else {
	    block->alloc = 1024;
	    block->data = realloc(block->data, block->alloc + 8);
	    if (!block->data)
		return -1;
	    block->data[0] = 0; // initialise first byte of buffer
	}
    }

    /* fits in current bit-field */
    if (nbits <= block->bit+1) {
	block->data[block->byte] |= (val << (block->bit+1-nbits));

	if ((block->bit-=nbits) == -1) {
	    block->bit = 7;
	    block->byte++;
	    block->data[block->byte] = 0;
	}
	return 0;
    }

    block->data[block->byte] |= (val >> (nbits -= block->bit+1));
    block->bit = 7;
    block->byte++;
    block->data[block->byte] = 0;
				 
    mask = 1<<(nbits-1);
    do {
	if (val & mask)
	    block->data[block->byte] |= (1 << block->bit);
	if (--block->bit == -1) {
	    block->bit = 7;
	    block->byte++;
	    block->data[block->byte] = 0;
	}
	mask >>= 1;
    } while(--nbits);

    return 0;
}

/*
 * Returns the next 'size' bytes from a block, or NULL if insufficient
 * data left.This is just a pointer into the block data and not an
 * allocated object, so do not free the result.
 */
static char *cram_extract_block(cram_block *b, int size) {
    char *cp = (char *)b->data + b->idx;
    b->idx += size;
    if (b->idx > b->uncomp_size)
	return NULL;

    return cp;
}

/*
 * ---------------------------------------------------------------------------
 * EXTERNAL
 */
// static int inline safe_itf8_geti(const char *cp, const char *endp, int32_t *val_p) {
//     const unsigned char *up = (unsigned char *)cp;
// 
//     if (endp - cp < 5 &&
// 	(cp >= endp || endp - cp < itf8_bytes[up[0]>>4])) {
//         *val_p = 0;
//         return 0;
//     }
// 
//     if (up[0] < 0x80) {
// 	*val_p =   up[0];
// 	return 1;
//     } else if (up[0] < 0xc0) {
// 	*val_p = ((up[0] <<8) |  up[1])                           & 0x3fff;
// 	return 2;
//     } else if (up[0] < 0xe0) {
// 	*val_p = ((up[0]<<16) | (up[1]<< 8) |  up[2])             & 0x1fffff;
// 	return 3;
//     } else if (up[0] < 0xf0) {
// 	*val_p = ((up[0]<<24) | (up[1]<<16) | (up[2]<<8) | up[3]) & 0x0fffffff;
// 	return 4;
//     } else {
// 	*val_p = ((up[0] & 0x0f)<<28) | (up[1]<<20) | (up[2]<<12) | (up[3]<<4) | (up[4] & 0x0f);
// 	return 5;
//     }
// }

int cram_external_decode_int(cram_slice *slice, cram_codec *c,
			     cram_block *in, char *out, int *out_size) {
    char *cp;
    cram_block *b;

    /* Find the external block */
    b = cram_get_block_by_id(slice, c->external.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = (char *)b->data + b->idx;
    // E_INT and E_LONG are guaranteed single item queries
    int err = 0;
    *(int32_t *)out = c->vv->varint_get32(&cp, (char *)b->data + b->uncomp_size, &err);
    b->idx = cp - (char *)b->data;
    *out_size = 1;

    return err ? -1 : 0;
}

int cram_external_decode_long(cram_slice *slice, cram_codec *c,
			      cram_block *in, char *out, int *out_size) {
    char *cp;
    cram_block *b;

    /* Find the external block */
    b = cram_get_block_by_id(slice, c->external.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = (char *)b->data + b->idx;
    // E_INT and E_LONG are guaranteed single item queries
    int err = 0;
    *(int64_t *)out = c->vv->varint_get64(&cp, (char *)b->data + b->uncomp_size, &err);
    b->idx = cp - (char *)b->data;
    *out_size = 1;

    return err ? -1 : 0;
}

int cram_external_decode_char(cram_slice *slice, cram_codec *c,
			      cram_block *in, char *out,
			      int *out_size) {
    char *cp;
    cram_block *b;

    /* Find the external block */
    b = cram_get_block_by_id(slice, c->external.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = cram_extract_block(b, *out_size);
    if (!cp)
	return -1;

    if (out)
	memcpy(out, cp, *out_size);
    return 0;
}

static int cram_external_decode_block(cram_slice *slice, cram_codec *c,
				      cram_block *in, char *out_,
				      int *out_size) {
    char *cp;
    cram_block *out = (cram_block *)out_;
    cram_block *b = NULL;

    /* Find the external block */
    b = cram_get_block_by_id(slice, c->external.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = cram_extract_block(b, *out_size);
    if (!cp)
	return -1;

    BLOCK_APPEND(out, cp, *out_size);
    return 0;
}

void cram_external_decode_free(cram_codec *c) {
    if (c)
	free(c);
}

int cram_external_decode_size(cram_slice *slice, cram_codec *c) {
    cram_block *b;

    /* Find the external block */
    b = cram_get_block_by_id(slice, c->external.content_id);
    if (!b)
        return -1;

    return b->uncomp_size;
}

cram_block *cram_external_get_block(cram_slice *slice, cram_codec *c) {
    return cram_get_block_by_id(slice, c->external.content_id);
}

cram_codec *cram_external_decode_init(cram_block_compression_hdr *hdr,
				      char *data, int size,
				      enum cram_external_type option,
				      int version, varint_vec *vv) {
    cram_codec *c;
    char *cp = data;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_EXTERNAL;
    if (option == E_INT)
	c->decode = cram_external_decode_int;
    else if (option == E_LONG)
	c->decode = cram_external_decode_long;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
	c->decode = cram_external_decode_char;
    else
	c->decode = cram_external_decode_block;
    c->free   = cram_external_decode_free;
    c->size   = cram_external_decode_size;
    c->get_block = cram_external_get_block;
    
    c->external.content_id = vv->varint_get32(&cp, NULL, NULL);

    if (cp - data != size) {
	fprintf(stderr, "Malformed external header stream\n");
	free(c);
	return NULL;
    }

    c->external.type = option;

    return c;
}

int cram_external_encode_int(cram_slice *slice, cram_codec *c,
			     char *in, int in_size) {
    uint32_t *i32 = (uint32_t *)in;

    c->vv->varint_put32_blk(c->out, *i32);
    return 0;
}

int cram_external_encode_long(cram_slice *slice, cram_codec *c,
			     char *in, int in_size) {
    uint64_t *i64 = (uint64_t *)in;

    c->vv->varint_put64_blk(c->out, *i64);
    return 0;
}

int cram_external_encode_char(cram_slice *slice, cram_codec *c,
			      char *in, int in_size) {
    BLOCK_APPEND(c->out, in, in_size);
    return 0;
}

void cram_external_encode_free(cram_codec *c) {
    if (!c)
	return;
    free(c);
}

int cram_external_encode_store(cram_codec *c, cram_block *b, char *prefix,
			       int version) {
    char tmp[99], *tp = tmp;
    int len = 0;

    if (prefix) {
	size_t l = strlen(prefix);
	BLOCK_APPEND(b, prefix, l);
	len += l;
    }

    tp += c->vv->varint_put32(tp, NULL, c->e_external.content_id);
    len += c->vv->varint_put32_blk(b, c->codec);
    len += c->vv->varint_put32_blk(b, tp-tmp);
    BLOCK_APPEND(b, tmp, tp-tmp);
    len += tp-tmp;

    return len;
}

cram_codec *cram_external_encode_init(cram_stats *st,
				      enum cram_external_type option,
				      void *dat,
				      int version, varint_vec *vv) {
    cram_codec *c;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec = E_EXTERNAL;
    c->free = cram_external_encode_free;
    if (option == E_INT)
	c->encode = cram_external_encode_int;
    else if (option == E_LONG)
	c->encode = cram_external_encode_long;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
	c->encode = cram_external_encode_char;
    else
	abort();
    c->store = cram_external_encode_store;
    c->flush = NULL;

    c->e_external.content_id = (size_t)dat;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * BETA
 */
void cram_nop_decode_reset(cram_codec *c) {}

int cram_beta_decode_long(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int64_t *out_i = (int64_t *)out;
    int i, n = *out_size;

    if (c->beta.nbits) {
        if (cram_not_enough_bits(in, c->beta.nbits * n))
	    return -1;

	for (i = 0; i < n; i++)
	    out_i[i] = get_bits_MSB(in, c->beta.nbits) - c->beta.offset;
    } else {
	for (i = 0; i < n; i++)
	    out_i[i] = -c->beta.offset;
    }

    return 0;
}

int cram_beta_decode_int(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n = *out_size;

    if (c->beta.nbits) {
        if (cram_not_enough_bits(in, c->beta.nbits * n))
	    return -1;

	for (i = 0; i < n; i++)
	    out_i[i] = get_bits_MSB(in, c->beta.nbits) - c->beta.offset;
    } else {
	for (i = 0; i < n; i++)
	    out_i[i] = -c->beta.offset;
    }

    return 0;
}

int cram_beta_decode_char(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int i, n = *out_size;

    if (c->beta.nbits) {
        if (cram_not_enough_bits(in, c->beta.nbits * n))
            return -1;

	if (out)
	    for (i = 0; i < n; i++)
		out[i] = get_bits_MSB(in, c->beta.nbits) - c->beta.offset;
	else
	    for (i = 0; i < n; i++)
		get_bits_MSB(in, c->beta.nbits);
    } else {
	if (out)
	    for (i = 0; i < n; i++)
		out[i] = -c->beta.offset;
    }

    return 0;
}

void cram_beta_decode_free(cram_codec *c) {
    if (c)
	free(c);
}

cram_codec *cram_beta_decode_init(cram_block_compression_hdr *hdr,
				  char *data, int size,
				  enum cram_external_type option,
				  int version, varint_vec *vv) {
    cram_codec *c;
    char *cp = data;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_BETA;
    if (option == E_LONG)
	c->decode = cram_beta_decode_long;
    else if (option == E_INT)
	c->decode = cram_beta_decode_int;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
	c->decode = cram_beta_decode_char;
    else {
	fprintf(stderr, "BYTE_ARRAYs not supported by this codec\n");
	return NULL;
    }
    c->free   = cram_beta_decode_free;
    
    c->beta.nbits = -1;
    c->beta.offset = vv->varint_get32(&cp, NULL, NULL);
    c->beta.nbits  = vv->varint_get32(&cp, NULL, NULL);

    if (cp - data != size
        || c->beta.nbits < 0 || c->beta.nbits > 8 * sizeof(int64_t)) {
	fprintf(stderr, "Malformed beta header stream\n");
	free(c);
	return NULL;
    }

    return c;
}

int cram_beta_encode_store(cram_codec *c, cram_block *b,
			   char *prefix, int version) {
    int len = 0;

    if (prefix) {
	size_t l = strlen(prefix);
	BLOCK_APPEND(b, prefix, l);
	len += l;
    }

    len += c->vv->varint_put32_blk(b, c->codec);
    // codec length
    len += c->vv->varint_put32_blk(b, c->vv->varint_size(c->e_beta.offset)
				   +  c->vv->varint_size(c->e_beta.nbits));
    len += c->vv->varint_put32_blk(b, c->e_beta.offset);
    len += c->vv->varint_put32_blk(b, c->e_beta.nbits);

    return len;
}

int cram_beta_encode_long(cram_slice *slice, cram_codec *c,
			  char *in, int in_size) {
    int64_t *syms = (int64_t *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
	r |= store_bits_MSB(c->out, syms[i] + c->e_beta.offset,
			    c->e_beta.nbits);

    return r;
}

int cram_beta_encode_int(cram_slice *slice, cram_codec *c,
			 char *in, int in_size) {
    int *syms = (int *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
	r |= store_bits_MSB(c->out, syms[i] + c->e_beta.offset,
			    c->e_beta.nbits);

    return r;
}

int cram_beta_encode_char(cram_slice *slice, cram_codec *c,
			 char *in, int in_size) {
    unsigned char *syms = (unsigned char *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
	r |= store_bits_MSB(c->out, syms[i] + c->e_beta.offset,
			    c->e_beta.nbits);

    return r;
}

void cram_beta_encode_free(cram_codec *c) {
    if (c) free(c);
}

cram_codec *cram_beta_encode_init(cram_stats *st,
				  enum cram_external_type option,
				  void *dat,
				  int version, varint_vec *vv) {
    cram_codec *c;
    int64_t min_val, max_val, len = 0;
    int64_t range;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec  = E_BETA;
    c->free   = cram_beta_encode_free;
    if (option == E_LONG)
	c->encode = cram_beta_encode_long;
    else if (option == E_INT)
	c->encode = cram_beta_encode_int;
    else
	c->encode = cram_beta_encode_char;
    c->store  = cram_beta_encode_store;
    c->flush = NULL;

    if (dat) {
	min_val = ((int64_t *)dat)[0];
	max_val = ((int64_t *)dat)[1];
    } else {
	min_val = INT64_MAX;
	max_val = INT64_MIN;
	int i;
	for (i = 0; i < MAX_STAT_VAL; i++) {
	    if (!st->freqs[i])
		continue;
	    if (min_val > i)
		min_val = i;
	    max_val = i;
	}
	if (st->h) {
	    HashItem *hi;
	    HashIter *iter = HashTableIterCreate();
	    while ((hi = HashTableIterNext(st->h, iter))) {
		i = (int64_t)(size_t)hi->key;
		if (min_val > i)
		    min_val = i;
		if (max_val < i)
		    max_val = i;
	    }
	    HashTableIterDestroy(iter);
	}
    }

    assert(max_val >= min_val);
    c->e_beta.offset = -min_val;
    range = (int64_t) max_val - min_val;
    while (range) {
	len++;
	range >>= 1;
    }
    c->e_beta.nbits = len;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * XPACK: BETA as a transform instead of writing to CORE block.
 *
 * This also has the additional requirement that the data series is not
 * interleaved with another, permitting efficient encoding and decoding
 * of all elements enmasse instead of needing to only extract the bits
 * necessary per item.
 */
int cram_xpack_decode_long(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int64_t *out_i = (int64_t *)out;
    int i, n = *out_size;

    if (c->xpack.nbits) {
	for (i = 0; i < n; i++)
	    out_i[i] = c->xpack.rmap[get_bits_MSB(in, c->xpack.nbits)];
    } else {
	for (i = 0; i < n; i++)
	    out_i[i] = c->xpack.rmap[0];
    }

    return 0;
}

int cram_xpack_decode_int(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n = *out_size;

    if (c->xpack.nbits) {
        if (cram_not_enough_bits(in, c->xpack.nbits * n))
	    return -1;

	for (i = 0; i < n; i++)
	    out_i[i] = c->xpack.rmap[get_bits_MSB(in, c->xpack.nbits)];
    } else {
	for (i = 0; i < n; i++)
	    out_i[i] = c->xpack.rmap[0];
    }

    return 0;
}

static int cram_xpack_decode_expand_char(cram_slice *slice, cram_codec *c) {
    cram_block *b = slice->block_by_id[512 + c->codec_id];
    if (b)
	return 0;

    // get sub-codec data.
    cram_block *sub_b = c->xpack.sub_codec->get_block(slice, c->xpack.sub_codec);

    // Allocate local block to expand into
    b = slice->block_by_id[512 + c->codec_id] = cram_new_block(0, 0);
    int sub_size = c->xpack.sub_codec->size(slice, c->xpack.sub_codec);
    int n = sub_b->uncomp_size * 8/c->xpack.nbits;
    BLOCK_GROW(b, n);
    b->uncomp_size = n;

    // Expand; FIXME: use hts_unpack
    unsigned char *dat = sub_b->data;
    unsigned char *out = b->data;
    int i, j = 0;

    switch (c->xpack.nbits) {
    case 8:
#pragma GCC unroll 8
	for (i = 0; i < n; i++)
	    out[i] = c->xpack.rmap[dat[j++]];
	break;

    case 4: {
	union {
	    uint16_t w;
	    uint8_t c[2];
	} map[256];
	int Y, Z, P=0;
	for (Y = 0; Y < 16; Y++) {
	    for (Z = 0; Z < 16; Z++, P++) {
		map[P].c[0] = c->xpack.rmap[Y];
		map[P].c[1] = c->xpack.rmap[Z];
	    }
	}

	uint16_t *out16 = (uint16_t *)out;
	for (j = 0; j < (n & ~1)/2; j++)
	    out16[j] = map[dat[j]].w;
	i = j*2;
	sub_b->byte += i;

	for (;i < n; i++)
	    out[i] = c->xpack.rmap[get_bits_MSB(sub_b, c->xpack.nbits)];
	break;
    }

    case 2: {
	union {
	    uint32_t w;
	    uint8_t c[4];
	} map[256];
	int W, X, Y, Z, P=0;
	for (W = 0; W < 4; W++)
	    for (X = 0; X < 4; X++)
		for (Y = 0; Y < 4; Y++)
		    for (Z = 0; Z < 4; Z++, P++) {
			map[P].c[0] = c->xpack.rmap[W];
			map[P].c[1] = c->xpack.rmap[X];
			map[P].c[2] = c->xpack.rmap[Y];
			map[P].c[3] = c->xpack.rmap[Z];
		    }

	uint32_t *out32 = (uint32_t *)out;
	for (j = 0; j < (n & ~3)/4; j++)
	    out32[j] = map[dat[j]].w;
	i = j*4;
	sub_b->byte += i;

	for (;i < n; i++)
	    out[i] = c->xpack.rmap[get_bits_MSB(sub_b, c->xpack.nbits)];
	break;
    }

    case 1:
	for (i = 0; i < (n & ~7); i += 8) {
	    out[i+0] = c->xpack.rmap[((sub_b->data[b->byte] >> 7)&1)];
	    out[i+1] = c->xpack.rmap[((sub_b->data[b->byte] >> 6)&1)];
	    out[i+2] = c->xpack.rmap[((sub_b->data[b->byte] >> 5)&1)];
	    out[i+3] = c->xpack.rmap[((sub_b->data[b->byte] >> 4)&1)];
	    out[i+4] = c->xpack.rmap[((sub_b->data[b->byte] >> 3)&1)];
	    out[i+5] = c->xpack.rmap[((sub_b->data[b->byte] >> 2)&1)];
	    out[i+6] = c->xpack.rmap[((sub_b->data[b->byte] >> 1)&1)];
	    out[i+7] = c->xpack.rmap[( sub_b->data[b->byte]      &1)];
	    sub_b->byte++;
	}
	for (;i < n; i++)
	    out[i] = c->xpack.rmap[get_bits_MSB(sub_b, c->xpack.nbits)];
	break;

    default:
	abort();
    }

    return 0;
}

int cram_xpack_decode_char(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int i, n = *out_size;

    // FIXME: we need to ban data-series interleaving in the spec for this to work.

    // Remember this may be called when threaded and multi-slice per container.
    // Hence one cram_codec instance, multiple slices, multiple blocks.
    // We therefore have to cache appropriate block info in slice and not codec.
    //    b = cram_get_block_by_id(slice, c->external.content_id);
    cram_xpack_decode_expand_char(slice, c);
    cram_block *b = slice->block_by_id[512 + c->codec_id];

    if (out)
	memcpy(out, b->data + b->byte, *out_size);
    b->byte += *out_size;

    return 0;
}

void cram_xpack_decode_free(cram_codec *c) {
    if (!c) return;

    if (c->xpack.sub_codec)
	c->xpack.sub_codec->free(c->xpack.sub_codec);

    //free(slice->block_by_id[512 + c->codec_id]);
    //slice->block_by_id[512 + c->codec_id] = 0;

    free(c);
}

int cram_xpack_decode_size(cram_slice *slice, cram_codec *c) {
    cram_xpack_decode_expand_char(slice, c);
    return slice->block_by_id[512 + c->codec_id]->uncomp_size;
}

cram_block *cram_xpack_get_block(cram_slice *slice, cram_codec *c) {
    cram_xpack_decode_expand_char(slice, c);
    return slice->block_by_id[512 + c->codec_id];
}

cram_codec *cram_xpack_decode_init(cram_block_compression_hdr *hdr,
				   char *data, int size,
				   enum cram_external_type option,
				   int version, varint_vec *vv) {
    cram_codec *c;
    char *cp = data;
    char *endp = data+size;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_XPACK;
    if (option == E_LONG)
	c->decode = cram_xpack_decode_long;
    else if (option == E_INT)
	c->decode = cram_xpack_decode_int;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
	c->decode = cram_xpack_decode_char;
    else {
	fprintf(stderr, "BYTE_ARRAYs not supported by this codec\n");
	return NULL;
    }
    c->free = cram_xpack_decode_free;
    c->size = cram_xpack_decode_size;
    c->get_block = cram_xpack_get_block;
    
    c->xpack.nbits = -1;
    c->xpack.nbits = vv->varint_get32(&cp, endp, NULL);
    c->xpack.nval  = vv->varint_get32(&cp, endp, NULL);
    int i;
    for (i = 0; i < c->xpack.nval; i++) {
	uint32_t v = vv->varint_get32(&cp, endp, NULL);
	if (v >= 256) return NULL;
	c->xpack.rmap[i] = v; // reverse map: e.g 0-3 to P,A,C,K
    }

    int encoding = vv->varint_get32(&cp, endp, NULL);
    int sub_size = vv->varint_get32(&cp, endp, NULL);
    if (sub_size < 0 || endp - cp < sub_size)
        goto malformed;
    c->xpack.sub_codec = cram_decoder_init(hdr, encoding, cp, sub_size,
					   option, version, vv);
    if (c->xpack.sub_codec == NULL)
        goto malformed;
    cp += sub_size;

    if (cp - data != size
        || c->xpack.nbits < 0 || c->xpack.nbits > 8 * sizeof(int64_t)) {
    malformed:
	fprintf(stderr, "Malformed xpack header stream\n");
	free(c);
	return NULL;
    }

    return c;
}

int cram_xpack_encode_flush(cram_codec *c) {
    store_bits_MSB(c->out, 0, 7); // ensure whole byte flushed

    // We've buffered up data to c->e_xpack->out.
    // We now need to pass this through the next layer of transform first.
    if (c->e_xpack.sub_codec->encode(/*slice*/NULL,
				     c->e_xpack.sub_codec,
				     (char *)BLOCK_DATA(c->out),
				     BLOCK_SIZE(c->out)))
	return -1;

    if (c->e_xpack.sub_codec->flush)
	return c->e_xpack.sub_codec->flush(c->e_xpack.sub_codec);
    
    return 0;
}

int cram_xpack_encode_store(cram_codec *c, cram_block *b,
			    char *prefix, int version) {
    int len = 0;

    if (prefix) {
	size_t l = strlen(prefix);
	BLOCK_APPEND(b, prefix, l);
	len += l;
    }

    // Store sub-codec
    cram_codec *tc = c->e_xpack.sub_codec;
    cram_block *tb = cram_new_block(0, 0);
    int len2 = tc->store(tc, tb, NULL, version);

    len += c->vv->varint_put32_blk(b, c->codec);

    // codec length
    int len1 = 0, i, n;
    for (i = 0; i < c->e_xpack.nval; i++)
	len1 += c->vv->varint_size(c->e_xpack.rmap[i]);
    len += c->vv->varint_put32_blk(b, c->vv->varint_size(c->e_xpack.nbits)
				   +  c->vv->varint_size(c->e_xpack.nval)
				   + len1 + len2);

    // The map and sub-codec
    len += c->vv->varint_put32_blk(b, c->e_xpack.nbits);
    len += c->vv->varint_put32_blk(b, c->e_xpack.nval);
    for (n = i = 0; i < c->e_xpack.nval; i++)
	len += c->vv->varint_put32_blk(b, c->e_xpack.rmap[i]);

    BLOCK_APPEND(b, BLOCK_DATA(tb), BLOCK_SIZE(tb));

    cram_free_block(tb);

    return len + len2;
}

// Same as cram_beta_encode_long
int cram_xpack_encode_long(cram_slice *slice, cram_codec *c,
			   char *in, int in_size) {
    int64_t *syms = (int64_t *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
	r |= store_bits_MSB(c->out, c->e_xpack.map[syms[i]], c->e_xpack.nbits);

    return r;
}

int cram_xpack_encode_int(cram_slice *slice, cram_codec *c,
			  char *in, int in_size) {
    int *syms = (int *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
	r |= store_bits_MSB(c->out, c->e_xpack.map[syms[i]], c->e_xpack.nbits);

    return r;
}

int cram_xpack_encode_char(cram_slice *slice, cram_codec *c,
			   char *in, int in_size) {
    unsigned char *syms = (unsigned char *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
	r |= store_bits_MSB(c->out, c->e_xpack.map[syms[i]], c->e_xpack.nbits);

    return r;
}

void cram_xpack_encode_free(cram_codec *c) {
    if (c) free(c);
}

cram_codec *cram_xpack_encode_init(cram_stats *st,
				   enum cram_external_type option,
				   void *dat,
				   int version, varint_vec *vv) {
    cram_codec *c;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec  = E_XPACK;
    c->free   = cram_xpack_encode_free;
    if (option == E_LONG)
	c->encode = cram_xpack_encode_long;
    else if (option == E_INT)
	c->encode = cram_xpack_encode_int;
    else
	c->encode = cram_xpack_encode_char;
    c->store  = cram_xpack_encode_store;
    c->flush  = cram_xpack_encode_flush;

    cram_xpack_encoder *e = (cram_xpack_encoder *)dat;
    c->e_xpack.nbits = e->nbits;
    c->e_xpack.nval = e->nval;
    c->e_xpack.sub_codec = cram_encoder_init(e->sub_encoding, NULL,
					     E_BYTE_ARRAY, e->sub_codec_dat,
					     version, vv);

    // Initialise fwd and rev maps
    memcpy(c->e_xpack.map, e->map, sizeof(e->map)); // P,A,C,K to 0,1,2,3
    int i, n;
    for (i = n = 0; i < 256; i++)
	if (e->map[i] != -1)
	    c->e_xpack.rmap[n++] = i;               // 0,1,2,3 to P,A,C,K
    if (n != e->nval) {
	fprintf(stderr, "Incorrectly specified number of map items in PACK\n");
	return NULL;
    }

    // Temporary buffer
    if (!(c->out = cram_new_block(0, 0)))
	return NULL;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * XRLE
 *
 * This also has the additional requirement that the data series is not
 * interleaved with another, permitting efficient encoding and decoding
 * of all elements enmasse instead of needing to only extract the bits
 * necessary per item.
 */
int cram_xrle_decode_long(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    // TODO
    return 0;
}

int cram_xrle_decode_int(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    // TODO
    return 0;
}

// Expands an XRLE transform and caches result in slice->block_by_id[]
static int cram_xrle_decode_expand_char(cram_slice *slice, cram_codec *c) {
    cram_block *b = slice->block_by_id[512 + c->codec_id];
    if (b)
	return 0;

    b = slice->block_by_id[512 + c->codec_id] = cram_new_block(0, 0);
    cram_block *lit_b = c->xrle.lit_codec->get_block(slice, c->xrle.lit_codec);
    unsigned char *lit_dat = lit_b->data;
    unsigned int lit_sz = lit_b->uncomp_size;
    unsigned int len_sz = c->xrle.len_codec->size(slice, c->xrle.len_codec);

    cram_block *len_b = c->xrle.len_codec->get_block(slice, c->xrle.len_codec);
    char *len_start = (char *)BLOCK_DATA(len_b);
    char *len_end   = len_start + len_b->uncomp_size;

    int guess; // guess expanded size
    BLOCK_GROW(b, guess = len_sz*5 + lit_sz);
    unsigned char *dat = b->data;
    unsigned char *dat_end = dat + guess;

    // Local faster copy of xrle.rep_score[]
    int i;
    uint8_t do_rle[256];
    for (i = 0; i < 256; i++)
	do_rle[i] = c->xrle.rep_score[i] > 0;

    int lit_idx;
    int cur_len = 0;
    unsigned char cur_lit = 0;
    
    unsigned char *lit_start = lit_dat;
    unsigned char *lit_end = lit_dat + lit_sz;
    while (lit_dat < lit_end) {
	// Get symbol and if needed length
	cur_lit = *lit_dat++;
	if (do_rle[cur_lit]) {
	    // Direct approach; marginal speed gain
	    int err = 0;
	    cur_len = c->vv->varint_get32(&len_start, len_end, &err) + 1;
	    if (err)
		return -1;

//	    int one;
//	    if (c->xrle.len_codec->decode(slice, c->xrle.len_codec, NULL,
//					  (char *)&cur_len, &one) < 0)
//		return -1;
//	    cur_len++;
	} else {
	    cur_len = 1;
	}

	// Expand memory guess if required
	if (dat + cur_len > dat_end) {
	    guess *= (lit_sz + 1.0) / (lit_dat - lit_start) * 1.05;
	    guess += cur_len;
	    ptrdiff_t off = dat - b->data;
	    BLOCK_RESIZE(b, guess);
	    dat = b->data + off;
	    dat_end = dat + guess;
	}

	switch(cur_len) {
	case 1: *dat++ = cur_lit;
	    break;

	default:
	    // Repeat symbol
	    memset(dat, cur_lit, cur_len);
	    dat += cur_len;
	}
    }
    b->uncomp_size = dat - b->data;

    return 0;
}

int cram_xrle_decode_size(cram_slice *slice, cram_codec *c) {
    cram_xrle_decode_expand_char(slice, c);
    return slice->block_by_id[512 + c->codec_id]->uncomp_size;
}

cram_block *cram_xrle_get_block(cram_slice *slice, cram_codec *c) {
    cram_xrle_decode_expand_char(slice, c);
    return slice->block_by_id[512 + c->codec_id];
}

int cram_xrle_decode_char(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int n = *out_size;

    cram_xrle_decode_expand_char(slice, c);
    cram_block *b = slice->block_by_id[512 + c->codec_id];

    memcpy(out, b->data + b->idx, n);
    b->idx += n;
    return 0;

    // Old code when not cached
    while (n > 0) {
	if (c->xrle.cur_len == 0) {
	    unsigned char lit;
	    int one = 1;
	    if (c->xrle.lit_codec->decode(slice, c->xrle.lit_codec, in,
					  (char *)&lit, &one) < 0)
		return -1;
	    c->xrle.cur_lit = lit;

	    if (c->xrle.rep_score[lit] > 0) {
		if (c->xrle.len_codec->decode(slice, c->xrle.len_codec, in,
					      (char *)&c->xrle.cur_len, &one) < 0)
		    return -1;
	    } // else cur_len still zero
	    //else fprintf(stderr, "%d\n", lit);

	    c->xrle.cur_len++;
	}

	if (n >= c->xrle.cur_len) {
	    memset(out, c->xrle.cur_lit, c->xrle.cur_len);
	    out += c->xrle.cur_len;
	    n -= c->xrle.cur_len;
	    c->xrle.cur_len = 0;
	} else {
	    memset(out, c->xrle.cur_lit, n);
	    out += n;
	    c->xrle.cur_len -= n;
	    n = 0;
	}
    }

    return 0;
}

void cram_xrle_decode_free(cram_codec *c) {
    if (!c) return;

    if (c->xrle.len_codec)
	c->xrle.len_codec->free(c->xrle.len_codec);

    if (c->xrle.lit_codec)
	c->xrle.lit_codec->free(c->xrle.lit_codec);

    free(c);
}

cram_codec *cram_xrle_decode_init(cram_block_compression_hdr *hdr,
				   char *data, int size,
				   enum cram_external_type option,
				   int version, varint_vec *vv) {
    cram_codec *c;
    char *cp = data;
    char *endp = data+size;
    int err = 0;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_XRLE;
    if (option == E_LONG)
	c->decode = cram_xrle_decode_long;
    else if (option == E_INT)
	c->decode = cram_xrle_decode_int;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
	c->decode = cram_xrle_decode_char;
    else {
	fprintf(stderr, "BYTE_ARRAYs not supported by this codec\n");
	return NULL;
    }
    c->free   = cram_xrle_decode_free;
    c->size   = cram_xrle_decode_size;
    c->get_block = cram_xrle_get_block;
    c->xrle.cur_len = 0;
    c->xrle.cur_lit = -1;

    // RLE map
    int i, j, nrle = vv->varint_get32(&cp, endp, &err);
    memset(c->xrle.rep_score, 0, 256*sizeof(*c->xrle.rep_score));
    for (i = 0; i < nrle && i < 256; i++) {
	j = vv->varint_get32(&cp, endp, &err);
	if (j >= 0 && j < 256)
	    c->xrle.rep_score[j] = 1;
    }
    
    // Length and literal sub encodings
    c->xrle.len_encoding = vv->varint_get32(&cp, endp, &err);
    int sub_size = vv->varint_get32(&cp, endp, &err);
    if (sub_size < 0 || endp - cp < sub_size)
        goto malformed;
    c->xrle.len_codec = cram_decoder_init(hdr, c->xrle.len_encoding,
					  cp, sub_size, E_INT, version, vv);
    if (c->xrle.len_codec == NULL)
        goto malformed;
    cp += sub_size;

    c->xrle.lit_encoding = vv->varint_get32(&cp, endp, &err);
    sub_size = vv->varint_get32(&cp, endp, &err);
    if (sub_size < 0 || endp - cp < sub_size)
        goto malformed;
    c->xrle.lit_codec = cram_decoder_init(hdr, c->xrle.lit_encoding,
					  cp, sub_size, option, version, vv);
    if (c->xrle.lit_codec == NULL)
        goto malformed;
    cp += sub_size;

    if (err)
	goto malformed;

    return c;

 malformed:
    fprintf(stderr, "Malformed xrle header stream\n");
    free(c);
    return NULL;
}

int cram_xrle_encode_flush(cram_codec *c) {
    if (c->e_xrle.rep_score[c->e_xrle.cur_lit] > 0) {
	// delayed last symbol
	if (c->e_xrle.len_codec->encode(/*slice*/NULL,
					c->e_xrle.len_codec,
					(char *)&c->e_xrle.cur_len, 1))
	    return -1;
//	unsigned char lit = c->e_xrle.cur_lit;
//	if (c->e_xrle.lit_codec->encode(/*slice*/NULL,
//					c->e_xrle.lit_codec,
//					(char *)&lit, 1))
//	    return -1;
//	fprintf(stderr, "Flush %d x %d\n", lit, c->e_xrle.cur_len);
    }

    if (c->e_xrle.len_codec->flush)
	return c->e_xrle.len_codec->flush(c->e_xrle.len_codec);
    if (c->e_xrle.lit_codec->flush)
	return c->e_xrle.lit_codec->flush(c->e_xrle.lit_codec);
    
    return 0;
}

int cram_xrle_encode_store(cram_codec *c, cram_block *b,
			    char *prefix, int version) {
    int len = 0;
    cram_codec *tc;
    cram_block *b_rle, *b_len, *b_lit;

    if (prefix) {
	size_t l = strlen(prefix);
	BLOCK_APPEND(b, prefix, l);
	len += l;
    }

    // List of symbols to RLE
    b_rle = cram_new_block(0, 0);
    int i, nrle = 0, len1 = 0;
    for (i = 0; i < 256; i++) {
	if (c->e_xrle.rep_score[i] > 0) {
	    nrle++;
	    len1 += c->vv->varint_put32_blk(b_rle,i);
	}
    }

    // Store length and literal sub-codecs to get encoded length
    tc = c->e_xrle.len_codec;
    b_len = cram_new_block(0, 0);
    int len2 = tc->store(tc, b_len, NULL, version);

    tc = c->e_xrle.lit_codec;
    b_lit = cram_new_block(0, 0);
    int len3 = tc->store(tc, b_lit, NULL, version);

    len += c->vv->varint_put32_blk(b, c->codec);
    len += c->vv->varint_put32_blk(b, len1 + len2 + len3
				   + c->vv->varint_size(nrle));
    len += c->vv->varint_put32_blk(b, nrle);
    // FIXME: plus symbols we want to do RLE on.  For now it's all
    BLOCK_APPEND(b, BLOCK_DATA(b_rle), BLOCK_SIZE(b_rle));
    BLOCK_APPEND(b, BLOCK_DATA(b_len), BLOCK_SIZE(b_len));
    BLOCK_APPEND(b, BLOCK_DATA(b_lit), BLOCK_SIZE(b_lit));

    cram_free_block(b_rle);
    cram_free_block(b_len);
    cram_free_block(b_lit);

    return len + len1 + len2 + len3;
}

int cram_xrle_encode_long(cram_slice *slice, cram_codec *c,
			   char *in, int in_size) {
    // FIXME
    return 0;
}

int cram_xrle_encode_int(cram_slice *slice, cram_codec *c,
			  char *in, int in_size) {
    // FIXME
    return 0;
}

int cram_xrle_encode_char(cram_slice *slice, cram_codec *c,
			  char *in, int in_size) {
    unsigned char *syms = (unsigned char *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++) {
	unsigned char lit = syms[i];

	//fprintf(stderr, "lit[%d]=%d\n", i, lit);

	// cur_lit is basically last symbol
	if (lit != c->xrle.cur_lit) {
	    // Diffrent sym
	    if (c->xrle.cur_lit >= 0 && c->e_xrle.rep_score[c->xrle.cur_lit] > 0) {
		c->e_xrle.len_codec->encode(slice, c->e_xrle.len_codec,
					    (char *)&c->xrle.cur_len, 1);
		//fprintf(stderr, "last %d x %d\n", c->xrle.cur_lit, c->xrle.cur_len);
	    }
	    //fprintf(stderr, "(sym) %d\n", lit);
	    c->e_xrle.lit_codec->encode(slice, c->e_xrle.lit_codec,
					(char *)&lit, 1);
	    c->e_xrle.cur_lit = syms[i];
	    c->xrle.cur_len = 0;

	} else if (c->xrle.cur_lit >= 0 && c->e_xrle.rep_score[c->xrle.cur_lit] <= 0) {
	    // Same sym, but cost of run len is too high
	    c->e_xrle.lit_codec->encode(slice, c->e_xrle.lit_codec,
					(char *)&lit, 1);
	    //fprintf(stderr, "sym %d\n", lit);
	    c->e_xrle.cur_lit = syms[i];
	    c->xrle.cur_len = 0;
	} else {
	    // Same sym and we're accumulating a run length
	    c->xrle.cur_len++;
	}
    }

    return r;
}

void cram_xrle_encode_free(cram_codec *c) {
    if (c) free(c);
}

cram_codec *cram_xrle_encode_init(cram_stats *st,
				   enum cram_external_type option,
				   void *dat,
				   int version, varint_vec *vv) {
    cram_codec *c;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec  = E_XRLE;
    c->free   = cram_xrle_encode_free;
    if (option == E_LONG)
	c->encode = cram_xrle_encode_long;
    else if (option == E_INT)
	c->encode = cram_xrle_encode_int;
    else
	c->encode = cram_xrle_encode_char;
    c->store  = cram_xrle_encode_store;
    c->flush  = cram_xrle_encode_flush;

    cram_xrle_encoder *e = (cram_xrle_encoder *)dat;

    c->e_xrle.len_codec = cram_encoder_init(e->len_encoding, NULL,
					    E_INT, e->len_dat,
					    version, vv);
    c->e_xrle.lit_codec = cram_encoder_init(e->lit_encoding, NULL,
					    E_BYTE, e->lit_dat,
					    version, vv);
    c->e_xrle.cur_lit = -1;
    c->e_xrle.cur_len = -1;

    memcpy(c->e_xrle.rep_score, e->rep_score, 256*sizeof(*c->e_xrle.rep_score));

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * SUBEXP
 */
int cram_subexp_decode(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int n, count;
    int k = c->subexp.k;

    for (count = 0, n = *out_size; count < n; count++) {
	int i = 0, tail;
	int val;

	/* Get number of 1s */
	//while (get_bit_MSB(in) == 1) i++;
	i = get_one_bits_MSB(in);
        if (i < 0 || cram_not_enough_bits(in, i > 0 ? i + k - 1 : k))
            return -1;
	/*
	 * Val is
	 * i > 0:  2^(k+i-1) + k+i-1 bits
	 * i = 0:  k bits
	 */
	if (i) {
	    tail = i + k-1;
	    val = 0;
	    while (tail) {
		//val = val<<1; val |= get_bit_MSB(in);
		GET_BIT_MSB(in, val);
		tail--;
	    }
	    val += 1 << (i + k-1);
	} else {
	    tail = k;
	    val = 0;
	    while (tail) {
		//val = val<<1; val |= get_bit_MSB(in);
		GET_BIT_MSB(in, val);
		tail--;
	    }
	}

	out_i[count] = val - c->subexp.offset;
    }

    return 0;
}

void cram_subexp_decode_free(cram_codec *c) {
    if (c)
	free(c);
}

cram_codec *cram_subexp_decode_init(cram_block_compression_hdr *hdr,
				    char *data, int size,
				    enum cram_external_type option,
				    int version, varint_vec *vv) {
    cram_codec *c;
    char *cp = data;

    if (option == E_BYTE_ARRAY_BLOCK) {
	fprintf(stderr, "BYTE_ARRAYs not supported by this codec\n");
	return NULL;
    }

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_SUBEXP;
    c->decode = cram_subexp_decode;
    c->free   = cram_subexp_decode_free;
    c->subexp.k = -1;

    c->subexp.offset = vv->varint_get32(&cp, data + size, NULL);
    c->subexp.k      = vv->varint_get32(&cp, data + size, NULL);

    if (cp - data != size || c->subexp.k < 0) {
	fprintf(stderr, "Malformed subexp header stream\n");
	free(c);
	return NULL;
    }

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * GAMMA
 */
int cram_gamma_decode(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n;

    for (i = 0, n = *out_size; i < n; i++) {
	int nz = 0;
	int val;
	//while (get_bit_MSB(in) == 0) nz++;
	nz = get_zero_bits_MSB(in);
        if (cram_not_enough_bits(in, nz))
            return -1;
	val = 1;
	while (nz > 0) {
	    //val <<= 1; val |= get_bit_MSB(in);
	    GET_BIT_MSB(in, val);
	    nz--;
	}

	out_i[i] = val - c->gamma.offset;
    }

    return 0;
}

void cram_gamma_decode_free(cram_codec *c) {
    if (c)
	free(c);
}

cram_codec *cram_gamma_decode_init(cram_block_compression_hdr *hdr,
				   char *data, int size,
				   enum cram_external_type option,
				   int version, varint_vec *vv) {
    cram_codec *c;
    char *cp = data;

    if (option == E_BYTE_ARRAY_BLOCK) {
	fprintf(stderr, "BYTE_ARRAYs not supported by this codec\n");
	return NULL;
    }

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_GAMMA;
    c->decode = cram_gamma_decode;
    c->free   = cram_gamma_decode_free;
    
    c->gamma.offset = vv->varint_get32(&cp, NULL, NULL);

    if (cp - data != size) {
	fprintf(stderr, "Malformed gamma header stream\n");
	free(c);
	return NULL;
    }

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * HUFFMAN
 */

static int code_sort(const void *vp1, const void *vp2) {
    const cram_huffman_code *c1 = (const cram_huffman_code *)vp1;
    const cram_huffman_code *c2 = (const cram_huffman_code *)vp2;

    if (c1->len != c2->len)
	return c1->len - c2->len;
    else
	return c1->symbol - c2->symbol;
}

void cram_huffman_decode_free(cram_codec *c) {
    if (!c)
	return;

    if (c->huffman.codes)
	free(c->huffman.codes);
    free(c);
}

int cram_huffman_decode_null(cram_slice *slice, cram_codec *c,
			     cram_block *in, char *out, int *out_size) {
    return -1;
}

int cram_huffman_decode_char0(cram_slice *slice, cram_codec *c,
			      cram_block *in, char *out, int *out_size) {
    int i, n;

    if (!out)
	return 0;

    /* Special case of 0 length codes */
    for (i = 0, n = *out_size; i < n; i++) {
	out[i] = c->huffman.codes[0].symbol;
    }
    return 0;
}

int cram_huffman_decode_char(cram_slice *slice, cram_codec *c,
			     cram_block *in, char *out, int *out_size) {
    int i, n, ncodes = c->huffman.ncodes;
    const cram_huffman_code * const codes = c->huffman.codes;

    for (i = 0, n = *out_size; i < n; i++) {
	int idx = 0;
	int val = 0, len = 0, last_len = 0;

	for (;;) {
	    int dlen = codes[idx].len - last_len;
	    if (cram_not_enough_bits(in, dlen))
		return -1;

	    //val <<= dlen;
	    //val  |= get_bits_MSB(in, dlen);
	    //last_len = (len += dlen);

	    last_len = (len += dlen);
	    for (; dlen; dlen--) GET_BIT_MSB(in, val);

	    idx = val - codes[idx].p;
	    if (idx >= ncodes || idx < 0)
		return -1;

	    if (codes[idx].code == val && codes[idx].len == len) {
		if (out)
		    out[i] = codes[idx].symbol;
		break;
	    }
	}
    }

    return 0;
}

int cram_huffman_decode_int0(cram_slice *slice, cram_codec *c,
			     cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n;
    const cram_huffman_code * const codes = c->huffman.codes;

    /* Special case of 0 length codes */
    for (i = 0, n = *out_size; i < n; i++) {
	out_i[i] = codes[0].symbol;
    }
    return 0;
}

int cram_huffman_decode_int(cram_slice *slice, cram_codec *c,
			    cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n, ncodes = c->huffman.ncodes;
    const cram_huffman_code * const codes = c->huffman.codes;

    for (i = 0, n = *out_size; i < n; i++) {
	int idx = 0;
	int val = 0, len = 0, last_len = 0;

	// Now one bit at a time for remaining checks
	for (;;) {
	    int dlen = codes[idx].len - last_len;
	    if (cram_not_enough_bits(in, dlen))
		return -1;
	    
	    //val <<= dlen;
	    //val  |= get_bits_MSB(in, dlen);
	    //last_len = (len += dlen);

	    last_len = (len += dlen);
	    for (; dlen; dlen--) GET_BIT_MSB(in, val);

	    idx = val - codes[idx].p;
	    if (idx >= ncodes || idx < 0)
		return -1;

	    if (codes[idx].code == val && codes[idx].len == len) {
		out_i[i] = codes[idx].symbol;
		break;
	    }
	}
    }

    return 0;
}

int cram_huffman_decode_long0(cram_slice *slice, cram_codec *c,
			      cram_block *in, char *out, int *out_size) {
    int64_t *out_i = (int64_t *)out;
    int i, n;
    const cram_huffman_code * const codes = c->huffman.codes;

    /* Special case of 0 length codes */
    for (i = 0, n = *out_size; i < n; i++) {
	out_i[i] = codes[0].symbol;
    }
    return 0;
}

int cram_huffman_decode_long(cram_slice *slice, cram_codec *c,
			     cram_block *in, char *out, int *out_size) {
    int64_t *out_i = (int64_t *)out;
    int i, n, ncodes = c->huffman.ncodes;
    const cram_huffman_code * const codes = c->huffman.codes;

    for (i = 0, n = *out_size; i < n; i++) {
	int idx = 0;
	int val = 0, len = 0, last_len = 0;

	// Now one bit at a time for remaining checks
	for (;;) {
	    int dlen = codes[idx].len - last_len;
	    if (cram_not_enough_bits(in, dlen))
		return -1;

	    //val <<= dlen;
	    //val  |= get_bits_MSB(in, dlen);
	    //last_len = (len += dlen);

	    last_len = (len += dlen);
	    for (; dlen; dlen--) GET_BIT_MSB(in, val);

	    idx = val - codes[idx].p;
	    if (idx >= ncodes || idx < 0)
		return -1;

	    if (codes[idx].code == val && codes[idx].len == len) {
		out_i[i] = codes[idx].symbol;
		break;
	    }
	}
    }

    return 0;
}

/*
 * Initialises a huffman decoder from an encoding data stream.
 */
cram_codec *cram_huffman_decode_init(cram_block_compression_hdr *hdr,
				     char *data, int size,
				     enum cram_external_type option,
				     int version, varint_vec *vv) {
    int32_t ncodes = 0, i, j;
    char *cp = data, *data_end = &data[size];
    cram_codec *h;
    cram_huffman_code *codes;
    int32_t val, last_len, max_len = 0;
    int err = 0;
    
    if (option == E_BYTE_ARRAY_BLOCK) {
	fprintf(stderr, "BYTE_ARRAYs not supported by this codec\n");
	return NULL;
    }

    ncodes = vv->varint_get32(&cp, data_end, &err);
    h = calloc(1, sizeof(*h));
    if (!h)
	return NULL;

    h->codec  = E_HUFFMAN;
    h->free   = cram_huffman_decode_free;

    h->huffman.ncodes = ncodes;
    codes = h->huffman.codes = malloc(ncodes * sizeof(*codes));
    if (!codes) {
	free(h);
	return NULL;
    }

    /* Read symbols and bit-lengths */
    if (option == E_LONG) {
	for (i = 0; i < ncodes; i++)
	    codes[i].symbol = vv->varint_get64(&cp, data_end, &err);
    } else {
	for (i = 0; i < ncodes; i++)
	    codes[i].symbol = vv->varint_get32(&cp, data_end, &err);
    }

    if (err) {
	fprintf(stderr, "Malformed huffman header stream\n");
	free(h);
	return NULL;
    }
    i = vv->varint_get32(&cp, data_end, &err);
    if (i != ncodes) {
	fprintf(stderr, "Malformed huffman header stream\n");
	free(h);
	return NULL;
    }

    if (ncodes == 0) {
	/* NULL huffman stream.  Ensure it returns an error if
           anything tries to use it. */
        h->decode = cram_huffman_decode_null;
        return h;
    }

    for (i = 0; i < ncodes; i++) {
        codes[i].len = vv->varint_get32(&cp, data_end, &err);
	if (max_len < codes[i].len)
	    max_len = codes[i].len;
    }
    if (err || cp - data != size || max_len >= ncodes) {
	fprintf(stderr, "Malformed huffman header stream\n");
	free(h);
	return NULL;
    }

    /* Sort by bit length and then by symbol value */
    qsort(codes, ncodes, sizeof(*codes), code_sort);

    /* Assign canonical codes */
    val = -1, last_len = 0;
    for (i = 0; i < ncodes; i++) {
	val++;
	if (codes[i].len > last_len) {
	    while (codes[i].len > last_len) {
		val <<= 1;
		last_len++;
	    }
	}
	codes[i].code = val;
    }

    /*
     * Compute the next starting point, offset by the i'th value.
     * For example if codes 10, 11, 12, 13 are 30, 31, 32, 33 then
     * codes[10..13].p = 30 - 10.
     */
    last_len = 0;
    for (i = j = 0; i < ncodes; i++) {
	if (codes[i].len > last_len) {
	    j = codes[i].code - i;
	    last_len = codes[i].len;
	}
	codes[i].p = j;
    }

//    puts("==HUFF LEN==");
//    for (i = 0; i <= last_len+1; i++) {
//	printf("len %d=%d prefix %d\n", i, h->huffman.lengths[i], h->huffman.prefix[i]); 
//    }
//    puts("===HUFFMAN CODES===");
//    for (i = 0; i < ncodes; i++) {
//	int j;
//	printf("%d: %d %d %d ", i, codes[i].symbol, codes[i].len, codes[i].code);
//	j = codes[i].len;
//	while (j) {
//	    putchar(codes[i].code & (1 << --j) ? '1' : '0');
//	}
//	printf(" %d\n", codes[i].code);
//    }

    if (option == E_BYTE || option == E_BYTE_ARRAY) {
	if (h->huffman.codes[0].len == 0)
	    h->decode = cram_huffman_decode_char0;
	else
	    h->decode = cram_huffman_decode_char;
    } else if (option == E_LONG) {
	if (h->huffman.codes[0].len == 0)
	    h->decode = cram_huffman_decode_long0;
	else
	    h->decode = cram_huffman_decode_long;
    } else if (option == E_INT) {
	if (h->huffman.codes[0].len == 0)
	    h->decode = cram_huffman_decode_int0;
	else
	    h->decode = cram_huffman_decode_int;
    } else {
	return NULL;
    }

    return (cram_codec *)h;
}

int cram_huffman_encode_char0(cram_slice *slice, cram_codec *c,
			      char *in, int in_size) {
    return 0;
}

int cram_huffman_encode_char(cram_slice *slice, cram_codec *c,
			     char *in, int in_size) {
    int i, code, len, r = 0;
    unsigned char *syms = (unsigned char *)in;

    while (in_size--) {
	int sym = *syms++;
	if (sym >= -1 && sym < MAX_HUFF) {
	    i = c->e_huffman.val2code[sym+1];
	    assert(c->e_huffman.codes[i].symbol == sym);
	    code = c->e_huffman.codes[i].code;
	    len  = c->e_huffman.codes[i].len;
	} else {
	    /* Slow - use a lookup table for when sym < MAX_HUFF? */
	    for (i = 0; i < c->e_huffman.nvals; i++) {
		if (c->e_huffman.codes[i].symbol == sym)
		    break;
	    }
	    if (i == c->e_huffman.nvals)
		return -1;
    
	    code = c->e_huffman.codes[i].code;
	    len  = c->e_huffman.codes[i].len;
	}

	r |= store_bits_MSB(c->out, code, len);
    }

    return r;
}

int cram_huffman_encode_long0(cram_slice *slice, cram_codec *c,
			      char *in, int in_size) {
    return 0;
}

int cram_huffman_encode_long(cram_slice *slice, cram_codec *c,
			     char *in, int in_size) {
    int i, code, len, r = 0;
    int64_t *syms = (int64_t *)in;

    while (in_size--) {
	int sym = *syms++;

	if (sym >= -1 && sym < MAX_HUFF) {
	    i = c->e_huffman.val2code[sym+1];
	    assert(c->e_huffman.codes[i].symbol == sym);
	    code = c->e_huffman.codes[i].code;
	    len  = c->e_huffman.codes[i].len;
	} else {
	    /* Slow - use a lookup table for when sym < MAX_HUFFMAN_SYM? */
	    for (i = 0; i < c->e_huffman.nvals; i++) {
		if (c->e_huffman.codes[i].symbol == sym)
		    break;
	    }
	    if (i == c->e_huffman.nvals)
		return -1;
    
	    code = c->e_huffman.codes[i].code;
	    len  = c->e_huffman.codes[i].len;
	}

	r |= store_bits_MSB(c->out, code, len);
    }

    return r;
}

int cram_huffman_encode_int0(cram_slice *slice, cram_codec *c,
			     char *in, int in_size) {
    return 0;
}

int cram_huffman_encode_int(cram_slice *slice, cram_codec *c,
			    char *in, int in_size) {
    int i, code, len, r = 0;
    int *syms = (int *)in;

    while (in_size--) {
	int sym = *syms++;

	if (sym >= -1 && sym < MAX_HUFF) {
	    i = c->e_huffman.val2code[sym+1];
	    assert(c->e_huffman.codes[i].symbol == sym);
	    code = c->e_huffman.codes[i].code;
	    len  = c->e_huffman.codes[i].len;
	} else {
	    /* Slow - use a lookup table for when sym < MAX_HUFFMAN_SYM? */
	    for (i = 0; i < c->e_huffman.nvals; i++) {
		if (c->e_huffman.codes[i].symbol == sym)
		    break;
	    }
	    if (i == c->e_huffman.nvals)
		return -1;
    
	    code = c->e_huffman.codes[i].code;
	    len  = c->e_huffman.codes[i].len;
	}

	r |= store_bits_MSB(c->out, code, len);
    }

    return r;
}

void cram_huffman_encode_free(cram_codec *c) {
    if (!c)
	return;

    if (c->e_huffman.codes)
	free(c->e_huffman.codes);
    free(c);
}

/*
 * Encodes a huffman tree.
 * Returns number of bytes written.
 */
int cram_huffman_encode_store(cram_codec *c, cram_block *b, char *prefix,
			      int version) {
    int i, len = 0;
    cram_huffman_code *codes = c->e_huffman.codes;
    /*
     * Up to code length 127 means 2.5e+26 bytes of data required (worst
     * case huffman tree needs symbols with freqs matching the Fibonacci
     * series). So guaranteed 1 byte per code.
     *
     * Symbols themselves could be 5 bytes (eg -1 is 5 bytes in itf8)
     * and 9 bytes in ltf8.
     *
     * Therefore 10*ncodes + 5 + 5 + 1 + 5 is max memory
     */
    char *tmp = malloc(12*c->e_huffman.nvals+16);
    char *tp = tmp;

    if (!tmp)
	return -1;

    if (prefix) {
	size_t l = strlen(prefix);
	BLOCK_APPEND(b, prefix, l);
	len += l;
    }

    tp += c->vv->varint_put32(tp, NULL, c->e_huffman.nvals);
    if (c->e_huffman.option == E_LONG) {
	for (i = 0; i < c->e_huffman.nvals; i++) {
	    tp += c->vv->varint_put64(tp, NULL, codes[i].symbol);
	}
    } else {
	for (i = 0; i < c->e_huffman.nvals; i++) {
	    tp += c->vv->varint_put32(tp, NULL, codes[i].symbol);
	}
    }

    tp += c->vv->varint_put32(tp, NULL, c->e_huffman.nvals);
    for (i = 0; i < c->e_huffman.nvals; i++) {
	tp += c->vv->varint_put32(tp, NULL, codes[i].len);
    }

    len += c->vv->varint_put32_blk(b, c->codec);
    len += c->vv->varint_put32_blk(b, tp-tmp);
    BLOCK_APPEND(b, tmp, tp-tmp);
    len += tp-tmp;

    free(tmp);

    return len;
}

cram_codec *cram_huffman_encode_init(cram_stats *st,
				     enum cram_external_type option,
				     void *dat,
				     int version, varint_vec *vv) {
    int *vals = NULL, *freqs = NULL, vals_alloc = 0, *lens, code, len;
    int nvals, i, ntot = 0, k;
    int64_t max_val = 0, min_val = INT64_MAX;
    cram_codec *c;
    cram_huffman_code *codes;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec = E_HUFFMAN;

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
		free(c);
		return NULL;
	    }
	}
	vals[nvals] = i;
	freqs[nvals] = st->freqs[i];
	assert(st->freqs[i] > 0);
	ntot += freqs[nvals];
	if (max_val < i) max_val = i;
	if (min_val > i) min_val = i;
	nvals++;
    }
    if (st->h) {
	HashIter *iter=  HashTableIterCreate();
	HashItem *hi;

	while ((hi = HashTableIterNext(st->h, iter))) {
	    if (nvals >= vals_alloc) {
		vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
		vals  = realloc(vals,  vals_alloc * sizeof(int));
		freqs = realloc(freqs, vals_alloc * sizeof(int));
		if (!vals || !freqs)
		    return NULL;
	    }
	    vals[nvals]=(size_t)hi->key;
	    freqs[nvals] = hi->data.i;
	    assert(hi->data.i > 0);
	    ntot += freqs[nvals];
	    if (max_val < i) max_val = i;
	    if (min_val > i) min_val = i;
	    nvals++;
	}
	HashTableIterDestroy(iter);
    }

    assert(nvals > 0);

    freqs = realloc(freqs, 2*nvals*sizeof(*freqs));
    lens = calloc(2*nvals, sizeof(*lens));
    if (!lens || !freqs)
	return NULL;

    /* Inefficient, use pointers to form chain so we can insert and maintain
     * a sorted list? This is currently O(nvals^2) complexity.
     */
    for (;;) {
	int64_t low1 = INT64_MAX, low2 = INT64_MAX;
	int ind1 = 0, ind2 = 0;
	for (i = 0; i < nvals; i++) {
	    if (freqs[i] < 0)
		continue;
	    if (low1 > freqs[i]) 
		low2 = low1, ind2 = ind1, low1 = freqs[i], ind1 = i;
	    else if (low2 > freqs[i])
		low2 = freqs[i], ind2 = i;
	}
	if (low2 == INT64_MAX)
	    break;

	freqs[nvals] = low1 + low2;
	lens[ind1] = nvals;
	lens[ind2] = nvals;
	freqs[ind1] *= -1;
	freqs[ind2] *= -1;
	nvals++;
    }
    nvals = nvals/2+1;

    /* Assign lengths */
    for (i = 0; i < nvals; i++) {
	int code_len = 0;
	for (k = lens[i]; k; k = lens[k])
	    code_len++;
	lens[i] = code_len;
	freqs[i] *= -1;
	//fprintf(stderr, "%d / %d => %d\n", vals[i], freqs[i], lens[i]);
    }


    /* Sort, need in a struct */
    if (!(codes = malloc(nvals * sizeof(*codes))))
	return NULL;
    for (i = 0; i < nvals; i++) {
	codes[i].symbol = vals[i];
	codes[i].len = lens[i];
    }
    qsort(codes, nvals, sizeof(*codes), code_sort);

    /*
     * Generate canonical codes from lengths.
     * Sort by length.
     * Start with 0.
     * Every new code of same length is +1.
     * Every new code of new length is +1 then <<1 per extra length.
     *
     * /\
     * a/\
     * /\/\
     * bcd/\
     *    ef
     * 
     * a 1  0
     * b 3  4 (0+1)<<2
     * c 3  5
     * d 3  6
     * e 4  14  (6+1)<<1
     * f 5  15     
     */
    code = 0; len = codes[0].len;
    for (i = 0; i < nvals; i++) {
	while (len != codes[i].len) {
	    code<<=1;
	    len++;
	}
	codes[i].code = code++;

	if (codes[i].symbol >= -1 && codes[i].symbol < MAX_HUFF)
	    c->e_huffman.val2code[codes[i].symbol+1] = i;

	//fprintf(stderr, "sym %d, code %d, len %d\n",
	//	codes[i].symbol, codes[i].code, codes[i].len);
    }

    free(lens);
    free(vals);
    free(freqs);

    c->e_huffman.codes = codes;
    c->e_huffman.nvals = nvals;
    c->e_huffman.option = option;

    c->free = cram_huffman_encode_free;
    if (option == E_BYTE || option == E_BYTE_ARRAY) {
	if (c->e_huffman.codes[0].len == 0)
	    c->encode = cram_huffman_encode_char0;
	else
	    c->encode = cram_huffman_encode_char;
    } else if (option == E_INT) {
	if (c->e_huffman.codes[0].len == 0)
	    c->encode = cram_huffman_encode_int0;
	else
	    c->encode = cram_huffman_encode_int;
    } else if (option == E_LONG) {
	if (c->e_huffman.codes[0].len == 0)
	    c->encode = cram_huffman_encode_long0;
	else
	    c->encode = cram_huffman_encode_long;
    } else {
	return NULL;
    }
    c->store = cram_huffman_encode_store;
    c->flush = NULL;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * BYTE_ARRAY_LEN
 */
int cram_byte_array_len_decode(cram_slice *slice, cram_codec *c,
			       cram_block *in, char *out,
			       int *out_size) {
    /* Fetch length */
    int32_t len = 0, one = 1;
    int r;

    r = c->byte_array_len.len_codec->decode(slice, c->byte_array_len.len_codec,
                                            in, (char *)&len, &one);
    //printf("ByteArray Len=%d\n", len);

    if (!r && c->byte_array_len.val_codec && len >= 0) {
	r = c->byte_array_len.val_codec->decode(slice,
						c->byte_array_len.val_codec,
						in, out, &len);
    } else {
	return -1;
    }

    *out_size = len;

    return r;
}

void cram_byte_array_len_decode_free(cram_codec *c) {
    if (!c) return;

    if (c->byte_array_len.len_codec)
	c->byte_array_len.len_codec->free(c->byte_array_len.len_codec);

    if (c->byte_array_len.val_codec)
	c->byte_array_len.val_codec->free(c->byte_array_len.val_codec);

    free(c);
}

cram_codec *cram_byte_array_len_decode_init(cram_block_compression_hdr *hdr,
					    char *data, int size,
					    enum cram_external_type option,
					    int version, varint_vec *vv) {
    cram_codec *c;
    char *cp   = data;
    char *endp = data + size;
    int32_t encoding = 0;
    int32_t sub_size = -1;
    int err = 0;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_BYTE_ARRAY_LEN;
    c->decode = cram_byte_array_len_decode;
    c->free   = cram_byte_array_len_decode_free;
    
    encoding = vv->varint_get32(&cp, endp, &err);
    sub_size = vv->varint_get32(&cp, endp, &err);
    if (sub_size < 0 || endp - cp < sub_size)
        goto malformed;
    c->byte_array_len.len_codec = cram_decoder_init(hdr, encoding, cp, sub_size,
						    E_INT, version, vv);
    if (c->byte_array_len.len_codec == NULL)
        goto no_codec;
    cp += sub_size;

    sub_size = -1;
    encoding = vv->varint_get32(&cp, endp, &err);
    sub_size = vv->varint_get32(&cp, endp, &err);
    if (sub_size < 0 || endp - cp < sub_size)
        goto malformed;
    c->byte_array_len.val_codec = cram_decoder_init(hdr, encoding, cp, sub_size,
						    option, version, vv);
    if (c->byte_array_len.val_codec == NULL)
        goto no_codec;
    cp += sub_size;

    if (cp - data != size)
        goto malformed;

    if (!err)
	return c;

 malformed:
    fprintf(stderr, "Malformed byte_array_len header stream\n");
 no_codec:
    free(c);
    return NULL;
}

int cram_byte_array_len_encode(cram_slice *slice, cram_codec *c,
			       char *in, int in_size) {
    int32_t i32 = in_size;
    int r = 0;

    r |= c->e_byte_array_len.len_codec->encode(slice,
					       c->e_byte_array_len.len_codec,
					       (char *)&i32, 1);
    r |= c->e_byte_array_len.val_codec->encode(slice,
					       c->e_byte_array_len.val_codec,
					       in, in_size);
    return r;
}

void cram_byte_array_len_encode_free(cram_codec *c) {
    if (!c)
	return;

    if (c->e_byte_array_len.len_codec)
	c->e_byte_array_len.len_codec->free(c->e_byte_array_len.len_codec);

    if (c->e_byte_array_len.val_codec)
	c->e_byte_array_len.val_codec->free(c->e_byte_array_len.val_codec);

    free(c);
}

int cram_byte_array_len_encode_store(cram_codec *c, cram_block *b,
				     char *prefix, int version) {
    int len = 0, len2, len3;
    cram_codec *tc;
    cram_block *b_len, *b_val;

    if (prefix) {
	size_t l = strlen(prefix);
	BLOCK_APPEND(b, prefix, l);
	len += l;
    }

    tc = c->e_byte_array_len.len_codec; 
    b_len = cram_new_block(0, 0);
    len2 = tc->store(tc, b_len, NULL, version);

    tc = c->e_byte_array_len.val_codec;
    b_val = cram_new_block(0, 0);
    len3 = tc->store(tc, b_val, NULL, version);

    len += c->vv->varint_put32_blk(b, c->codec);
    len += c->vv->varint_put32_blk(b, len2+len3);
    BLOCK_APPEND(b, BLOCK_DATA(b_len), BLOCK_SIZE(b_len));
    BLOCK_APPEND(b, BLOCK_DATA(b_val), BLOCK_SIZE(b_val));

    cram_free_block(b_len);
    cram_free_block(b_val);

    return len + len2 + len3;
}

cram_codec *cram_byte_array_len_encode_init(cram_stats *st,
					    enum cram_external_type option,
					    void *dat,
					    int version, varint_vec *vv) {
    cram_codec *c;
    cram_byte_array_len_encoder *e = (cram_byte_array_len_encoder *)dat;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec = E_BYTE_ARRAY_LEN;
    c->free = cram_byte_array_len_encode_free;
    c->encode = cram_byte_array_len_encode;
    c->store = cram_byte_array_len_encode_store;
    c->flush = NULL;

    c->e_byte_array_len.len_codec = cram_encoder_init(e->len_encoding,
						      st, E_INT, 
						      e->len_dat,
						      version, vv);
    c->e_byte_array_len.val_codec = cram_encoder_init(e->val_encoding,
						      NULL, E_BYTE_ARRAY, 
						      e->val_dat,
						      version, vv);

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * BYTE_ARRAY_STOP
 */
static int cram_byte_array_stop_decode_char(cram_slice *slice, cram_codec *c,
					    cram_block *in, char *out,
					    int *out_size) {
    char *cp, ch;
    cram_block *b = NULL;

    b = cram_get_block_by_id(slice, c->byte_array_stop.content_id);
    if (!b)
        return *out_size?-1:0;

    if (b->idx >= b->uncomp_size)
	return -1;

    cp = (char *)b->data + b->idx;
    if (out) {
	while ((ch = *cp) != (char)c->byte_array_stop.stop) {
	    if (cp - (char *)b->data >= b->uncomp_size)
		return -1;
	    *out++ = ch;
	    cp++;
	}
    } else {
	// Consume input, but produce no output
	while ((ch = *cp) != (char)c->byte_array_stop.stop) {
	    if (cp - (char *)b->data >= b->uncomp_size)
		return -1;
	    cp++;
	}
    }

    *out_size = cp - (char *)(b->data + b->idx);
    b->idx = cp - (char *)b->data + 1;

    return 0;
}

int cram_byte_array_stop_decode_block(cram_slice *slice, cram_codec *c,
				      cram_block *in, char *out_,
				      int *out_size) {
    cram_block *b;
    cram_block *out = (cram_block *)out_;
    char *cp, *out_cp, *cp_end;
    char stop;

    b = cram_get_block_by_id(slice, c->byte_array_stop.content_id);
    if (!b)
        return *out_size?-1:0;

    if (b->idx >= b->uncomp_size)
	return -1;
    cp = (char *)b->data + b->idx;
    cp_end = (char *)b->data + b->uncomp_size;
    out_cp = (char *)BLOCK_END(out);

    stop = c->byte_array_stop.stop;
    if (cp_end - cp < out->alloc - out->byte) {
	while (cp != cp_end && *cp != stop)
	    *out_cp++ = *cp++;
	BLOCK_SIZE(out) = out_cp - (char *)BLOCK_DATA(out);
    } else {
	char *cp_start;
	for (cp_start = cp; cp != cp_end && *cp != stop; cp++)
	    ;
	BLOCK_APPEND(out, cp_start, cp - cp_start);
	BLOCK_GROW(out, cp - cp_start);
    }

    *out_size = cp - (char *)(b->data + b->idx);
    b->idx = cp - (char *)b->data + 1;

    return 0;
}

void cram_byte_array_stop_decode_free(cram_codec *c) {
    if (!c) return;

    free(c);
}

cram_codec *cram_byte_array_stop_decode_init(cram_block_compression_hdr *hdr,
					     char *data, int size,
					     enum cram_external_type option,
					     int version, varint_vec *vv) {
    cram_codec *c;
    unsigned char *cp = (unsigned char *)data;
    int err = 0;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_BYTE_ARRAY_STOP;
    switch (option) {
    case E_BYTE_ARRAY_BLOCK:
        c->decode = cram_byte_array_stop_decode_block;
	break;
    case E_BYTE_ARRAY:
        c->decode = cram_byte_array_stop_decode_char;
	break;
    default:
      fprintf(stderr, "byte_array_stop codec only supports BYTE_ARRAYs.\n");
      free(c);
      return NULL;
    }
    c->free   = cram_byte_array_stop_decode_free;
    
    c->byte_array_stop.stop = *cp++;
    if (CRAM_MAJOR_VERS(version) == 1) {
	c->byte_array_stop.content_id = cp[0] + (cp[1]<<8) + (cp[2]<<16)
	    + (cp[3]<<24);
	cp += 4;
    } else {
	c->byte_array_stop.content_id = vv->varint_get32((char **)&cp, NULL, &err);
    }

    if ((char *)cp - data != size || err) {
	fprintf(stderr, "Malformed byte_array_stop header stream\n");
	free(c);
	return NULL;
    }

    return c;
}

int cram_byte_array_stop_encode(cram_slice *slice, cram_codec *c,
				char *in, int in_size) {
    BLOCK_APPEND(c->out, in, in_size);
    BLOCK_APPEND_CHAR(c->out, c->e_byte_array_stop.stop);
    return 0;
}

void cram_byte_array_stop_encode_free(cram_codec *c) {
    if (!c)
	return;
    free(c);
}

int cram_byte_array_stop_encode_store(cram_codec *c, cram_block *b,
				      char *prefix, int version) {
    int len = 0;
    char buf[20], *cp = buf;

    if (prefix) {
	size_t l = strlen(prefix);
	BLOCK_APPEND(b, prefix, l);
	len += l;
    }

    cp += c->vv->varint_put32(cp, NULL, c->codec);

    if (CRAM_MAJOR_VERS(version) == 1) {
	cp += c->vv->varint_put32(cp, NULL, 5);
	*cp++ = c->e_byte_array_stop.stop;
	*cp++ = (c->e_byte_array_stop.content_id >>  0) & 0xff;
	*cp++ = (c->e_byte_array_stop.content_id >>  8) & 0xff;
	*cp++ = (c->e_byte_array_stop.content_id >> 16) & 0xff;
	*cp++ = (c->e_byte_array_stop.content_id >> 24) & 0xff;
    } else {
	cp += c->vv->varint_put32(cp, NULL, 1 +
				  c->vv->varint_size(c->e_byte_array_stop.content_id));
	*cp++ = c->e_byte_array_stop.stop;
	cp += c->vv->varint_put32(cp, NULL, c->e_byte_array_stop.content_id);
    }

    BLOCK_APPEND(b, buf, cp-buf);
    len += cp-buf;

    return len;
}

cram_codec *cram_byte_array_stop_encode_init(cram_stats *st,
					     enum cram_external_type option,
					     void *dat,
					     int version, varint_vec *vv) {
    cram_codec *c;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec = E_BYTE_ARRAY_STOP;
    c->free = cram_byte_array_stop_encode_free;
    c->encode = cram_byte_array_stop_encode;
    c->store = cram_byte_array_stop_encode_store;
    c->flush = NULL;

    c->e_byte_array_stop.stop = ((int *)dat)[0];
    c->e_byte_array_stop.content_id = ((int *)dat)[1];

    return c;
}

/*
 * ---------------------------------------------------------------------------
 */

const char *cram_encoding2str(enum cram_encoding t) {
    switch (t) {
    case E_NULL:            return "NULL";
    case E_EXTERNAL:        return "EXTERNAL";
    case E_GOLOMB:          return "GOLOMB";
    case E_HUFFMAN:         return "HUFFMAN";
    case E_BYTE_ARRAY_LEN:  return "BYTE_ARRAY_LEN";
    case E_BYTE_ARRAY_STOP: return "BYTE_ARRAY_STOP";
    case E_BETA:            return "BETA";
    case E_SUBEXP:          return "SUBEXP";
    case E_GOLOMB_RICE:     return "GOLOMB_RICE";
    case E_GAMMA:           return "GAMMA";
    case E_XPACK:           return "XPACK";
    case E_XRLE:            return "XRLE";
    case E_NUM_CODECS:
    default:                return "?";
    }
}

static cram_codec *(*decode_init[])(cram_block_compression_hdr *hdr,
				    char *data,
				    int size,
				    enum cram_external_type option,
				    int version, varint_vec *vv) = {
    NULL,
    cram_external_decode_init,
    NULL,
    cram_huffman_decode_init,
    cram_byte_array_len_decode_init,
    cram_byte_array_stop_decode_init,
    cram_beta_decode_init,
    cram_subexp_decode_init,
    NULL,
    cram_gamma_decode_init,
    NULL,
    cram_xpack_decode_init,
    cram_xrle_decode_init,
};

cram_codec *cram_decoder_init(cram_block_compression_hdr *hdr,
			      enum cram_encoding codec,
			      char *data, int size,
			      enum cram_external_type option,
			      int version, varint_vec *vv) {
    if (codec >= E_NULL && codec < E_NUM_CODECS && decode_init[codec]) {
	cram_codec *r = decode_init[codec](hdr, data, size, option, version, vv);
	if (r) {
	    r->vv = vv;
	    r->codec_id = hdr->ncodecs++;
	}
	return r;
    } else {
	fprintf(stderr, "Unimplemented codec of type %s\n", cram_encoding2str(codec));
	return NULL;
    }
}

static cram_codec *(*encode_init[])(cram_stats *stx,
				    enum cram_external_type option,
				    void *opt,
				    int version, varint_vec *vv) = {
    NULL,
    cram_external_encode_init,
    NULL,
    cram_huffman_encode_init,
    cram_byte_array_len_encode_init,
    cram_byte_array_stop_encode_init,
    cram_beta_encode_init,
    NULL, //cram_subexp_encode_init,
    NULL,
    NULL, //cram_gamma_encode_init,
    NULL,
    cram_xpack_encode_init,
    cram_xrle_encode_init,
};

cram_codec *cram_encoder_init(enum cram_encoding codec,
			      cram_stats *st,
			      enum cram_external_type option,
			      void *dat,
			      int version, varint_vec *vv) {
    if (st && !st->nvals)
	return NULL;

    if (encode_init[codec]) {
	cram_codec *r;
	if ((r = encode_init[codec](st, option, dat, version, vv)))
	    r->out = NULL;
	r->vv = vv;
	return r;
    } else {
	fprintf(stderr, "Unimplemented codec of type %s\n", cram_encoding2str(codec));
	abort();
    }
}

/*
 * Returns the content_id used by this codec, also in id2 if byte_array_len.
 * Returns -1 for the CORE block and -2 for unneeded.
 * id2 is only filled out for BYTE_ARRAY_LEN which uses 2 codecs.
 */
int cram_codec_to_id(cram_codec *c, int *id2) {
    int bnum1, bnum2 = -2;

    switch (c->codec) {
    case E_HUFFMAN:
	bnum1 = c->huffman.ncodes == 1 ? -2 : -1;
	break;
    case E_GOLOMB:
    case E_BETA:
    case E_SUBEXP:
    case E_GOLOMB_RICE:
    case E_GAMMA:
	bnum1 = -1;
	break;
    case E_EXTERNAL:
	bnum1 = c->external.content_id;
	break;
    case E_BYTE_ARRAY_LEN:
	bnum1 = cram_codec_to_id(c->byte_array_len.len_codec, NULL);
	bnum2 = cram_codec_to_id(c->byte_array_len.val_codec, NULL);
	break;
    case E_XPACK:
	bnum1 = cram_codec_to_id(c->xpack.sub_codec, NULL);
	break;
    case E_XRLE:
	bnum1 = cram_codec_to_id(c->xrle.len_codec, NULL);
	bnum2 = cram_codec_to_id(c->xrle.lit_codec, NULL);
	break;
    case E_BYTE_ARRAY_STOP:
	bnum1 = c->byte_array_stop.content_id;
	break;
    case E_NULL:
	bnum1 = -2;
	break;
    default:
	fprintf(stderr, "Unknown codec type %d\n", c->codec);
	bnum1 = -1;
    }

    if (id2)
	*id2 = bnum2;
    return bnum1;
}

/*
 * cram_codec structures are specialised for decoding or encoding.
 * Unfortunately this makes turning a decoder into an encoder (such as
 * when transcoding files) problematic.
 *
 * This function converts a cram decoder codec into an encoder version
 * in-place (ie it modifiers the codec itself).
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int cram_codec_decoder2encoder(cram_fd *fd, cram_codec *c) {
    int j;

    switch (c->codec) {
    case E_EXTERNAL:
	// shares struct with decode
	c->free = cram_external_encode_free;
	c->store = cram_external_encode_store;
	if (c->decode == cram_external_decode_int)
	    c->encode = cram_external_encode_int;
	else if (c->decode == cram_external_decode_char)
	    c->encode = cram_external_encode_char;
	else if (c->decode == cram_external_decode_long)
	    c->encode = cram_external_encode_long;
	else if (c->decode == cram_external_decode_block)
	    c->encode = cram_external_encode_char;
	else
	    return -1;
	break;

    case E_HUFFMAN: {
	// New structure, so switch.
	// FIXME: we huffman and e_huffman structs amended, we could
	// unify this.
	cram_codec *t = malloc(sizeof(*t));
	t->codec = E_HUFFMAN;
	t->free = cram_huffman_encode_free;
	t->store = cram_huffman_encode_store;
	t->e_huffman.codes = c->huffman.codes;
	t->e_huffman.nvals = c->huffman.ncodes;
	for (j = 0; j < t->e_huffman.nvals; j++) {
	    int32_t sym = t->e_huffman.codes[j].symbol;
	    if (sym >= -1 && sym < MAX_HUFF)
		t->e_huffman.val2code[sym+1] = j;
	}

	if (c->decode == cram_huffman_decode_char0)
	    t->encode = cram_huffman_encode_char0;
	else if (c->decode == cram_huffman_decode_char)
	    t->encode = cram_huffman_encode_char;
	else if (c->decode == cram_huffman_decode_int0)
	    t->encode = cram_huffman_encode_int0;
	else if (c->decode == cram_huffman_decode_int)
	    t->encode = cram_huffman_encode_int;
	else if (c->decode == cram_huffman_decode_long0)
	    t->encode = cram_huffman_encode_long0;
	else if (c->decode == cram_huffman_decode_long)
	    t->encode = cram_huffman_encode_long;
	else {
	    free(t);
	    return -1;
	}
	*c = *t;
	free(t);
	break;
    }

    case E_BETA:
	// shares struct with decode
	c->free = cram_beta_encode_free;
	c->store = cram_beta_encode_store;
	if (c->decode == cram_beta_decode_long)
	    c->encode = cram_beta_encode_long;
	else if (c->decode == cram_beta_decode_int)
	    c->encode = cram_beta_encode_int;
	else if (c->decode == cram_beta_decode_char)
	    c->encode = cram_beta_encode_char;
	else
	    return -1;
	break;

    case E_XPACK: {
	// shares struct with decode
	cram_codec t = *c;
	t.free = cram_xpack_encode_free;
	t.store = cram_xpack_encode_store;
	if (t.decode == cram_xpack_decode_long)
	    t.encode = cram_xpack_encode_long;
	else if (t.decode == cram_xpack_decode_int)
	    t.encode = cram_xpack_encode_int;
	else if (t.decode == cram_xpack_decode_char)
	    t.encode = cram_xpack_encode_char;
	else
	    return -1;
	t.e_xpack.sub_codec = t.xpack.sub_codec;
	if (cram_codec_decoder2encoder(fd, t.e_xpack.sub_codec) == -1)
	    return -1;
	*c = t;
	break;
    }

    case E_BYTE_ARRAY_LEN: {
	cram_codec *t = malloc(sizeof(*t));
	t->codec = E_BYTE_ARRAY_LEN;
	t->free   = cram_byte_array_len_encode_free;
	t->store  = cram_byte_array_len_encode_store;
	t->encode = cram_byte_array_len_encode;
	t->e_byte_array_len.len_codec = c->byte_array_len.len_codec;
	t->e_byte_array_len.val_codec = c->byte_array_len.val_codec;
	if (cram_codec_decoder2encoder(fd, t->e_byte_array_len.len_codec) == -1 ||
	    cram_codec_decoder2encoder(fd, t->e_byte_array_len.val_codec) == -1) {
	    t->free(t);
	    return -1;
	}

	// {len,val}_{encoding,dat} are undefined, but unused.
	// Leaving them unset here means we can test that assertion.
	*c = *t;
	free(t);
	break;
    }

    case E_BYTE_ARRAY_STOP:
	// shares struct with decode
	c->free   = cram_byte_array_stop_encode_free;
	c->store  = cram_byte_array_stop_encode_store;
	c->encode = cram_byte_array_stop_encode;
	break;

    default:
	return -1;
    }

    return 0;
}

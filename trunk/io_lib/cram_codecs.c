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

#include "io_lib/cram.h"

static char *codec2str(enum cram_encoding codec) {
    switch (codec) {
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
    }

    return "(unknown)";
}

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
    do {
	b = block->data[block->byte] >> block->bit;
	if (--block->bit == -1) {
	    block->bit = 7;
	    block->byte++;
	}
	n++;
    } while (b&1);

    return n-1;
}

static int get_zero_bits_MSB(cram_block *block) {
    int n = 0, b;
    do {
	b = block->data[block->byte] >> block->bit;
	if (--block->bit == -1) {
	    block->bit = 7;
	    block->byte++;
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
static signed int get_bits_MSB(cram_block *block, int nbits) {
    unsigned int val = 0;
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

    for (i = 0; i < nbits; i++)
	//val = (val << 1) | get_bit_MSB(block);
	GET_BIT_MSB(block, val);
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
static void store_bits_MSB(cram_block *block, unsigned int val, int nbits) {
    /* fprintf(stderr, " store_bits: %02x %d\n", val, nbits); */

    /*
     * Use slow mode until we tweak the huffman generator to never generate
     * codes longer than 24-bits.
     */
    unsigned int mask;

    if (block->byte+4 >= block->alloc) {
	if (block->byte) {
	    block->alloc *= 2;
	    block->data = realloc(block->data, block->alloc + 4);
	} else {
	    block->alloc = 1024;
	    block->data = realloc(block->data, block->alloc + 4);
	    block->data[0] = 0; // initialise first byte of buffer
	}
    }

    
    
    if (nbits <= block->bit+1) {
	block->data[block->byte] |= (val << (block->bit+1-nbits));
	if ((block->bit-=nbits) == -1) {
	    block->bit = 7;
	    block->byte++;
	    block->data[block->byte] = 0;
	}
	return;
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
int cram_external_decode_int(cram_slice *slice, cram_codec *c,
			     cram_block *in, char *out, int *out_size) {
    int i;
    char *cp;
    cram_block *b = NULL;

    /* Find the external block: FIXME replace with a lookup table */
    if (slice->block_by_id) {
	if (!(b = slice->block_by_id[c->external.content_id]))
	    return -1;
    } else {
	for (i = 0; i < slice->hdr->num_blocks; i++) {
	    b = slice->block[i];
	    if (b->content_type == EXTERNAL &&
		b->content_id == c->external.content_id) {
		break;
	    }
	}
	if (i == slice->hdr->num_blocks)
	    return -1;
    }

    cp = (char *)b->data + b->idx;
    // E_INT and E_LONG are guaranteed single item queries
    b->idx += itf8_get(cp, (int32_t *)out);
    *out_size = 1;

    return 0;
}

int cram_external_decode_char(cram_slice *slice, cram_codec *c,
			      cram_block *in, char *out,
			      int *out_size) {
    int i;
    char *cp;
    cram_block *b = NULL;

    /* Find the external block: FIXME replace with a lookup table */
    if (slice->block_by_id) {
	if (!(b = slice->block_by_id[c->external.content_id]))
	    return -1;
    } else {
	for (i = 0; i < slice->hdr->num_blocks; i++) {
	    b = slice->block[i];
	    if (b->content_type == EXTERNAL &&
		b->content_id == c->external.content_id) {
		break;
	    }
	}
	if (i == slice->hdr->num_blocks)
	    return -1;
    }

    cp = cram_extract_block(b, *out_size);
    if (!cp)
	return -1;

    memcpy(out, cp, *out_size);
    return 0;
}

int cram_external_decode_block(cram_slice *slice, cram_codec *c,
			      cram_block *in, char *out_,
			      int *out_size) {
    int i;
    char *cp;
    cram_block *b = NULL;
    cram_block *out = (cram_block *)out_;

    /* Find the external block: FIXME replace with a lookup table */
    if (slice->block_by_id) {
	if (!(b = slice->block_by_id[c->external.content_id]))
	    return -1;
    } else {
	for (i = 0; i < slice->hdr->num_blocks; i++) {
	    b = slice->block[i];
	    if (b->content_type == EXTERNAL &&
		b->content_id == c->external.content_id) {
		break;
	    }
	}
	if (i == slice->hdr->num_blocks)
	    return -1;
    }

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

cram_codec *cram_external_decode_init(char *data, int size, enum cram_external_type option) {
    cram_codec *c;
    char *cp = data;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_EXTERNAL;
    if (option == E_INT || option == E_LONG)
	c->decode = cram_external_decode_int;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
	c->decode = cram_external_decode_char;
    else
	c->decode = cram_external_decode_block;
    c->free   = cram_external_decode_free;
    
    cp += itf8_get(cp, &c->external.content_id);

    if (cp - data != size) {
	fprintf(stderr, "Malformed external header stream\n");
	free(c);
	return NULL;
    }

    c->external.type = option;

    return c;
}

int cram_external_encode(cram_slice *slice, cram_codec *c,
			cram_block *out, char *in, int in_size) {
    return -1; // not imp.
}

void cram_external_encode_free(cram_codec *c) {
    if (!c)
	return;
    free(c);
}

int cram_external_encode_store(cram_codec *c, char *buf, char *prefix) {
    char tmp[8192], *cp = buf, *tp = tmp;

    if (prefix) {
	while ((*cp++ = *prefix++))
	    ;
	cp--; // skip nul
    }

    tp += itf8_put(tp, c->e_external.content_id);
    cp += itf8_put(cp, c->codec);
    cp += itf8_put(cp, tp-tmp);
    memcpy(cp, tmp, tp-tmp);
    cp += tp-tmp;

    return cp - buf;
}

cram_codec *cram_external_encode_init(cram_stats *st,
				      enum cram_external_type option,
				      void *dat) {
    cram_codec *c;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec = E_EXTERNAL;
    c->free = cram_external_encode_free;
    c->encode = cram_external_encode;
    c->store = cram_external_encode_store;

    c->e_external.content_id = (int)dat;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * BETA
 */
int cram_beta_decode(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n;

    if (c->beta.nbits) {
	for (i = 0, n = *out_size; i < n; i++)
	    out_i[i] = get_bits_MSB(in, c->beta.nbits) - c->beta.offset;
    } else {
	for (i = 0, n = *out_size; i < n; i++)
	    out_i[i] = 0;
    }

    return 0;
}

void cram_beta_decode_free(cram_codec *c) {
    if (c)
	free(c);
}

cram_codec *cram_beta_decode_init(char *data, int size, enum cram_external_type option) {
    cram_codec *c;
    char *cp = data;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_BETA;
    c->decode = cram_beta_decode;
    c->free   = cram_beta_decode_free;
    
    cp += itf8_get(cp, &c->beta.offset);
    cp += itf8_get(cp, &c->beta.nbits);

    if (cp - data != size) {
	fprintf(stderr, "Malformed beta header stream\n");
	free(c);
	return NULL;
    }

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

cram_codec *cram_subexp_decode_init(char *data, int size, enum cram_external_type option) {
    cram_codec *c;
    char *cp = data;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_SUBEXP;
    c->decode = cram_subexp_decode;
    c->free   = cram_subexp_decode_free;
    
    cp += itf8_get(cp, &c->subexp.offset);
    cp += itf8_get(cp, &c->subexp.k);

    if (cp - data != size) {
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

cram_codec *cram_gamma_decode_init(char *data, int size, enum cram_external_type option) {
    cram_codec *c;
    char *cp = data;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_GAMMA;
    c->decode = cram_gamma_decode;
    c->free   = cram_gamma_decode_free;
    
    cp += itf8_get(cp, &c->gamma.offset);

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
    if (c->huffman.plen)
	free(c->huffman.plen);
    free(c);
}

int cram_huffman_decode_char(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int i, n;

    /* Special case of 0 length codes */
    if (c->huffman.codes[0].len == 0) {
	for (i = 0, n = *out_size; i < n; i++) {
	    out[i] = c->huffman.codes[0].symbol;
	}
	return 0;
    }

    for (i = 0, n = *out_size; i < n; i++) {
	int idx = 0;
	int val = 0, len = 0, last_len = 0;

	for (;;) {
#if 0
	    GET_BIT_MSB(in, val);
	    len++;
#else    
	    int dlen = c->huffman.codes[idx].len - last_len;
	    val <<= dlen;
	    val  |= get_bits_MSB(in, dlen);
	    last_len = (len  += dlen);
#endif

	    //idx = c->huffman.lengths[len] + val-c->huffman.prefix[len];
	    idx = val-c->huffman.plen[len].p + c->huffman.plen[len].l;
	    if (c->huffman.codes[idx].code == val && 
		c->huffman.codes[idx].len == len) {
		out[i] = c->huffman.codes[idx].symbol;
		break;
	    }
	}
    }

    return 0;
}

int cram_huffman_decode_int(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n;

    /* Special case of 0 length codes */
    if (c->huffman.codes[0].len == 0) {
	for (i = 0, n = *out_size; i < n; i++) {
	    out_i[i] = c->huffman.codes[0].symbol;
	}
	return 0;
    }

    for (i = 0, n = *out_size; i < n; i++) {
	int idx = 0;
	int val = 0, len = 0, last_len = 0;

	// Now one bit at a time for remaining checks
	for (;;) {
#if 0
	    GET_BIT_MSB(in, val);
	    len++;
#else    
	    int dlen = c->huffman.codes[idx].len - last_len;
	    val <<= dlen;
	    val  |= get_bits_MSB(in, dlen);
	    last_len = (len  += dlen);
#endif

	    //idx = c->huffman.lengths[len] + val-c->huffman.prefix[len];
	    idx = val-c->huffman.plen[len].p + c->huffman.plen[len].l;
	    if (c->huffman.codes[idx].code == val && 
		c->huffman.codes[idx].len == len) {
		out_i[i] = c->huffman.codes[idx].symbol;
		break;
	    }
	}
    }

    return 0;
}

/*
 * Initialises a huffman decoder from an encoding data stream.
 */
cram_codec *cram_huffman_decode_init(char *data, int size, enum cram_external_type option) {
    int32_t ncodes, i, j;
    char *cp = data;
    cram_codec *h;
    cram_huffman_code *codes;
    int32_t val, last_len, max_len = 0;
    
    cp += itf8_get(cp, &ncodes);
    h = malloc(sizeof(*h));
    if (!h)
	return NULL;

    codes = h->huffman.codes = malloc(ncodes * sizeof(*codes));

    /* Read symbols and bit-lengths */
    for (i = 0; i < ncodes; i++) {
	cp += itf8_get(cp, &codes[i].symbol);
    }

    cp += itf8_get(cp, &i);
    if (i != ncodes) {
	fprintf(stderr, "Malformed huffman header stream\n");
	free(h);
	return NULL;
    }

    for (i = 0; i < ncodes; i++) {
	cp += itf8_get(cp, &codes[i].len);
	if (max_len < codes[i].len)
	    max_len = codes[i].len;
    }
    if (cp - data != size) {
	fprintf(stderr, "Malformed huffman header stream\n");
	free(h);
	return NULL;
    }

    /* Sort by bit length and then by symbol value */
    qsort(codes, ncodes, sizeof(*codes), code_sort);

    /* Assign canonical codes */
    h->huffman.plen = calloc(max_len+2, sizeof(prefix_len));
    val = -1, last_len = 0;
    for (i = 0; i < ncodes; i++) {
	val++;
	if (codes[i].len > last_len) {
	    while (codes[i].len > last_len) {
		val <<= 1;
		last_len++;
	    }
	    h->huffman.plen[last_len].p = val;
	}
	codes[i].code = val;
    }
    h->huffman.plen[last_len+1].p = -1;

    /* Identify length starting points */
    h->huffman.plen[last_len+1].l = -1;
    last_len = -1;
    for (i = j = 0; i < ncodes; i++) {
	if (codes[i].len > last_len) {
	    last_len = codes[i].len;
	    while (j <= last_len)
		h->huffman.plen[j++].l = i;
	}
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

    h->codec  = E_HUFFMAN;
    if (option == E_BYTE || option == E_BYTE_ARRAY)
	h->decode = cram_huffman_decode_char;
    else if (option == E_BYTE_ARRAY_BLOCK)
	abort();
    else
	h->decode = cram_huffman_decode_int;
    h->free   = cram_huffman_decode_free;

    return (cram_codec *)h;
}

int cram_huffman_encode_char(cram_slice *slice, cram_codec *c,
			     cram_block *out, char *in, int in_size) {
    int i, code, len;
    unsigned char *syms = (unsigned char *)in;

    /* Special case of 0 length codes */
    if (c->e_huffman.codes[0].len == 0)
	return 0;

    do {
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

	store_bits_MSB(out, code, len);
    } while (--in_size);

    return 0;
}

int cram_huffman_encode_int(cram_slice *slice, cram_codec *c,
			    cram_block *out, char *in, int in_size) {
    int i, code, len;
    int *syms = (int *)in;

    /* Special case of 0 length codes */
    if (c->e_huffman.codes[0].len == 0)
	return 0;

    do {
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

	store_bits_MSB(out, code, len);
    } while (--in_size);

    return 0;
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
int cram_huffman_encode_store(cram_codec *c, char *buf, char *prefix) {
    int i;
    cram_huffman_code *codes = c->e_huffman.codes;
    /*
     * Up to code length 127 means 2.5e+26 bytes of data required (worst
     * case huffman tree needs symbols with freqs matching the Fibonacci
     * series). So guaranteed 1 byte per code.
     *
     * Symbols themselves could be 5 bytes (eg -1 is 5 bytes in itf8).
     *
     * Therefore 6*ncodes + 5 + 5 + 1 + 5 is max memory
     */
    char *tmp = malloc(6*c->e_huffman.nvals+16);
    char *cp = buf, *tp = tmp;

    if (prefix) {
	while ((*cp++ = *prefix++))
	    ;
	cp--; // skip nul
    }

    tp += itf8_put(tp, c->e_huffman.nvals);
    for (i = 0; i < c->e_huffman.nvals; i++) {
	tp += itf8_put(tp, codes[i].symbol);
    }

    tp += itf8_put(tp, c->e_huffman.nvals);
    for (i = 0; i < c->e_huffman.nvals; i++) {
	tp += itf8_put(tp, codes[i].len);
    }

    cp += itf8_put(cp, c->codec);
    cp += itf8_put(cp, tp-tmp);
    memcpy(cp, tmp, tp-tmp);
    cp += tp-tmp;

    //assert(cp-buf < 8192);
    free(tmp);

    return cp - buf;
}

cram_codec *cram_huffman_encode_init(cram_stats *st,
				     enum cram_external_type option,
				     void *dat) {
    int *vals = NULL, *freqs = NULL, vals_alloc = 0, *lens, code, len;
    int nvals, i, ntot = 0, max_val = 0, min_val = INT_MAX, k;
    cram_codec *c;
    cram_huffman_code *codes;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec = E_HUFFMAN;
    c->free = cram_huffman_encode_free;
    if (option == E_BYTE || option == E_BYTE_ARRAY)
	c->encode = cram_huffman_encode_char;
    else
	c->encode = cram_huffman_encode_int;
    c->store = cram_huffman_encode_store;

    /* Count number of unique symbols */
    for (nvals = i = 0; i < MAX_STAT_VAL; i++) {
	if (!st->freqs[i])
	    continue;
	if (nvals >= vals_alloc) {
	    vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
	    vals  = realloc(vals,  vals_alloc * sizeof(int));
	    freqs = realloc(freqs, vals_alloc * sizeof(int));
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
	    }
	    vals[nvals]=(int)hi->key;
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

    /* Inefficient, use pointers to form chain so we can insert and maintain
     * a sorted list? This is currently O(nvals^2) complexity.
     */
    for (;;) {
	int low1 = INT_MAX, low2 = INT_MAX;
	int ind1 = 0, ind2 = 0;
	for (i = 0; i < nvals; i++) {
	    if (freqs[i] < 0)
		continue;
	    if (low1 > freqs[i]) 
		low2 = low1, ind2 = ind1, low1 = freqs[i], ind1 = i;
	    else if (low2 > freqs[i])
		low2 = freqs[i], ind2 = i;
	}
	if (low2 == INT_MAX)
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
    codes = malloc(nvals * sizeof(*codes));
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
    int32_t len, one = 1;

    c->byte_array_len.len_codec->decode(slice, c->byte_array_len.len_codec, in, (char *)&len, &one);
    //printf("ByteArray Len=%d\n", len);

    if (c->byte_array_len.value_codec) {
	c->byte_array_len.value_codec->decode(slice,
					      c->byte_array_len.value_codec,
					      in, out, &len);
    } else {
	return -1;
    }

    *out_size = len;

    return 0;
}

void cram_byte_array_len_decode_free(cram_codec *c) {
    if (!c) return;

    if (c->byte_array_len.len_codec)
	c->byte_array_len.len_codec->free(c->byte_array_len.len_codec);

    if (c->byte_array_len.value_codec)
	c->byte_array_len.value_codec->free(c->byte_array_len.value_codec);

    free(c);
}

cram_codec *cram_byte_array_len_decode_init(char *data, int size, enum cram_external_type option) {
    cram_codec *c;
    char *cp = data;
    int32_t encoding;
    int32_t sub_size;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_BYTE_ARRAY_LEN;
    c->decode = cram_byte_array_len_decode;
    c->free   = cram_byte_array_len_decode_free;
    
    cp += itf8_get(cp, &encoding);
    cp += itf8_get(cp, &sub_size);
    c->byte_array_len.len_codec = cram_decoder_init(encoding, cp, sub_size, E_INT);
    cp += sub_size;

    cp += itf8_get(cp, &encoding);
    cp += itf8_get(cp, &sub_size);
    c->byte_array_len.value_codec = cram_decoder_init(encoding, cp, sub_size, option);
    cp += sub_size;

    if (cp - data != size) {
	fprintf(stderr, "Malformed byte_array_len header stream\n");
	free(c);
	return NULL;
    }

    return c;
}

int cram_byte_array_len_encode(cram_slice *slice, cram_codec *c,
			cram_block *out, char *in, int in_size) {
    return -1; // not imp.
}

void cram_byte_array_len_encode_free(cram_codec *c) {
    if (!c)
	return;
    free(c);
}

int cram_byte_array_len_encode_store(cram_codec *c, char *buf, char *prefix) {
    char *cp = buf;

    if (prefix) {
	while ((*cp++ = *prefix++))
	    ;
	cp--; // skip nul
    }

    cp += itf8_put(cp, c->codec);
    cp += itf8_put(cp, c->e_byte_array_len.len_len +
		       c->e_byte_array_len.val_len);
    memcpy(cp, c->e_byte_array_len.len_dat, c->e_byte_array_len.len_len);
    cp += c->e_byte_array_len.len_len;
    memcpy(cp, c->e_byte_array_len.val_dat, c->e_byte_array_len.val_len);
    cp += c->e_byte_array_len.val_len;

    return cp - buf;
}

cram_codec *cram_byte_array_len_encode_init(cram_stats *st,
					    enum cram_external_type option,
					    void *dat) {
    cram_codec *c;
    cram_byte_array_len_encoder *e = (cram_byte_array_len_encoder *)dat;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec = E_BYTE_ARRAY_LEN;
    c->free = cram_byte_array_len_encode_free;
    c->encode = cram_byte_array_len_encode;
    c->store = cram_byte_array_len_encode_store;

    c->e_byte_array_len.len_len = e->len_len;
    c->e_byte_array_len.len_dat = e->len_dat;
    c->e_byte_array_len.val_len = e->val_len;
    c->e_byte_array_len.val_dat = e->val_dat;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * BYTE_ARRAY_STOP
 */
int cram_byte_array_stop_decode_char(cram_slice *slice, cram_codec *c,
				     cram_block *in, char *out,
				     int *out_size) {
    int i;
    cram_block *b = NULL;
    char *cp, ch;

    if (slice->block_by_id) {
	if (!(b = slice->block_by_id[c->byte_array_stop.content_id]))
	    return -1;
    } else {
	for (i = 0; i < slice->hdr->num_blocks; i++) {
	    b = slice->block[i];
	    if (b->content_type == EXTERNAL &&
		b->content_id == c->byte_array_stop.content_id) {
		break;
	    }
	}
	if (i == slice->hdr->num_blocks)
	    return -1;
    }

    assert(b->idx < b->uncomp_size);
    cp = (char *)b->data + b->idx;
    while ((ch = *cp) != (char)c->byte_array_stop.stop) {
	assert(cp - (char *)b->data < b->uncomp_size);
	*out++ = ch;
	cp++;
    }

    *out_size = cp - (char *)(b->data + b->idx);
    b->idx = cp - (char *)b->data + 1;

    return 0;
}

int cram_byte_array_stop_decode_block(cram_slice *slice, cram_codec *c,
				      cram_block *in, char *out_,
				      int *out_size) {
    int i, space = 256;
    cram_block *b = NULL;
    cram_block *out = (cram_block *)out_;
    char *cp, ch, *out_cp;

    if (slice->block_by_id) {
	if (!(b = slice->block_by_id[c->byte_array_stop.content_id]))
	    return -1;
    } else {
	for (i = 0; i < slice->hdr->num_blocks; i++) {
	    b = slice->block[i];
	    if (b->content_type == EXTERNAL &&
		b->content_id == c->byte_array_stop.content_id) {
		break;
	    }
	}
	if (i == slice->hdr->num_blocks)
	    return -1;
    }

    assert(b->idx < b->uncomp_size);
    cp = (char *)b->data + b->idx;
    BLOCK_GROW(out, space);
    i = 0;
    out_cp = BLOCK_END(out);
    while ((ch = *cp) != (char)c->byte_array_stop.stop) {
	assert(cp - (char *)b->data < b->uncomp_size);
	*out_cp++ = ch;
	cp++;

	if (++i == space) {
	    BLOCK_SIZE(out) = out_cp - (char *)BLOCK_DATA(out);
	    space *= 2;
	    BLOCK_GROW(out, space);
	    i = 0;
	    out_cp = BLOCK_END(out);
	}
    }
    BLOCK_SIZE(out) = out_cp - (char *)BLOCK_DATA(out);

    *out_size = cp - (char *)(b->data + b->idx);
    b->idx = cp - (char *)b->data + 1;

    return 0;
}

void cram_byte_array_stop_decode_free(cram_codec *c) {
    if (!c) return;

    free(c);
}

cram_codec *cram_byte_array_stop_decode_init(char *data, int size, enum cram_external_type option) {
    cram_codec *c;
    unsigned char *cp = (unsigned char *)data;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_BYTE_ARRAY_STOP;
    c->decode = (option == E_BYTE_ARRAY_BLOCK)
	? cram_byte_array_stop_decode_block
	: cram_byte_array_stop_decode_char;
    c->free   = cram_byte_array_stop_decode_free;
    
    c->byte_array_stop.stop = *cp++;
    c->byte_array_stop.content_id = cp[0] + (cp[1]<<8) + (cp[2]<<16) + (cp[3]<<24);
    cp += 4;

    if ((char *)cp - data != size) {
	fprintf(stderr, "Malformed byte_array_stop header stream\n");
	free(c);
	return NULL;
    }

    return c;
}

int cram_byte_array_stop_encode(cram_slice *slice, cram_codec *c,
			cram_block *out, char *in, int in_size) {
    return -1; // not imp.
}

void cram_byte_array_stop_encode_free(cram_codec *c) {
    if (!c)
	return;
    free(c);
}

int cram_byte_array_stop_encode_store(cram_codec *c, char *buf, char *prefix) {
    char *cp = buf;

    if (prefix) {
	while ((*cp++ = *prefix++))
	    ;
	cp--; // skip nul
    }

    cp += itf8_put(cp, c->codec);
    cp += itf8_put(cp, 5);
    *cp++ = c->e_byte_array_stop.stop;
    *cp++ = (c->e_byte_array_stop.content_id >>  0) & 0xff;
    *cp++ = (c->e_byte_array_stop.content_id >>  8) & 0xff;
    *cp++ = (c->e_byte_array_stop.content_id >> 16) & 0xff;
    *cp++ = (c->e_byte_array_stop.content_id >> 23) & 0xff;

    return cp - buf;
}

cram_codec *cram_byte_array_stop_encode_init(cram_stats *st,
					     enum cram_external_type option,
					     void *dat) {
    cram_codec *c;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec = E_BYTE_ARRAY_STOP;
    c->free = cram_byte_array_stop_encode_free;
    c->encode = cram_byte_array_stop_encode;
    c->store = cram_byte_array_stop_encode_store;

    c->e_byte_array_stop.stop = ((int *)dat)[0];
    c->e_byte_array_stop.content_id = ((int *)dat)[1];

    return c;
}

/*
 * ---------------------------------------------------------------------------
 */

char *cram_encoding2str(enum cram_encoding t) {
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
    }
    return "?";
}

static cram_codec *(*decode_init[])(char *data,
				    int size,
				    enum cram_external_type option) = {
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
};

cram_codec *cram_decoder_init(enum cram_encoding codec,
			      char *data, int size,
			      enum cram_external_type option) {
    if (decode_init[codec]) {
	return decode_init[codec](data, size, option);
    } else {
	fprintf(stderr, "Unimplemented codec of type %s\n", codec2str(codec));
	return NULL;
    }
}

static cram_codec *(*encode_init[])(cram_stats *stx,
				    enum cram_external_type option,
				    void *opt) = {
    NULL,
    cram_external_encode_init,
    NULL,
    cram_huffman_encode_init,
    cram_byte_array_len_encode_init,
    cram_byte_array_stop_encode_init,
    NULL, //cram_beta_encode_init,
    NULL, //cram_subexp_encode_init,
    NULL,
    NULL, //cram_gamma_encode_init,
};

cram_codec *cram_encoder_init(enum cram_encoding codec,
			      cram_stats *st,
			      enum cram_external_type option,
			      void *dat) {
    if (st && !st->nvals)
	return NULL;

    if (encode_init[codec]) {
	return encode_init[codec](st, option, dat);
    } else {
	fprintf(stderr, "Unimplemented codec of type %s\n", codec2str(codec));
	return NULL;
    }
}

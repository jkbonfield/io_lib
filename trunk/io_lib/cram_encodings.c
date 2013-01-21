#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdlib.h>
#include <string.h>

#include "io_lib/cram.h"

/*
 * ---------------------------------------------------------------------------
 * EXTERNAL
 */
int cram_external_decode(cram_slice *slice, cram_codec *c, block_t *in, char *out, int *out_size) {
    int i;
    char *cp;
    cram_block *b = NULL;

    /* Find the external block: FIXME replace with a lookup table */
    for (i = 0; i < slice->hdr->num_blocks; i++) {
	b = slice->block[i];
	if (b->content_type == EXTERNAL &&
	    b->content_id == c->external.content_id) {
	    break;
	}
    }
    if (i == slice->hdr->num_blocks)
	return -1;

    /* FIXME: how to tell string externals from integer externals? */
    if (c->external.type == E_INT || c->external.type == E_LONG) {
	int32_t i32;
	int n;

	cp = b->data + b->idx;

	for (i = 0, n = *out_size; i < n; i+=4) {
	    cp += itf8_get(cp, &i32);
	    ((int32_t *)out)[i] = i32;
	}

	b->idx = cp - b->data;
    } else {
	cp = cram_extract_block(b, *out_size);
	if (!cp)
	    return -1;

	memcpy(out, cp, *out_size);
    }

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
    c->decode = cram_external_decode;
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

/*
 * ---------------------------------------------------------------------------
 * BETA
 */
int cram_beta_decode(cram_slice *slice, cram_codec *c, block_t *in, char *out, int *out_size) {
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
int cram_subexp_decode(cram_slice *slice, cram_codec *c, block_t *in, char *out, int *out_size) {
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

	out_i[count] = val;
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
int cram_gamma_decode(cram_slice *slice, cram_codec *c, block_t *in, char *out, int *out_size) {
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
    if (c->huffman.lengths)
	free(c->huffman.lengths);
    free(c);
}

int cram_huffman_decode(cram_slice *slice, cram_codec *c, block_t *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n;
    int val = 0, len = 0;

    /* Special case of 0 length codes */
    if (c->huffman.codes[0].len == 0) {
	for (i = 0, n = *out_size; i < n; i++) {
	    out_i[i] = c->huffman.codes[0].symbol;
	}
	*out_size = n;

	return 0;
    }

    for (i = 0, n = *out_size; i < n; i++) {
	int idx = 0;

	for (;;) {
	    //val = (val<<1) | get_bit_MSB(in);
	    GET_BIT_MSB(in, val);
	    len++;

	    //printf("val=%d, len=%d\n", val, len);
	    idx = c->huffman.lengths[len];
	    while (c->huffman.codes[idx].len == len) {
		//printf("   check val %d\n", c->huffman.codes[idx].code);
		if (c->huffman.codes[idx].code == val)
		    break;
		idx++;
	    }

	    if (c->huffman.codes[idx].code == val &&
		c->huffman.codes[idx].len  == len) {
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
    int32_t val, last_len;
    
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
    }
    if (cp - data != size) {
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
	while (codes[i].len > last_len) {
	    val <<= 1;
	    last_len++;
	}
	codes[i].code = val;
    }

    /* Identify length starting points */
    h->huffman.lengths = malloc((last_len+1) * sizeof(int));
    last_len = -1;
    for (i = j = 0; i < ncodes; i++) {
	if (codes[i].len > last_len) {
	    last_len = codes[i].len;
	    while (j <= last_len)
		h->huffman.lengths[j++] = i;
	}
    }

//    for (i = 0; i <= last_len; i++) {
//	printf("len %d=%d\n", i, h->huffman.lengths[i]); 
//    }
//    for (i = 0; i < ncodes; i++) {
//	int j;
//	printf("%d: %d %d %d ", i, codes[i].symbol, codes[i].len, codes[i].code);
//	j = codes[i].len;
//	while (j) {
//	    putchar(codes[i].code & (1 << --j) ? '1' : '0');
//	}
//	putchar('\n');
//    }

    h->codec  = E_HUFFMAN;
    h->decode = cram_huffman_decode;
    h->free   = cram_huffman_decode_free;

    return (cram_codec *)h;
}

/*
 * ---------------------------------------------------------------------------
 * BYTE_ARRAY_LEN
 */
int cram_byte_array_len_decode(cram_slice *slice, cram_codec *c, block_t *in, char *out, int *out_size) {
    /* Fetch length */
    int32_t len, one = 1;

    c->byte_array_len.len_codec->decode(slice, c->byte_array_len.len_codec, in, (char *)&len, &one);
    //printf("ByteArray Len=%d\n", len);

    if (c->byte_array_len.value_codec) {
	c->byte_array_len.value_codec->decode(slice, c->byte_array_len.value_codec, in, out, &len);
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
    c->byte_array_len.len_codec = cram_decoder_init(encoding, cp, sub_size, option);
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

/*
 * ---------------------------------------------------------------------------
 * BYTE_ARRAY_STOP
 */
int cram_byte_array_stop_decode(cram_slice *slice, cram_codec *c, block_t *in, char *out, int *out_size) {
    int i;
    cram_block *b = NULL;
    char *cp, ch;

    for (i = 0; i < slice->hdr->num_blocks; i++) {
	b = slice->block[i];
	if (b->content_type == EXTERNAL &&
	    b->content_id == c->byte_array_stop.content_id) {
	    break;
	}
    }
    if (i == slice->hdr->num_blocks)
	return -1;

    cp = b->data + b->idx;
    while ((ch = *cp) != c->byte_array_stop.stop) {
	*out++ = ch;
	cp++;
    }

    *out_size = cp - (b->data + b->idx);
    b->idx = cp - b->data + 1;

    return 0;
}

void cram_byte_array_stop_decode_free(cram_codec *c) {
    if (!c) return;

    free(c);
}

cram_codec *cram_byte_array_stop_decode_init(char *data, int size, enum cram_external_type option) {
    cram_codec *c;
    unsigned char *cp = (unsigned char *)data;
    int32_t encoding;
    int32_t sub_size;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_BYTE_ARRAY_STOP;
    c->decode = cram_byte_array_stop_decode;
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

static cram_codec *(*codec_init[])(char *data, int size, enum cram_external_type option) = {
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

cram_codec *cram_decoder_init(enum cram_encoding codec, char *data, int size, enum cram_external_type option) {
    if (codec_init[codec]) {
	return codec_init[codec](data, size, option);
    } else {
	return NULL;
    }
}

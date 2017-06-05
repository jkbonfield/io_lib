#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>
#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <io_lib/hash_table.h>
#include <io_lib/string_alloc.h>
#include <io_lib/pooled_alloc.h>
#include <io_lib/cram_block_compression.h>
#include <io_lib/rANS_static4x16.h>
#include <bzlib.h>
#include <assert.h>


//=============================================================================
// Name encoder
//=============================================================================

typedef enum {
    C_CAT, C_RLE, C_RANS0, C_RANS1,
    C_PACK0, C_PACK1, C_RLE0, C_RLE1,
    C_PACK_RLE0, C_PACK_RLE1, C_X4
} codec_t;

static int compress_t3(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len, int no_X4);
static int uncompress_t3(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len);

//#define DEBUG

//-----------------------------------------------------------------------------
// Simple variable sized unsigned integers
static int i7put(uint8_t *buf, uint64_t val) {
    uint8_t *b = buf;
    do {
	*b++ = (val & 0x7f) | ((val >= 0x80) << 7);
	val >>= 7;
    } while (val);

    return b-buf;
}

static int i7get(uint8_t *buf, uint64_t *val) {
    uint64_t v = 0;
    uint8_t *b = buf;
    int s = 0;
    uint8_t c;

    do {
	c = *b++;
	v |= (c & 0x7f) << s;
	s += 7;
    } while (c & 0x80);

    *val = v;
    return b - buf;
}

//-----------------------------------------------------------------------------
// Cat: the nul codec, ideal for small data

static int cat_encode(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len) {
    uint8_t *o = out;
    *o++ = C_CAT;
    o += i7put(o, in_len);
    memcpy(o, in, in_len);

    *out_len = o+in_len - out;
    return 0;
}

// Returns number of bytes read from 'in' on success,
//        -1 on failure.
static int64_t cat_decode(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len) {
    uint8_t *i = in;
    uint64_t ulen;
    assert(*i == C_CAT);
    i++;
    i += i7get(i, &ulen);
    assert(ulen <= *out_len);

    memcpy(out, i, ulen);
    *out_len = ulen;

    return (i-in) + ulen;
}

//-----------------------------------------------------------------------------
// Run length encoding.
#ifndef GUARD
#  define GUARD 233
#endif

#ifndef RUN_LEN
#  define RUN_LEN 4
#endif

static int rle_encode(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len) {
    uint64_t i, k = 0;
    int last = -1;
    int run_len = 0;

    out[k++] = C_RLE;
    //out[k++] = GUARD;
    k += i7put(&out[k], in_len);
    
    for (i = 0; i < in_len; i++) {
	if (in[i] == last) {
	    run_len++;
	} else {
	    if (++run_len >= RUN_LEN) {
		k -= run_len * (1 + (last==GUARD));
		out[k++] = GUARD;
		out[k++] = (run_len & 0x7f) + (run_len>=128 ?128 : 0), run_len>>=7;
		if (run_len) out[k++] = (run_len & 0x7f) + (run_len>=128 ?128 : 0), run_len>>=7;
		if (run_len) out[k++] = (run_len & 0x7f) + (run_len>=128 ?128 : 0), run_len>>=7;
		if (run_len) out[k++] = (run_len & 0x7f) + (run_len>=128 ?128 : 0), run_len>>=7;
		if (run_len) out[k++] = (run_len & 0x7f) + (run_len>=128 ?128 : 0), run_len>>=7;
		out[k++] = last;
	    }
	    run_len = 0;
	}

	if (in[i] == GUARD) {
	    out[k++] = GUARD;
	    out[k++] = 0;
	} else {
	    out[k++] = in[i];
	}
	last = in[i];
    }

    // Trailing run
    if (++run_len >= RUN_LEN) {
	k -= run_len * (1 + (last==GUARD));
	out[k++] = GUARD;
	out[k++] = (run_len & 0x7f) + (run_len>=128 ?128 : 0), run_len>>=7;
	if (run_len) out[k++] = (run_len & 0x7f) + (run_len>=128 ?128 : 0), run_len>>=7;
	if (run_len) out[k++] = (run_len & 0x7f) + (run_len>=128 ?128 : 0), run_len>>=7;
	if (run_len) out[k++] = (run_len & 0x7f) + (run_len>=128 ?128 : 0), run_len>>=7;
	if (run_len) out[k++] = (run_len & 0x7f) + (run_len>=128 ?128 : 0), run_len>>=7;
	out[k++] = last;
    }

    *out_len = k;
    return 0;
}

// Returns number of bytes read from 'in' on success,
//        -1 on failure.
static int64_t rle_decode(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len) {
    uint64_t ulen;
    uint64_t i, j;

    assert(in[0] == C_RLE);
    i = 1+i7get(&in[1], &ulen);
    assert(ulen <= *out_len);

    for (j = 0; i < in_len && j < ulen; i++) {
	if (in[i] == GUARD) {
	    if (in[++i] == 0) {
		assert(j+1 <= ulen);
	        out[j++] = GUARD;
	    } else {
		uint32_t run_len = 0;
		unsigned char c, s = 0;
		do {
		    c = in[i++];
		    run_len |= (c & 0x7f) << s;
		    s += 7;
		} while (c & 0x80);
		assert(j+run_len <= ulen);
		memset(&out[j], in[i], run_len);
		j += run_len;
	    }
	} else {
	    assert(j+1 <= ulen);
	    out[j++] = in[i];
	}
    }

    *out_len = ulen;

    return i; // number of compressed bytes consumed.
}

//-----------------------------------------------------------------------------
// rANS codec
static int rans_encode(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len, int method) {
    unsigned int olen = *out_len-10;
    *out = C_RANS0;

    if (rans_compress_to_4x16(in, in_len, out+10, &olen, method) == NULL)
	return -1;

    int nb = i7put(out+1, in_len);
    nb += i7put(out+1+nb, olen);

    assert(nb <= 9);
    if (nb < 9)
	memmove(out+1+nb, out+10, olen);

    //fprintf(stderr, "ENC %d %d\n", (int)in_len, (int)olen);

    *out_len = 1+olen+nb;
    return 0;
}

// Returns number of bytes read from 'in' on success,
//        -1 on failure.
static int64_t rans_decode(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len, int method) {
    unsigned int olen = *out_len;
    assert(*in == C_RANS0);

    uint64_t ulen, clen;
    int nb = i7get(in+1, &ulen);
    nb += i7get(in+1+nb, &clen);

    olen = ulen;
    if (rans_uncompress_to_4x16(in+1+nb, in_len, out, &olen, method) == NULL)
	return -1;

    //fprintf(stderr, "DEC %d %d %d %d\n", (int)in_len, (int)olen, (int)*out_len, (int)olen2);
    *out_len = olen;
    return 1 + nb + clen;
}

static int rans0_encode(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len) {
    unsigned int olen = *out_len-5;
    *out = C_RANS0;

    if (rans_compress_to_4x16(in, in_len, out+5, &olen, 0) == NULL)
	return -1;
    *(uint32_t *)(out+1) = olen;

    *out_len = olen+5;
    return 0;
}

// Returns number of bytes read from 'in' on success,
//        -1 on failure.
static int64_t rans0_decode(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len) {
    unsigned int olen = *out_len;
    assert(*in == C_RANS0);

    if (rans_uncompress_to_4x16(in+5, in_len, out, &olen, 0) == NULL)
	return -1;

    *out_len = olen;
    return 5 + *(uint32_t *)(in+1);
}

static int rans1_encode(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len) {
    unsigned int olen = *out_len-5;
    *out = C_RANS1;

    if (rans_compress_to_4x16(in, in_len, out+5, &olen, 1) == NULL)
	return -1;
    *(uint32_t *)(out+1) = olen;

    *out_len = olen+5;
    return 0;
}

// Returns number of bytes read from 'in' on success,
//        -1 on failure.
static int64_t rans1_decode(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len) {
    unsigned int olen = *out_len;
    assert(*in == C_RANS1);

    if (rans_uncompress_to_4x16(in+5, in_len, out, &olen, 1) == NULL)
	return -1;

    *out_len = olen;
    return 5 + *(uint32_t *)(in+1);
}

//-----------------------------------------------------------------------------
// X4: splitting 32-bit data into 4 8-bit data streams and encoding
// separately.
static int x4_encode(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len) {
    uint64_t j, i4[4];
    uint64_t len4 = (in_len+3)&~3, olen4, olen4_space;
    uint8_t *in4 = malloc(len4);
    uint8_t *out4;
    if (!in4) return -1;

    // Split the 4 interleave signals into their respective arrays
    len4 /= 4;
    for (j = 0; j < 4; j++)
	i4[j] = j*len4;

    for (j = 0; i4[0] < len4; i4[0]++, i4[1]++, i4[2]++, i4[3]++, j+=4) {
	in4[i4[0]] = in[j+0];
	in4[i4[1]] = in[j+1];
	in4[i4[2]] = in[j+2];
	in4[i4[3]] = in[j+3];
    }

    // Encode each signal using the best method per portion.
    out4 = out;
    *out4++ = C_X4;
    out4 += i7put(out4, in_len);
    olen4 = out4-out;

    for (j = 0; j < 4; j++) {
	olen4_space = *out_len - olen4;
	if (compress_t3(in4 + j*len4, len4, out4, &olen4_space, 1) < 0) return -1;
	olen4 += olen4_space;
	out4  += olen4_space;
    }

    *out_len = olen4;
    free(in4);

    return 0;
}

static int x4_decode(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len) {
    uint8_t *i = in, *o, *o_orig;
    uint64_t ulen, olen, i4[4], j;
    assert(*i == C_X4);
    i++;
    i += i7get(i, &ulen);

    o = o_orig = malloc(ulen);
    if (!o)
	return -1;

    in_len -= (i-in);

    // Uncompress
    uint64_t len4 = (ulen+3)&~3;
    len4 /= 4;
    for (j = 0; j < 4; j++) {
	i4[j] = j*len4;
	olen = *out_len - (o-o_orig);
	int64_t clen = uncompress_t3(i, in_len, o, &olen);
	if (clen < 0) {
	    free(o_orig);
	    return -1;
	}
	i += clen;
	o += olen;
    }

    // Reorder
    j = 0;
    while (j < ulen) {
	out[j++] = o_orig[i4[0]++];
	out[j++] = o_orig[i4[1]++];
	out[j++] = o_orig[i4[2]++];
	out[j++] = o_orig[i4[3]++];
    }

    free(o_orig);

    *out_len = j;
    return i-in;
}

static int compress_t3(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len, int no_X4) {
    uint64_t best_sz = UINT64_MAX;
    codec_t best = C_CAT;
    uint64_t olen = *out_len;

    *out_len = olen;
    if (cat_encode(in, in_len, out, out_len) < 0) return -1;
#ifdef DEBUG
    fprintf(stderr, "CAT   -> %ld\n", (long)*out_len);
#endif
    if (best_sz > *out_len) {
	best_sz = *out_len;
	best = C_CAT;
    }

    *out_len = olen;
    if (rle_encode(in, in_len, out, out_len) < 0) return -1;
#ifdef DEBUG
    fprintf(stderr, "RLE   -> %ld\n", (long)*out_len);
#endif
    if (best_sz > *out_len) {
	best_sz = *out_len;
	best = C_RLE;
    }

    int rmethods[] = {0,1,128,129,64,65,192,193}, m;
    for (m = 0; m < 2; m++) {
	//for (m = 0; m < 8; m++) {
	*out_len = olen;
	if (rans_encode(in, in_len, out, out_len, rmethods[m]) < 0) return -1;
#ifdef DEBUG
	fprintf(stderr, "RANS0+%d -> %ld\n", m, (long)*out_len);
#endif
	if (best_sz > *out_len) {
	    best_sz = *out_len;
	    best = C_RANS0+m;
	}
    }

//    *out_len = olen;
//    if (rans0_encode(in, in_len, out, out_len) < 0) return -1;
//#ifdef DEBUG
//    fprintf(stderr, "RANS0 -> %ld\n", (long)*out_len);
//#endif
//    if (best_sz > *out_len) {
//	best_sz = *out_len;
//	best = RANS0;
//    }
//
//    if (in_len >= 4) {
//	*out_len = olen;
//	if (rans1_encode(in, in_len, out, out_len) < 0) return -1;
//#ifdef DEBUG
//	fprintf(stderr, "RANS1 -> %ld\n", (long)*out_len);
//#endif
//	if (best_sz > *out_len) {
//	    best_sz = *out_len;
//	    best = RANS1;
//	}
//    }

    if (!no_X4 && in_len%4 == 0 && in_len >= 32) {
#ifdef DEBUG
	fprintf(stderr, "\n");
#endif
	*out_len = olen;
	if (x4_encode(in, in_len, out, out_len) < 0) return -1;
#ifdef DEBUG
	fprintf(stderr, "X4    -> %ld\n", (long)*out_len);
#endif
	if (best_sz > *out_len) {
	    best_sz = *out_len;
	    best = C_X4;
	}
    }

#ifdef DEBUG
    fprintf(stderr, "Best method = %d, %ld -> %ld\n", best, (long)in_len, (long)best_sz);
#endif

    switch (best) {
    case C_CAT:
	if (cat_encode(in, in_len, out, out_len) < 0) return -1;
	break;

    case C_RLE:
	if (rle_encode(in, in_len, out, out_len) < 0) return -1;
	break;

    case C_RANS0:
    case C_RANS1:
    case C_PACK0:
    case C_PACK1:
    case C_RLE0:
    case C_RLE1:
    case C_PACK_RLE0:
    case C_PACK_RLE1:
	*out_len = olen;
	int rmethods[] = {0,1,128,129,64,65,192,193}, m;
	if (rans_encode(in, in_len, out, out_len, rmethods[best-C_RANS0]) < 0) return -1;
	break;

//    case RANS0:
//	*out_len = olen;
//	if (rans0_encode(in, in_len, out, out_len) < 0) return -1;
//	break;
//
//    case RANS1:
//	*out_len = olen;
//	if (rans1_encode(in, in_len, out, out_len) < 0) return -1;
//	// last method and already in this format.

    default:
	break;
    }

    return 0;
}

static uint64_t uncompressed_size(uint8_t *in, uint64_t in_len) {
    uint64_t ulen;

    switch(*in) {
    case C_CAT:
    case C_RLE:
    case C_X4:
    case C_RANS0:
    case C_RANS1:
    case C_PACK0:
    case C_PACK1:
    case C_RLE0:
    case C_RLE1:
    case C_PACK_RLE0:
    case C_PACK_RLE1:
	i7get(in+1, &ulen);
	break;

    default:
	return -1;
    }

    return ulen;
}

static int uncompress_t3(uint8_t *in, uint64_t in_len, uint8_t *out, uint64_t *out_len) {
    switch (*in) {
    case C_CAT:
	return cat_decode(in, in_len, out, out_len);

    case C_RLE:
	return rle_decode(in, in_len, out, out_len);

    case C_RANS0:
    case C_RANS1:
    case C_PACK0:
    case C_PACK1:
    case C_RLE0:
    case C_RLE1:
    case C_PACK_RLE0:
    case C_PACK_RLE1: {
	int rmethods[] = {0,1,128,129,64,65,192,193}, m;
	return rans_decode(in, in_len, out, out_len, rmethods[*in - C_RANS0]);
    }

//    case RANS0:
//	return rans0_decode(in, in_len, out, out_len);
//
//    case RANS1:
//	return rans1_decode(in, in_len, out, out_len);

    case C_X4:
	return x4_decode(in, in_len, out, out_len);

    default:
	return -1;
    }
}

//-----------------------------------------------------------------------------
// tokenise_name3 code

// FIXME
#define MAX_TOKENS 129
#define MAX_DESCRIPTORS (MAX_TOKENS<<4)

// Number of names per block
// FIXME
#define MAX_NAMES 100000

enum name_type {N_ERR = -1, N_TYPE = 0, N_ALPHA, N_CHAR, N_DZLEN, N_DIGITS0, N_DUP, N_DIFF, 
		N_DIGITS, N_D1, N_D2, N_D3, N_DDELTA, N_DDELTA0, N_MATCH, N_END};

char *types[]={"TYPE", "ALPHA", "CHAR", "DZLEN", "DIG0", "DUP", "DIFF",
	       "DIGITS", "", "", "", "DDELTA", "DDELTA0", "MATCH", "END"};

typedef struct trie trie_t;

typedef struct {
    char *last_name;
    int last_ntok;
    enum name_type last_token_type[MAX_TOKENS];
    int last_token_int[MAX_TOKENS];
    int last_token_str[MAX_TOKENS];
    //int last_token_delta[MAX_TOKENS];
} last_context;


typedef struct {
    uint8_t *buf;
    size_t buf_a, buf_l; // alloc and used length.
    int tnum, ttype;
    int dup_from;
} descriptor;

typedef struct {
    last_context *lc;

    // For finding entire line dups
    int counter;

    // Trie used in encoder only
    trie_t *t_head;
    pool_alloc_t *pool;

    // FIXME: inefficient to have a large MAX size
    //descriptor *desc;
    descriptor desc[MAX_DESCRIPTORS];
} name_context;

name_context *create_context(int max_names) {
    name_context *ctx = malloc(sizeof(*ctx) + max_names*sizeof(*ctx->lc));
    if (!ctx) return NULL;

    ctx->counter = 0;
    ctx->t_head = NULL;

    ctx->lc = (last_context *)(((char *)ctx) + sizeof(*ctx));
    ctx->pool = NULL;

    memset(ctx->desc, 0, MAX_DESCRIPTORS * sizeof(*ctx->desc));
    //ctx->desc = desc;

    return ctx;
}

void free_trie(trie_t *t);
void free_context(name_context *ctx) {
    if (!ctx)
	return;

//    if (ctx->t_head)
//	free_trie(ctx->t_head);
    if (ctx->pool)
	pool_destroy(ctx->pool);

    free(ctx);
}

//-----------------------------------------------------------------------------
// Fast unsigned integer printing code.
// Returns number of bytes written.
static int append_uint32_fixed(char *cp, uint32_t i, uint8_t l) {
    switch (l) {
    case 9:*cp++ = i / 100000000 + '0', i %= 100000000;
    case 8:*cp++ = i / 10000000  + '0', i %= 10000000;
    case 7:*cp++ = i / 1000000   + '0', i %= 1000000;
    case 6:*cp++ = i / 100000    + '0', i %= 100000;
    case 5:*cp++ = i / 10000     + '0', i %= 10000;
    case 4:*cp++ = i / 1000      + '0', i %= 1000;
    case 3:*cp++ = i / 100       + '0', i %= 100;
    case 2:*cp++ = i / 10        + '0', i %= 10;
    case 1:*cp++ = i             + '0';
    case 0:break;
    }
    return l;
}

static int append_uint32_var(char *cp, uint32_t i) {
    char *op = cp;
    uint32_t j;

    //if (i < 10)         goto b0;
    if (i < 100)        goto b1;
    //if (i < 1000)       goto b2;
    if (i < 10000)      goto b3;
    //if (i < 100000)     goto b4;
    if (i < 1000000)    goto b5;
    //if (i < 10000000)   goto b6;
    if (i < 100000000)  goto b7;

    if ((j = i / 1000000000)) {*cp++ = j + '0'; i -= j*1000000000; goto x8;}
    if ((j = i / 100000000))  {*cp++ = j + '0'; i -= j*100000000;  goto x7;}
 b7:if ((j = i / 10000000))   {*cp++ = j + '0'; i -= j*10000000;   goto x6;}
    if ((j = i / 1000000))    {*cp++ = j + '0', i -= j*1000000;    goto x5;}
 b5:if ((j = i / 100000))     {*cp++ = j + '0', i -= j*100000;     goto x4;}
    if ((j = i / 10000))      {*cp++ = j + '0', i -= j*10000;      goto x3;}
 b3:if ((j = i / 1000))       {*cp++ = j + '0', i -= j*1000;       goto x2;}
    if ((j = i / 100))        {*cp++ = j + '0', i -= j*100;        goto x1;}
 b1:if ((j = i / 10))         {*cp++ = j + '0', i -= j*10;         goto x0;}
    if (i)                     *cp++ = i + '0';
    return cp-op;

 x8:*cp++ = i / 100000000 + '0', i %= 100000000;
 x7:*cp++ = i / 10000000  + '0', i %= 10000000;
 x6:*cp++ = i / 1000000   + '0', i %= 1000000;
 x5:*cp++ = i / 100000    + '0', i %= 100000;
 x4:*cp++ = i / 10000     + '0', i %= 10000;
 x3:*cp++ = i / 1000      + '0', i %= 1000;
 x2:*cp++ = i / 100       + '0', i %= 100;
 x1:*cp++ = i / 10        + '0', i %= 10;
 x0:*cp++ = i             + '0';

    return cp-op;
}

//-----------------------------------------------------------------------------
// Example descriptor encoding and IO.
//
// Here we just append to a buffer so we can dump out the results.
// These could then be passed through a static entropy encoder that
// encodes the entire buffer.
//
// Alternatively an adaptive entropy encoder could be place inline
// here to encode as it goes using additional knowledge from the
// supplied context.

// Ensure room for sz more bytes.
static int descriptor_grow(descriptor *fd, uint32_t sz) {
    while (fd->buf_l + sz > fd->buf_a) {
	size_t buf_a = fd->buf_a ? fd->buf_a*2 : 65536;
	unsigned char *buf = realloc(fd->buf, buf_a);
	if (!buf)
	    return -1;
	fd->buf = buf;
	fd->buf_a = buf_a;
    }

    return 0;
}

static int encode_token_type(name_context *ctx, int ntok,
			     enum name_type type) {
    int id = ntok<<4;

    if (descriptor_grow(&ctx->desc[id], 1) < 0) return -1;

    ctx->desc[id].buf[ctx->desc[id].buf_l++] = type;

    return 0;
}

static int encode_token_match(name_context *ctx, int ntok) {
    return encode_token_type(ctx, ntok, N_MATCH);
}

static int encode_token_end(name_context *ctx, int ntok) {
    return encode_token_type(ctx, ntok, N_END);
}

static enum name_type decode_token_type(name_context *ctx, int ntok) {
    int id = ntok<<4;
    if (ctx->desc[id].buf_l >= ctx->desc[id].buf_a) return -1;
    return ctx->desc[id].buf[ctx->desc[id].buf_l++];
}

// int stored as 32-bit quantities
static int encode_token_int(name_context *ctx, int ntok,
			    enum name_type type, uint32_t val) {
    int id = (ntok<<4) | type;

    if (encode_token_type(ctx, ntok, type) < 0) return -1;
    if (descriptor_grow(&ctx->desc[id], 4) < 0)	return -1;

    // Assumes little endian and unalign access OK.
    *(uint32_t *)(ctx->desc[id].buf + ctx->desc[id].buf_l) = val;
    ctx->desc[id].buf_l += 4;

    return 0;
}

// Return 0 on success, -1 on failure;
static int decode_token_int(name_context *ctx, int ntok,
			    enum name_type type, uint32_t *val) {
    int id = (ntok<<4) | type;
    // FIXME: add checks

    // Assumes little endian and unalign access OK.
    *val = *(uint32_t *)(ctx->desc[id].buf + ctx->desc[id].buf_l);
    ctx->desc[id].buf_l += 4;

    return 0;
}

// 8 bit integer quantity
static int encode_token_int1(name_context *ctx, int ntok,
			     enum name_type type, uint32_t val) {
    int id = (ntok<<4) | type;

    if (encode_token_type(ctx, ntok, type) < 0) return -1;
    if (descriptor_grow(&ctx->desc[id], 1) < 0)	return -1;

    ctx->desc[id].buf[ctx->desc[id].buf_l++] = val;

    return 0;
}

static int encode_token_int1_(name_context *ctx, int ntok,
			      enum name_type type, uint32_t val) {
    int id = (ntok<<4) | type;

    if (descriptor_grow(&ctx->desc[id], 1) < 0)	return -1;

    ctx->desc[id].buf[ctx->desc[id].buf_l++] = val;

    return 0;
}

// Return 0 on success, -1 on failure;
static int decode_token_int1(name_context *ctx, int ntok,
			     enum name_type type, uint32_t *val) {
    int id = (ntok<<4) | type;
    // FIXME: add checks

    *val = ctx->desc[id].buf[ctx->desc[id].buf_l++];

    return 0;
}

// Int stored in 4 data series as 4x8 bit quantities
static int encode_token_int4(name_context *ctx, int ntok,
			    enum name_type type, uint32_t val) {
    int id = (ntok<<4) | type;

    if (encode_token_type(ctx, ntok, type) < 0) return -1;
    if (descriptor_grow(&ctx->desc[id  ], 1) < 0)	return -1;
    if (descriptor_grow(&ctx->desc[id+1], 1) < 0)	return -1;
    if (descriptor_grow(&ctx->desc[id+2], 1) < 0)	return -1;
    if (descriptor_grow(&ctx->desc[id+3], 1) < 0)	return -1;

    ctx->desc[id  ].buf[ctx->desc[id  ].buf_l++] = val>>0;
    ctx->desc[id+1].buf[ctx->desc[id+1].buf_l++] = val>>8;
    ctx->desc[id+2].buf[ctx->desc[id+2].buf_l++] = val>>16;
    ctx->desc[id+3].buf[ctx->desc[id+3].buf_l++] = val>>24; 

    return 0;
}

// Return 0 on success, -1 on failure;
static int decode_token_int4(name_context *ctx, int ntok,
			     enum name_type type, uint32_t *val) {
    int id = (ntok<<4) | type;
    // FIXME: add checks

    *val = 
	(ctx->desc[id  ].buf[ctx->desc[id  ].buf_l++] << 0 ) |
	(ctx->desc[id+1].buf[ctx->desc[id+1].buf_l++] << 8 ) |
	(ctx->desc[id+2].buf[ctx->desc[id+2].buf_l++] << 16) |
	(ctx->desc[id+3].buf[ctx->desc[id+3].buf_l++] << 24);

    return 0;
}

// 7 bits at a time with variable size.
static int encode_token_int7(name_context *ctx, int ntok,
			     enum name_type type, uint32_t val) {
    int id = (ntok<<4) | type;

    if (encode_token_type(ctx, ntok, type) < 0) return -1;
    if (descriptor_grow(&ctx->desc[id], 5) < 0)	return -1;

    do {
	ctx->desc[id].buf[ctx->desc[id].buf_l++] = (val & 0x7f) | ((val >= 0x80)<<7);
	val >>= 7;
    } while (val);

    return 0;
}

// Return 0 on success, -1 on failure;
static int decode_token_int7(name_context *ctx, int ntok,
			     enum name_type type, uint32_t *val) {
    int id = (ntok<<4) | type;
    uint32_t v = 0, s = 0;
    uint8_t c;

    // FIXME: add checks
    do {
	c = ctx->desc[id].buf[ctx->desc[id].buf_l++];
	v |= (c & 0x7f) << s;
	s += 7;
    } while (c & 0x80);

    *val = v;
    return 0;
}

//#define encode_token_int encode_token_int7
//#define decode_token_int decode_token_int7

//#define encode_token_int encode_token_int4
//#define decode_token_int decode_token_int4



// Basic C-string style for now.
//
// Maybe XOR with previous string as context?
// This permits partial match to be encoded efficiently.
static int encode_token_alpha(name_context *ctx, int ntok,
			    char *str, int len) {
    int id = (ntok<<4) | N_ALPHA;

    if (encode_token_type(ctx, ntok, N_ALPHA) < 0)  return -1;
    if (descriptor_grow(&ctx->desc[id], len+1) < 0) return -1;
    memcpy(&ctx->desc[id].buf[ctx->desc[id].buf_l], str, len);
    ctx->desc[id].buf[ctx->desc[id].buf_l+len] = 0;
    ctx->desc[id].buf_l += len+1;

    return 0;
}

//// Strings using a separate length model
//static int encode_token_alpha_len(name_context *ctx, int ntok,
//				  char *str, int len) {
//    assert(len < 256); // FIXME
//    int id = (ntok<<4) | N_ALPHA;
//
//    if (encode_token_type(ctx, ntok, N_ALPHA) < 0)  return -1;
//    if (descriptor_grow(&ctx->desc[id],   len) < 0) return -1;
//    if (descriptor_grow(&ctx->desc[id+1], 1) < 0) return -1;
//    memcpy(&ctx->desc[id].buf[ctx->desc[id].buf_l], str, len);
//    ctx->desc[id].buf[ctx->desc[id].buf_l+len] = 0;
//    ctx->desc[id].buf_l += len;
//    ctx->desc[id+1].buf[ctx->desc[id+1].buf_l++] = len;
//
//    return 0;
//}
//#define encode_token_alpha encode_token_alpha_len

// FIXME: need limit on string length for security
// Return length on success, -1 on failure;
static int decode_token_alpha(name_context *ctx, int ntok, char *str) {
    int id = (ntok<<4) | N_ALPHA;
    char c;
    int len = 0;
    do {
	// FIXME: add checks
	c = ctx->desc[id].buf[ctx->desc[id].buf_l++];
	str[len++] = c;
    } while(c);

    return len-1;
}

static int encode_token_char(name_context *ctx, int ntok, char c) {
    int id = (ntok<<4) | N_CHAR;

    if (encode_token_type(ctx, ntok, N_CHAR) < 0) return -1;
    if (descriptor_grow(&ctx->desc[id], 1) < 0)    return -1;
    ctx->desc[id].buf[ctx->desc[id].buf_l++] = c;

    return 0;
}

// FIXME: need limit on string length for security
// Return length on success, -1 on failure;
static int decode_token_char(name_context *ctx, int ntok, char *str) {
    int id = (ntok<<4) | N_CHAR;

    // FIXME: add checks
    *str = ctx->desc[id].buf[ctx->desc[id].buf_l++];

    return 1;
}


// A duplicated name
static int encode_token_dup(name_context *ctx, uint32_t val) {
    return encode_token_int(ctx, 0, N_DUP, val);
}

// Which read to delta against
static int encode_token_diff(name_context *ctx, uint32_t val) {
    return encode_token_int(ctx, 0, N_DIFF, val);
}


//-----------------------------------------------------------------------------
// Trie implementation for tracking common name prefixes.
typedef struct trie {
    char c;
    int count;
    //struct trie *next[128];
    struct trie *next, *sibling;
    int n; // Nth line
} trie_t;

//static trie_t *t_head = NULL;

void free_trie(trie_t *t) {
    trie_t *x, *n;
    for (x = t->next; x; x = n) {
	n = x->sibling;
	free_trie(x);
    }
    free(t);
}

int build_trie(name_context *ctx, char *data, size_t len, int n) {
    int nlines = 0;
    size_t i;
    trie_t *t;

    if (!ctx->t_head)
	ctx->t_head = calloc(1, sizeof(*ctx->t_head));

    // Build our trie, also counting input lines
    for (nlines = i = 0; i < len; i++, nlines++) {
	t = ctx->t_head;
	t->count++;
	while (i < len && data[i] > '\n') {
	    unsigned char c = data[i++];
	    if (c & 0x80)
		//fprintf(stderr, "8-bit ASCII is unsupported\n");
		abort();
	    c &= 127;

	    trie_t *x = t->next, *l = NULL;
	    while (x && x->c != c) {
		l = x; x = x->sibling;
	    }
	    if (!x) {
		if (!ctx->pool)
		    ctx->pool = pool_create(sizeof(trie_t));
		x = (trie_t *)pool_alloc(ctx->pool);
		memset(x, 0, sizeof(*x));
		if (!l)
		    x = t->next    = x;
		else
		    x = l->sibling = x;
		x->n = n;
		x->c = c;
	    }
	    t = x;
	    t->c = c;
	    t->count++;
	}
    }

    return 0;
}

void dump_trie(trie_t *t, int depth) {
    if (depth == 0) {
	printf("graph x_%p {\n    splines = ortho\n    ranksep=2\n", t);
	printf("    p_%p [label=\"\"];\n", t);
	dump_trie(t, 1);
	printf("}\n");
    } else {
	int j, k, count;//, cj;
	char label[100], *cp;
	trie_t *tp = t;

//    patricia:
//	for (count = j = 0; j < 128; j++)
//	    if (t->next[j])
//		count++, cj=j;
//
//	if (count == 1) {
//	    t = t->next[cj];
//	    *cp++ = cj;
//	    goto patricia;
//	}

	trie_t *x;
	for (x = t->next; x; x = x->sibling) {
	    printf("    p_%p [label=\"%c\"];\n", x, x->c);
	    printf("    p_%p -- p_%p [label=\"%d\", penwidth=\"%f\"];\n", tp, x, x->count, MAX((log(x->count)-3)*2,1));
	    //if (depth <= 11)
		dump_trie(x, depth+1);
	}

#if 0	    
	for (j = 0; j < 128; j++) {
	    trie_t *tn;

	    if (!t->next[j])
		continue;

	    cp = label;
	    tn = t->next[j];
	    *cp++ = j;
//	patricia:

	    for (count = k = 0; k < 128; k++)
		if (tn->next[k])
		    count++;//, cj=k;

//	    if (count == 1) {
//		tn = tn->next[cj];
//		*cp++ = cj;
//		goto patricia;
//	    }
	    *cp++ = 0;

	    printf("    p_%p [label=\"%s\"];\n", tn, label);
	    printf("    p_%p -- p_%p [label=\"%d\", penwidth=\"%f\"];\n", tp, tn, tn->count, MAX((log(tn->count)-3)*2,1));
	    if (depth <= 11)
		dump_trie(tn, depth+1);
	}
#endif
    }
}

int search_trie(name_context *ctx, char *data, size_t len, int n, int *exact, int *is_fixed, int *fixed_len) {
    int nlines = 0;
    size_t i, j = -1;
    trie_t *t;
    int from = -1, p3 = -1;

    // Horrid hack for the encoder only.
    // We optimise per known name format here.
    int prefix_len;
    char *d = *data == '@' ? data+1 : data;
    int l   = *data == '@' ? len-1  : len;
    int f = (*data == '>') ? 1 : 0;
    if (l > 70 && d[f+0] == 'm' && d[7] == '_' && d[f+14] == '_' && d[f+61] == '/') {
	prefix_len = 60; // PacBio
	*is_fixed = 0;
    } else if (l == 17 && d[f+5] == ':' && d[f+11] == ':') {
	prefix_len = 7;  // IonTorrent
	*fixed_len = 7;
	*is_fixed = 1;
    } else if (l > 37 && d[f+8] == '-' && d[f+13] == '-' && d[f+18] == '-' && d[f+23] == '-' &&
	       ((d[f+0] >= '0' && d[f+0] <='9') || (d[f+0] >= 'a' && d[f+0] <= 'f')) &&
	       ((d[f+35] >= '0' && d[f+35] <='9') || (d[f+35] >= 'a' && d[f+35] <= 'f'))) {
	// ONT: f33d30d5-6eb8-4115-8f46-154c2620a5da_Basecall_1D_template...
	prefix_len = 37;
	*fixed_len = 37;
	*is_fixed = 1;
    } else {
	// Anything else we give up on the trie method, but we still want to search
	// for exact matches;
	prefix_len = INT_MAX;
	*is_fixed = 0;
    }
    //prefix_len = INT_MAX;

    if (!ctx->t_head)
	ctx->t_head = calloc(1, sizeof(*ctx->t_head));

    // Find an item in the trie
    for (nlines = i = 0; i < len; i++, nlines++) {
	t = ctx->t_head;
	while (i < len && data[i] > '\n') {
	    unsigned char c = data[i++];
	    if (c & 0x80)
		//fprintf(stderr, "8-bit ASCII is unsupported\n");
		abort();
	    c &= 127;

	    trie_t *x = t->next;
	    while (x && x->c != c)
		x = x->sibling;
	    t = x;

//	    t = t->next[c];

	    from = t->n;
	    if (i == prefix_len) p3 = t->n;
	    //if (t->count >= .0035*ctx->t_head->count && t->n != n) p3 = t->n; // pacbio
	    //if (i == 60) p3 = t->n; // pacbio
	    //if (i == 7) p3 = t->n; // iontorrent
	    t->n = n;
	}
    }

    //printf("Looked for %d, found %d, prefix %d\n", n, from, p3);

    *exact = (n != from);
    return *exact ? from : p3;
}


//=============================================================================
// Name encoder
//=============================================================================

/*
 * Tokenises a read name using ctx as context as the previous
 * tokenisation.
 *
 * Parsed elements are then emitted for encoding by calling the
 * encode_token() function with the context, token number (Nth token
 * in line), token type and token value.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
static int encode_name(name_context *ctx, char *name, int len) {
    int i, is_fixed, fixed_len;

    int exact;
    int cnum = ctx->counter++;
    int pnum = search_trie(ctx, name, len, cnum, &exact, &is_fixed, &fixed_len);
    if (pnum < 0) pnum = cnum ? cnum-1 : 0;
    //pnum = pnum & (MAX_NAMES-1);
    //cnum = cnum & (MAX_NAMES-1);
    //if (pnum == cnum) {pnum = cnum ? cnum-1 : 0;}
#ifdef ENC_DEBUG
    fprintf(stderr, "%d: pnum=%d (%d), exact=%d\n%s\n%s\n",
	    ctx->counter, pnum, cnum-pnum, exact, ctx->lc[pnum].last_name, name);
#endif

    // Return DUP or DIFF switch, plus the distance.
    if (exact && len == strlen(ctx->lc[pnum].last_name)) {
	encode_token_dup(ctx, cnum-pnum);
	ctx->lc[cnum].last_name = name;
	ctx->lc[cnum].last_ntok = ctx->lc[pnum].last_ntok;
	// FIXME: optimise this
	int nc = ctx->lc[cnum].last_ntok ? ctx->lc[cnum].last_ntok : MAX_TOKENS;
	memcpy(ctx->lc[cnum].last_token_type, ctx->lc[pnum].last_token_type, nc * sizeof(int));
	memcpy(ctx->lc[cnum].last_token_int , ctx->lc[pnum].last_token_int , nc * sizeof(int));
	memcpy(ctx->lc[cnum].last_token_str , ctx->lc[pnum].last_token_str , nc * sizeof(int));
	return 0;
    }

    encode_token_diff(ctx, cnum-pnum);

    int ntok = 1;
    i = 0;
    if (is_fixed) {
	if (pnum < cnum && ntok < ctx->lc[pnum].last_ntok && ctx->lc[pnum].last_token_type[ntok] == N_ALPHA) {
	    if (ctx->lc[pnum].last_token_int[ntok] == fixed_len && memcmp(name, ctx->lc[pnum].last_name, fixed_len) == 0) {
		encode_token_match(ctx, ntok);
	    } else {
		encode_token_alpha(ctx, ntok, name, fixed_len);
	    }
	} else {
	    encode_token_alpha(ctx, ntok, name, fixed_len);
	}
	ctx->lc[cnum].last_token_int[ntok] = fixed_len;
	ctx->lc[cnum].last_token_str[ntok] = 0;
	ctx->lc[cnum].last_token_type[ntok++] = N_ALPHA;
	i = fixed_len;
    }

    for (; i < len; i++) {
	/* Determine data type of this segment */
	if (isalpha(name[i])) {
	    int s = i+1;

	    // FIXME: try which of these is best.  alnum is good sometimes.
	    //putchar(name[i]);
	    while (s < len && isalpha(name[s]))
	    //while (s < len && isalnum(name[s]))
		//putchar(name[s]),
		s++;

	    // Single byte strings are better encoded as chars.
	    if (s-i == 1) goto n_char;

	    if (pnum < cnum && ntok < ctx->lc[pnum].last_ntok && ctx->lc[pnum].last_token_type[ntok] == N_ALPHA) {
		if (s-i == ctx->lc[pnum].last_token_int[ntok] &&
		    memcmp(&name[i], 
			   &ctx->lc[pnum].last_name[ctx->lc[pnum].last_token_str[ntok]],
			   s-i) == 0) {
#ifdef ENC_DEBUG
		    fprintf(stderr, "Tok %d (alpha-mat, %.*s)\n", N_MATCH, s-i, &name[i]);
#endif
		    if (encode_token_match(ctx, ntok) < 0) return -1;
		} else {
#ifdef ENC_DEBUG
		    fprintf(stderr, "Tok %d (alpha, %.*s / %.*s)\n", N_ALPHA,
		    	    s-i, &ctx->lc[pnum].last_name[ctx->lc[pnum].last_token_str[ntok]], s-i, &name[i]);
#endif
		    // same token/length, but mismatches
		    if (encode_token_alpha(ctx, ntok, &name[i], s-i) < 0) return -1;
		}
	    } else {
#ifdef ENC_DEBUG
		fprintf(stderr, "Tok %d (new alpha, %.*s)\n", N_ALPHA, s-i, &name[i]);
#endif
		if (encode_token_alpha(ctx, ntok, &name[i], s-i) < 0) return -1;
	    }

	    ctx->lc[cnum].last_token_int[ntok] = s-i;
	    ctx->lc[cnum].last_token_str[ntok] = i;
	    ctx->lc[cnum].last_token_type[ntok] = N_ALPHA;

	    i = s-1;
	} else if (name[i] == '0') digits0: {
	    // Digits starting with zero; encode length + value
	    uint32_t s = i;
	    uint32_t v = 0;
	    int d = 0;

	    while (s < len && isdigit(name[s]) && s-i < 8) {
		v = v*10 + name[s] - '0';
		//putchar(name[s]);
		s++;
	    }

	    // TODO: optimise choice over whether to switch from DIGITS to DELTA
	    // regularly vs all DIGITS, also MATCH vs DELTA 0.
	    if (pnum < cnum && ntok < ctx->lc[pnum].last_ntok && ctx->lc[pnum].last_token_type[ntok] == N_DIGITS0) {
		d = v - ctx->lc[pnum].last_token_int[ntok];
		if (d == 0 && ctx->lc[pnum].last_token_str[ntok] == s-i) {
#ifdef ENC_DEBUG
		    fprintf(stderr, "Tok %d (dig-mat, %d)\n", N_MATCH, v);
#endif
		    if (encode_token_match(ctx, ntok) < 0) return -1;
		    //ctx->lc[pnum].last_token_delta[ntok]=0;
		} else if (d < 256 && d >= 0 && ctx->lc[pnum].last_token_str[ntok] == s-i) {
#ifdef ENC_DEBUG
		    fprintf(stderr, "Tok %d (dig-delta, %d / %d)\n", N_DDELTA, ctx->lc[pnum].last_token_int[ntok], v);
#endif
		    //if (encode_token_int1_(ctx, ntok, N_DZLEN, s-i) < 0) return -1;
		    if (encode_token_int1(ctx, ntok, N_DDELTA0, d) < 0) return -1;
		    //ctx->lc[pnum].last_token_delta[ntok]=1;
		} else {
#ifdef ENC_DEBUG
		    fprintf(stderr, "Tok %d (dig, %d / %d)\n", N_DIGITS, ctx->lc[pnum].last_token_int[ntok], v);
#endif
		    if (encode_token_int1_(ctx, ntok, N_DZLEN, s-i) < 0) return -1;
		    if (encode_token_int(ctx, ntok, N_DIGITS0, v) < 0) return -1;
		    //ctx->lc[pnum].last_token_delta[ntok]=0;
		}
	    } else {
#ifdef ENC_DEBUG
		fprintf(stderr, "Tok %d (new dig, %d)\n", N_DIGITS, v);
#endif
		if (encode_token_int1_(ctx, ntok, N_DZLEN, s-i) < 0) return -1;
		if (encode_token_int(ctx, ntok, N_DIGITS0, v) < 0) return -1;
		//ctx->lc[pnum].last_token_delta[ntok]=0;
	    }

	    ctx->lc[cnum].last_token_str[ntok] = s-i; // length
	    ctx->lc[cnum].last_token_int[ntok] = v;
	    ctx->lc[cnum].last_token_type[ntok] = N_DIGITS0;

	    i = s-1;
	} else if (isdigit(name[i])) {
	    // digits starting 1-9; encode value
	    uint32_t s = i;
	    uint32_t v = 0;
	    int d = 0;

	    while (s < len && isdigit(name[s]) && s-i < 8) {
		v = v*10 + name[s] - '0';
		//putchar(name[s]);
		s++;
	    }

	    // If the last token was DIGITS0 and we are the same length, then encode
	    // using that method instead as it seems likely the entire column is fixed
	    // width, sometimes with leading zeros.
	    if (pnum < cnum && ntok < ctx->lc[pnum].last_ntok &&
		ctx->lc[pnum].last_token_type[ntok] == N_DIGITS0 &&
		ctx->lc[pnum].last_token_str[ntok] == s-i)
		goto digits0;
	    
	    // TODO: optimise choice over whether to switch from DIGITS to DELTA
	    // regularly vs all DIGITS, also MATCH vs DELTA 0.
	    if (pnum < cnum && ntok < ctx->lc[pnum].last_ntok && ctx->lc[pnum].last_token_type[ntok] == N_DIGITS) {
		d = v - ctx->lc[pnum].last_token_int[ntok];
		if (d == 0) {
#ifdef ENC_DEBUG
		    fprintf(stderr, "Tok %d (dig-mat, %d)\n", N_MATCH, v);
#endif
		    if (encode_token_match(ctx, ntok) < 0) return -1;
		    //ctx->lc[pnum].last_token_delta[ntok]=0;
		} else if (d < 256 && d >= 0) {
#ifdef ENC_DEBUG
		    fprintf(stderr, "Tok %d (dig-delta, %d / %d)\n", N_DDELTA, ctx->lc[pnum].last_token_int[ntok], v);
#endif
		    if (encode_token_int1(ctx, ntok, N_DDELTA, d) < 0) return -1;
		    //ctx->lc[pnum].last_token_delta[ntok]=1;
		} else {
#ifdef ENC_DEBUG
		    fprintf(stderr, "Tok %d (dig, %d / %d)\n", N_DIGITS, ctx->lc[pnum].last_token_int[ntok], v);
#endif
		    if (encode_token_int(ctx, ntok, N_DIGITS, v) < 0) return -1;
		    //ctx->lc[pnum].last_token_delta[ntok]=0;
		}
	    } else {
#ifdef ENC_DEBUG
		fprintf(stderr, "Tok %d (new dig, %d)\n", N_DIGITS, v);
#endif
		if (encode_token_int(ctx, ntok, N_DIGITS, v) < 0) return -1;
		//ctx->lc[pnum].last_token_delta[ntok]=0;
	    }

	    ctx->lc[cnum].last_token_int[ntok] = v;
	    ctx->lc[cnum].last_token_type[ntok] = N_DIGITS;

	    i = s-1;
	} else {
	n_char:
	    //if (!isalpha(name[i])) putchar(name[i]);
	    if (pnum < cnum && ntok < ctx->lc[pnum].last_ntok && ctx->lc[pnum].last_token_type[ntok] == N_CHAR) {
		if (name[i] == ctx->lc[pnum].last_token_int[ntok]) {
#ifdef ENC_DEBUG
		    fprintf(stderr, "Tok %d (chr-mat, %c)\n", N_MATCH, name[i]);
#endif
		    if (encode_token_match(ctx, ntok) < 0) return -1;
		} else {
#ifdef ENC_DEBUG
		    fprintf(stderr, "Tok %d (chr, %c / %c)\n", N_CHAR, ctx->lc[pnum].last_token_int[ntok], name[i]);
#endif
		    if (encode_token_char(ctx, ntok, name[i]) < 0) return -1;
		}
	    } else {
#ifdef ENC_DEBUG
		fprintf(stderr, "Tok %d (new chr, %c)\n", N_CHAR, name[i]);
#endif
		if (encode_token_char(ctx, ntok, name[i]) < 0) return -1;
	    }

	    ctx->lc[cnum].last_token_int[ntok] = name[i];
	    ctx->lc[cnum].last_token_type[ntok] = N_CHAR;
	}

	ntok++;
	//putchar(' ');
    }

#ifdef ENC_DEBUG
    fprintf(stderr, "Tok %d (end)\n", N_END);
#endif
    if (encode_token_end(ctx, ntok) < 0) return -1;

    //printf("Encoded %.*s with %d tokens\n", len, name, ntok);
    
    ctx->lc[cnum].last_name = name;
    ctx->lc[cnum].last_ntok = ntok;

    return 0;
}


//-----------------------------------------------------------------------------
// Name decoder

// FIXME: should know the maximum name length for safety.
static int decode_name(name_context *ctx, char *name) {
    int t0 = decode_token_type(ctx, 0);
    uint32_t dist;
    int pnum, cnum = ctx->counter++;

    if (t0 < 0)
	return 0;

    decode_token_int(ctx, 0, t0, &dist);
    if ((pnum = cnum - dist) < 0) pnum = 0;

    //fprintf(stderr, "t0=%d, dist=%d, pnum=%d, cnum=%d\n", t0, dist, pnum, cnum);

    if (t0 == N_DUP) {
	strcpy(name, ctx->lc[pnum].last_name);
	// FIXME: optimise this
	ctx->lc[cnum].last_name = name;
	ctx->lc[cnum].last_ntok = ctx->lc[pnum].last_ntok;
	int nc = ctx->lc[cnum].last_ntok ? ctx->lc[cnum].last_ntok : MAX_TOKENS;
	memcpy(ctx->lc[cnum].last_token_type, ctx->lc[pnum].last_token_type, nc * sizeof(int));
	memcpy(ctx->lc[cnum].last_token_int , ctx->lc[pnum].last_token_int , nc * sizeof(int));
	memcpy(ctx->lc[cnum].last_token_str , ctx->lc[pnum].last_token_str , nc * sizeof(int));

	return strlen(name);
    }

    *name = 0;
    int ntok, len = 0, len2;

    for (ntok = 1; ntok < MAX_TOKENS; ntok++) {
	uint32_t v, vl;
	enum name_type tok;
	tok = decode_token_type(ctx, ntok);
	//fprintf(stderr, "Tok %d = %d\n", ntok, tok);

	switch (tok) {
	case N_CHAR:
	    decode_token_char(ctx, ntok, &name[len]);
	    //fprintf(stderr, "Tok %d CHAR %c\n", ntok, name[len]);
	    ctx->lc[cnum].last_token_type[ntok] = N_CHAR;
	    ctx->lc[cnum].last_token_int [ntok] = name[len++];
	    break;

	case N_ALPHA:
	    len2 = decode_token_alpha(ctx, ntok, &name[len]);
	    //fprintf(stderr, "Tok %d ALPHA %.*s\n", ntok, len2, &name[len]);
	    ctx->lc[cnum].last_token_type[ntok] = N_ALPHA;
	    ctx->lc[cnum].last_token_str [ntok] = len;
	    ctx->lc[cnum].last_token_int [ntok] = len2;
	    len += len2;
	    break;

	case N_DIGITS0: // [0-9]*
	    decode_token_int1(ctx, ntok, N_DZLEN, &vl);
	    decode_token_int(ctx, ntok, N_DIGITS0, &v);
	    len += append_uint32_fixed(&name[len], v, vl);
	    //fprintf(stderr, "Tok %d DIGITS0 %0*d\n", ntok, vl, v);
	    ctx->lc[cnum].last_token_type[ntok] = N_DIGITS0;
	    ctx->lc[cnum].last_token_int [ntok] = v;
	    ctx->lc[cnum].last_token_str [ntok] = vl;
	    break;

	case N_DDELTA0:
	    decode_token_int1(ctx, ntok, N_DDELTA0, &v);
	    v += ctx->lc[pnum].last_token_int[ntok];
	    len += append_uint32_fixed(&name[len], v, ctx->lc[pnum].last_token_str[ntok]);
	    //fprintf(stderr, "Tok %d DELTA0 %0*d\n", ntok, ctx->lc[pnum].last_token_str[ntok], v);
	    ctx->lc[cnum].last_token_type[ntok] = N_DIGITS0;
	    ctx->lc[cnum].last_token_int [ntok] = v;
	    ctx->lc[cnum].last_token_str [ntok] = ctx->lc[pnum].last_token_str[ntok];
	    break;

	case N_DIGITS: // [1-9][0-9]*
	    decode_token_int(ctx, ntok, N_DIGITS, &v);
	    len += append_uint32_var(&name[len], v);
	    //fprintf(stderr, "Tok %d DIGITS %d\n", ntok, v);
	    ctx->lc[cnum].last_token_type[ntok] = N_DIGITS;
	    ctx->lc[cnum].last_token_int [ntok] = v;
	    break;

	case N_DDELTA:
	    decode_token_int1(ctx, ntok, N_DDELTA, &v);
	    v += ctx->lc[pnum].last_token_int[ntok];
	    len += append_uint32_var(&name[len], v);
	    //fprintf(stderr, "Tok %d DELTA %d\n", ntok, v);
	    ctx->lc[cnum].last_token_type[ntok] = N_DIGITS;
	    ctx->lc[cnum].last_token_int [ntok] = v;
	    break;

	case N_MATCH:
	    switch (ctx->lc[pnum].last_token_type[ntok]) {
	    case N_CHAR:
		name[len++] = ctx->lc[pnum].last_token_int[ntok];
		//fprintf(stderr, "Tok %d MATCH CHAR %c\n", ntok, ctx->lc[pnum].last_token_int[ntok]);
		ctx->lc[cnum].last_token_type[ntok] = N_CHAR;
		ctx->lc[cnum].last_token_int [ntok] = ctx->lc[pnum].last_token_int[ntok];
		break;

	    case N_ALPHA:
		memcpy(&name[len],
		       &ctx->lc[pnum].last_name[ctx->lc[pnum].last_token_str[ntok]],
		       ctx->lc[pnum].last_token_int[ntok]);
		//fprintf(stderr, "Tok %d MATCH ALPHA %.*s\n", ntok, ctx->lc[pnum].last_token_int[ntok], &name[len]);
		ctx->lc[cnum].last_token_type[ntok] = N_ALPHA;
		ctx->lc[cnum].last_token_str [ntok] = len;
		ctx->lc[cnum].last_token_int [ntok] = ctx->lc[pnum].last_token_int[ntok];
		len += ctx->lc[pnum].last_token_int[ntok];
		break;

	    case N_DIGITS:
		len += append_uint32_var(&name[len], ctx->lc[pnum].last_token_int[ntok]);
		//fprintf(stderr, "Tok %d MATCH DIGITS %d\n", ntok, ctx->lc[pnum].last_token_int[ntok]);
		ctx->lc[cnum].last_token_type[ntok] = N_DIGITS;
		ctx->lc[cnum].last_token_int [ntok] = ctx->lc[pnum].last_token_int[ntok];
		break;

	    case N_DIGITS0:
		len += append_uint32_fixed(&name[len], ctx->lc[pnum].last_token_int[ntok], ctx->lc[pnum].last_token_str[ntok]);
		//fprintf(stderr, "Tok %d MATCH DIGITS %0*d\n", ntok, ctx->lc[pnum].last_token_str[ntok], ctx->lc[pnum].last_token_int[ntok]);
		ctx->lc[cnum].last_token_type[ntok] = N_DIGITS0;
		ctx->lc[cnum].last_token_int [ntok] = ctx->lc[pnum].last_token_int[ntok];
		ctx->lc[cnum].last_token_str [ntok] = ctx->lc[pnum].last_token_str[ntok];
		break;

	    default:
		abort();
	    }
	    break;

	case N_END:
	    name[len++] = 0;
	    ctx->lc[cnum].last_token_type[ntok] = N_END;
	    // FIXME: avoid using memcpy, just keep pointer into buffer?
	    ctx->lc[cnum].last_name = name;
	    ctx->lc[cnum].last_ntok = ntok;
	    
	    return len;

	default:
	    return 0; // eof we hope!
	}
    }


    return -1;
}

// Large enough for whole file for now.
#ifndef BLK_SIZE
#  define BLK_SIZE 10*1024*1024
#endif

static char *decode_block(unsigned char *dat, size_t dat_size, size_t *out_size) {
    dstring_t *ds = dstring_create(NULL);
    char *blk = malloc(BLK_SIZE*2); // FIXME

    if (!blk)
	return NULL;

    uint32_t sz;
    while (dat_size > 0) {
	sz = *(uint32_t *)dat;
	dat += 4; dat_size -= 4;

	if (sz > dat_size) {
	    free(blk);
	    return NULL;
	}

	uint8_t *in = dat;
	dat += sz;
	dat_size -= sz;
	
	// 2.0-2.3 / 1.1-1.2
	name_context *ctx = create_context(MAX_NAMES);
	char *line;
	int i, c, o = 0;

	// Unpack descriptors
	int tnum = -1;
	while (o < sz) {
	    uint8_t ttype = in[o++];
	    if (ttype == 255) {
		uint16_t j = *(uint16_t *)&in[o];
		o += 2;
		ttype = in[o++];
		if (ttype == 0)
		    tnum++;
		i = (tnum<<4) | ttype;

		ctx->desc[i].buf_l = 0;
		ctx->desc[i].buf_a = ctx->desc[j].buf_a;
		ctx->desc[i].buf = malloc(ctx->desc[i].buf_a);
		memcpy(ctx->desc[i].buf, ctx->desc[j].buf, ctx->desc[i].buf_a);
		//fprintf(stderr, "Copy ttype %d, i=%d,j=%d, size %d\n", ttype, i, j, (int)ctx->desc[i].buf_a);
		continue;
	    }

	    if (ttype == 0)
		tnum++;

	    //fprintf(stderr, "Read %02x\n", c);

	    // Load compressed block
	    uint8_t ctype[10];
	    int nb;
	    uint64_t clen, ulen = uncompressed_size(&in[o], sz-o);
	    if (ulen < 0) {
		free(blk);
		return NULL;
	    }
	    i = (tnum<<4) | ttype;

	    ctx->desc[i].buf_l = 0;
	    ctx->desc[i].buf = malloc(ulen);

	    ctx->desc[i].buf_a = ulen;
	    clen = uncompress_t3(&in[o], sz-o, ctx->desc[i].buf, &ctx->desc[i].buf_a);
	    assert(ctx->desc[i].buf_a == ulen);

//	    fprintf(stderr, "%d: Decode tnum %d type %d clen %d ulen %d via %d\n",
//		    o, tnum, ttype, (int)clen, (int)ctx->desc[i].buf_a, ctx->desc[i].buf[0]);

	    o += clen;

	    // Encode tnum 0 type 0 ulen 100000 clen 12530 via 2
	    // Encode tnum 0 type 6 ulen 196800 clen 43928 via 3
	    // Encode tnum 0 type 7 ulen 203200 clen 17531 via 3
	    // Encode tnum 1 type 0 ulen 50800 clen 10 via 1
	    // Encode tnum 1 type 1 ulen 3 clen 5 via 0
	    // Encode tnum 2 type 0 ulen 50800 clen 10 via 1
	    // 	
	}

	int ret;

	line = blk;
	while ((ret = decode_name(ctx, line)) > 0) {
	    dstring_nappend(ds, line, strlen(line)+1);
	    line += ret+1;
	}

	for (i = 0; i < MAX_DESCRIPTORS; i++) {
	    if (ctx->desc[i].buf) {
		free(ctx->desc[i].buf);
		ctx->desc[i].buf = 0;
	    }
	}

	free_context(ctx);

	if (ret < 0) {
	    free(blk);
	    return NULL;
	}
    }

    char *uncomp = DSTRING_STR(ds);
    DSTRING_STR(ds) = NULL;

    dstring_destroy(ds);
    free(blk);

    return uncomp;
}

static char *encode_block(char *blk, int len, size_t *out_size) {
    // Construct trie
    int i, j, ctr = 0;
    //name_context *ctx = calloc(1, sizeof(*ctx));
    name_context *ctx = create_context(MAX_NAMES);
    char *out_blk = malloc(2*len), *out = out_blk; // FIXME

    for (i = j = 0; i < len; j=++i) {
	while (i < len && blk[i] > '\n')
	    i++;
	if (blk[i] >= '\n')
	    break;

	//blk[i] = '\0';
	build_trie(ctx, &blk[j], i-j, ctr++);
    }

    memset(&ctx->desc[0], 0, MAX_DESCRIPTORS * sizeof(ctx->desc[0]));

    //fprintf(stderr, "Processed %d of %d in block, line %d\n", last_start, len, ctr);

    // Encode name
    for (i = j = 0; i < len; j=++i) {
	while (i < len && blk[i] > '\n')
	    i++;
	if (blk[i] > '\n')
	    break;

	blk[i] = '\0';
	if (encode_name(ctx, &blk[j], i-j) < 0) {
	    free_context(ctx);
	    return NULL;
	}
    }

    //dump_trie(ctx->t_head, 0);

    // Serialise descriptors
    int last_tnum = -1;
    uint32_t tot_size = 0;
    int ndesc = 0;
    for (i = 0; i < MAX_DESCRIPTORS; i++) {
	if (!ctx->desc[i].buf_l) continue;

	ndesc++;

	int tnum = i>>4;
	int ttype = i&15;

	if (ttype == 0) {
	    assert(tnum == last_tnum+1);
	    last_tnum = tnum;
	}

	uint64_t out_len = 1.5 * rans_compress_bound_4x16(ctx->desc[i].buf_l, 1); // guesswork
	uint8_t *out = malloc(out_len);
	assert(out);

	//uint8_t ttype8 = ttype;
	//write(1, &ttype8, 1);

	if (compress_t3(ctx->desc[i].buf, ctx->desc[i].buf_l, out, &out_len, 0) < 0)
	    abort();

	free(ctx->desc[i].buf);
	ctx->desc[i].buf = out;
	ctx->desc[i].buf_l = out_len;
	ctx->desc[i].tnum = tnum;
	ctx->desc[i].ttype = ttype;

	// Find dups
	int j;
	for (j = 0; j < i; j++) {
	    if (!ctx->desc[j].buf)
		continue;
	    if (ctx->desc[i].buf_l != ctx->desc[j].buf_l)
		continue;
	    if (memcmp(ctx->desc[i].buf, ctx->desc[j].buf, ctx->desc[i].buf_l) == 0)
		break;
	}
	if (j < i) {
	    //fprintf(stderr, "Dup %d %d size %d\n", i, j, (int)ctx->desc[i].buf_l);
	    ctx->desc[i].dup_from = j;
	    tot_size += 4; // flag, dup_from, ttype
	    //fprintf(stderr, "Desc %d %d/%d => DUP %d\n", i, tnum, ttype, j);
	} else {
	    ctx->desc[i].dup_from = 0;
	    tot_size += out_len + 1; // ttype
	    //fprintf(stderr, "Desc %d %d/%d => %d\n", i, tnum, ttype, (int)ctx->desc[i].buf_l);
	}
	    
	//	    fprintf(stderr, "Encode tnum %d type %d ulen %d clen %d via %d\n",
	//		    tnum, ttype, (int)ctx->desc[i].buf_l, (int)out_len, *out);

	//if (out_len != write(1, out, out_len))
	//    abort();

	//free(out);
	//free(ctx->desc[i].buf);
    }
    //fprintf(stderr, "Serialised %d descriptors\n", ndesc);

    // Write
    *(uint32_t *)out = tot_size;
    out += 4;
    //write(1, &tot_size, 4);
    for (i = 0; i < MAX_DESCRIPTORS; i++) {
	if (!ctx->desc[i].buf_l) continue;
	uint8_t ttype8 = ctx->desc[i].ttype;
	if (ctx->desc[i].dup_from) {
	    *out++ = 255;
	    *(uint16_t *)out = ctx->desc[i].dup_from;
	    out+=2;
	    *out++ = ttype8;
	    //uint8_t x = 255;
	    //write(1, &x, 1);
	    //uint16_t y = ctx->desc[i].dup_from;
	    //write(1, &y, 2);
	    //write(1, &ttype8, 1);
	} else {
	    *out++ = ttype8;
	    memcpy(out, ctx->desc[i].buf, ctx->desc[i].buf_l);
	    out += ctx->desc[i].buf_l;
	    //write(1, &ttype8, 1);
	    //write(1, ctx->desc[i].buf, ctx->desc[i].buf_l);
	}
    }

    // Tidy up memory
    for (i = 0; i < MAX_DESCRIPTORS; i++) {
	if (!ctx->desc[i].buf_l) continue;
	free(ctx->desc[i].buf);
    }

    free_context(ctx);
    
    *out_size = out - out_blk;
    return out_blk;
}


//=============================================================================
// Plugin definitions
//=============================================================================


static const char *name(void) {
    return "Names3 compression";
}

unsigned char *compress_block(int method,
			      int level,
			      cram_slice *s,
			      unsigned char *data,
			      size_t len,
			      size_t *comp_len) {

    return encode_block(data, len, comp_len);
}

unsigned char *uncompress_block(cram_slice *s,
				unsigned char *in,
				size_t in_size,
				size_t *out_size) {
    return decode_block(in, in_size, out_size);
}

static cram_compressor c = {
    'T',//FOUR_CC("\0BSC"),
    1<<DS_RN, // Names only
    1.5,
    name,
    compress_block,
    uncompress_block,
};

cram_compressor *cram_compressor_init(void) {
    return &c;
}

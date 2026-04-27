/*
 * Copyright (c) 2013 Genome Research Ltd.
 * Author(s): James Bonfield, Rob Davies
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
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2010-3
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <fcntl.h>
#include <zlib.h>
#ifdef HAVE_LIBDEFLATE
#include <libdeflate.h>
#endif
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stddef.h>

#include <pthread.h>

#include "io_lib/bam.h"
#include "io_lib/os.h"
#include "io_lib/thread_pool.h"
#include "io_lib/crc32.h"
#include "io_lib/bgzip.h"
#include "io_lib/sam_header.h"

// On later gcc releases the ALLOW_UAC code causes the vectorizor to
// use aligned SIMD instructions on unaligned memory access.  This is due
// to our own abuse of char to int aliasing, but doing things the legal
// way is still slower overall (5-14% depending on system and compiler).
// However by explicitly altering the alignment of the integer types
// we can persuade the compiler to generate the unaligned SIMD
// instructions instead.
#if defined(__GNUC__) && !(defined(__clang__) || defined(__ICC))
typedef  int16_t  int16_u __attribute__ ((aligned (1)));
typedef uint16_t uint16_u __attribute__ ((aligned (1)));
typedef  int32_t  int32_u __attribute__ ((aligned (1)));
typedef uint32_t uint32_u __attribute__ ((aligned (1)));
#else
typedef  int16_t  int16_u;
typedef uint16_t uint16_u;
typedef  int32_t  int32_u;
typedef uint32_t uint32_u;
#endif

// If using Cloudflare's zlib, consider just switching to the zlib crc.
// However this doesn't work on older machines.

//#define iolib_crc32 crc32

#define STORE_UINT16(ucp, val)			\
    *(ucp)++ = ((uint16_t) val)      & 0xff;	\
    *(ucp)++ = ((uint16_t) val >> 8) & 0xff;

#define STORE_UINT32(ucp, val)			\
    *(ucp)++ = ((uint32_t) (val))       & 0xff;	\
    *(ucp)++ = ((uint32_t) (val) >>  8) & 0xff;	\
    *(ucp)++ = ((uint32_t) (val) >> 16) & 0xff;	\
    *(ucp)++ = ((uint32_t) (val) >> 24) & 0xff;

#define STORE_UINT64(ucp, val)			\
    *(ucp)++ = ((uint64_t) (val))       & 0xff;	\
    *(ucp)++ = ((uint64_t) (val) >>  8) & 0xff;	\
    *(ucp)++ = ((uint64_t) (val) >> 16) & 0xff;	\
    *(ucp)++ = ((uint64_t) (val) >> 24) & 0xff; \
    *(ucp)++ = ((uint64_t) (val) >> 32) & 0xff; \
    *(ucp)++ = ((uint64_t) (val) >> 40) & 0xff; \
    *(ucp)++ = ((uint64_t) (val) >> 48) & 0xff; \
    *(ucp)++ = ((uint64_t) (val) >> 56) & 0xff;

/* Custom strtol for aux tags, always base 10 */
static int64_t inline STRTOL64(const char *v, const char **rv, int b) {
    int64_t n = 0;
    int neg = 1;
    switch(*v) {
    case '-':
	neg=-1;
	break;
    case '+':
	break;
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':
	n = *v - '0';
	break;
    default:
	*rv = v;
	return 0;
    }
    v++;

    while (isdigit(*v))
	n = n*10 + *v++ - '0';
    *rv = v;
    return neg*n;
}

static int8_t aux_type_size[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 1, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 8, 0, 4, 0, 0, 4, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

/*
 * Skips to the next tag start.
 * Returns next tag pointer (or end of buffer) on success,
 *         NULL on failure
 */
char *bam_aux_skip(const char *s_) {
    const uint8_t *s = (const uint8_t *)s_;
    int sz;
    if ((sz = aux_type_size[s[2]])) {
	s += sz + 3;
    } else {
	switch(s[2]) {
	case 'Z':
	case 'H': {  /* Variable length, null terminated */
	    s += 3;
	    while (*s++);
	    ;
	    break;
	}
	case 'B': {  /* Array types */
	    uint32_t count;
	    if ((sz = aux_type_size[s[3]])) {
		count = (   (uint32_t) s[4]
			    + ((uint32_t) s[5] << 8)
			    + ((uint32_t) s[6] << 16)
			    + ((uint32_t) s[7] << 24));
		s += 8 + count * sz;
	    } else {
		return NULL;
	    }
	    break;
	}
	default:
	    return NULL;
	}
    }

    return (char *)s;
}

/*
 * Looks for aux field 'key' and returns the type + value.
 * The type is the first char and the value is the 2nd character onwards.
 *
 * Returns NULL if not found.
 */
char *bam_aux_find(bam_seq_t *b, const char *key) {
    char *cp = bam_aux(b);

    while (*cp) {
	if (cp[0] == key[0] && cp[1] == key[1])
	    return cp+2;

	if (!(cp = bam_aux_skip(cp)))
	    return NULL;
    }

    return NULL;
}

int32_t bam_aux_i(const uint8_t *dat) {
    switch(dat[0]) {
    case 'i':
	return (int32_t)(dat[1] + (dat[2]<<8) + (dat[3]<<16) + (dat[4]<<24));
    case 'I':
	return (uint32_t)(dat[1] + (dat[2]<<8) + (dat[3]<<16) + (dat[4]<<24));
	break;
    case 's':
	return (int16_t)(dat[1] + (dat[2]<<8));
    case 'S':
	return (uint16_t)(dat[1] + (dat[2]<<8));
    case 'c':
	return (int8_t)dat[1];
    case 'C':
	return (uint8_t)dat[1];
    }

    abort();
}

float bam_aux_f(const uint8_t *dat) {
    assert(dat[0] == 'f');
    union {
	uint32_t i;
	float    f;
    } f;
    f.i = (dat[1]<<0) + (dat[2]<<8) + (dat[3]<<16) + (dat[4]<<24);
    return f.f;
}

double bam_aux_d(const uint8_t *dat) { 
    assert(dat[0] == 'd');
    union {
	uint64_t i;
	double   d;
    } d;
    d.i = (((uint64_t)dat[1])<<0)+
	  (((uint64_t)dat[2])<<8)+
	  (((uint64_t)dat[3])<<16)+
	  (((uint64_t)dat[4])<<24)+
	  (((uint64_t)dat[5])<<32)+
	  (((uint64_t)dat[6])<<40)+
	  (((uint64_t)dat[7])<<48)+
	  (((uint64_t)dat[8])<<54);
    return d.d;
}

char bam_aux_A(const uint8_t *dat) {
    assert(dat[0] == 'A');
    return dat[1];
}

char *bam_aux_Z(const uint8_t *dat) {
    assert(dat[0] == 'Z' || dat[0] == 'H');
    return (char *)(dat+1);
}

/*
 * An iterator on bam aux fields. NB: This code is not reentrant or multi-
 * thread capable. The values returned are valid until the next call to
 * this function.
 * key:  points to an array of 3 characters (eg "RGi", "NMC")
 * type: points to an address of 1 character (eg 'Z', 'i') for the SAM type
 * val:  points to an address of a bam_aux_t union.
 *
 * The first two bytes of key are the real key and the next byte is the
 * BAM type field. Note that this may differ to the returned SAM type
 * field. For example key[2] == 'S' for unsigned short while *type is
 * set to 'i'.
 *
 * Pass in *iter_handle as NULL to initialise the search and then
 * pass in the modified value on each subsequent call to continue the search.
 *
 * Returns 0 if the next value is valid, setting key, type and val.
 *        -1 when no more found.
 */
int bam_aux_iter_full(bam_seq_t *b, char **iter_handle,
		      char *key, char *type, bam_aux_t *val) {
    char *s;

    if (!iter_handle || !*iter_handle) {
	s = (char *)bam_aux(b);
    } else {
	s = *iter_handle;
    }

    /* We null terminate our aux list for ease */
    if (s[0] == 0)
	return -1;

    key[0] = s[0];
    key[1] = s[1];
    key[2] = s[2];
    
    switch (s[2]) {
    case 'A':
	if (type) *type = 'A';
	if (val) val->i = *(s+3);
	s+=4;
	break;

    case 'C':
	if (type) *type = 'i';
	if (val) val->i = *(uint8_t *)(s+3);
	s+=4;
	break;

    case 'c':
	if (type) *type = 'i';
	if (val) val->i = *(int8_t *)(s+3);
	s+=4;
	break;

    case 'S':
	if (type) *type = 'i';
	if (val)
	    val->i = (uint16_t)((((unsigned char *)s)[3]<< 0) +
				(((unsigned char *)s)[4]<< 8));
	s+=5;
	break;

    case 's':
	if (type) *type = 'i';
	if (val)
	    val->i = (int16_t)((((unsigned char *)s)[3]<< 0) +
			       (((unsigned char *)s)[4]<< 8));
	s+=5;
	break;

    case 'I':
	if (type) *type = 'i';
	if (val)
	    val->i = (uint32_t)((((unsigned char *)s)[3]<< 0) +
				(((unsigned char *)s)[4]<< 8) +
				(((unsigned char *)s)[5]<<16) +
				(((unsigned char *)s)[6]<<24));
	s+=7;
	break;

    case 'i':
	if (type) *type = 'i';
	if (val)
	    val->i = (int32_t)((((unsigned char *)s)[3]<< 0) +
			       (((unsigned char *)s)[4]<< 8) +
			       (((unsigned char *)s)[5]<<16) +
			       (((unsigned char *)s)[6]<<24));
	s+=7;
	break;

    case 'f':
	if (type) *type = 'f';
	if (val) /* Assume same endianness as integer */
	    val->i = (int32_t)((((unsigned char *)s)[3]<< 0) +
			       (((unsigned char *)s)[4]<< 8) +
			       (((unsigned char *)s)[5]<<16) +
			       (((unsigned char *)s)[6]<<24));
	s+=7;
	break;

    case 'd':
	if (type) *type = 'd';
	if (val) /* Assume same endianness as integer */
	    val->i64 = (uint64_t)(((uint64_t)(((unsigned char *)s)[ 3])<< 0) +
				  ((uint64_t)(((unsigned char *)s)[ 4])<< 8) +
				  ((uint64_t)(((unsigned char *)s)[ 5])<<16) +
				  ((uint64_t)(((unsigned char *)s)[ 6])<<24) +
				  ((uint64_t)(((unsigned char *)s)[ 7])<<32) +
				  ((uint64_t)(((unsigned char *)s)[ 8])<<40) +
				  ((uint64_t)(((unsigned char *)s)[ 9])<<48) +
				  ((uint64_t)(((unsigned char *)s)[10])<<54));
	s+=11;
	break;

    case 'Z': case 'H':
	if (type) *type = s[2];
	s+=3;
	if (val) val->s = s;
	while (*s++);
	break;

    case 'B': {
	uint32_t count;
	if (type) *type = 'B';
	count = (unsigned int)((((unsigned char *)s)[4]<< 0) +
			       (((unsigned char *)s)[5]<< 8) +
			       (((unsigned char *)s)[6]<<16) +
			       (((unsigned char *)s)[7]<<24));

	if (val) {
	    val->B.n = count;
	    val->B.t = s[3];
	    val->B.s = (unsigned char *)s+8;
	}
	s+=8;

	switch(val->B.t) {
	case 'c': case 'C': s +=   count; break;
	case 's': case 'S': s += 2*count; break;
	case 'i': case 'I': s += 4*count; break;
	case 'f':           s += 4*count; break;
	default:
	    fprintf(stderr, "Unknown sub-type '%c' for aux type 'B'\n",
		    val->B.t);
	    return -1;
	}
	break;
    }

    default:
	fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	return -1;
    }

    if (iter_handle)
	*iter_handle = s;

    return 0;
}

/*
 * As above, but only 2 characters of the key are returned so the
 * original BAM type is not visible.
 *
 * Note this can cause ambiguities if you wish to distinguish between
 * -1 billion and +3 billion.
 */
int bam_aux_iter(bam_seq_t *b, char **iter_handle,
		 char *key, char *type, bam_aux_t *val) {
    char k3[3];
    int r = bam_aux_iter_full(b, iter_handle, k3, type, val);

    if (r == 0) {
	key[0] = k3[0];
	key[1] = k3[1];
    }

    return r;
}

static int reg2bin(int start, int end) {
    if (end>start) end--;
    if ((start>>14) == (end>>14)) return ((1<<15)-1)/7 + (start>>14);
    if ((start>>17) == (end>>17)) return ((1<<12)-1)/7 + (start>>17);
    if ((start>>20) == (end>>20)) return ((1<<9 )-1)/7 + (start>>20);
    if ((start>>23) == (end>>23)) return ((1<<6 )-1)/7 + (start>>23);
    if ((start>>26) == (end>>26)) return ((1<<3 )-1)/7 + (start>>26);
    return 0;
}

/*
 * Constructs a bam_seq_t from components.
 * Ignores auxiliary tags for now.
 *
 * Returns -1 on error
 *          number of bytes written to bam_seq_t on success (ie tag offset)
 */
int bam_construct_seq(bam_seq_t **b, size_t extra_len,
		      const char *qname, size_t qname_len,
		      int flag,
		      int rname,      // Ref ID
		      int64_t pos, // first aligned base (1-based)
		      int64_t end, // last aligned base (to calculate bin)
		      int mapq,
		      uint32_t ncigar, const uint32_t *cigar,
		      int mrnm,       // Mate Ref ID
		      int64_t mpos,
		      int64_t isize,
		      int len,
		      const char *seq,
		      const char *qual) {
    size_t required;
    char *cp;
    int i;
    uint32_t *ip;

    /*
     * cp = "=ACMGRSVTWYHKDBN";
     * memset(L, 15, 256);
     * for (i = 0; i < 16; i++) {
     *     L[cp[i]] = L[tolower(cp[i])] = i;
     * }
     */
    static const char L[256] = {
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15, 0,15,15,
	15, 1,14, 2,13,15,15, 4,11,15,15,12,15, 3,15,15,
	15,15, 5, 6, 8,15, 7, 9,15,10,15,15,15,15,15,15,
	15, 1,14, 2,13,15,15, 4,11,15,15,12,15, 3,15,15,
	15,15, 5, 6, 8,15, 7, 9,15,10,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15
    };

    /* Sanity checks */
    if (NULL == b) return -1;
    if (len < 0) return -1;  /* not sure why the spec has it as an int */
    if (qname_len > 0 && NULL == qname) return -1;
    if (ncigar > 0 && NULL == cigar) return -1;
    if (len > 0 && NULL == seq) return -1;

    /* Reallocate if needed */
    required = (sizeof(**b)             /* the struct itself */
		+ round4(qname_len + 1) /* query name (aligned) */
		+ 4 * ncigar            /* CIGAR string */
		+ (len + 1) / 2         /* Sequence, 2 bases per byte */
		+ len                   /* Quality */
		+ extra_len + 1);       /* Extra for optional tags */

    if (NULL == *b || (*b)->alloc < required) {
	bam_seq_t *new_bam = realloc(*b, required);
	if (NULL == new_bam) return -1;
	*b = new_bam;
	(*b)->alloc = required;
    }

    (*b)->ref = rname;
    (*b)->pos = pos-1;
    bam_set_map_qual(*b, mapq);
    bam_set_name_len(*b, qname_len+1);
    bam_set_flag(*b, flag);
    bam_set_cigar_len(*b, ncigar);
    (*b)->len = len;
    (*b)->mate_ref = mrnm;
    (*b)->mate_pos = mpos-1;
    (*b)->ins_size = isize;

    cp = bam_name(*b);
    memcpy(cp, qname, qname_len);
    cp[qname_len] = 0;

    /* Cigar */
    cp = (char *)bam_cigar(*b);
    ip = (uint32_t *)cp;
    for (i = 0; i < ncigar; i++) {
	ip[i] = cigar[i];
    }
    cp += ncigar*4;

    /* Bin */
    if (!((*b)->flag & BAM_CIGAR32)) {
	if (0 == end) { /* Calculate end from pos and cigar */
	    end = pos;
	    for (i = 0; i < ncigar; i++) {
		if (BAM_CONSUME_REF(cigar[i] & BAM_CIGAR_MASK))
		    end += cigar[i] >> BAM_CIGAR_SHIFT;
	    }
	}

	// range is [beg,end) and zero based.
	bam_set_bin(*b, reg2bin(pos-1,end));
    }

    /* Seq */
    for (i = 0; i < len-1; i += 2) {
	*cp++ = (L[(uc)seq[i]]<<4) + L[(uc)seq[i+1]];
    }
    if (i < len)
	*cp++ = L[(uc)seq[i]]<<4;

    /* Qual */
    if (qual) {
	memcpy(cp, qual, len);
	cp += len;
    } else {
	for (i = 0; i < len; i++) {
	    *cp++ = '\xff';
	}
    }

    *cp = 0; /* terminate aux list, for ease of parsing later */

    /* cp now points to the auxiliary tags if required */
    (*b)->blk_size = (int)(cp-(char *)&(*b)->ref);
    return (int)(cp-(char *)(*b));
}

int bam_aux_add(bam_seq_t **b, const char tag[2], char type,
		uint32_t array_len, const void *data) {
    int tlen;
    size_t len;
    size_t used;
    uint8_t *cp;
#ifndef SP_LITTLE_ENDIAN
    uint32_t i;
#endif

    if (NULL == b || NULL == *b) return -1;

    /* Find size  of data type and how much space is needed */
    if (0 == (tlen = aux_type_size[(uint8_t) type])) {
	if (type == 'H' || type == 'Z') { /* Variable length types */
	    if (array_len != 0) return -1; /* No arrays for these allowed */
	    tlen = strlen((const char *) data) + 1;
	} else {
	    /* unknown type */
	    return -1;
	}
    }

    len = array_len > 0 ? 8 + tlen * array_len : 3 + tlen;
    
    /* Find the end of the existing tags and ensure there is enough space */
    cp = (uint8_t *)&(*b)->ref + (*b)->blk_size;
    used = cp - (uint8_t *)(*b);

    if ((*b)->alloc < used + len + 1) { /* + 1 for NUL terminator */
	size_t required = used + len + 1;
	bam_seq_t *new_bam = realloc((*b), required);
	if (NULL == new_bam) return -1;
	*b = new_bam;
	(*b)->alloc = required;
	cp = (uint8_t *)new_bam + used;
    }

    /* Append the data */
    *cp++ = tag[0];
    *cp++ = tag[1];
    if (array_len > 0) {  /* Array type */
	*cp++ = 'B';
	*cp++ = type;
	STORE_UINT32(cp, array_len);
    } else {
	*cp++ = type;
    }

    if (array_len == 0) array_len = 1;
#ifdef SP_LITTLE_ENDIAN
    memcpy(cp, data, array_len * tlen);
    cp += array_len * tlen;
#else
    switch (type) {
    case 'A': case 'c': case 'C':
	memcpy(cp, data, array_len);
	cp += array_len;
	break;
    case 's':
    case 'S': {
	uint16_t *sdata = (uint16_t *) data;
	for (i = 0; i < array_len; i++) {
	    STORE_UINT16(cp, sdata[i]);
	}
	break;
    }
    case 'i':
    case 'I':
    case 'f': {
	uint32_t *idata = (uint32_t *) data;
	for (i = 0; i < array_len; i++) {
	    STORE_UINT32(cp, idata[i]);
	}
	break;
    }
    case 'd': {
	uint64_t *ddata = (uint64_t *) data;
	for (i = 0; i < array_len; i++) {
	    STORE_UINT64(cp, ddata[i]);
	}
	break;
    }
    case 'H': case 'Z':
	memcpy(cp, data, tlen);
	cp += tlen;
	break;
    }
#endif
    
    /* Put a NUL at the end for bam_aux_iter */
    *cp = 0;

    /* Update block_size */
    (*b)->blk_size = (uint32_t)(cp - (uint8_t *)&(*b)->ref);
    
    return 0;
}

/* Calculate space needed to store tags */

ssize_t bam_aux_size_vec(uint32_t count, bam_aux_tag_t *tags) {
    uint32_t i;
    ssize_t len = 0;
    int sz;

    if (NULL == tags) return -1;
    
    for (i = 0; i < count; i++) {
	switch (tags[i].type) {
	case 'C': case 'S': case 'I':
	    if (tags[i].value.ui < 256) {
		sz = 1;
	    } else if (tags[i].value.ui < 65536) {
		sz = 2;
	    } else {
		sz = 4;
	    }
	    break;
	case 'c': case 's': case 'i':
	    if (tags[i].value.i >= -128 && tags[i].value.i < 128) {
		sz = 1;
	    } else if (tags[i].value.i >= -32768 && tags[i].value.i < 32768) {
		sz = 2;
	    } else {
		sz = 4;
	    }
	    break;
	case 'A':
	    sz = 1;
	    break;
	case 'f':
	    sz = 4;
	    break;
	case 'd':
	    sz = 8;
	    break;
	case 'H': case 'Z':
	    if (tags[i].array_len != 0) return -1;
	    sz = strlen(tags[i].value.z) + 1;
	    break;
	default:
	    return -1; /* bad data type */
	}
	len += tags[i].array_len == 0 ? sz + 3 : sz * tags[i].array_len + 8;
    }

    return len + 1;  /* + 1 for NUL byte at end */
}

int bam_aux_add_vec(bam_seq_t **b, uint32_t count, bam_aux_tag_t *tags) {
    ssize_t required = bam_aux_size_vec(count, tags);
    uint32_t i;
    size_t used;

    if (required < 0) return -1;
    if (NULL == b || NULL == *b) return -1;

    /* Find the end of the existing tags and ensure there is enough space.
       Do this once for the entire vector so we don't keep reallocing */
    used = (uint8_t *)&(*b)->ref + (*b)->blk_size - (uint8_t *)(*b);

    if ((*b)->alloc < used + required) {
	bam_seq_t *new_bam = realloc((*b), used + required);
	if (NULL == new_bam) return -1;
	*b = new_bam;
	(*b)->alloc = used + required;
    }

    /* Add the tags, storing integers in the most appropriate size */
    for (i = 0; i < count; i++) {
	if (tags[i].array_len > 0) {
	    /* Deal with array tags */
	    if (bam_aux_add(b, tags[i].tag, tags[i].type,
			    tags[i].array_len, tags[i].value.array)) return -1;
	    continue;
	}

	/* Non-array tags, storing integers as the most appropriate size */
	switch (tags[i].type) {
	case 'C': case 'S': case 'I':
	    if (tags[i].value.ui < 256) {
		uint8_t byte = tags[i].value.ui;
		if (bam_aux_add(b, tags[i].tag, 'C', 0, &byte)) return -1;
	    } else if (tags[i].value.ui < 65536) {
		uint16_t word = tags[i].value.ui;
		if (bam_aux_add(b, tags[i].tag, 'S', 0, &word)) return -1;
	    } else {
		if (bam_aux_add(b, tags[i].tag, 'I', 0, &tags[i].value.ui)) {
		    return -1;
		}
	    }
	    break;
	case 'c': case 's': case 'i':
	    if (tags[i].value.i >= -128 && tags[i].value.i < 128) {
		int8_t byte = tags[i].value.i;
		if (bam_aux_add(b, tags[i].tag, 'c', 0, &byte)) return -1;
	    } else if (tags[i].value.i >= -32768 && tags[i].value.i < 32768) {
		int16_t word = tags[i].value.i;
		if (bam_aux_add(b, tags[i].tag, 's', 0, &word)) return -1;
	    } else {
		if (bam_aux_add(b, tags[i].tag, 'i', 0, &tags[i].value.i))
		    return -1;
	    }
	    break;
	case 'A':
	    if (bam_aux_add(b, tags[i].tag, 'A', 0, &tags[i].value.a))
		return -1;
	    break;
	case 'f': case 'd':
	    if (bam_aux_add(b, tags[i].tag, tags[i].type, 0, &tags[i].value.f))
		return -1;
	    break;
	case 'H': case 'Z':
	    if (bam_aux_add(b, tags[i].tag, tags[i].type, 0, tags[i].value.z))
		return -1;
	    break;
	default:
	    return -1; /* unknown type */
	}	    
    }
    return 0;
}

/* Add SAM-formatted aux tags to a bam_seq_t struct.
   This is basically a copy of the code in sam_next_seq.  Unfortunately
   trying to get them to use a common version slows sam_next_seq down
   rather a lot, even when inlined.  Hence this extra copy. */

int bam_aux_add_from_sam(bam_seq_t **bsp, char *sam) {
    unsigned char *cpf = (unsigned char *) sam;
    unsigned char *cpt = (unsigned char *)&(*bsp)->ref + (*bsp)->blk_size;
    unsigned char *end = (unsigned char *)(*bsp) + (*bsp)->alloc;

    while (*cpf) {
	unsigned char *key = cpf, *value;
	size_t max_len;

	if (!(key[0] && key[1] && key[2] == ':' && key[3] && key[4] == ':'))
	    return -1;
	cpf += 5;

	value = cpf;
	while (*cpf && *cpf != '\t')
	    cpf++;

	if (aux_type_size[key[3]]) {
            max_len = aux_type_size[key[3]] + 3;
        } else if (key[3] != 'B') {
            max_len = cpf - value + 4;
        } else {
	    /* Worst case */
            max_len = (cpf - value) * 4 + 8;
        }

	/* ensure we have enough room */
        if (end - cpt < max_len) {
            size_t used = cpt - (unsigned char *)(*bsp);
            bam_seq_t *new_bam = realloc(*bsp, used + max_len);
            if (NULL == new_bam) return -1;
            *bsp = new_bam;
            (*bsp)->alloc += used + max_len;
            cpt = (unsigned char *)(*bsp) + used;
            end = (unsigned char *)(*bsp) + (*bsp)->alloc;
        }

	*cpt++ = key[0];
	*cpt++ = key[1];

	switch(key[3]) {
	    int64_t n;

	case 'A':
	    *cpt++ = 'A';
	    *cpt++ = *value;
	    break;

	case 'i':
	    //n = atoi((char *)value);
	    n = STRTOL64((char *)value, (const char **)&value, 10);
	    if (n >= 0) {
		if (n < 256) {
		    *cpt++ = 'C';
		    *cpt++ = n;
		} else if (n < 65536) {
		    *cpt++ = 'S';
		    STORE_UINT16(cpt, n);
		} else {
		    *cpt++ = 'I';
		    STORE_UINT32(cpt, n);
		}
	    } else {
		if (n >= -128 && n < 128) {
		    *cpt++ = 'c';
		    *cpt++ = n;
		} else if (n >= -32768 && n < 32768) {
		    *cpt++ = 's';
		    STORE_UINT16(cpt, n);
		} else {
		    *cpt++ = 'i';
		    STORE_UINT32(cpt, n);
		}
	    }
	    break;

	case 'f': {
	    union {
		float f;
		int i;
	    } u;
	    u.f = atof((char *)value);
	    *cpt++ = 'f';
	    STORE_UINT32(cpt, u.i);
	    break;
	}

	case 'Z':
	    *cpt++ = 'Z';
	    while (value != cpf)
		*cpt++=*value++;
	    *cpt++ = 0;
	    break;

	case 'H':
	    *cpt++ = 'H';
	    while (value != cpf)
		*cpt++=*value++;
	    *cpt++ = 0;
	    break;

	case 'B': {
	    char subtype = *value++;
	    unsigned char *sz;
	    int count = 0;

	    *cpt++ = 'B';
	    *cpt++ = subtype;
	    sz = cpt; cpt += 4; /* Fill out later */

	    while (*value == ',') {
		value++;
		switch (subtype) {
		case 'c': case 'C':
		    *cpt++ = strtol((char *)value, (char **)&value, 10);
		    break;

		case 's': case 'S':
		    n = strtol((char *)value, (char **)&value, 10);
		    STORE_UINT16(cpt, n);
		    break;
		    
		case 'i': case 'I':
		    n = strtoll((char *)value, (char **)&value, 10);
		    STORE_UINT32(cpt, n);
		    break;
		    
		case 'f': {
		    union {
			float f;
			int i;
		    } u;

		    u.f = strtod((char *)value, (char **)&value);
		    STORE_UINT32(cpt, u.i);
		    break;
		}
		}
		count++;
	    }
	    if (value != cpf) {
		fprintf(stderr, "Malformed %c%c:B:... auxiliary field\n",
			key[0], key[1]);
		value = cpf;
	    }
	    STORE_UINT32(sz, count);
	    break;
	}

	default:
	    cpt -= 2;
	    fprintf(stderr, "Unknown aux format code '%c'\n", key[3]);
	    break;
	}

	if (*cpf == '\t')
	    cpf++;
    }

    if (cpt == end) { /* Hopefully very unlikely */
        size_t used = cpt - (unsigned char *)(*bsp);
        bam_seq_t *new_bam = realloc(*bsp, used + 1);
        if (NULL == new_bam) return -1;
        *bsp = new_bam;
        (*bsp)->alloc += used + 1;
        cpt = (unsigned char *)(*bsp) + used;
    }

    *cpt = 0;
    (*bsp)->blk_size = cpt - (unsigned char *)&(*bsp)->ref;
    return 0;
}

/*! Add preformated raw aux data to the bam_seq.
 *
 * Consider using bam_aux_add instead if you have information in a more
 * integer or string form.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
int bam_aux_add_data(bam_seq_t **b, const char tag[2], char type,
		     size_t len, const uint8_t *data) {
    uint8_t *cp;
    size_t used;

    if (NULL == b || NULL == data) return -1;

    /* Find the end of the existing tags and ensure there is enough space */
    cp = (uint8_t *)&(*b)->ref + (*b)->blk_size;
    used = cp - (uint8_t *)(*b);

    if ((*b)->alloc < used + len + 4) {
	size_t required = used + len + 4;
	bam_seq_t *new_bam = realloc((*b), required);
	if (NULL == new_bam) return -1;
	*b = new_bam;
	(*b)->alloc = required;
	cp = (uint8_t *) new_bam + used;
    }

    *cp++ = tag[0];
    *cp++ = tag[1];
    *cp++ = type;
    memcpy(cp, data, len);
    cp += len;
    *cp = 0;

    (*b)->blk_size = (uint32_t)(cp - (uint8_t *)&(*b)->ref);

    return 0;
}

/*! Add raw data to a bam structure.
 *
 * This could be useful if you wish to manually construct your own bam
 * entries or if you need to append an entire block of preformatting
 * aux data.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
int bam_add_raw(bam_seq_t **b, size_t len, const uint8_t *data) {
    uint8_t *cp;
    size_t used;

    if (NULL == b || NULL == data) return -1;

    /* Find the end of the existing tags and ensure there is enough space */
    cp = (uint8_t *)&(*b)->ref + (*b)->blk_size;
    used = cp - (uint8_t *)(*b);

    if ((*b)->alloc < used + len + 1) {
	size_t required = used + len + 1;
	bam_seq_t *new_bam = realloc((*b), required);
	if (NULL == new_bam) return -1;
	*b = new_bam;
	(*b)->alloc = required;
	cp = (uint8_t *) new_bam + used;
    }

    memcpy(cp, data, len);
    cp += len;
    *cp = 0;

    (*b)->blk_size = (uint32_t)(cp - (uint8_t *)&(*b)->ref);

    return 0;
}


/*! Duplicates a bam_seq_t structure.
 *
 * @return
 * Returns the new bam_seq_t pointer on success;
 *         NULL on failure.
 */
bam_seq_t *bam_dup(bam_seq_t *b) {
    bam_seq_t *d;
    int a = ((int)((b->alloc+15)/16))*16;

    if (!b)
	return NULL;

    if (!(d = malloc(a)))
	return NULL;

    memcpy(d, b, b->alloc);
    d->alloc = a;
    return d;
}


unsigned char *append_int(unsigned char *cp, int32_t i) {
    int32_t j;

    if (i < 0) {
	*cp++ = '-';
	if (i == INT_MIN) {
	    *cp++ = '2'; *cp++ = '1'; *cp++ = '4'; *cp++ = '7';
	    *cp++ = '4'; *cp++ = '8'; *cp++ = '3'; *cp++ = '6';
	    *cp++ = '4'; *cp++ = '8';
	    return cp;
	}

	i = -i;
    } else if (i == 0) {
	*cp++ = '0';
	return cp;
    }

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
 b7: if ((j = i / 10000000))   {*cp++ = j + '0'; i -= j*10000000;   goto x6;}
    if ((j = i / 1000000))    {*cp++ = j + '0', i -= j*1000000;    goto x5;}
 b5: if ((j = i / 100000))     {*cp++ = j + '0', i -= j*100000;     goto x4;}
    if ((j = i / 10000))      {*cp++ = j + '0', i -= j*10000;      goto x3;}
 b3: if ((j = i / 1000))       {*cp++ = j + '0', i -= j*1000;       goto x2;}
    if ((j = i / 100))        {*cp++ = j + '0', i -= j*100;        goto x1;}
 b1: if ((j = i / 10))         {*cp++ = j + '0', i -= j*10;         goto x0;}
    if (i)                     *cp++ = i + '0';
    return cp;

 x8: *cp++ = i / 100000000 + '0', i %= 100000000;
 x7: *cp++ = i / 10000000  + '0', i %= 10000000;
 x6: *cp++ = i / 1000000   + '0', i %= 1000000;
 x5: *cp++ = i / 100000    + '0', i %= 100000;
 x4: *cp++ = i / 10000     + '0', i %= 10000;
 x3: *cp++ = i / 1000      + '0', i %= 1000;
 x2: *cp++ = i / 100       + '0', i %= 100;
 x1: *cp++ = i / 10        + '0', i %= 10;
 x0: *cp++ = i             + '0';

    return cp;
}

/*
 * Unsigned version of above.
 * Only differs when the int has the top bit set (~2.15 billion and above).
 */
unsigned char *append_uint(unsigned char *cp, uint32_t i) {
    uint32_t j;

    if (i == 0) {
	*cp++ = '0';
	return cp;
    }

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
 b7: if ((j = i / 10000000))   {*cp++ = j + '0'; i -= j*10000000;   goto x6;}
    if ((j = i / 1000000))    {*cp++ = j + '0', i -= j*1000000;    goto x5;}
 b5: if ((j = i / 100000))     {*cp++ = j + '0', i -= j*100000;     goto x4;}
    if ((j = i / 10000))      {*cp++ = j + '0', i -= j*10000;      goto x3;}
 b3: if ((j = i / 1000))       {*cp++ = j + '0', i -= j*1000;       goto x2;}
    if ((j = i / 100))        {*cp++ = j + '0', i -= j*100;        goto x1;}
 b1: if ((j = i / 10))         {*cp++ = j + '0', i -= j*10;         goto x0;}
    if (i)                     *cp++ = i + '0';
    return cp;

 x8: *cp++ = i / 100000000 + '0', i %= 100000000;
 x7: *cp++ = i / 10000000  + '0', i %= 10000000;
 x6: *cp++ = i / 1000000   + '0', i %= 1000000;
 x5: *cp++ = i / 100000    + '0', i %= 100000;
 x4: *cp++ = i / 10000     + '0', i %= 10000;
 x3: *cp++ = i / 1000      + '0', i %= 1000;
 x2: *cp++ = i / 100       + '0', i %= 100;
 x1: *cp++ = i / 10        + '0', i %= 10;
 x0: *cp++ = i             + '0';

    return cp;
}

unsigned char *append_int64(unsigned char *cp, int64_t i) {
    int64_t j;

    if (i < 0) {
	*cp++ = '-';
	if (i == INT_MIN) {
	    *cp++ = '2'; *cp++ = '1'; *cp++ = '4'; *cp++ = '7';
	    *cp++ = '4'; *cp++ = '8'; *cp++ = '3'; *cp++ = '6';
	    *cp++ = '4'; *cp++ = '8';
	    return cp;
	}

	i = -i;
    } else if (i == 0) {
	*cp++ = '0';
	return cp;
    }

    //if (i < 10)         goto b0;
    if (i < 100)        goto b1;
    //if (i < 1000)       goto b2;
    if (i < 10000)      goto b3;
    //if (i < 100000)     goto b4;
    if (i < 1000000)    goto b5;
    //if (i < 10000000)   goto b6;
    if (i < 100000000)  goto b7;

    if ((j = i / 1000000000000000000))  {*cp++ = j + '0'; i -= j*1000000000000000000; goto xh;}
    if ((j = i / 100000000000000000))   {*cp++ = j + '0'; i -= j*100000000000000000; goto xg;}
    if ((j = i / 10000000000000000))    {*cp++ = j + '0'; i -= j*10000000000000000; goto xf;}
    if ((j = i / 1000000000000000))     {*cp++ = j + '0'; i -= j*1000000000000000; goto xe;}
    if ((j = i / 100000000000000))      {*cp++ = j + '0'; i -= j*100000000000000; goto xd;}
    if ((j = i / 10000000000000))       {*cp++ = j + '0'; i -= j*10000000000000; goto xc;}
    if ((j = i / 1000000000000))        {*cp++ = j + '0'; i -= j*1000000000000; goto xb;}
    if ((j = i / 100000000000))         {*cp++ = j + '0'; i -= j*100000000000; goto xa;}
    if ((j = i / 10000000000))          {*cp++ = j + '0'; i -= j*10000000000; goto x9;}
    if ((j = i / 1000000000)) {*cp++ = j + '0'; i -= j*1000000000; goto x8;}
    if ((j = i / 100000000))  {*cp++ = j + '0'; i -= j*100000000;  goto x7;}
 b7: if ((j = i / 10000000))   {*cp++ = j + '0'; i -= j*10000000;   goto x6;}
    if ((j = i / 1000000))    {*cp++ = j + '0', i -= j*1000000;    goto x5;}
 b5: if ((j = i / 100000))     {*cp++ = j + '0', i -= j*100000;     goto x4;}
    if ((j = i / 10000))      {*cp++ = j + '0', i -= j*10000;      goto x3;}
 b3: if ((j = i / 1000))       {*cp++ = j + '0', i -= j*1000;       goto x2;}
    if ((j = i / 100))        {*cp++ = j + '0', i -= j*100;        goto x1;}
 b1: if ((j = i / 10))         {*cp++ = j + '0', i -= j*10;         goto x0;}
    if (i)                     *cp++ = i + '0';
    return cp;

 xh: *cp++ = i / 100000000000000000  + '0', i %= 100000000000000000;
 xg: *cp++ = i / 10000000000000000   + '0', i %= 10000000000000000;
 xf: *cp++ = i / 1000000000000000    + '0', i %= 1000000000000000;
 xe: *cp++ = i / 100000000000000     + '0', i %= 100000000000000;
 xd: *cp++ = i / 10000000000000	     + '0', i %= 10000000000000;
 xc: *cp++ = i / 1000000000000       + '0', i %= 1000000000000;
 xb: *cp++ = i / 100000000000  	     + '0', i %= 100000000000;
 xa: *cp++ = i / 10000000000         + '0', i %= 10000000000;
 x9: *cp++ = i / 1000000000          + '0', i %= 1000000000;
 x8: *cp++ = i / 100000000           + '0', i %= 100000000;
 x7: *cp++ = i / 10000000            + '0', i %= 10000000;
 x6: *cp++ = i / 1000000             + '0', i %= 1000000;
 x5: *cp++ = i / 100000              + '0', i %= 100000;
 x4: *cp++ = i / 10000               + '0', i %= 10000;
 x3: *cp++ = i / 1000                + '0', i %= 1000;
 x2: *cp++ = i / 100                 + '0', i %= 100;
 x1: *cp++ = i / 10                  + '0', i %= 10;
 x0: *cp++ = i                       + '0';

    return cp;
}

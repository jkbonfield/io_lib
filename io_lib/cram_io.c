/*
 * Low level I/O primitives.
 *
 * - ITF8 encoding and decoding.
 * - Block based I/O
 * - Zlib inflating and deflating (memory)
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
#include "io_lib/deflate_interlaced.h" /* For block_create() */

/*
 * Reads an integer in ITF-8 encoding from 'cp' and stores it in
 * *val.
 *
 * Returns the number of bytes read on success
 *        -1 on failure
 */
int itf8_decode(cram_fd *fd, int32_t *val_p) {
    static int nbytes[16] = {
	0,0,0,0, 0,0,0,0,                               // 0000xxxx - 0111xxxx
	1,1,1,1,                                        // 1000xxxx - 1011xxxx
	2,2,                                            // 1100xxxx - 1101xxxx
	3,                                              // 1110xxxx
	4,                                              // 1111xxxx
    };

    static int nbits[16] = {
	0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, // 0000xxxx - 0111xxxx
	0x3f, 0x3f, 0x3f, 0x3f,                         // 1000xxxx - 1011xxxx
	0x1f, 0x1f,                                     // 1100xxxx - 1101xxxx
	0x0f,                                           // 1110xxxx
	0x0f,                                           // 1111xxxx
    };

    int32_t val = getc(fd->fp);
    if (val == -1)
	return -1;

    int i = nbytes[val>>4];
    val &= nbits[val>>4];

    switch(i) {
    case 0:
	*val_p = val;
	return 1;

    case 1:
	val = (val<<8) | (unsigned char)getc(fd->fp);
	*val_p = val;
	return 2;

    case 2:
	val = (val<<8) | (unsigned char)getc(fd->fp);
	val = (val<<8) | (unsigned char)getc(fd->fp);
	*val_p = val;
	return 3;

    case 3:
	val = (val<<8) | (unsigned char)getc(fd->fp);
	val = (val<<8) | (unsigned char)getc(fd->fp);
	val = (val<<8) | (unsigned char)getc(fd->fp);
	*val_p = val;
	return 4;

    case 4: // really 3.5 more, why make it different?
	val = (val<<8) | (unsigned char)getc(fd->fp);
	val = (val<<8) | (unsigned char)getc(fd->fp);
	val = (val<<8) | (unsigned char)getc(fd->fp);
	val = (val<<4) | (((unsigned char)getc(fd->fp)) & 0x0f);
	*val_p = val;
    }

    return 5;
}

/*
 * Encodes and writes a single integer in ITF-8 format.
 * Returns 0 on success
 *        -1 on failure
 */
int itf8_encode(cram_fd *fd, int32_t val) {
    char buf[5];
    int len = itf8_put(buf, val);
    return fwrite(buf, 1, len, fd->fp) == len ? 0 : -1;
}

#ifndef ITF8_MACROS
/*
 * As above, but decoding from memory
 */
int itf8_get(char *cp, int32_t *val_p) {
    unsigned char *up = (unsigned char *)cp;
    
    if (up[0] < 0x80) {
	*val_p =   up[0];
	return 1;
    } else if (up[0] < 0xc0) {
	*val_p = ((up[0] <<8) |  up[1])                           & 0x3fff;
	return 2;
    } else if (up[0] < 0xe0) {
	*val_p = ((up[0]<<16) | (up[1]<< 8) |  up[2])             & 0x1fffff;
	return 3;
    } else if (up[0] < 0xf0) {
	*val_p = ((up[0]<<24) | (up[1]<<16) | (up[2]<<8) | up[3]) & 0x0fffffff;
	return 4;
    } else {
	*val_p = ((up[0] & 0x0f)<<28) | (up[1]<<20) | (up[2]<<12) | (up[3]<<4) | (up[4] & 0x0f);
	return 5;
    }
}

/*
 * Stores a value to memory in ITF-8 format.
 *
 * Returns the number of bytes required to store the number.
 * This is a maximum of 5 bytes.
 */
int itf8_put(char *cp, int32_t val) {
    if        (!(val & ~0x00000007f)) { // 1 byte
	*cp = val;
	return 1;
    } else if (!(val & ~0x00003fff)) { // 2 byte
	*cp++ = (val >> 8 ) | 0x80;
	*cp   = val & 0xff;
	return 2;
    } else if (!(val & ~0x01fffff)) { // 3 byte
	*cp++ = (val >> 16) | 0xc0;
	*cp++ = (val >> 8 ) & 0xff;
	*cp   = val & 0xff;
	return 3;
    } else if (!(val & ~0x0fffffff)) { // 4 byte
	*cp++ = (val >> 24) | 0xe0;
	*cp++ = (val >> 16) & 0xff;
	*cp++ = (val >> 8 ) & 0xff;
	*cp   = val & 0xff;
	return 4;
    } else {                           // 5 byte
	*cp++ = 0xf0 | ((val>>28) & 0xff);
	*cp++ = (val >> 20) & 0xff;
	*cp++ = (val >> 12) & 0xff;
	*cp++ = (val >> 4 ) & 0xff;
	*cp = val & 0x0f;
	return 5;
    }
}
#endif

/* ----------------------------------------------------------------------- */

/* zlib compression code - from Gap5's tg_iface_g.c */
static char *zlib_mem_inflate(char *cdata, size_t csize, size_t *size) {
    z_stream s;
    unsigned char *data = NULL; /* Uncompressed output */
    int data_alloc = 0;
    int err;

    /* Starting point at uncompressed size, 4x compressed */
    data = malloc(data_alloc = csize*4+10);

    /* Initialise zlib stream */
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)cdata;
    s.avail_in = csize;
    s.total_in = 0;
    s.next_out  = data;
    s.avail_out = data_alloc;
    s.total_out = 0;

    //err = inflateInit(&s);
    err = inflateInit2(&s, 15 + 32);
    if (err != Z_OK) {
	fprintf(stderr, "zlib inflateInit error: %s\n", s.msg);
	return NULL;
    }

    /* Decode to 'data' array */
    for (;s.avail_in;) {
	s.next_out = &data[s.total_out];
	err = inflate(&s, Z_NO_FLUSH);
	if (err == Z_STREAM_END)
	    break;

	if (err != Z_OK) {
	    fprintf(stderr, "zlib inflate error: %s\n", s.msg);
	    break;
	}

	/* More to come, so realloc */
	data = realloc(data, data_alloc += s.avail_in*4 + 10);
	s.avail_out += s.avail_in*4+10;
    }
    inflateEnd(&s);

    *size = s.total_out;
    return (char *)data;
}

static char *zlib_mem_deflate(char *data, size_t size, size_t *cdata_size,
			      int level, int strat) {
    z_stream s;
    unsigned char *cdata = NULL; /* Compressed output */
    int cdata_alloc = 0;
    int cdata_pos = 0;
    int err;

    cdata = malloc(cdata_alloc = size*1.05+100);
    cdata_pos = 0;

    /* Initialise zlib stream */
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)data;
    s.avail_in = size;
    s.total_in = 0;
    s.next_out  = cdata;
    s.avail_out = cdata_alloc;
    s.total_out = 0;
    s.data_type = Z_BINARY;

    err = deflateInit2(&s, level, Z_DEFLATED, 15|16, 9, strat);
    if (err != Z_OK) {
	fprintf(stderr, "zlib deflateInit2 error: %s\n", s.msg);
	return NULL;
    }

    /* Encode to 'cdata' array */
    for (;s.avail_in;) {
	s.next_out = &cdata[cdata_pos];
	s.avail_out = cdata_alloc - cdata_pos;
	if (cdata_alloc - cdata_pos <= 0) {
	    fprintf(stderr, "Deflate produced larger output than expected. Abort\n"); 
	    return NULL;
	}
	err = deflate(&s, Z_NO_FLUSH);
	cdata_pos = cdata_alloc - s.avail_out;
	if (err != Z_OK) {
	    fprintf(stderr, "zlib deflate error: %s\n", s.msg);
	    break;
	}
    }
    if (deflate(&s, Z_FINISH) != Z_STREAM_END) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    *cdata_size = s.total_out;

    if (deflateEnd(&s) != Z_OK) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    return (char *)cdata;
}

/* ----------------------------------------------------------------------- */

void cram_uncompress_block(cram_block *b) {
    char *uncomp;
    size_t uncomp_size = 0;

    switch (b->method) {
    case RAW:
	b->uncomp_size = b->comp_size;
	return;

    case GZIP:
	uncomp = zlib_mem_inflate((char *)b->data, b->comp_size, &uncomp_size);
	assert((int)uncomp_size == b->uncomp_size);
	free(b->data);
	b->data = (unsigned char *)uncomp;
	b->method = RAW;
	break;

    case BZIP2:
	fprintf(stderr, "Bzip2 compression not yet implemented\n");
	abort();
	break;
    }
}

/*
 * Compresses a block using one of two different zlib strategies. If we only
 * want one choice set strat2 to be -1.
 *
 * The logic here is that sometimes Z_RLE does a better job than Z_FILTERED
 * or Z_DEFAULT_STRATEGY on quality data. If so, we'd rather use it as it is
 * significantly faster.
 */
void cram_compress_block(cram_fd *fd, cram_block *b, cram_metrics *metrics,
			 int level,  int strat,
			 int level2, int strat2) {
    char *comp;
    size_t comp_size = 0;

    if (level == 0) {
	b->method = RAW;
	b->comp_size = b->uncomp_size;
	return;
    }

    if (b->method != RAW) {
	fprintf(stderr, "Attempt to compress an already compressed block.\n");
	return;
    }

    if (strat2 >= 0)
	if (fd->verbose > 1)
	    fprintf(stderr, "metrics trial %d, next_trial %d, m1 %d, m2 %d\n",
		    metrics->trial, metrics->next_trial,
		    metrics->m1, metrics->m2);

    if (strat2 >= 0 && (metrics->trial || --metrics->next_trial == 0)) {
	char *c1, *c2;
	size_t s1, s2;

	if (metrics->next_trial == 0) {
	    metrics->next_trial = 100;
	    metrics->trial = 2;
	    metrics->m1 = metrics->m2 = 0;
	} else {
	    metrics->trial--;
	}

	c1 = zlib_mem_deflate((char *)b->data, b->uncomp_size,
			      &s1, level, strat);
	c2 = zlib_mem_deflate((char *)b->data, b->uncomp_size,
			      &s2, level2, strat2);
	if (s1 < s2) {
	    if (fd->verbose > 1)
		fprintf(stderr, "M1 wins\n");
	    comp = c1; comp_size = s1;
	    free(c2);
	    metrics->m1++;
	} else {
	    if (fd->verbose > 1)
		fprintf(stderr, "M2 wins\n");
	    comp = c2; comp_size = s2;
	    free(c1);
	    metrics->m2++;
	}
    } else if (strat2 >= 0) {
	comp = zlib_mem_deflate((char *)b->data, b->uncomp_size, &comp_size,
				metrics->m1 > metrics->m2 ? level : level2,
				metrics->m1 > metrics->m2 ? strat : strat2);
    } else {
	comp = zlib_mem_deflate((char *)b->data, b->uncomp_size, &comp_size,
				level, strat);
    }

    free(b->data);
    b->data = (unsigned char *)comp;
    b->method = GZIP;
    b->comp_size = comp_size;

    if (fd->verbose)
	fprintf(stderr, "Compressed block ID %d from %d to %d\n",
		b->content_id, b->uncomp_size, b->comp_size);
}

cram_metrics *cram_new_metrics(void) {
    cram_metrics *m = malloc(sizeof(*m));
    m->m1 = m->m2 = 0;
    m->trial = 2;
    m->next_trial = 100;
    return m;
}


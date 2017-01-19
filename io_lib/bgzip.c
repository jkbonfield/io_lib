/*
 * Copyright (c) 2016 Genome Research Ltd.
 * Author(s): James Bonfield, Rob Davies.
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
 * Minimal bgzip .gzi index support for gzipped references.
 *
 * The .gzi format is undocumented, but used extensively!
 * It consists of a series of 64-bit little endian integers starting
 * with N = number_of_pairs and N (compressed_offset,
 * uncompressed_offset) pairs.
 *
 * The user is expected to do (eg) a binary search to convert
 * uncompressed offsets to compressed offsets, and then start reading
 * from that point onwards.
 */

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include <zlib.h>
#include "io_lib/bgzip.h"
#include "io_lib/os.h"

/* ----------------------------------------------------------------------
 * bgzip .gzi index support
 */

/* Loads an bgzip index and returns it.
 * Returns NULL on failure.
 */
gzi *gzi_index_load(const char *fn) {
    gzi *idx = malloc(sizeof(*idx));
    FILE *fp;

    if (strlen(fn) >= 4 && strcmp(fn+strlen(fn)-4, ".gzi") == 0) {
	// We were given the .gzi filename itself
	fp = fopen(fn, "rb");
	if (!fp) perror(fn);
    } else {
	// Append .gzi suffix and hope it exists
	char fn2[8192];
	snprintf(fn2, 8192, "%s.gzi", fn);

	fp = fopen(fn2, "rb");
    }

    if (!fp)
	goto err;

    uint64_t n, i;

    if (8 != fread(&n, 1, 8, fp))
	goto err;
    n = le_int8(n);

    if (n < 0 || n >= INT_MAX/8 - 1)
	goto err;

    idx->n = n;
    idx->c_off = malloc(8*n+8);
    idx->u_off = malloc(8*n+8);
    if (!idx->c_off || !idx->u_off)
	goto err;

    idx->c_off[0] = idx->u_off[0] = 0;
    for (i = 1; i <= n; i++) {
	if (8 != fread(&idx->c_off[i], 1, 8, fp) || 
	    8 != fread(&idx->u_off[i], 1, 8, fp))
	    goto err;
	idx->c_off[i] = le_int8(idx->c_off[i]);
	idx->u_off[i] = le_int8(idx->u_off[i]);
    }

    return idx;

 err:
    if (fp)
	fclose(fp);

    if (idx)
	free(idx);

    return NULL;
}

void gzi_index_free(gzi *idx) {
    if (idx) {
	free(idx->c_off);
	free(idx->u_off);
	free(idx);
    }
}

gzi *gzi_index_init() {
    gzi *idx = calloc(1, sizeof(*idx));
    return idx;
}

/*
 * Adds a compressed / uncompressed map to the gzi structure.
 * Returns 0 on success;
 *        -1 on failure.
 */
int gzi_index_add_block(gzi *idx, uint64_t c_off, uint64_t u_off) {
    uint64_t n;
    idx->n++;
    n = idx->n;
    idx->c_off = realloc(idx->c_off, 8*n+8);
    idx->u_off = realloc(idx->u_off, 8*n+8);
    if (!idx->c_off || !idx->u_off)
	return -1;
    if (n == 1){
	idx->c_off[n-1] = c_off;
	idx->u_off[n-1] = u_off;
    } else {
	idx->c_off[n-1] = c_off+idx->c_off[n-2];
	idx->u_off[n-1] = u_off+idx->u_off[n-2];
    }
    return 0;
}


/*
 * Writes a gzi file to a file with basename 'bname' and optional suffix
 * (or NULL if unused).
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int gzi_index_dump(gzi *idx, const char *bname, const char *suffix) {
    char *tmp = (char *)bname;
    if (!idx)
        return -1;

    if (suffix) {
        int blen = strlen(bname);
        int slen = strlen(suffix);
        if (!(tmp = (char*) malloc(blen + slen + 1)))
	    return -1;
        memcpy(tmp, bname, blen);
        memcpy(tmp+blen, suffix, slen+1);
    }

    FILE *idx_f = fopen(tmp, "wb");
    if (!idx_f) {
	perror(tmp);
	if (tmp != bname)
	    free(tmp);
	return -1;
    }
    if (tmp != bname)
	free(tmp);

    int i;
    uint64_t n = idx->n;
    if (fwrite(le_int8(&n), sizeof(n), 1, idx_f) < 0)
	goto fail;
    for (i=0; i<idx->n; i++) {
	if (fwrite(le_int8(&idx->c_off[i]), sizeof idx->c_off[i], 1, idx_f) < 0)
	    goto fail;
	if (fwrite(le_int8(&idx->u_off[i]), sizeof idx->u_off[i], 1, idx_f) < 0)
	    goto fail;
    }

    if (fclose(idx_f) < 0)
        return -1;

    return 0;

 fail:
    fclose(idx_f);
    return -1;
}

/*
 * Uncompressed offset to virtual offset.
 * A virtual offset is a compressed offset << 16 ORed with
 * the uncompressed relative offset since the statr of that
 * compressed offset.
 *
 * Eg 12345 uncompressed offset may map to a block at
 * 12000 uncompressed, 9000 compressed, which then becomes
 * (9000<<16)|345 virtual offset.
 *
 * *sz is returned as the size of the compressed block containig
 * uoff, or 0 if unknown (determine from EOF instead).
 */
static int64_t gzi_uoff_to_voff(gzi *idx, uint64_t uoff, int *sz) {
    /* Binary search */
    int lo = 0, hi = idx->n, x;

    while (hi - lo > 1) {
	x = (hi + lo)/2;
	if (idx->u_off[x] > uoff)
	    hi = x;
	else
	    lo = x;
    }

    x = (hi > lo && idx->u_off[hi] > uoff) ? lo : hi;

    if (uoff - idx->u_off[x] >= 65536)
	return -1;

    if (sz) {
	if (x < idx->n)
	    *sz = idx->c_off[x+1] - idx->c_off[x];
	
	else
	    *sz = 0;
    }

    return (idx->c_off[x]<<16) | (uoff - idx->u_off[x]);
}

uint64_t gzi_load(FILE *fp, gzi *idx, uint64_t ustart, uint64_t uend, char *out) {
    int csz = 0, err;
    int64_t vstart = gzi_uoff_to_voff(idx, ustart, 0);
    int64_t vend   = gzi_uoff_to_voff(idx, uend, &csz);

    off_t cstart = vstart >> 16;
    off_t cend   = vend   >> 16;

    uint64_t out_sz = 0;

    if (!csz) {
	// go to EOF to find size of last blockx
	fseeko(fp, 0, SEEK_END);
	csz = ftello(fp) - cstart;
    } else {
	csz += cend - cstart;
    }


    // Load the compressed blocks
    char *comp = malloc(csz);
    if (!comp)
	return 0;
    
    if (fseeko(fp, cstart, SEEK_SET) < 0)
	return 0;

    if (csz != fread(comp, 1, csz, fp))
	return 0;


    z_stream z;
    z.zalloc = 0;
    z.zfree = 0;
    if (inflateInit2(&z, 31) != Z_OK) {
	fprintf(stderr, "Zlib err: %s\n", z.msg);
	free(comp);
	return 0;
    }

    z.next_in = (unsigned char *)comp;
    z.avail_in = csz;

    // Discard initial portion
    unsigned char buf[65536];
    z.next_out = buf;
    z.avail_out = vstart & 0xffff;
    if (z.avail_out) {
	int err = inflate(&z, Z_FINISH);
	if (err != Z_OK && err != Z_BUF_ERROR) {
	    fprintf(stderr, "Zlib err: %s\n", z.msg);
	    free(comp);
	    return 0;
	}
    }

    // Decode remainder, in a loop as we have concatenated zib streams.
    z.total_out = 0;
    z.next_out = (unsigned char *)out;
    z.avail_out = uend-ustart+1;

    do {
	err = inflate(&z, Z_FINISH);
	out_sz += z.total_out;
	if (err == Z_STREAM_END && z.avail_out && z.avail_in)
	    inflateReset(&z);
    } while ((err == Z_STREAM_END || err == Z_OK) && z.avail_out != 0);

    inflateEnd(&z);
    free(comp);
    return (err == Z_STREAM_END || err == Z_OK || err == Z_BUF_ERROR) ? out_sz : 0;
}


/* ----------------------------------------------------------------------
 * A FILE* wrapper that can read and seek either into uncompressed or
 * bgzip compressed files.
 *
 * Note, this is crude and not at all good at handling small reads efficiently
 * due to no cachine and pointless seeks!  It got bolted on without
 * the necessary redesigns.
 */
typedef struct bzi_FILE {
    FILE *fp;
    gzi  *idx;
    uint64_t pos;
} bzi_FILE;

void bzi_close(bzi_FILE *zp) {
    if (!zp)
	return;

    if (zp->fp) fclose(zp->fp);
    gzi_index_free(zp->idx);
    free(zp);
}

bzi_FILE *bzi_open(const char *path, const char *mode) {
    if (*mode != 'r')
	return NULL;

    bzi_FILE *zp = calloc(1, sizeof(*zp));
    if (!zp) goto err;
    if (!(zp->fp = fopen(path, mode))) goto err;

    // Try loading the index, but assume failure means it's a normal file.
    zp->idx = gzi_index_load(path);

    return zp;

 err:
    bzi_close(zp);
    return NULL;
}

// NOTE: every read is new seek + load.  Not intended for use on 
// lots of small reads.
size_t bzi_read(void *ptr, size_t size, size_t nmemb, bzi_FILE *zp) {
    if (!zp->idx) {
	return fread(ptr, size, nmemb, zp->fp);
    } else {
	uint64_t n = gzi_load(zp->fp, zp->idx,
			      zp->pos, zp->pos + size*nmemb -1, ptr);
	zp->pos += n;
	return n;
    }
}

int bzi_seek(bzi_FILE *zp, off_t offset, int whence) {
    if (!zp->idx) {
	return fseeko(zp->fp, offset, whence);
    } else {
	switch (whence) {
	case SEEK_SET:
	    zp->pos = offset;
	    break;
	    
	case SEEK_CUR:
	    zp->pos += offset;

	default:
	    // SEEK_END not supported
	    return -1;
	}

	return 0;
    }
}

/* ----------------------------------------------------------------------
 */

#ifdef TEST_MAIN
int main(int argc, char **argv) {
    if (argc != 4) {
	fprintf(stderr, "Usage: %s input.gz start end\n", argv[0]);
	return 1;
    }

    gzi *idx = gzi_index_load(argv[1]);
    uint64_t ustart = atoll(argv[2]), uend = atoll(argv[3]);

    if (!idx) {
	fprintf(stderr, "Unable to open index: %s\n", argv[1]);
	return 1;
    }

    FILE *fp = fopen(argv[1], "rb");
    if (!fp) {
	perror(argv[1]);
	return 1;
    }
    char *buf = malloc(uend - ustart + 1);
    if (!buf)
	return 1;

    uint64_t sz = gzi_load(fp, idx, ustart, uend, buf);
    if (sz != write(1, buf, sz))
	return 1;

    free(buf);
    gzi_index_free(idx);

    return 0;
}
#endif

#ifdef TEST_MAIN2
int main(int argc, char **argv) {
    if (argc != 4) {
	fprintf(stderr, "Usage: %s input.gz start end\n", argv[0]);
	return 1;
    }

    bzi_FILE *zp = bzi_open(argv[1], "rb");
    uint64_t ustart = atoll(argv[2]), uend = atoll(argv[3]);

    if (!zp) {
	perror(argv[1]);
	return 1;
    }

    char *buf = malloc(uend - ustart + 1);
    if (!buf)
	return 1;

    bzi_seek(zp, ustart, SEEK_SET);
    uint64_t sz = bzi_read(buf, 1, uend-ustart+1, zp);
    if (sz != write(1, buf, sz))
	return 1;

    free(buf);
    bzi_close(zp);

    return 0;
}
#endif

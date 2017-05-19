/*
 * Copyright (c) 2013-2015 Genome Research Ltd.
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
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013-2015
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
#ifdef HAVE_LIBBZ2
#include <bzlib.h>
#endif
#ifdef HAVE_LIBLZMA
#include <lzma.h>
#endif

#include <sys/types.h>
#include <dirent.h>
#include <dlfcn.h>
#include <unistd.h>
#include <time.h>

#include "io_lib/cram_block_compression.h"
#include "io_lib/cram.h"
#include "io_lib/os.h"
#include "io_lib/rANS_static.h"
#include "io_lib/rANS_static4x16.h"
#include "io_lib/hash_table.h"
#include "io_lib/crc32.h"


static cram_compressor **codecs = NULL;
static int ncodecs = 0;

static HashTable *codec_hash = NULL;

#ifdef DEBUG_TIME
int64_t btm[1<<24] = {0};

void dump_decode_time(void) {
    int i;
    for (i = 0; i < 1<<24; i++) {
	if (!btm[i]) continue;
	if (i < 256) {
	    fprintf(stderr, "%6d: %f seconds\n", i, btm[i] / 1000000000.0);
	} else {
	    char key[3] = {i>>16, i>>8, i};
	    fprintf(stderr, "   %.3s: %f seconds\n", key, btm[i] / 1000000000.0);
	}
    }
}
#endif

/*
 * Uncompresses a CRAM block, if compressed.
 */
int cram_uncompress_block(cram_slice *s, cram_block *b) {
    char *uncomp;
    size_t uncomp_size = 0;
    uint32_t method;

    if (b->crc32_checked == 0) {
	uint32_t crc = iolib_crc32(b->crc_part, b->data ? b->data : (uc *)"", b->alloc);
	b->crc32_checked = 1;
	if (crc != b->crc32) {
	    fprintf(stderr, "Block CRC32 failure\n");
	    return -1;
	}
    }

    if (b->uncomp_size == 0) {
	// blank block
	b->method = RAW;
	return 0;
    }

#ifdef DEBUG_TIME
    struct timespec ts1, ts2;
    clock_gettime(CLOCK_REALTIME, &ts1);
#endif

    switch ((method = b->method)) {
    case RAW:
	//b->uncomp_size = b->comp_size;
	return 0;

    default: {
	HashItem *hi;

	// Match FourCC to codec[] array entry
	if (!(hi = HashTableSearch(codec_hash, (char *)&method, 4))) {
	    fprintf(stderr, "Unknown codec id %d\n", method);
	    return -1;
	}

	method = hi->data.i;
	if (method >= ncodecs)
	    return -1;
	if (!codecs[method])
	    return -1;
    }

	// Deliberate fall through
    case GZIP:
    case BZIP2:
    case LZMA:
    case RANS0:
    case RANS1:
    case RANS_PR0:
	uncomp_size = b->uncomp_size;
	uncomp = (char *)codecs[method]->uncompress_block(s, b->data,b->comp_size,
							  &uncomp_size);
	if (!uncomp)
	    return -1;
	if ((int)uncomp_size != b->uncomp_size) {
	    free(uncomp);
	    return -1;
	}
	free(b->data);
	b->data = (unsigned char *)uncomp;
	b->alloc = uncomp_size;
	b->method = RAW;
	break;
    }

#ifdef DEBUG_TIME
    clock_gettime(CLOCK_REALTIME, &ts2);
    btm[b->content_id] += (ts2.tv_sec - ts1.tv_sec) * 1000000000 + ts2.tv_nsec - ts1.tv_nsec;
#endif

    return 0;
}

#define EBASE 65536
//static double entropy16(unsigned short *data, int len) {
//    double E[EBASE];
//    double P[EBASE];
//    double e = 0;
//    int i;
//    
//    for (i = 0; i < EBASE; i++)
//        P[i] = 0;
//
//    for (i = 0; i < len; i++)
//        P[data[i]]++;
//
//    for (i = 0; i < EBASE; i++) {
//        if (P[i]) {
//            P[i] /= len;
//            E[i] = -(log(P[i])/log(EBASE));
//        } else {
//            E[i] = 0;
//        }
//    }
//
//    for (e = i = 0; i < len; i++)
//        e += E[data[i]];
//
//    return e * log(EBASE)/log(256);
//}
//
//#define EBASE2 256
//static double entropy8(unsigned char *data, int len) {
//    int F[EBASE2];
//    double e = 0;
//    int i;
//    
//    for (i = 0; i < EBASE2; i++)
//        F[i] = 0;
//
//    for (i = 0; i < len; i++)
//        F[data[i]]++;
//
//    for (i = 0; i < EBASE2; i++) {
//        if (F[i]) {
//	    e += -log((double)F[i]/len) * F[i];
//        }
//    }
//
//    return e / log(EBASE2);
//}

static char *cram_compress_by_method(char *in, size_t in_size,
				     size_t *out_size,
				     cram_slice *s,
				     enum cram_block_method method,
				     int level, int strat) {
    if (method == GZIP && strat == Z_RLE)
	method = GZIP_RLE;

    switch (method) {
    default:
	if (method >= ncodecs)
	    return NULL;
	if (!codecs[method])
	    return NULL;

	// Deliberate fall through

    case GZIP:
    case GZIP_RLE:
    case BZIP2:
    case LZMA:
    case RANS0:
    case RANS1:
    case RANS_PR0:
    case RANS_PR1:
    case RANS_PR64:
    case RANS_PR65:
    case RANS_PR128:
    case RANS_PR129:
    case RANS_PR192:
    case RANS_PR193:
	return (char *)codecs[method]->compress_block(method,
						      level, s,
						      (uint8_t *)in, in_size,
						      out_size);

    case RAW:
	break;
    }

    return NULL;
}

/*
 * Compresses a block using one of two different zlib strategies. If we only
 * want one choice set strat2 to be -1.
 *
 * The logic here is that sometimes Z_RLE does a better job than Z_FILTERED
 * or Z_DEFAULT_STRATEGY on quality data. If so, we'd rather use it as it is
 * significantly faster.
 */
int cram_compress_block(cram_fd *fd, cram_slice *s, cram_block *b, cram_metrics *metrics,
			int method, int level) {

    char *comp = NULL;
    size_t comp_size = 0;
    int strat;
    int m;

#ifdef DEBUG_TIME
    struct timespec ts1, ts2;
    clock_gettime(CLOCK_REALTIME, &ts1);
#endif

    if (b->method != RAW) {
	// Maybe already compressed if s->block[0] was compressed and
	// we have e.g. s->block[DS_BA] set to s->block[0] due to only
	// one base type present and hence using E_HUFFMAN on block 0.
	// A second explicit attempt to compress the same block then
	// occurs.
	return 0;
    }

    // Clear unavailable methods
    for (m = 0; m < 32; m++) {
	if (m >= ncodecs || !codecs[m])
	    method &= ~(1<<m);
    }

    //fprintf(stderr, "IN: block %d, sz %d\n", b->content_id, b->uncomp_size);

    if (method == RAW || level == 0 || b->uncomp_size == 0) {
	b->method = RAW;
	b->comp_size = b->uncomp_size;
	//fprintf(stderr, "Skip block id %d\n", b->content_id);
	return 0;
    }

    if (metrics) {
	if (fd->metrics_lock) pthread_mutex_lock(fd->metrics_lock);
	if (metrics->trial > 0 || --metrics->next_trial <= 0) {
	    int m;
	    size_t sz_best = INT_MAX;
	    size_t sz[CRAM_MAX_METHOD] = {0};
	    int method_best = 0;
	    char *c_best = NULL, *c = NULL;

	    if (metrics->revised_method)
		method = metrics->revised_method;
	    else
		metrics->revised_method = method;

	    if (metrics->next_trial == 0) {
		int i;
		metrics->next_trial = TRIAL_SPAN;
		metrics->trial = NTRIALS;
		for (i = 0; i < CRAM_MAX_METHOD; i++)
		    metrics->sz[i] /= 2;
	    }

	    if (fd->metrics_lock) pthread_mutex_unlock(fd->metrics_lock);

	    // Compress this block using the best method
	    for (m = 0; m < 32; m++) {
		if (method & (1<<m)) {
		    int strat;
		    if (m == GZIP)
			strat = Z_FILTERED;
		    else if (m == GZIP_RLE)
			strat = Z_RLE;
		    else
			strat = 0;

		    if (!codecs[m])
			continue;

		    if (codecs[m]->content_ids &&
			!(codecs[m]->content_ids & b->id_bits))
			continue;
		    c = cram_compress_by_method((char *)b->data, b->uncomp_size,
						&sz[m], s, m, level, strat);
		    if (c && sz_best > sz[m]) {
			sz_best = sz[m];
			method_best = m;
			if (c_best)
			    free(c_best);
			c_best = c;
		    } else if (c) {
			free(c);
		    } else {
			sz[m] = b->uncomp_size*2+1000; // ie. worse than raw
		    }
		    
		    //fprintf(stderr, "Method %d; Block %d; %d->%d\n", m, b->content_id, b->uncomp_size, sz[m]);
		}
	    }

	    //fprintf(stderr, "sz_best = %d\n", sz_best);

	    free(b->data);
	    b->data = (uint8_t *)c_best;
	    //printf("method_best = %s\n", cram_block_method2str(method_best));
	    //b->method = method_best == GZIP_RLE ? GZIP : method_best;
	    b->method = codecs[method_best]->code;
	    b->comp_size = sz_best;

	    // Accumulate stats for all methods tried
	    if (fd->metrics_lock) pthread_mutex_lock(fd->metrics_lock);
	    for (m = 0; m < CRAM_MAX_METHOD; m++)
		metrics->sz[m] += sz[m];

	    // When enough trials performed, find the best on average
	    if (--metrics->trial == 0) {
		int best_method = RAW;
		int best_sz = INT_MAX;

		// Scale methods by cost
		if (fd->level <= 3) {
		    for (m = 0; m < CRAM_MAX_METHOD && m < ncodecs; m++) {
			if (!codecs[m])
			    continue;
			metrics->sz[m] *= codecs[m]->cost;
		    }
		} else if (fd->level <= 6) {
		    for (m = 0; m < CRAM_MAX_METHOD && m < ncodecs; m++) {
			if (!codecs[m])
			    continue;
			metrics->sz[m] *= 1+(codecs[m]->cost-1)/2;
		    }
		} // else cost is ignored

		for (m = 0; m < CRAM_MAX_METHOD; m++) {
		    if (!metrics->sz[m])
			continue;

		    if (method & (1<<m) && best_sz > metrics->sz[m])
			best_sz = metrics->sz[m], best_method = m;
		}

		if (best_method == GZIP_RLE) {
		    metrics->strat  = Z_RLE;
		} else {
		    metrics->strat  = Z_FILTERED;
		}
		metrics->method = best_method;

		// If we see at least MAXFAIL trials in a row for a specific
		// compression method with more than MAXDELTA aggregate
		// size then we drop this from the list of methods used
		// for this block type.
#define MAXDELTA 0.20
#define MAXFAILS 4
		for (m = 0; m < CRAM_MAX_METHOD; m++) {
		    if (best_method == m) {
			metrics->cnt[m] = 0;
			metrics->extra[m] = 0;
		    } else if (best_sz < metrics->sz[m]) {
			double r = (double)metrics->sz[m] / best_sz - 1;
			if (++metrics->cnt[m] >= MAXFAILS && 
			    (metrics->extra[m] += r) >= MAXDELTA)
			    method &= ~(1<<m);
		    }
		}

		//if (method != metrics->revised_method)
		//    fprintf(stderr, "%d: method from %x to %x\n",
		//	    b->content_id, metrics->revised_method, method);
		metrics->revised_method = method;
	    }
	    if (fd->metrics_lock) pthread_mutex_unlock(fd->metrics_lock);
	} else {
	    strat = metrics->strat;
	    method = metrics->method;

	    if (fd->metrics_lock) pthread_mutex_unlock(fd->metrics_lock);
	    comp = cram_compress_by_method((char *)b->data, b->uncomp_size,
					   &comp_size, s, method,
					   level, strat);
	    if (!comp)
		return -1;
	    free(b->data);
	    b->data = (unsigned char *)comp;
	    b->comp_size = comp_size;
	    b->method = codecs[method]->code;
	}

    } else {
	// no cached metrics, so just do zlib?
	comp = cram_compress_by_method((char *)b->data, b->uncomp_size,
				       &comp_size, s, GZIP, level, Z_FILTERED);
	if (!comp) {
	    fprintf(stderr, "Compression failed!\n");
	    return -1;
	}
	free(b->data);
	b->data = (unsigned char *)comp;
	b->comp_size = comp_size;
	b->method = GZIP;
    }

    if (fd->verbose)
	fprintf(stderr, "Compressed block ID %d from %d to %d by method %s\n",
		b->content_id, b->uncomp_size, b->comp_size,
		cram_block_method2str(b->method));

    if (b->method == RANS1)
	b->method = RANS0; // Spec just has RANS (not 0/1) with auto-sensing

#ifdef DEBUG_TIME
    clock_gettime(CLOCK_REALTIME, &ts2);
    if (metrics)
	metrics->tm += (ts2.tv_sec - ts1.tv_sec) * 1000000000 + ts2.tv_nsec - ts1.tv_nsec;
#endif

    return 0;
}


/* -----------------------------------------------------------------------------
 * GZIP (default and RLE).
 */
/* ----------------------------------------------------------------------
 * zlib compression code - from Gap5's tg_iface_g.c
 * They're static here as they're only used within the cram_compress_block
 * and cram_uncompress_block functions, which are the external interface.
 */
static char *zlib_mem_inflate(char *cdata, size_t csize, size_t *size) {
    z_stream s;
    unsigned char *data = NULL; /* Uncompressed output */
    int data_alloc = 0;
    int err;

    /* Starting point at uncompressed size, and scale after that */
    data = malloc(data_alloc = csize*1.2+100);
    if (!data)
	return NULL;

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
	free(data);
	return NULL;
    }

    /* Decode to 'data' array */
    for (;s.avail_in;) {
	unsigned char *data_tmp;
	int alloc_inc;

	s.next_out = &data[s.total_out];
	err = inflate(&s, Z_NO_FLUSH);
	if (err == Z_STREAM_END)
	    break;

	if (err != Z_OK) {
	    fprintf(stderr, "zlib inflate error: %s\n", s.msg);
	    break;
	}

	/* More to come, so realloc based on growth so far */
	alloc_inc = (double)s.avail_in/s.total_in * s.total_out + 100;
	data = realloc((data_tmp = data), data_alloc += alloc_inc);
	if (!data) {
	    free(data_tmp);
	    return NULL;
	}
	s.avail_out += alloc_inc;
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
    if (!cdata)
	return NULL;
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

static const char *gzip_codec_name(void) {
    return "Gzip (default)";
}

static const char *gzip_rle_codec_name(void) {
    return "Gzip (RLE)";
}

unsigned char *gzip_codec_compress(int method,
				   int level,
				   cram_slice *s,
				   unsigned char *in,
				   size_t in_size,
				   size_t *out_size) {
    unsigned char *comp;

    comp = (uint8_t *)zlib_mem_deflate((char *)in, in_size,
				       out_size, level, Z_FILTERED);
    if (!comp)
        return NULL;

    return comp;
}

unsigned char *gzip_rle_codec_compress(int method,
				       int level,
				       cram_slice *s,
				       unsigned char *in,
				       size_t in_size,
				       size_t *out_size) {
    unsigned char *comp;

    comp = (uint8_t *)zlib_mem_deflate((char *)in, in_size,
				       out_size, 1, Z_RLE);
    if (!comp)
        return NULL;

    return comp;
}

unsigned char *gzip_codec_uncompress(cram_slice *s,
				     unsigned char *in,
				     size_t in_size,
				     size_t *out_size) {
    unsigned char *uncomp;

    if (!(uncomp = (uint8_t *)zlib_mem_inflate((char *)in, in_size, out_size)))
        return NULL;

    return uncomp;
}

static cram_compressor gzip_codec = {
    GZIP, //"GZIP",
    0,
    1.04,
    gzip_codec_name,
    gzip_codec_compress,
    gzip_codec_uncompress,
};

static cram_compressor gzip_rle_codec = {
    GZIP, //"GZIP",
    0,
    1.0,
    gzip_rle_codec_name,
    gzip_rle_codec_compress,
    gzip_codec_uncompress,
};


/* -----------------------------------------------------------------------------
 * Bzip2
 */
static const char *bzip2_codec_name(void) {
    return "Bzip2";
}

unsigned char *bzip2_codec_compress(int method,
				    int level,
				    cram_slice *s,
				    unsigned char *in,
				    size_t in_size,
				    size_t *out_size) {
    unsigned int comp_size = 1.01*in_size+600;
    unsigned char *comp;

    if (!(comp = (unsigned char *)malloc(comp_size)))
        return NULL;

    if (BZ_OK != BZ2_bzBuffToBuffCompress((char *)comp, &comp_size,
                                          (char *)in, in_size,
                                          level, 0, 30)) {
        free(comp);
        return NULL;
    }
    *out_size = comp_size;

    return comp;
}

unsigned char *bzip2_codec_uncompress(cram_slice *s,
				      unsigned char *in,
				      size_t in_size,
				      size_t *out_size) {
    char *uncomp;
    unsigned int out_sz = *out_size;

    if (!(uncomp = (char *)malloc(out_sz)))
        return NULL;

    if (BZ_OK != BZ2_bzBuffToBuffDecompress(uncomp, &out_sz,
                                            (char *)in, in_size,
                                            0, 0)) {
        free(uncomp);
        return NULL;
    }

    *out_size = out_sz;
    return (unsigned char *)uncomp;
}

static cram_compressor bzip2_codec = {
    BZIP2, // "BZIP",
    0,
    1.08,
    bzip2_codec_name,
    bzip2_codec_compress,
    bzip2_codec_uncompress,
};


/* -----------------------------------------------------------------------------
 * LZMA
 */
#ifdef HAVE_LIBLZMA
/* ------------------------------------------------------------------------ */
/*
 * Data compression routines using liblzma (xz)
 *
 * On a test set this shrunk the main db from 136157104 bytes to 114796168, but
 * caused tg_index to grow from 2m43.707s to 15m3.961s. Exporting as bfastq
 * went from 18.3s to 36.3s. So decompression suffers too, but not as bad
 * as compression times.
 *
 * For now we disable this functionality. If it's to be reenabled make sure you
 * improve the mem_inflate implementation as it's just a test hack at the
 * moment.
 */

static char *lzma_mem_deflate(char *data, size_t size, size_t *cdata_size,
			      int level) {
    char *out;
    size_t out_size = lzma_stream_buffer_bound(size);
    *cdata_size = 0;

    out = malloc(out_size);

    /* Single call compression */
    if (LZMA_OK != lzma_easy_buffer_encode(level, LZMA_CHECK_CRC32, NULL,
					   (uint8_t *)data, size,
					   (uint8_t *)out, cdata_size,
					   out_size))
    	return NULL;

    return out;
}

static char *lzma_mem_inflate(char *cdata, size_t csize, size_t *size) {
    lzma_stream strm = LZMA_STREAM_INIT;
    size_t out_size = 0, out_pos = 0;
    char *out = NULL;
    int r;

    /* Initiate the decoder */
    if (LZMA_OK != lzma_stream_decoder(&strm, lzma_easy_decoder_memusage(9), 0))
	return NULL;

    /* Decode loop */
    strm.avail_in = csize;
    strm.next_in = (uint8_t *)cdata;

    for (;strm.avail_in;) {
	if (strm.avail_in > out_size - out_pos) {
	    out_size += strm.avail_in * 4 + 32768;
	    out = realloc(out, out_size);
	}
	strm.avail_out = out_size - out_pos;
	strm.next_out = (uint8_t *)&out[out_pos];

	r = lzma_code(&strm, LZMA_RUN);
	if (LZMA_OK != r && LZMA_STREAM_END != r) {
	    //fprintf(stderr, "r=%d\n", r);
	    //fprintf(stderr, "mem=%"PRId64"d\n", (int64_t)lzma_memusage(&strm));
	    return NULL;
	}

	out_pos = strm.total_out;

	if (r == LZMA_STREAM_END)
	    break;
    }

    /* finish up any unflushed data; necessary? */
    r = lzma_code(&strm, LZMA_FINISH);
    if (r != LZMA_OK && r != LZMA_STREAM_END) {
	fprintf(stderr, "r=%d\n", r);
	return NULL;
    }

    out = realloc(out, strm.total_out);
    *size = strm.total_out;

    lzma_end(&strm);

    return out;
}

static const char *lzma_codec_name(void) {
    return "LZMA";
}

unsigned char *lzma_codec_compress(int method,
				   int level,
				   cram_slice *s,
				   unsigned char *in,
				   size_t in_size,
				   size_t *out_size) {
    return (uint8_t *)lzma_mem_deflate((char *)in, in_size, out_size, level);
}

unsigned char *lzma_codec_uncompress(cram_slice *s,
				     unsigned char *in,
				     size_t in_size,
				     size_t *out_size) {
    return (uint8_t *)lzma_mem_inflate((char *)in, in_size, out_size);
}

static cram_compressor lzma_codec = {
    LZMA, // "LZMA",
    0,
    1.10,
    lzma_codec_name,
    lzma_codec_compress,
    lzma_codec_uncompress,
};
#endif


/* -----------------------------------------------------------------------------
 * RANS0/RANS1
 */
static const char *rans0_codec_name(void) {
    return "RANS0";
}

static const char *rans1_codec_name(void) {
    return "RANS1";
}

unsigned char *rans0_codec_compress(int method,
				    int level,
				    cram_slice *s,
				    unsigned char *in,
				    size_t in_size,
				    size_t *out_size) {
    unsigned int i_out_size = *out_size;
    unsigned char *comp;

    comp = (unsigned char *)rans_compress(in, in_size, &i_out_size, 0);
    *out_size = i_out_size;

    return comp;
}

unsigned char *rans1_codec_compress(int method,
				    int level,
				    cram_slice *s,
				    unsigned char *in,
				    size_t in_size,
				    size_t *out_size) {
    unsigned int i_out_size = *out_size;
    unsigned char *comp;

    comp = rans_compress(in, in_size, &i_out_size, 1);
    *out_size = i_out_size;

    return comp;
}

unsigned char *rans_codec_uncompress(cram_slice *s,
				     unsigned char *in,
				     size_t in_size,
				     size_t *out_size) {
    unsigned int i_out_size = *out_size;
    unsigned char *uncomp;
    
    uncomp = rans_uncompress(in, in_size, &i_out_size, 0/*FIXME: unused*/);
    *out_size = i_out_size;

    return uncomp;
}

static cram_compressor rans0_codec = {
    RANS0, //FOUR_CC("RANS"),
    0,
    1.0,
    rans0_codec_name,
    rans0_codec_compress,
    rans_codec_uncompress,
};

static cram_compressor rans1_codec = {
    RANS0, //FOUR_CC("RANS"),
    0,
    1.02,
    rans1_codec_name,
    rans1_codec_compress,
    rans_codec_uncompress,
};

/* -----------------------------------------------------------------------------
 * RANS 4x16 "PR".
 * Order 0 / order 1
 * With / without bit "P"acking
 * With / without "R"LE.
 *
 * Hence 2x2x2 = 8 variants.  We use the same functions for each with different
 * parameters.
 */
static const char *rans_pr_codec_name(void) {
    return "RANS_PR";
}

unsigned char *rans_pr_codec_compress(int method,
				      int level,
				      cram_slice *s,
				      unsigned char *in,
				      size_t in_size,
				      size_t *out_size) {
    unsigned int i_out_size = *out_size;
    unsigned char *comp;
    static int rmethod[] = {0,1,64,65,128,129,192,193};

    comp = (unsigned char *)rans_compress_4x16(in, in_size, &i_out_size,
					       rmethod[method-RANS_PR0]);
    *out_size = i_out_size;

    return comp;
}

unsigned char *rans_pr_codec_uncompress(cram_slice *s,
					unsigned char *in,
					size_t in_size,
					size_t *out_size) {
    unsigned int i_out_size = *out_size;
    unsigned char *uncomp;
    
    uncomp = rans_uncompress_4x16(in, in_size, &i_out_size, 0/*FIXME: unused*/);
    *out_size = i_out_size;

    return uncomp;
}

static cram_compressor rans4x16pr_codec = {
    RANS_PR0, //FOUR_CC("RANS"),
    0,
    1.0,
    rans_pr_codec_name,
    rans_pr_codec_compress,
    rans_pr_codec_uncompress,
};


/* -----------------------------------------------------------------------------
 * Codec initialisation
 */

/*
 * Searches $CRAM_CODEC_DIR for dynamic loadable libraries (.so/.dll).
 * Any containing a cram_compressor_init() symbol are loaded up and
 * initialised.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int cram_compression_codec_init(void) {
    DIR *dir;
    struct dirent *ent;
    char *plugins = getenv("CRAM_CODEC_DIR");
    static int init_done = 0;

    if (init_done)
	return 0;

    /* FIXME: Add thread blocking here to force execution by one thread only */

    codec_hash = HashTableCreate(16, HASH_DYNAMIC_SIZE | HASH_NONVOLATILE_KEYS);
    if (!codec_hash)
	return -1;

    /* Add the inbuilt ones; cram_block_method enum */
    codecs = (cram_compressor **)calloc(BLOCK_METHOD_END, sizeof(*codecs));
    if (!codecs)
	return -1;
    //ncodecs = BLOCK_METHOD_END;
    ncodecs = 21; // matches 0xfffc0000 in cram_encode.c

    codecs[GZIP]     = &gzip_codec;
#ifdef HAVE_LIBBZ2
    codecs[BZIP2]    = &bzip2_codec;
#endif
#ifdef HAVE_LIBLZMA
    codecs[LZMA]     = &lzma_codec;
#endif
    codecs[RANS0]       = &rans0_codec;
    codecs[RANS1]       = &rans1_codec;
    codecs[RANS_PR0]    = &rans4x16pr_codec;
    codecs[RANS_PR1]    = &rans4x16pr_codec;
    codecs[RANS_PR64]   = &rans4x16pr_codec;
    codecs[RANS_PR65]   = &rans4x16pr_codec;
    codecs[RANS_PR128]  = &rans4x16pr_codec;
    codecs[RANS_PR129]  = &rans4x16pr_codec;
    codecs[RANS_PR192]  = &rans4x16pr_codec;
    codecs[RANS_PR193]  = &rans4x16pr_codec;
    codecs[GZIP_RLE] = &gzip_rle_codec;

    if (!plugins || !*plugins)
	return 0;

    if (!(dir = opendir(plugins)))
	return -1;

    while ((ent = readdir(dir))) {
        size_t len;
        void *dl_handle;
        char path[PATH_MAX];
        cram_compressor *(*codec_init)(void);
	HashData hd;

        if (ent->d_type != DT_REG && ent->d_type != DT_UNKNOWN)
            continue;

        len = strlen(ent->d_name);
        if (len < 3 || strcmp(&ent->d_name[len-3], ".so") != 0)
            continue;

        snprintf(path, PATH_MAX, "%s/%s", plugins, ent->d_name);
        if (!(dl_handle = dlopen(path, RTLD_NOW))) {
            perror(path);
            puts(dlerror());
            continue;
        }

        codec_init = (cram_compressor *(*)(void))dlsym(dl_handle, "cram_compressor_init");

        if (!codec_init) {
            dlclose(dl_handle);
            continue;
        }

        fprintf(stderr, "Adding codec %2d %.256s\n", ncodecs, ent->d_name);

        codecs = (cram_compressor **)realloc(codecs, (ncodecs+1) * sizeof(*codecs));
        codecs[ncodecs] = (*codec_init)();


	hd.i = ncodecs;
	codecs[ncodecs]->code &= 255; // CRAM has 1 byte
	HashTableAdd(codec_hash, (char *)&codecs[ncodecs]->code, 4, hd, NULL);

	ncodecs++;
    }

    closedir(dir);

    init_done = 1;

    return 0;
}

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
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 */

/*! \file
 * Include cram.h instead.
 *
 * This is an internal part of the CRAM system and is automatically included
 * when you #include cram.h.
 *
 * Implements the low level block compression methods.
 */

#ifndef _CRAM_BLOCK_COMPRESSION_H_
#define _CRAM_BLOCK_COMPRESSION_H_

#include "cram.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Frequency of metrics tests */
#define TRIAL_SPAN 50
#define NTRIALS 3

#define FOUR_CC(s) (((s)[0]<<24)+((s)[1]<<16)+((s)[2]<<8)+((s)[3]))

/*
 * A structure detailing how a cram_compressor looks.
 */
typedef struct {
    // 4 character code;
    // NB only bottom 8-bits are permitted in file format currently.
    uint32_t code;

    // bit-wise OR of DS_* content ids
    int64_t content_ids;

    // Weightings for compression and uncompression speeds.
    // Treat 1.0 as gzip equiv.
//    int compress_weight;
//    int uncompress_weight;

    // How many times more costly than gzip
    double cost;

    // String identifier, for debugging
    const char *(*name)(void);

    // Compresses a block.
    //
    // The compressed block size is returned in *out_size
    // Returns malloced buffer on success;
    //         NULL on failure
    unsigned char *(*compress_block)(int level,
                                     unsigned char *in,
                                     size_t in_size,
                                     size_t *out_size);

    // Uncompresses a block.
    //
    // The expected uncompressed size is supplied in *out_size, but the
    // function should also return the exact size.
    //
    // Returns malloced buffer on success;
    //         NULL on failure
    unsigned char *(*uncompress_block)(unsigned char *in,
                                       size_t in_size,
                                       size_t *out_size);
} cram_compressor;


/*! Uncompresses a CRAM block, if compressed.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_uncompress_block(cram_block *b);


/*! Compresses a block.
 *
 * Compresses a block using one of two different zlib strategies. If we only
 * want one choice set strat2 to be -1.
 *
 * The logic here is that sometimes Z_RLE does a better job than Z_FILTERED
 * or Z_DEFAULT_STRATEGY on quality data. If so, we'd rather use it as it is
 * significantly faster.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_compress_block(cram_fd *fd, cram_block *b, cram_metrics *metrics,
			int method, int level);

/*
 *! Initialises the inbuilt and external block compression codecs.
 *
 * Searches $CRAM_CODEC_DIR for dynamic loadable libraries (.so/.dll).
 * Any containing a cram_compressor_init() symbol are loaded up and
 * initialised.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure.
 */
int cram_compression_codec_init(void);

#ifdef __cplusplus
}
#endif

#endif /* _CRAM_BLOCK_COMPRESSION_H_ */

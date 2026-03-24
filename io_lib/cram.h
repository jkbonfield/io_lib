/*
Copyright (c) 2012-2013, 2015, 2018 Genome Research Ltd.
Author: James Bonfield <jkb@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*! \file
 * CRAM interface.
 *
 * Consider using the higher level hts_*() API for programs that wish to
 * be file format agnostic (see htslib/hts.h).
 *
 * This API should be used for CRAM specific code. The specifics of the
 * public API are implemented in cram_io.h, cram_encode.h and cram_decode.h
 * although these should not be included directly (use this file instead).
 */

#ifndef CRAM_ALL_H
#define CRAM_ALL_H

#include "bam.h"
#include "sam_header.h"

#ifdef WITH_CRAM
#else
typedef struct {
    void *header;
    void *refs;
} cram_fd;
typedef void refs_t;
#endif

#define SEQS_PER_SLICE 10000
#define BASES_PER_SLICE (SEQS_PER_SLICE*500)
#define SLICE_PER_CNT  1

enum cram_option {
    CRAM_OPT_DECODE_MD,
    CRAM_OPT_PREFIX,
    CRAM_OPT_VERBOSITY,
    CRAM_OPT_SEQS_PER_SLICE,
    CRAM_OPT_SLICES_PER_CONTAINER,
    CRAM_OPT_RANGE,
    CRAM_OPT_VERSION,
    CRAM_OPT_EMBED_REF,
    CRAM_OPT_IGNORE_MD5,
    CRAM_OPT_REFERENCE,
    CRAM_OPT_MULTI_SEQ_PER_SLICE,
    CRAM_OPT_NO_REF,
    CRAM_OPT_USE_BZIP2,
    CRAM_OPT_SHARED_REF,
    CRAM_OPT_NTHREADS,
    CRAM_OPT_THREAD_POOL,
    CRAM_OPT_BINNING,
    CRAM_OPT_USE_ARITH,
    CRAM_OPT_USE_LZMA,
    CRAM_OPT_REQUIRED_FIELDS,
    CRAM_OPT_USE_RANS,
    CRAM_OPT_IGNORE_CHKSUM,
    CRAM_OPT_BASES_PER_SLICE,
    CRAM_OPT_LOSSY_READ_NAMES,
    CRAM_OPT_PRESERVE_AUX_ORDER,
    CRAM_OPT_PRESERVE_AUX_SIZE,
    CRAM_OPT_WITH_BGZIP_INDEX,
    CRAM_OPT_OUTPUT_BGZIP_IDX,
    CRAM_OPT_USE_BSC,
    CRAM_OPT_USE_ZSTD,
    CRAM_OPT_USE_FQZ,
    CRAM_OPT_EMBED_CONS,
    CRAM_OPT_USE_TOK,
    CRAM_OPT_PROFILE
};

// REQUIRED_FIELDS
enum sam_fields {
    SAM_QNAME = 0x00000001,
    SAM_FLAG  = 0x00000002,
    SAM_RNAME = 0x00000004,
    SAM_POS   = 0x00000008,
    SAM_MAPQ  = 0x00000010,
    SAM_CIGAR = 0x00000020,
    SAM_RNEXT = 0x00000040,
    SAM_PNEXT = 0x00000080,
    SAM_TLEN  = 0x00000100,
    SAM_SEQ   = 0x00000200,
    SAM_QUAL  = 0x00000400,
    SAM_AUX   = 0x00000800,
    SAM_RGAUX = 0x00001000,
};

//cram_fd *cram_open(const char *filename, const char *mode);
//int cram_close(cram_fd *fd);
//int cram_flush(cram_fd *fd);
//int cram_write_eof_block(cram_fd *fd);
//int cram_eof(cram_fd *fd);
//int cram_set_option(cram_fd *fd, enum cram_option opt, ...);
//int cram_set_voption(cram_fd *fd, enum cram_option opt, va_list args);
//int cram_load_reference(cram_fd *fd, char *fn);


#endif

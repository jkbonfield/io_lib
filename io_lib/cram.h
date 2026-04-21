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

#define HTS_NO_SAM_HDR
#include <htslib/cram.h>

// Minimal implementation for gap5's export_snps.c to build, but not to run.
// (That code is internal debugging and not exported to the users.)
typedef struct {
    char *name;
    char *fn;
    int64_t length;
    int64_t offset;
    int bases_per_line;
    int line_length;
    int64_t count;         // for shared references so we know to dealloc seq
    char *seq;
    //mFILE *mf;
    int is_md5;            // Reference comes from a raw seq found by MD5
} ref_entry;

typedef struct refs_t {
    // big enough to swallow htslib's copy (112)
    uint8_t htslib_private[1024];

    // io_lib bits accessed by Gap5.  This ensures any modification doesn't
    // overwrite htslib's data.
    void *fp;
    ref_entry **ref_id;
    int nref;
} refs_t;

static inline
char *load_ref_portion(void *fp, ref_entry *e, int start, int end) {
    return NULL;
}

static inline
refs_t *refs_load_fai(refs_t *r_orig, char *fn, int is_err) {
    return NULL;
}

static inline
void refs_free(refs_t *r) {}

#define SEQS_PER_SLICE 10000
#define BASES_PER_SLICE (SEQS_PER_SLICE*500)
#define SLICE_PER_CNT  1
#define CRAM_OPT_PROFILE HTS_OPT_PROFILE

/* ---------------------------------------------------------------------------
 * CRAM_IO_CUSTOM_BUFFERING mode, used by libmaus
 */

typedef size_t (*cram_io_C_FILE_fread_t)(void *ptr, size_t size, size_t nmemb, void *stream);
typedef size_t (*cram_io_C_FILE_fwrite_t)(void *ptr, size_t size, size_t nmemb, void *stream);
typedef int    (*cram_io_C_FILE_fseek_t)(void * fd, off_t offset, int whence);
typedef off_t  (*cram_io_C_FILE_ftell_t)(void * fd);

typedef struct {
    void                   *user_data;
    cram_io_C_FILE_fread_t  fread_callback;
    cram_io_C_FILE_fseek_t  fseek_callback;
    cram_io_C_FILE_ftell_t  ftell_callback;
} cram_io_input_t;

typedef struct {
    void                   *user_data;
    cram_io_C_FILE_fwrite_t fwrite_callback;
    cram_io_C_FILE_ftell_t  ftell_callback;
} cram_io_output_t;

typedef cram_io_input_t * (*cram_io_allocate_read_input_t)(char const * filename, int const decompress);
typedef cram_io_input_t * (*cram_io_deallocate_read_input_t)(cram_io_input_t * obj);

typedef cram_io_output_t * (*cram_io_allocate_write_output_t)(char const * filename);
typedef cram_io_output_t * (*cram_io_deallocate_write_output_t)(cram_io_output_t * obj);

typedef struct {
    int refid;
    int64_t start;
    int64_t end;
} cram_range;


//typedef struct {
//    /* input buffer size */
//    size_t         fp_in_buf_size;
//    /* input buffer base pointer */
//    char          *fp_in_buffer;
//    /* position of buffer start in file */
//    uint64_t       fp_in_buf_start;
//    /* start of window pointer; same as fp_in_buffer */
//    char          *fp_in_buf_pa;
//    /* window current pointer */
//    char          *fp_in_buf_pc;
//    /* window end pointer;  same as fp_in_buffer + fp_in_buf_size (no seeks) */
//    char          *fp_in_buf_pe;    
//} cram_fd_input_buffer;
//
//typedef struct {
//    /* output buffer size */
//    size_t         fp_out_buf_size;
//    /* output buffer base pointer */
//    char          *fp_out_buffer;
//    /* position of buffer start in file */
//    uint64_t       fp_out_buf_start;
//    /* start of window pointer; same as fp_out_buffer */
//    char          *fp_out_buf_pa;
//    /* window current pointer */
//    char          *fp_out_buf_pc;
//    /* window end pointer */
//    char          *fp_out_buf_pe;    
//} cram_fd_output_buffer;

#endif

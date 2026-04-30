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
    int is_md5;            // Reference comes from a raw seq found by MD5
} ref_entry;

typedef struct refs_t {
    // io_lib bits accessed by Gap5.  This ensures any modification doesn't
    // overwrite htslib's data.
    // Please use sam_hdr_* API in sam_header.h in preference.
    void *fp;
    ref_entry **ref_id;
    int nref;
} refs_t;

#define SEQS_PER_SLICE 10000
#define BASES_PER_SLICE (SEQS_PER_SLICE*500)
#define SLICE_PER_CNT  1
#define CRAM_OPT_PROFILE HTS_OPT_PROFILE

typedef struct {
    int refid;
    int64_t start;
    int64_t end;
} cram_range;

#endif

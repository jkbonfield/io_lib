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
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013, 2026
 */

#ifndef _SAM_HDR_H_
#define _SAM_HDR_H_

#ifdef __cplusplus
extern "C" {
#endif


/*
 * This is a shim over htslib's header implementation which has a lot of
 * identical function names that would otherwise cause clashes.
 *
 * However the reason for this is htslib's implementation was derived from
 * io_lib's anyway and then updated, so we may as well switch to it with a
 * bit of typedefing and shim functions.
 */
#include <stdarg.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include "io_lib/dstring.h"

typedef sam_hdr_t       SAM_hdr;

#define sam_hdr_add      sam_hdr_add_line
//#define sam_hdr_free     sam_hdr_destroy
#define sam_hdr_name2ref sam_hdr_name2tid
#define sam_hdr_add_PG   sam_hdr_add_pg

// These are private in htslib.  I'm not sure why
enum sam_sort_order {
    ORDER_UNKNOWN  =-1,
    ORDER_UNSORTED = 0,
    ORDER_NAME     = 1,
    ORDER_COORD    = 2
};

enum sam_group_order {
    ORDER_NONE      =-1,
    ORDER_QUERY     = 0,
    ORDER_REFERENCE = 1
};

static inline enum sam_sort_order sam_hrecs_sort_order(sam_hdr_t *hdr) {
    kstring_t str = KS_INITIALIZE;
    if (sam_hdr_find_tag_hd(hdr, "SO", &str) < 0)
        return ORDER_UNKNOWN;

    int ret;
    if (strcmp(str.s, "coordinate") == 0)
        ret = ORDER_COORD;
    else if (strcmp(str.s, "name") == 0)
        ret = ORDER_NAME;
    else if (strcmp(str.s, "unsorted") == 0)
        ret = ORDER_UNSORTED;
    else
        ret = ORDER_UNKNOWN;

    ks_free(&str);
    return ret;
}

#ifdef __cplusplus
}
#endif

#endif /* _SAM_HDR_H_ */

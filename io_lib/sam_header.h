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

/*! Parsed \@SQ lines */
typedef struct {
    char *name;
    uint32_t len;
} SAM_SQ;


// A container for htslib's header API instead (which is derived from this
// code originally).
typedef struct {
    dstring_t *text;          //!< concatenated text, indexed by SAM_hdr_tag
    sam_hdr_t *hdr;           //!<htslib header struct
    int nref;                 //!< Number of \@SQ lines
    SAM_SQ *ref;              //!< Array of parsed \@SQ lines
} SAM_hdr;

//typedef sam_hdr_t       SAM_hdr;

#define sam_hdr_add(h,t,...) sam_hdr_add_line((h),(t),__VA_ARGS__)
#define sam_hdr_add_PG(h,n,...) sam_hdr_add_pg((h),(n),__VA_ARGS__)


int sam_hdr_name2ref(SAM_hdr *h, const char *name);

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

SAM_hdr *sam_hdr_convert(sam_hdr_t *hdr);
void sam_hdr_free(SAM_hdr *hdr);

// io_lib's header.
static inline sam_hdr_t *sam_hdr_parse_htslib(const char *hdr, int len) {
    return sam_hdr_parse(len, hdr);
}

static inline SAM_hdr *sam_hdr_parse_iolib(const char *hdr, int len) {
    SAM_hdr *h = sam_hdr_convert(sam_hdr_parse(len, hdr));
    if (!h)
        return NULL;
    return h;
}

// Map io_lib header calls to sam_hdr_parse_ which converts from htslib's
// identically named sam_hdr_parse function.
#define sam_hdr_parse(h,l) sam_hdr_parse_iolib(h,l)


static inline int SAM_hdr_nref(SAM_hdr *h) {
    return h->nref;
}

static inline const char *SAM_hdr_tid2name(SAM_hdr *h, int id) {
    return sam_hdr_tid2name(h->hdr, id);
}

static inline hts_pos_t SAM_hdr_tid2len(SAM_hdr *h, int id) {
    return sam_hdr_tid2len(h->hdr, id);
}

static inline SAM_hdr *SAM_hdr_dup(SAM_hdr *h) {
    return sam_hdr_convert(sam_hdr_dup(h->hdr));
}

#ifdef __cplusplus
}
#endif

#endif /* _SAM_HDR_H_ */

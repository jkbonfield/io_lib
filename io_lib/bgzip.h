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

#ifndef _BGZIP_H_
#define _BGZIP_H_

typedef struct gzi {
    uint64_t n;
    uint64_t *c_off;
    uint64_t *u_off;
} gzi;

int gzi_index_add_block(gzi *idx, uint64_t c_off, uint64_t u_off);
int gzi_index_dump(gzi *idx, const char *bname, const char *suffix);
gzi *gzi_index_init();
gzi *gzi_index_load(const char *fn);
void gzi_index_free(gzi *idx);
uint64_t gzi_load(FILE *fp, gzi *idx, uint64_t ustart, uint64_t uend, char *out);

struct bzi_FILE;
typedef struct bzi_FILE bzi_FILE;

bzi_FILE *bzi_open(const char *path, const char *mode);
void bzi_close(bzi_FILE *zp);
size_t bzi_read(void *ptr, size_t size, size_t nmemb, bzi_FILE *zp);
int bzi_seek(bzi_FILE *zp, off_t offset, int whence);

#endif /* _BGZIP_H_ */

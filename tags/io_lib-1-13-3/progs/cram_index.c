/*
 * Copyright (c) 2013 Genome Research Ltd.
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
 * A debugging program to dump out information on the layout of a CRAM file.
 * It's an abomination frankly, but isn't intended for production use.
 */

#include "io_lib_config.h"

#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include <io_lib/cram.h>
#include <io_lib/zfio.h>

int main(int argc, char **argv) {
    cram_fd *fd;
    cram_container *c;
    off_t cpos, spos, hpos;
    zfp *fp;
    char fn[PATH_MAX];

    if (argc != 2) {
	fprintf(stderr, "Usage: cram_index filename.cram\n");
	return 1;
    }

    if (NULL == (fd = cram_open(argv[1], "rb"))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", argv[1]);
	return 1;
    }

    sprintf(fn, "%s.crai", argv[1]);
    if (!(fp = zfopen(fn, "wz"))) {
	perror(fn);
	return 1;
    }

    cpos = ftello(fd->fp);
    while ((c = cram_read_container(fd))) {
	int j;

	if (fd->err) {
	    perror("Cram container read");
	    return 1;
	}

	hpos = ftello(fd->fp);

	if (!(c->comp_hdr_block = cram_read_block(fd)))
	    return 1;
	assert(c->comp_hdr_block->content_type == COMPRESSION_HEADER);

	c->comp_hdr = cram_decode_compression_header(fd, c->comp_hdr_block);
	if (!c->comp_hdr)
	    return 1;

	// 1.1 format
	for (j = 0; j < c->num_landmarks; j++) {
	    char buf[1024];
	    cram_slice *s;
	    int sz;
	    
	    spos = ftello(fd->fp);
	    assert(spos - cpos - c->offset == c->landmark[j]);

	    s = cram_read_slice(fd);

	    sz = (int)(ftello(fd->fp) - spos);

	    sprintf(buf, "%d\t%d\t%d\t%"PRId64"\t%d\t%d\n",
		    s->hdr->ref_seq_id, s->hdr->ref_seq_start,
		    s->hdr->ref_seq_span, (int64_t)cpos,
		    c->landmark[j], sz);
	    zfputs(buf, fp);

	    cram_free_slice(s);
	}

	cpos = ftello(fd->fp);
	assert(cpos == hpos + c->length);

	cram_free_container(c);
    }

    cram_close(fd);
    zfclose(fp);

    return 0;
}

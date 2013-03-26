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
    off_t first_c, cpos, spos, hpos;
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

    first_c = cpos = ftello(fd->fp);
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

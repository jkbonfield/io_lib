/*
 * -----------------------------------------------------------------------------
 * A debugging program to dump out information on the layout of a CRAM file.
 *
 * TODO:
 * 1) Much of the content of main() needs to be broken into separate functions
 *    and to make its way into the library side of things.
 *
 *    Ideally we should implement a cram sequence iterator.
 *
 * 2) Remove the fixed length buffers for things like name[1000].
 *
 * 3) Allow writing directly into the dstring rather than to a
 *    temporary buffer and copying in afterwards.
 *    This means implementing a dstring_grow() function too.
 *
 *    Even better is if we determine that the EXTERNAL block can be
 *    used "as is" without needing to copy as it has no other encoding
 *    layer applied. This is doable by checking the codec. If the slice
 *    codec for names is EXTERNAL then we can use the block->data directly,
 *    otherwise we need to decode and store in dstring.
 */

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

#include <io_lib/cram.h>


int main(int argc, char **argv) {
    cram_fd *fd;
    cram_container *c;
    size_t pos, pos2;
    bam_file_t *bfd;
    bam_seq_t *bam = NULL;
    refs *refs;
    size_t bam_alloc = 0;
    char mode[4] = {'w', '\0', '\0', '\0'};
    char *prefix = NULL;

    if (argc >= 2 && strcmp(argv[1], "-b") == 0) {
	mode[1] = 'b';
	argc--;
	argv++;
    }

    if (argc >= 2 && argv[1][0] == '-' && argv[1][1] >= '0' && argv[1][1] <= '9') {
	mode[2] = argv[1][1];
	argc--;
	argv++;
    }

    if (argc >= 2 && strcmp(argv[1], "-u") == 0) {
	mode[2] = '0';
	argc--;
	argv++;
    }

    if (argc >= 3 && strcmp(argv[1], "-p") == 0) {
	prefix = argv[2];
	argc-=2;
	argv+=2;
    }

    bfd = bam_open("-", mode);

    if (argc != 2 && argc != 3) {
	fprintf(stderr, "Usage: cram_dump [-b] [-0..9] [-u] filename.cram [ref.fa]\n");
	return 1;
    }


    if (argc == 3) {
	refs = load_reference(argv[2]);
    } else {
	refs = NULL;
    }

    if (NULL == (fd = cram_open(argv[1], "rb"))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", argv[1]);
	return 1;
    }
    if (prefix)
	cram_set_prefix(fd, prefix);

    bfd->header_len = fd->SAM_hdr->length;
    bfd->header = malloc(fd->SAM_hdr->length);
    memcpy(bfd->header, fd->SAM_hdr->header, fd->SAM_hdr->length);

    if (-1 == bam_parse_header(bfd))
        return -1;

    bam_write_header(bfd);
    if (refs)
	refs2id(refs, bfd);

    pos = ftello(fd->fp);
    while ((c = cram_read_container(fd))) {
	int j;

	if (fd->err) {
	    perror("Cram container read");
	    return 1;
	}

	for (j = 0; j < c->num_landmarks; j++) {
	    cram_slice *s;
	    int id, rec;
	    
	    pos2 = ftello(fd->fp);
	    assert(pos2 - pos - c->offset == c->landmark[j]);

	    s = cram_read_slice(fd);

	    // FIXME: this should be correct, but we have a bug in the
	    // Java CRAM implementation?
	    s->last_apos = s->hdr->ref_seq_start;
	    
	    for (id = 0; id < s->hdr->num_blocks; id++)
		cram_uncompress_block(s->block[id]);

	    /* Test decoding of 1st seq */
	    if (cram_decode_slice(c, s, bfd, refs) != 0) {
		fprintf(stderr, "Failure to decode slice\n");
		return 1;
	    }

	    /* Convert to SAM */
	    for (rec = 0; rec < s->hdr->num_records; rec++) {
		cram_to_bam(bfd, fd, s, &s->crecs[rec], rec, &bam, &bam_alloc);
		bam_put_seq(bfd, bam);
	    }
	    
	    cram_free_slice(s);
	}

	cram_free_container(c);
	pos = ftello(fd->fp);
    }

    cram_close(fd);
    bam_close(bfd);

    if (refs)
	free_refs(refs);

    free(bam);

    return 0;
}

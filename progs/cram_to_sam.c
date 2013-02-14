/*
 * Author: James Bonfield, Sanger Institute, 2013.
 *
 * Converts a CRAM file to a SAM or BAM file.
 */

#include "io_lib_config.h"

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

#include <io_lib/cram.h>

void usage(FILE *fp) {
    fprintf(fp, "Usage: cram_to_sam [-m] [-b] [-0..9] [-u] "
	    "filename.cram ref.fa\n\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "    -m             Generate MD and NM tags:\n");
    fprintf(fp, "    -b             Output in BAM (defaults to SAM)\n");
    fprintf(fp, "    -1 to -9       Set zlib compression level for BAM\n");
    fprintf(fp, "    -0 or -u       Output uncompressed, if BAM.\n");
    fprintf(fp, "    -p str         Set the prefix for auto-generated seq. names\n");
}

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
    int decode_md = 0;
    cram_opt opt;
    int C;


    while ((C = getopt(argc, argv, "bu0123456789mp:h")) != -1) {
	switch (C) {
	case 'b':
	    mode[1] = 'b';
	    break;

	case 'u':
	    mode[2] = '0';
	    break;

	case '0': case '1': case '2': case '3': case '4':
	case '5': case '6': case '7': case '8': case '9':
	    mode[2] = C;
	    break;

	case 'm':
	    decode_md = 1;
	    break;

	case 'p':
	    prefix = optarg;
	    break;

	case 'h':
	    usage(stdout);
	    return 0;

	case '?':
	    fprintf(stderr, "Unrecognised option: -%c\n", optopt);
	    usage(stderr);
	    return 1;
	}
    }

    bfd = bam_open("-", mode);

    if (argc - optind != 2) {
	usage(stderr);
	return 1;
    }

    if (NULL == (fd = cram_open(argv[optind], "rb"))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", argv[optind]);
	return 1;
    }

    if (prefix)
	opt.s = prefix, cram_set_option(fd, CRAM_OPT_PREFIX, &opt);

    if (decode_md)
	opt.i = decode_md, cram_set_option(fd, CRAM_OPT_DECODE_MD, &opt);

    cram_load_reference(fd, argv[optind+1]);

    bfd->header_len = fd->SAM_hdr->header_len;
    bfd->header = malloc(fd->SAM_hdr->header_len);
    memcpy(bfd->header, fd->SAM_hdr->header, fd->SAM_hdr->header_len);

    if (-1 == bam_parse_header(bfd))
        return -1;

    bam_write_header(bfd);
    if (fd->refs)
	refs2id(fd->refs, bfd);

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
	    if (cram_decode_slice(fd, c, s, bfd) != 0) {
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

    free(bam);

    return 0;
}

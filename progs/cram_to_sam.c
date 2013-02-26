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
	    "filename.cram ref.fa [output_filename]\n\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "    -m             Generate MD and NM tags:\n");
    fprintf(fp, "    -b             Output in BAM (defaults to SAM)\n");
    fprintf(fp, "    -1 to -9       Set zlib compression level for BAM\n");
    fprintf(fp, "    -0 or -u       Output uncompressed, if BAM.\n");
    fprintf(fp, "    -p str         Set the prefix for auto-generated seq. names\n");
}

int main(int argc, char **argv) {
    cram_fd *fd;
    bam_file_t *bfd;
    bam_seq_t *bam = NULL;
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

    if (argc - optind != 2 && argc - optind != 3) {
	usage(stderr);
	return 1;
    }

    if (argc - optind == 2) {
	if (NULL == (bfd = bam_open("-", mode))) {
	    fprintf(stderr, "Failed to open SAM/BAM output\n.");
	    return 1;
	}
    } else {
	if (NULL == (bfd = bam_open(argv[optind+2], mode))) {
	    fprintf(stderr, "Failed to open SAM/BAM output\n.");
	    perror(argv[optind+2]);
	    return 1;
	}
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
        return 1;

    bam_write_header(bfd);
    if (fd->refs)
	refs2id(fd->refs, bfd);

    while (cram_get_bam_seq(fd, &bam, &bam_alloc) == 0) {
	bam_put_seq(bfd, bam);
    }

    cram_close(fd);
    bam_close(bfd);

    free(bam);

    return 0;
}

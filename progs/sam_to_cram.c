/*
 * Author: James Bonfield, Sanger Institute, 2013.
 *
 * Converts a SAM or BAM file into a CRAM file.
 *
 * Usage:
 *     sam_to_cram [-level] input.sam reference.fasta [output.cram]
 */

#include "io_lib_config.h"

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

#include <io_lib/cram.h>

void usage(FILE *fp) {
    fprintf(fp, "Usage: sam_to_cram [-0..9] [-u] [-v] [-s] [-S] in.sam/bam ref.fa [output.cram]\n\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "    -1 to -9       Set zlib compression level for CRAM\n");
    fprintf(fp, "    -0 or -u       No zlib compression.\n");
    fprintf(fp, "    -v             Verbose output.\n");
    fprintf(fp, "    -s integer     Sequences per slice, default %d.\n",
	    SEQS_PER_SLICE);
    fprintf(fp, "    -S integer     Slices per container, default %d.\n",
	    SLICE_PER_CNT);
    fprintf(fp, "    -V version     Specify the CRAM format version to write (eg 1.1, 2.0)\n");
}

int main(int argc, char **argv) {
    cram_fd *out;
    bam_file_t *in;
    bam_seq_t *s = NULL;
    char *out_fn;
    int level = '\0'; // nul terminate string => auto level
    char out_mode[4];
    int c, verbose = 0;
    cram_opt opt;
    int s_opt = 0, S_opt = 0;
    char *arg_list;

    while ((c = getopt(argc, argv, "u0123456789hvs:S:V:")) != -1) {
	switch (c) {
	case '0': case '1': case '2': case '3': case '4':
	case '5': case '6': case '7': case '8': case '9':
	    level = c;
	    break;
	    
	case 'u':
	    level = '0';
	    break;

	case 'h':
	    usage(stdout);
	    return 0;

	case 'v':
	    verbose++;
	    break;

	case 's':
	    s_opt = atoi(optarg);
	    break;

	case 'S':
	    S_opt = atoi(optarg);
	    break;

	case 'V':
	    cram_set_option(NULL, CRAM_OPT_VERSION, (cram_opt *)&optarg);
	    break;

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

    /* opening */
    if (NULL == (in = bam_open(argv[optind], "rb"))) {
	perror(argv[optind]);
	return 1;
    }

    out_fn = argc - optind == 3 ? argv[optind+2] : "-";
    sprintf(out_mode, "wb%c", level);
    if (NULL == (out = cram_open(out_fn, out_mode))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", out_fn);
	return 1;
    }

    /* SAM Header */
    if (!(arg_list = stringify_argv(argc, argv)))
	return 1;
    sam_header_add_PG(in->header, "sam_to_cram",
		      "VN", PACKAGE_VERSION,
		      "CL", arg_list, NULL);
    free(arg_list);

    if (-1 == cram_write_SAM_hdr(out, in->header))
	return 1;

    out->SAM_hdr = in->header;
    
    cram_load_reference(out, argv[optind+1]);
    if (!out->refs)
	return 1;
    refs2id(out->refs, out->SAM_hdr);

    opt.i = verbose;
    cram_set_option(out, CRAM_OPT_VERBOSITY, &opt);
    if (s_opt) {
	opt.i = s_opt;
	cram_set_option(out, CRAM_OPT_SEQS_PER_SLICE, &opt);
    }
    if (S_opt) {
	opt.i = S_opt;
	cram_set_option(out, CRAM_OPT_SLICES_PER_CONTAINER, &opt);
    }

    /* Sequence iterators */
    while (bam_next_seq(in, &s) > 0) {
	if (-1 == cram_put_bam_seq(out, s)) {
	    fprintf(stderr, "Failed in cram_put_bam_seq()\n");
	    return 1;
	}
    }

    bam_close(in);
    out->SAM_hdr = NULL; // freed by bam_close()
    if (-1 == cram_close(out)) {
	fprintf(stderr, "Failed in cram_close()\n");
	return 1;
    }

    if (s)
	free(s);

    return 0;
}

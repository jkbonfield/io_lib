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
    fprintf(fp, "Usage: sam_to_cram [-r ref.fa] [-0..9] [-u] [-v] [-s int] "
	    "[-S int] in.sam/bam [output.cram]\n\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "    -r ref.fa      Specifies the reference file.\n");
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
    char *arg_list, *ref_fn = NULL;

    while ((c = getopt(argc, argv, "u0123456789hvs:S:V:r:")) != -1) {
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

	case 'r':
	    ref_fn = optarg;
	    break;

	case '?':
	    fprintf(stderr, "Unrecognised option: -%c\n", optopt);
	    usage(stderr);
	    return 1;
	}
    }

    if (argc - optind != 1 && argc - optind != 2) {
	usage(stderr);
	return 1;
    }

    /* opening */
    if (NULL == (in = bam_open(argv[optind], "rb"))) {
	perror(argv[optind]);
	return 1;
    }

    out_fn = argc - optind == 2 ? argv[optind+1] : "-";
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

    /* Find and load reference */
    if (!ref_fn) {
	SAM_hdr_type *ty = sam_header_find(in->header, "SQ", NULL, NULL);
	if (ty) {
	    int len;
	    ref_fn = sam_header_find_key2(in->header, ty, "UR", &len);
	    if (ref_fn) {
		ref_fn[len] = 0;
		ref_fn += 3;
		if (strncmp(ref_fn, "file:", 5) == 0)
		    ref_fn += 5;
	    }
	}
    }

    out->SAM_hdr = in->header;
    if (ref_fn)
	cram_load_reference(out, ref_fn);

    if (!out->refs) {
	fprintf(stderr, "Unable to open reference.\n"
		"Please specify a valid reference with -r ref.fa option.\n");
	return 1;
    }
    refs2id(out->refs, out->SAM_hdr);

    if (-1 == cram_write_SAM_hdr(out, in->header))
	return 1;

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

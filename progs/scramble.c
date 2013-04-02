#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <fcntl.h>
#include <zlib.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>

#include <io_lib/scram.h>
#include <io_lib/os.h>

static char *parse_format(char *str) {
    if (strcmp(str, "sam") == 0 || strcmp(str, "SAM") == 0)
	return "";

    if (strcmp(str, "bam") == 0 || strcmp(str, "BAM") == 0)
	return "b";

    if (strcmp(str, "cram") == 0 || strcmp(str, "CRAM") == 0)
	return "c";

    fprintf(stderr, "Unrecognised file format '%s'\n", str);
    exit(1);
}

static char *detect_format(char *fn) {
    char *cp = strrchr(fn, '.');

    if (!cp)
	return "";

    if (strcmp(cp, ".sam") == 0 || strcmp(cp, ".SAM") == 0)
	return "";
    if (strcmp(cp, ".bam") == 0 || strcmp(cp, ".BAM") == 0)
	return "b";
    if (strcmp(cp, ".cram") == 0 || strcmp(cp, ".CRAM") == 0)
	return "c";

    return "";
}

static void usage(FILE *fp) {
    fprintf(fp, "  -=- sCRAMble -=-     version %s\n", PACKAGE_VERSION);
    fprintf(fp, "Author: James Bonfield, Wellcome Trust Sanger Institute. 2013\n\n");

    fprintf(fp, "Usage:    scramble [options] [input_file [output_file]]\n");

    fprintf(fp, "Options:\n");
    fprintf(fp, "    -I format      Set input format:  \"bam\", \"sam\" or \"cram\".\n");
    fprintf(fp, "    -O format      Set output format: \"bam\", \"sam\" or \"cram\".\n");
    fprintf(fp, "    -1 to -9       Set zlib compression level.\n");
    fprintf(fp, "    -0 or -u       No zlib compression.\n");
    //fprintf(fp, "    -v             Verbose output.\n");
    fprintf(fp, "    -R range       [Cram] Specifies the refseq:start-end range\n");
    fprintf(fp, "    -r ref.fa      [Cram] Specifies the reference file.\n");
    fprintf(fp, "    -s integer     [Cram] Sequences per slice, default %d.\n",
	    SEQS_PER_SLICE);
    fprintf(fp, "    -S integer     [Cram] Slices per container, default %d.\n",
	    SLICE_PER_CNT);
    fprintf(fp, "    -V version     [Cram] Specify the file format version to write (eg 1.1, 2.0)\n");
    fprintf(fp, "    -X             [Cram] Embed reference sequence.\n");
}

int main(int argc, char **argv) {
    scram_fd *in, *out;
    bam_seq_t *s;
    char imode[10], *in_f = "", omode[10], *out_f = "";
    int level = '\0'; // nul terminate string => auto level
    int c, verbose = 0;
    int s_opt = 0, S_opt = 0, embed_ref = 0, ignore_md5 = 0;
    char *ref_fn = NULL;
    int start, end, multi_seq = 0;
    char ref_name[1024] = {0};
    refs *refs;

    /* Parse command line arguments */
    while ((c = getopt(argc, argv, "u0123456789hvs:S:V:r:XI:O:R:!M")) != -1) {
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
	    cram_set_option(NULL, CRAM_OPT_VERSION, optarg);
	    break;

	case 'r':
	    ref_fn = optarg;
	    break;

	case 'X':
	    embed_ref = 1;
	    break;

	case 'I':
	    in_f = parse_format(optarg);
	    break;

	case 'O':
	    out_f = parse_format(optarg);
	    break;

	case 'R': {
	    char *cp = strchr(optarg, ':');
	    if (cp) {
		*cp = 0;
		switch (sscanf(cp+1, "%d-%d", &start, &end)) {
		case 1:
		    end = start;
		    break;
		case 2:
		    break;
		default:
		    fprintf(stderr, "Malformed range format\n");
		    return 1;
		}
	    } else {
		start = INT_MIN;
		end   = INT_MAX;
	    }
	    strncpy(ref_name, optarg, 1023);
	    break;
	}

	case '!':
	    ignore_md5 = 1;
	    break;

	case 'M':
	    multi_seq = 1;
	    break;

	case '?':
	    fprintf(stderr, "Unrecognised option: -%c\n", optopt);
	    usage(stderr);
	    return 1;
	}
    }    

    if (argc - optind > 2) {
	fprintf(stderr, "Usage: scramble [input_file [output_file]]\n");
	return 1;
    }
    

    /* Open up input and output files */
    sprintf(imode, "r%s%c", in_f, level);
    if (argc - optind > 0) {
	if (*in_f == 0)
	    sprintf(imode, "r%s%c", detect_format(argv[optind]), level);
	if (!(in = scram_open(argv[optind], imode))) {
	    fprintf(stderr, "Failed to open bam file %s\n", argv[optind]);
	    return 1;
	}
    } else {
	if (!(in = scram_open("-", imode))) {
	    fprintf(stderr, "Failed to open bam file %s\n", argv[optind]);
	    return 1;
	}
    }
    if (!in->is_bam && ref_fn) {
	cram_load_reference(in->c, ref_fn);
	if (!in->c->refs && !embed_ref) {
	    fprintf(stderr, "Unable to find an appropriate reference.\n"
		    "Please specify a valid reference with "
		    "-r ref.fa option.\n");
	    return 1;
	}
    }

    sprintf(omode, "w%s%c", out_f, level);
    if (argc - optind > 1) {
	if (*out_f == 0)
	    sprintf(omode, "w%s%c", detect_format(argv[optind+1]), level);
	if (!(out = scram_open(argv[optind+1], omode))) {
	    fprintf(stderr, "Failed to open bam file %s\n", argv[optind+1]);
	    return 1;
	}
    } else {
	if (!(out = scram_open("-", omode))) {
	    fprintf(stderr, "Failed to open bam file %s\n", argv[optind+1]);
	    return 1;
	}
    }


    /* Set any format specific options */
    scram_set_refs(out, refs = scram_get_refs(in));

    scram_set_option(out, CRAM_OPT_VERBOSITY, verbose);
    if (s_opt)
	scram_set_option(out, CRAM_OPT_SEQS_PER_SLICE, s_opt);

    if (S_opt)
	scram_set_option(out, CRAM_OPT_SLICES_PER_CONTAINER, S_opt);

    if (embed_ref)
	scram_set_option(out, CRAM_OPT_EMBED_REF, embed_ref);

    if (multi_seq)
	scram_set_option(out, CRAM_OPT_MULTI_SEQ_PER_SLICE, multi_seq);

    if (ignore_md5)
	scram_set_option(in, CRAM_OPT_IGNORE_MD5, ignore_md5);
    

    /* Copy header and refs from in to out, for writing purposes */
    scram_set_header(out, scram_get_header(in));

    // Needs doing after loading the header.
    if (ref_fn) 
	scram_set_option(out, CRAM_OPT_REFERENCE, ref_fn);

    if (scram_get_header(in)) {
	if (scram_write_header(out))
	    return 1;
    }


    /* Support for sub-range queries, currently implemented for CRAM only */
    if (*ref_name != 0) {
	cram_range r;
	int refid;

	if (in->is_bam) {
	    fprintf(stderr, "Currently the -R option is only implemented for CRAM indices\n");
	    return 1;
	}
	    
	cram_index_load(in->c, argv[optind]);

	refid = sam_header_name2ref(in->c->SAM_hdr, ref_name);

	if (refid == -1 && *ref_name != '*') {
	    fprintf(stderr, "Unknown reference name '%s'\n", ref_name);
	    return 1;
	}
	r.refid = refid;
	r.start = start;
	r.end = end;
	scram_set_option(in, CRAM_OPT_RANGE, &r);
    }

    /* Do the actual file format conversion */
    s = NULL;
    while (scram_get_seq(in, &s) >= 0) {
	if (-1 == scram_put_seq(out, s))
	    return 1;
    }
    if (!scram_eof(in))
	return 1;

    /* Finally tidy up and close files */
    scram_set_header(out, NULL);

    if (refs == scram_get_refs(out)) {
	scram_set_refs(out, NULL);
    }

    scram_close(in);
    scram_close(out);

    if (s)
	free(s);

    return 0;
}

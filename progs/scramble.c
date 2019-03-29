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
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

// Enable if we want V3.1 support.  TODO: add a configure param for this
#define HAVE_FQZ

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <fcntl.h>
#include <zlib.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>

#if defined(__MINGW32__) || defined(__FreeBSD__) || defined(__APPLE__)
#   include <getopt.h>
#endif

#include <io_lib/scram.h>
#include <io_lib/os.h>

static char *parse_format(char *str) {
    if (strcmp(str, "sam") == 0 || strcmp(str, "SAM") == 0)
	return "s";

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
	return "s";
    if (strcmp(cp, ".bam") == 0 || strcmp(cp, ".BAM") == 0)
	return "b";
    if (strcmp(cp, ".cram") == 0 || strcmp(cp, ".CRAM") == 0)
	return "c";

    return "";
}

static void usage(FILE *fp) {
    fprintf(fp, "  -=- sCRAMble -=-     version %s\n", PACKAGE_VERSION);
    fprintf(fp, "Author: James Bonfield, Wellcome Trust Sanger Institute. 2013-2018\n\n");

    fprintf(fp, "Usage:    scramble [options] [input_file [output_file]]\n");

    fprintf(fp, "Options:\n");
    fprintf(fp, "    -I format      Set input format:  \"bam\", \"sam\" or \"cram\".\n");
    fprintf(fp, "    -O format      Set output format: \"bam\", \"sam\" or \"cram\".\n");
    fprintf(fp, "    -1 to -9       Set compression level.\n");
    fprintf(fp, "    -0 or -u       No compression.\n");
    //fprintf(fp, "    -v             Verbose output.\n");
    fprintf(fp, "    -H             [SAM] Do not print header\n");
    fprintf(fp, "    -R range       [Cram] Specifies the refseq:start-end range\n");
    fprintf(fp, "    -r ref.fa      [Cram] Specifies the reference file.\n");
    fprintf(fp, "    -b integer     [Cram] Max. bases per slice, default %d.\n",
	    BASES_PER_SLICE);
    fprintf(fp, "    -s integer     [Cram] Sequences per slice, default %d.\n",
	    SEQS_PER_SLICE);
    fprintf(fp, "    -S integer     [Cram] Slices per container, default %d.\n",
	    SLICE_PER_CNT);
    fprintf(fp, "    -V version     [Cram] Specify the file format version to write (eg 1.1, 2.0)\n");
    fprintf(fp, "    -e             [Cram] Embed reference sequence.\n");
    fprintf(fp, "    -x             [Cram] Non-reference based encoding.\n");
    fprintf(fp, "    -M             [Cram] Use multiple references per slice.\n");
    fprintf(fp, "    -m             [Cram] Generate MD and NM tags.\n");
    fprintf(fp, "    -a             [Cram] Also compress using arithmetic coder (V3.1+).\n");
#ifdef HAVE_LIBBZ2
    fprintf(fp, "    -j             [Cram] Also compress using bzip2.\n");
#endif
#ifdef HAVE_LIBLZMA
    fprintf(fp, "    -Z             [Cram] Also compress using lzma.\n");
#endif
#ifdef HAVE_LIBBSC
    fprintf(fp, "    -J             [Cram] Also compression using libbsc (V3.1+)\n");
#endif
    fprintf(fp, "    -f             [Cram] Also compression using fqzcomp (V3.1+)\n");
    fprintf(fp, "    -T             [Cram] Also compression using name tokeniser (V3.1+)\n");
    fprintf(fp, "    -n             [Cram] Discard read names where possible.\n");
    fprintf(fp, "    -P             Preserve all aux tags (incl RG,NM,MD)\n");
    fprintf(fp, "    -p             Preserve aux tag sizes ('i', 's', 'c')\n");
    fprintf(fp, "    -q             Don't add scramble @PG header line\n");
    fprintf(fp, "    -N integer     Stop decoding after 'integer' sequences\n");
    fprintf(fp, "    -t N           Use N threads (availability varies by format)\n");
    fprintf(fp, "    -B             Enable Illumina 8 quality-binning system (lossy)\n");
    fprintf(fp, "    -!             Disable all checking of checksums\n");
    fprintf(fp, "    -g FILE        Convert to Bam using index (file.gzi)\n");
    fprintf(fp, "    -G FILE        Output Bam index when bam input(file.gzi)\n");
    fprintf(fp, "    -X mode        [Cram] Mode is fast, default, small or archive.\n");
}

int main(int argc, char **argv) {
    scram_fd *in, *out;
    bam_seq_t *s;
    char imode[10], *in_f = "", omode[10], *out_f = "", *index_fn = NULL, *index_out_fn = NULL;
    int level = '\0'; // nul terminate string => auto level
    int c, verbose = 0;
    int s_opt = 0, S_opt = 0, embed_ref = 0, embed_cons = 0, ignore_md5 = 0, decode_md = 0;
    char *ref_fn = NULL;
    int start, end, multi_seq = -1, no_ref = 0;
    int use_bz2 = 0, use_bsc = 0, use_lzma = 0, use_fqz = 0, use_tok = 0, use_arith = 0;
    double vers = 3.0; // 3.0, 3.1, 4.0, etc
    char ref_name[1024] = {0};
    refs_t *refs;
    int nthreads = 1;
    t_pool *p = NULL;
    gzi *idx =NULL;
    int max_reads = -1;
    enum quality_binning binning = BINNING_NONE;
    int sam_fields = 0; // all
    int header = 1;
    int bases_per_slice = 0;
    int lossy_read_names = 0;
    int preserve_aux_order = 0;
    int preserve_aux_size = 0;
    int add_pg = 1;

    scram_init();

    /* Parse command line arguments */
    while ((c = getopt(argc, argv, "u0123456789hvs:S:V:r:xeEI:O:R:!MmajJZt:BN:F:Hb:nPpqg:G:fTX:")) != -1) {
	switch (c) {
	case 'X':
	    if (strcmp(optarg, "default") == 0) {
		// nothing
	    } else if (strcmp(optarg, "fast") == 0) {
		level = '1';
		s_opt = 1000;
	    } else if (strcmp(optarg, "small") == 0) {
		level = '7';
		if (vers >= 3.099)
		    use_arith = use_bz2 = use_tok = 1;
		else
		    use_bz2 = 1;
	    } else if (strcmp(optarg, "archive") == 0) {
		level = '7';
		use_arith = 1;
		if (vers >= 3.099)
		    use_bz2 = use_fqz = use_tok = 1;
		else
		    use_bz2 = use_lzma = 1;
		s_opt = 100000;
	    } else {
		fprintf(stderr, "Unknown parameter set: choose 'fast', 'default' or 'archive'\n");
		fprintf(stderr, "Assuming default\n");
	    }
	    break;

	case 'F':
	    sam_fields = strtol(optarg, NULL, 0); // undocumented for testing
	    break;

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

	case 'H':
	    header = 0;
	    break;

	case 'v':
	    verbose++;
	    break;

	case 's':
	    s_opt = atoi(optarg);
	    bases_per_slice = s_opt * 500; // guesswork...
	    break;

	case 'b':
	    bases_per_slice = atoi(optarg);
	    break;

	case 'S':
	    S_opt = atoi(optarg);
	    break;

	case 'm':
	    decode_md = 1;
	    break;

	case 'V':
	    vers = atof(optarg);
	    if (cram_set_option(NULL, CRAM_OPT_VERSION, optarg))
		return 1;
	    break;

	case 'r':
	    ref_fn = optarg;
	    break;

	case 'e':
	    embed_ref = 1;
	    break;
	case 'E':
	    embed_cons = 1;
	    break;

	case 'x':
	    no_ref = 1;
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

	case 'n':
	    lossy_read_names = 1;
	    break;

	case 'M':
	    multi_seq = 1;
	    break;

	case 'a':
	    use_arith = 1;
	    break;

	case 'j':
#ifdef HAVE_LIBBZ2
	    use_bz2 = 1;
#else
	    fprintf(stderr, "Warning: bzip2 support is not compiled into this"
		    " version.\nPlease recompile.\n");
#endif
	    break;

#ifdef HAVE_LIBBSC
	case 'J':
	    use_bsc = 1;
	    break;
#else
	    fprintf(stderr, "Warning: bsc support is not compiled into this"
		    " version.\nPlease recompile.\n");
#endif

	case 'Z':
#ifdef HAVE_LIBLZMA
	    use_lzma = 1;
#else
	    fprintf(stderr, "Warning: lzma support is not compiled into this"
		    " version.\nPlease recompile.\n");
#endif
	    break;

	case 'f':
#ifdef HAVE_FQZ
	    use_fqz = 1;
#else
	    fprintf(stderr, "Warning: FQZ support is not compiled into this"
		    " version.\nPlease recompile.\n");
#endif
	    break;

	case 'T':
	    use_tok = 1;
	    break;

	case 't':
	    nthreads = atoi(optarg);
	    if (nthreads < 1) {
		fprintf(stderr, "Number of threads needs to be >= 1\n");
		return 1;
	    }
	    break;

	case 'B':
	    binning = BINNING_ILLUMINA;
	    break;

	case 'P':
	    preserve_aux_order = 1;
	    break;

	case 'p':
	    preserve_aux_size = 1;
	    break;

	case 'q':
	    add_pg = 0;
	    break;

	case 'N':
	    max_reads = atoi(optarg);
	    break;

	case 'g':
	    index_fn = optarg;
	    break;

	case 'G':
	    index_out_fn = optarg;
	    break;

	case '?':
	    fprintf(stderr, "Unrecognised option: -%c\n", optopt);
	    usage(stderr);
	    return 1;
	}
    }    

    if (cram_default_version() <= 300 && (use_bsc || use_fqz)) {
	fprintf(stderr, "Libbsc and/or fqzcomp codecs are only permitted in CRAM v3.1 and 4.0.\n"
		"Note these CRAM versions are a technology demonstration only.\n"
		"Future versions of Scramble may not be able to read these files.\n");
	return 1;
    }

    if (cram_default_version() > 300) {
	fprintf(stderr, "\nWARNING: this version of CRAM is not a recognised GA4GH standard.\n"
		"Note this CRAM version is a technology demonstration only.\n"
		"Future versions of Scramble may not be able to read these files.\n\n");
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
	    fprintf(stderr, "Failed to open file %s\n", argv[optind]);
	    return 1;
	}
    } else {
        if (isatty(0)) {
	    usage(stdout);
	    return 0;
	}

	if (!(in = scram_open("-", imode))) {
	    fprintf(stderr, "Failed to open file %s\n", argv[optind]);
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
	    fprintf(stderr, "Failed to open file %s\n", argv[optind+1]);
	    return 1;
	}
    } else {
	if (!(out = scram_open("-", omode))) {
	    fprintf(stderr, "Failed to open file %s\n", argv[optind+1]);
	    return 1;
	}
    }


    /* Set any format specific options */
    scram_set_refs(out, refs = scram_get_refs(in));

    scram_set_option(out, CRAM_OPT_VERBOSITY, verbose);
    if (s_opt)
	if (scram_set_option(out, CRAM_OPT_SEQS_PER_SLICE, s_opt))
	    return 1;

    if (S_opt)
	if (scram_set_option(out, CRAM_OPT_SLICES_PER_CONTAINER, S_opt))
	    return 1;

    if (bases_per_slice)
	if (scram_set_option(out, CRAM_OPT_BASES_PER_SLICE, bases_per_slice))
	    return 1;

    if (embed_ref) {
	if (scram_get_header(in)->sort_order == ORDER_NAME ||
	    scram_get_header(in)->sort_order == ORDER_UNSORTED) {
	    fprintf(stderr, "Embeded reference with non-coordinate sorted data is "
		    "not supported.\nUsing -x for no-ref instead.\n");
	    if (scram_set_option(out, CRAM_OPT_NO_REF, 1))
		return 1;
	} else {
	    if (scram_set_option(out, CRAM_OPT_EMBED_REF, embed_ref))
		return 1;
	}
    }

    if (embed_cons) {
	if (scram_get_header(in)->sort_order == ORDER_NAME ||
	    scram_get_header(in)->sort_order == ORDER_UNSORTED) {
	    fprintf(stderr, "Embeded consensus with non-coordinate sorted data is "
		    "not supported.\n");
	} else {
	    if (scram_set_option(out, CRAM_OPT_EMBED_CONS, embed_cons))
		return 1;
	}
    }

    if (use_bz2)
	if (scram_set_option(out, CRAM_OPT_USE_BZIP2, use_bz2))
	    return 1;

    if (use_bsc)
	if (scram_set_option(out, CRAM_OPT_USE_BSC, use_bsc))
	    return 1;

    if (use_lzma)
	if (scram_set_option(out, CRAM_OPT_USE_LZMA, use_lzma))
	    return 1;

    if (use_fqz)
	if (scram_set_option(out, CRAM_OPT_USE_FQZ, use_fqz))
	    return 1;

    if (use_tok)
	if (scram_set_option(out, CRAM_OPT_USE_TOK, use_tok))
	    return 1;

    if (use_arith)
	if (scram_set_option(out, CRAM_OPT_USE_ARITH, use_arith))
	    return 1;

    if (binning != BINNING_NONE)
	if (scram_set_option(out, CRAM_OPT_BINNING, binning))
	    return 1;

    if (no_ref)
	if (scram_set_option(out, CRAM_OPT_NO_REF, no_ref))
	    return 1;

    if (multi_seq)
	if (scram_set_option(out, CRAM_OPT_MULTI_SEQ_PER_SLICE, multi_seq))
	    return 1;

    if (decode_md) {
	if (no_ref) {
	    fprintf(stderr, "Cannot use -m in conjunction with -x.\n");
	    return 1;
	}
	if (scram_set_option(in, CRAM_OPT_DECODE_MD, decode_md))
	    return 1;
    }

    if (index_fn) {
	if (NULL == (idx = gzi_index_load(index_fn))) {
	    fprintf(stderr, "Cannot open index file.\n");
	    return 1;
	}
	if (scram_set_option(out, CRAM_OPT_WITH_BGZIP_INDEX, idx))
	    return 1;
    }

    if (index_out_fn) {
	if (scram_set_option(in, CRAM_OPT_OUTPUT_BGZIP_IDX, index_out_fn))
	    return 1;
    }

    if (nthreads > 1) {
	if (NULL == (p = t_pool_init(nthreads*2, nthreads)))
	    return 1;

	if (scram_set_option(in,  CRAM_OPT_THREAD_POOL, p))
	    return 1;
	if (scram_set_option(out, CRAM_OPT_THREAD_POOL, p))
	    return 1;
    }

    if (ignore_md5) {
	if (scram_set_option(in, CRAM_OPT_IGNORE_MD5, ignore_md5))
	    return 1;
	if (scram_set_option(in, CRAM_OPT_IGNORE_CHKSUM, ignore_md5))
	    return 1;
	if (scram_set_option(out, CRAM_OPT_IGNORE_CHKSUM, ignore_md5))
	    return 1;
    }
    
    if (lossy_read_names) {
	if (scram_set_option(out, CRAM_OPT_LOSSY_READ_NAMES, lossy_read_names))
	    return 1;
    }

    if (preserve_aux_order)
	if (scram_set_option(out, CRAM_OPT_PRESERVE_AUX_ORDER, preserve_aux_order))
	    return 1;

    if (preserve_aux_size)
	if (scram_set_option(out, CRAM_OPT_PRESERVE_AUX_SIZE, preserve_aux_size))
	    return 1;

    if (sam_fields)
	scram_set_option(in, CRAM_OPT_REQUIRED_FIELDS, sam_fields);

    /* Copy header and refs from in to out, for writing purposes */
    scram_set_header(out, scram_get_header(in));

    // Needs doing after loading the header.
    if (ref_fn) {
	if (scram_set_option(out, CRAM_OPT_REFERENCE, ref_fn))
	    return 1;
    } else {
	// Attempt to fill out a cram->refs[] array from @SQ headers
	scram_set_option(out, CRAM_OPT_REFERENCE, NULL);
    }

    if (scram_get_header(out)) {
        if (add_pg) {
	    char *arg_list = stringify_argv(argc, argv);

	    if (!arg_list)
		return 1;

	
	    if (sam_hdr_add_PG(scram_get_header(out), "scramble",
			       "VN", PACKAGE_VERSION,
			       "CL", arg_list, NULL))
	        return 1;

	    free(arg_list);
	}

	if ((header || (omode[1] != 's' && omode[1] != '\0')) && scram_write_header(out) != 0)
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

	refid = sam_hdr_name2ref(in->c->header, ref_name);

	if (refid == -1 && *ref_name != '*') {
	    fprintf(stderr, "Unknown reference name '%s'\n", ref_name);
	    return 1;
	}
	r.refid = refid;
	r.start = start;
	r.end = end;
	if (scram_set_option(in, CRAM_OPT_RANGE, &r))
	    return 1;
    }

    /* Do the actual file format conversion */
    s = NULL;

    while (scram_get_seq(in, &s) >= 0) {
	if (-1 == scram_put_seq(out, s)) {
	    fprintf(stderr, "Failed to encode sequence\n");
	    return 1;
	}
	if (max_reads >= 0)
	    if (--max_reads == 0)
		break;
    }

    switch(scram_eof(in)) {
    case -1:
	fprintf(stderr, "Failed to decode sequence\n");
	return 1;
    case 0:
	if (max_reads == -1) {
	    fprintf(stderr, "Failed to decode sequence\n");
	    return 1;
	} else {
	    break;
	}
    case 2:
	fprintf(stderr, "Warning: no end-of-file block identified. "
		"File may be truncated.\n");
	break;
    case 1: default:
	// expected case
	break;
    }

    /* Finally tidy up and close files */
    if (scram_close(in)) {
	fprintf(stderr, "Failed in scram_close(in)\n");
	return 1;
    }
    if (scram_close(out)) {
	fprintf(stderr, "Failed in scram_close(out)\n");
	return 1;
    }

    if (p)
	t_pool_destroy(p, 0);

    if (s)
	free(s);

    return 0;
}

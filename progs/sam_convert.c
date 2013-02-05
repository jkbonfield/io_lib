#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <fcntl.h>
#include <zlib.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>

#include <io_lib/bam.h>
#include <io_lib/os.h>

int main(int argc, char **argv) {
    bam_file_t *in, *out;
    bam_seq_t *s;
    int nr = 0;
    int nm[2];
    char *mode = "rb";


    if (argc != 3) {
	fprintf(stderr, "Usage: sam_convert input_file output_file\n");
	fprintf(stderr, "Input and output files should end in .sam or .bam.\n");
	return 1;
    }
    
    {
	size_t l = strlen(argv[1]);
	if (l >= 4 && argv[1][l-3] == 's')
	    mode = "r";
    }

    in = bam_open(argv[1], mode);
    if (!in) {
	fprintf(stderr, "Failed to open bam file %s\n", argv[1]);
	return 1;
    }

    {
	size_t l = strlen(argv[2]);
	if (l >= 4 && argv[2][l-3] == 's')
	    mode = "w";
	else
	    mode = "wb";
    }

    out = bam_open(argv[2], mode);
    if (!out) {
	fprintf(stderr, "Failed to open bam file %s\n", argv[2]);
	return 1;
    }

    /* Copy header and refs from in to out, for writing purposes */
    out->header     = in->header;
    out->header_len = in->header_len;
    out->ref        = in->ref;
    out->nref       = in->nref;

    if (in->header) {
	if (bam_write_header(out))
	    return 1;
    }

    nm[0] = nm[1] = 0;
    s = NULL;
    while (bam_next_seq(in, &s) > 0) {
	out->ref = in->ref; // tmp hack
	if (-1 == bam_put_seq(out, s))
	    return 1;
	nr++;
	nm[(bam_flag(s) & BAM_FUNMAP) > 0]++;
    }

    printf("Mapped      = %d\n",nm[0]);
    printf("Unmpped     = %d\n",nm[1]);
    printf("Total reads = %d\n",nr);

    out->header = NULL;
    out->ref    = NULL;

    bam_close(in);
    bam_close(out);

    return 0;
}

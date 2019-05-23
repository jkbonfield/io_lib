#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <ctype.h>

#include "htscodecs/fqzcomp_qual.h"

#ifndef MAX_REC
#define MAX_REC 1000000
#endif

static fqz_slice fixed_slice = {0};

fqz_slice *fake_slice(size_t buf_len, int *len, int *r2, int *sel, int nlen) {
    fixed_slice.num_records = (nlen == 1) ? (buf_len+len[0]-1) / len[0] : nlen;
    assert(fixed_slice.num_records <= MAX_REC);
    int i, tlen = 0;
    if (!fixed_slice.crecs)
	fixed_slice.crecs = malloc(MAX_REC * sizeof(*fixed_slice.crecs));
    for (i = 0; i < fixed_slice.num_records; i++) {
	int idx = i < nlen ? i : nlen-1;
	fixed_slice.crecs[i].len = len[idx];
	fixed_slice.crecs[i].qual = tlen;
	fixed_slice.crecs[i].flags = r2 ? r2[idx]*FQZ_FREAD2 : 0;
	fixed_slice.crecs[i].flags |= sel ? (sel[idx]<<16) : 0;
	tlen += len[idx];
    }

    return &fixed_slice;
}

#define BS 1024*1024
static unsigned char *load(char *fn, size_t *lenp) {
    unsigned char *data = NULL;
    uint64_t dsize = 0;
    uint64_t dcurr = 0;
    signed int len;

    //build_rcp_freq();

    int fd = open(fn, O_RDONLY);
    if (!fd) {
	perror(fn);
	return NULL;
    }

    do {
	if (dsize - dcurr < BS) {
	    dsize = dsize ? dsize * 2 : BS;
	    data = realloc(data, dsize);
	}

	len = read(fd, data + dcurr, BS);
	if (len > 0)
	    dcurr += len;
    } while (len > 0);

    if (len == -1) {
	perror("read");
    }
    close(fd);

    *lenp = dcurr;
    return data;
}

#define BLK_SIZE 300*1000000
//#define BLK_SIZE 100*100000

int count_lines(unsigned char *in, size_t len) {
    size_t i;
    int lines = 0;

    for (i = 0; i < len; i++)
	if (in[i] == '\n')
	    lines++;

    return lines;
}

// QUAL [is_read2 [selector]]
void parse_lines(unsigned char *in, size_t len,
		 int *rec_len, int *rec_r2, int *rec_sel, size_t *new_len) {
    size_t i, j, start;
    int rec = 0;

    for (start = i = j = 0; i < len; i++) {
	if (in[i] == '\n' || in[i] == ' ' || in[i] == '\t') {
	    rec_len[rec] = i-start;

	    // Read2 marker
	    while (i < len && in[i] != '\n' && isspace(in[i]))
		i++;

	    if (in[i] != '\n')
		rec_r2[rec] = atoi((char *)&in[i]);
	    else
		rec_r2[rec] = 0;

	    while (i < len && !isspace(in[i]))
		i++;

	    // selector
	    while (i < len && in[i] != '\n' && isspace(in[i]))
		i++;

	    if (in[i] != '\n')
		rec_sel[rec] = atoi((char *)&in[i]);
	    else
		rec_sel[rec] = 0;

	    while (i < len && in[i] != '\n')
		i++;

	    start = i+1;
	    rec++;
	} else {
	    in[j++] = in[i]-33; // ASCII phred to qual
	}
    }
    *new_len = j;
}

int main(int argc, char **argv) {
    unsigned char *in, *out;
    size_t in_len, out_len;
    int decomp = 0, vers = 4;

    while (argc > 1 && argv[1][0] == '-') {
	if (argc > 1 && strcmp(argv[1], "-d") == 0) {
	    decomp = 1;
	    argv++;
	    argc--;
	}
	if (argc > 2 && strcmp(argv[1], "-s") == 0) {
	    vers += atoi(argv[2])*256;
	    argv+=2;
	    argc-=2;
	}
    }
    in = load(argc > 1 ? argv[1] : "/dev/stdin", &in_len);
    if (!in)
	exit(1);

    int blk_size = BLK_SIZE; // MAX
    if (argc > 3)
	blk_size = atoi(argv[3]);
    if (blk_size > BLK_SIZE)
	blk_size = BLK_SIZE;

    if (decomp) {
	unsigned char *in2 = in;
	while (in_len > 0) {
	    // Read sizes as 32-bit
	    size_t out_len = *(uint32_t *)in2;  in2 += 4;
	    size_t in2_len = *(uint32_t *)in2;  in2 += 4;

	    fprintf(stderr, "out_len %ld, in_len %ld\n", out_len, in2_len);

	    int *lengths = malloc(MAX_REC * sizeof(int));
	    out = (unsigned char *)fqz_decompress((char *)in2, in_len-8, &out_len, lengths);

	    // Convert from binary back to ASCII with newlines
	    int i = 0, j = 0;
	    while (j < out_len) {
		int k;
		char seq[MAX_SEQ];
		for (k = 0; k < lengths[i]; k++)
		    seq[k] = out[j+k]+33;
		seq[k] = 0;
		puts(seq);
		j += lengths[i++];
	    }
	    free(out);
	    in2 += in2_len;
	    in_len -= in2_len+8;

	    free(lengths);

	    break; // One cycle only until we fix blocking to be \n based
	}
    } else {
	// Convert from ASCII newline separated file to binary block.
	// We return an array of line lengths and optionally param selectors.
	int nlines = count_lines(in, in_len);
	fprintf(stderr, "nlines=%d\n", nlines);
	int *rec_len = calloc(nlines, sizeof(*rec_len));
	int *rec_r2  = calloc(nlines, sizeof(*rec_r2));
	int *rec_sel = calloc(nlines, sizeof(*rec_sel));
	parse_lines(in, in_len, rec_len, rec_r2, rec_sel, &in_len);

	unsigned char *in2 = in;
	long t_out = 0;
	out = NULL;
	while (in_len > 0) {
	    // FIXME: blk_size no longer working in test.  One cycle only!
	    size_t in2_len = in_len <= blk_size ? in_len : blk_size;
	    fqz_slice *s = fake_slice(in2_len, rec_len, rec_r2, rec_sel, nlines);
	    out = (unsigned char *)fqz_compress(vers, s, (char *)in2, in2_len, &out_len, 0);

	    // Write out 32-bit sizes.
	    uint32_t u32;
	    u32 = in2_len; if (write(1, &u32, 4) != 4) return 1;
	    u32 = out_len; if (write(1, &u32, 4) != 4) return 1;
	    if (write(1, out, out_len) < 0) return 1;
	    in_len -= in2_len;
	    in2 += in2_len;
	    t_out += out_len+16;

	    break; // One cycle only until we fix blocking to be \n based
	}
	free(out);
	free(rec_len);
	free(rec_r2);
	free(rec_sel);
	fprintf(stderr, "Total output = %ld\n", t_out);
    }

    free(in);

    return 0;
}

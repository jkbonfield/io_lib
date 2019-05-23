#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>

#include "htscodecs/rANS_static.h"

#ifndef BLK_SIZE
#  define BLK_SIZE 1024*1024
#endif

// Room to allow for expanded BLK_SIZE on worst case compression.
#define BLK_SIZE2 ((int)(1.05*BLK_SIZE))

/*-----------------------------------------------------------------------------
 * Main.
 *
 * This is a simple command line tool for testing order-0 and order-1
 * compression using the rANS codec. Simply compile with
 *
 * gcc -DTEST_MAIN -O3 -I. cram/rANS_static.c -o cram/rANS_static
 *
 * Usage: cram/rANS_static -o0 < file    > file.o0
 *        cram/rANS_static -d  < file.o0 > file2
 *
 *        cram/rANS_static -o1 < file    > file.o1
 *        cram/rANS_static -d  < file.o1 > file2
 */
int main(int argc, char **argv) {
    int opt, order = 0;
    unsigned char in_buf[BLK_SIZE2+257*257*3];
    int decode = 0, test = 0;
    FILE *infp = stdin, *outfp = stdout;
    struct timeval tv1, tv2, tv3;
    size_t bytes = 0;

    extern char *optarg;
    extern int optind;

    while ((opt = getopt(argc, argv, "o:dt")) != -1) {
	switch (opt) {
	case 'o':
	    order = atoi(optarg);
	    break;

	case 'd':
	    decode = 1;
	    break;
	    
	case 't':
	    test = 1;
	    break;
	}
    }

    order = order ? 1 : 0; // Only support O(0) and O(1)

    if (optind < argc) {
	if (!(infp = fopen(argv[optind], "rb"))) {
	    perror(argv[optind]);
	    return 1;
	}
	optind++;
    }

    if (optind < argc) {
	if (!(outfp = fopen(argv[optind], "wb"))) {
	    perror(argv[optind]);
	    return 1;
	}
	optind++;
    }

    gettimeofday(&tv1, NULL);

    if (test) {
	size_t len, in_sz = 0, out_sz = 0;
	typedef struct {
	    unsigned char *blk;
	    uint32_t sz;
	} blocks;
	blocks *b = NULL, *bc, *bu;
	int nb = 0, i;
	
	while ((len = fread(in_buf, 1, BLK_SIZE, infp)) != 0) {
	    // inefficient, but it'll do for testing
	    b = realloc(b, (nb+1)*sizeof(*b));
	    b[nb].blk = malloc(len);
	    b[nb].sz = len;
	    memcpy(b[nb].blk, in_buf, len);
	    nb++;
	    in_sz += len;
	}

	int trials = 10;
	while (trials--) {
	    bc = malloc(nb*sizeof(*bc));
	    bu = malloc(nb*sizeof(*bu));

	    gettimeofday(&tv1, NULL);

	    out_sz = 0;
	    for (i = 0; i < nb; i++) {
		bc[i].blk = rans_compress(b[i].blk, b[i].sz, &bc[i].sz, order);
		out_sz += 5 + bc[i].sz;
		bc[i].blk = realloc(bc[i].blk, bc[i].sz);
	    }
	
	    gettimeofday(&tv2, NULL);

	    for (i = 0; i < nb; i++) {
		bu[i].blk = rans_uncompress(bc[i].blk, bc[i].sz, &bu[i].sz);
	    }

	    gettimeofday(&tv3, NULL);

	    for (i = 0; i < nb; i++) {
		if (b[i].sz != bu[i].sz || memcmp(b[i].blk, bu[i].blk, b[i].sz))
		    fprintf(stderr, "Mismatch in block %d\n", i);
		free(bc[i].blk);
		free(bu[i].blk);
	    }
	    free(bc);
	    free(bu);

	    fprintf(stderr, "%5.1f MB/s enc, %5.1f MB/s dec\t %lld bytes -> %lld bytes\n",
		    (double)in_sz / ((long)(tv2.tv_sec - tv1.tv_sec)*1000000 +
				     tv2.tv_usec - tv1.tv_usec),
		    (double)in_sz / ((long)(tv3.tv_sec - tv2.tv_sec)*1000000 +
				     tv3.tv_usec - tv2.tv_usec),
		    (long long)in_sz, (long long)out_sz);
	}

	exit(0);
	
    }

    if (decode) {
	// Only used in some test implementations of RC_GetFreq()
	//RC_init();
	//RC_init2();

	for (;;) {
	    uint32_t in_size, out_size;
	    unsigned char *out;

	    order = fgetc(infp);
	    if (4 != fread(&in_size, 1, 4, infp))
		break;
	    if (in_size != fread(in_buf, 1, in_size, infp)) {
		fprintf(stderr, "Truncated input\n");
		exit(1);
	    }
	    out = rans_uncompress(in_buf, in_size, &out_size);
	    if (!out)
		abort();

	    fwrite(out, 1, out_size, outfp);
	    free(out);

	    bytes += out_size;
	}
    } else {
	for (;;) {
	    uint32_t in_size, out_size;
	    unsigned char *out;

	    in_size = fread(in_buf, 1, BLK_SIZE, infp);
	    if (in_size <= 0)
		break;

	    out = rans_compress(in_buf, in_size, &out_size,
				order && in_size >= 4);

	    fputc(order && in_size >= 4, outfp);
	    fwrite(&out_size, 1, 4, outfp);
	    fwrite(out, 1, out_size, outfp);
	    free(out);

	    bytes += in_size;
	}
    }

    gettimeofday(&tv2, NULL);

    fprintf(stderr, "Took %ld microseconds, %5.1f MB/s\n",
	    (long)(tv2.tv_sec - tv1.tv_sec)*1000000 +
	    tv2.tv_usec - tv1.tv_usec,
	    (double)bytes / ((long)(tv2.tv_sec - tv1.tv_sec)*1000000 +
			     tv2.tv_usec - tv1.tv_usec));
    return 0;
}

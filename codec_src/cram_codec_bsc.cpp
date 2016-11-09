//g++ -g -fPIC -shared -fopenmp cram_compress_bsc.cpp  -o libcram_compress_bsc.so -I.. -I. -I$HOME/ftp/compression/libbsc/ -L$HOME/ftp/compression/libbsc/ -lbsc

#include <stdio.h>
#include <stdlib.h>
#include <cstdint>
#include <limits>
#include "io_lib/cram_block_compression.h"

#include <libbsc/libbsc.h>

#define FEATURES LIBBSC_FEATURE_FASTMODE

static const char *name(void) {
    return "LibBSC compression";
}

unsigned char *compress_block(int level,
			      cram_slice *s,
			      unsigned char *in,
			      size_t in_size,
			      size_t *out_size) {
    unsigned char *comp = (unsigned char *)malloc(in_size+LIBBSC_HEADER_SIZE);

    *out_size = bsc_compress(in, comp, in_size,
			     0, // LZP hash size, 0 or 10..28
			     0, // LZP min match, 0 or 4..255
			     LIBBSC_BLOCKSORTER_BWT,
			     LIBBSC_CODER_QLFC_STATIC,
			     //LIBBSC_CODER_QLFC_ADAPTIVE, // maybe 50% slower?
			     FEATURES);

    if (*out_size <= 0) {
	free(comp);
	return NULL;
    }

    return comp;
}

unsigned char *uncompress_block(cram_slice *s,
				unsigned char *in,
				size_t in_size,
				size_t *out_size) {
    int block_size, data_size;
    unsigned char *uncomp;

    if (bsc_block_info(in, LIBBSC_HEADER_SIZE,
		       &block_size, &data_size,
		       FEATURES) != LIBBSC_NO_ERROR)
	return NULL;

    if (!(uncomp = (unsigned char *)malloc(data_size)))
	return NULL;

    if (bsc_decompress(in, block_size, uncomp, data_size,
		       FEATURES) != LIBBSC_NO_ERROR) {
	free(uncomp);
	return NULL;
    }

    *out_size = data_size;

    return uncomp;
}

static cram_compressor c = {
    'b',//FOUR_CC("\0BSC"),
    0, // all data series
    1.1,
    name,
    compress_block,
    uncompress_block,
};

/*
#define LIBBSC_DEFAULT_LZPHASHSIZE     16
#define LIBBSC_DEFAULT_LZPMINLEN       128
#define LIBBSC_DEFAULT_BLOCKSORTER     LIBBSC_BLOCKSORTER_BWT
#define LIBBSC_DEFAULT_FEATURES        LIBBSC_FEATURE_MULTITHREADING
*/

extern "C" {
cram_compressor *cram_compressor_init(void) {
    bsc_init(FEATURES);
    return &c;
}
}

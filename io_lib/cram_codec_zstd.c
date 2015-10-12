//gcc -g -fPIC -shared cram_codec_zstd.c -o libcram_codec_zstd.so -I.. -I. -I/nfs/users/nfs_j/jkb/ftp/compression/zstd/lib -L/nfs/users/nfs_j/jkb/ftp/compression/zstd/lib -lzstd -Wl,-rpath,/nfs/users/nfs_j/jkb/ftp/compression/zstd/lib

#include <stdio.h>
#include <stdlib.h>
#include "io_lib/cram_block_compression.h"

#include <zstd.h>

static const char *name(void) {
    return "ZSTD compression";
}

unsigned char *compress_block(int level,
			      cram_slice *s,
			      unsigned char *in,
			      size_t in_size,
			      size_t *out_size) {
    unsigned char *comp = (unsigned char *)malloc(ZSTD_compressBound(in_size));

    *out_size = ZSTD_compress(comp, ZSTD_compressBound(in_size),
			      in, in_size);
    if (*out_size <= 0) {
	free(comp);
	return NULL;
    }

    return comp;
}

unsigned char *uncompress_block(unsigned char *in,
				// cram_slice *s,
				size_t in_size,
				size_t *out_size) {
    int block_size, data_size;
    unsigned char *uncomp = malloc(*out_size);

    *out_size = ZSTD_decompress(uncomp, *out_size,
				in, in_size);

    return uncomp;
}

static cram_compressor c = {
    'x', //FOUR_CC("ZSTD"),
    0, // all data series
    1.0,
    name,
    compress_block,
    uncompress_block,
};

cram_compressor *cram_compressor_init(void) {
    return &c;
}



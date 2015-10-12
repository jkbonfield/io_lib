//gcc -g -fPIC -shared cram_codec_squash.c -o libcram_codec_squash.so -I.. -I. -I/nfs/sam_scratch/jkb/opt/squash/include/squash-0.7 -L/nfs/sam_scratch/jkb/opt/squash/lib/ -lsquash0.7 -Wl,-rpath,/nfs/sam_scratch/jkb/opt/squash/lib/

#include <stdio.h>
#include <stdlib.h>
#include "io_lib/cram_block_compression.h"

#include <squash/squash.h>

static const char *name(void) {
  return "Squash generic compression";
}

#ifndef CODEC_NAME
#  define CODEC_NAME "brotli"
#endif

#ifndef CODEC_LEVEL
#  define CODEC_LEVEL "11"
#endif

unsigned char *compress_block(int level,
			      cram_slice *s,
			      unsigned char *in,
			      size_t in_size,
			      size_t *out_size) {
  uint8_t *compressed;

  *out_size = squash_get_max_compressed_size(CODEC_NAME, in_size);
  compressed = (uint8_t *)malloc(*out_size);

  SquashStatus res =
    squash_compress (CODEC_NAME,
		     out_size, compressed,
		     in_size, (const uint8_t*)in,
		     "level", CODEC_LEVEL,
		     NULL);
    if (res != SQUASH_OK) {
      fprintf (stderr, "Unable to decompress data [%d]: %s\n",
	       res, squash_status_to_string (res));
      free(compressed);
      return NULL;
    }

  return compressed;
}

unsigned char *uncompress_block(unsigned char *in,
				// cram_slice *s,
				size_t in_size,
				size_t *out_size) {
  unsigned char *uncomp = malloc(*out_size);
  SquashStatus res =
    squash_decompress(CODEC_NAME, out_size, uncomp, in_size, in, NULL);
  if (res != SQUASH_OK) {
    fprintf (stderr, "Unable to decompress data [%d]: %s\n",
	     res, squash_status_to_string (res));
    exit (res);
  }

  return uncomp;
}

static cram_compressor c = {
    'Q', //FOUR_CC("SQSH"),
    0, // all data series
    1.0,
    name,
    compress_block,
    uncompress_block,
};

cram_compressor *cram_compressor_init(void) {
    return &c;
}



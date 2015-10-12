//gcc -g -fPIC -shared cram_codec_zstd.c -o libcram_codec_zstd.so -I.. -I. -I/nfs/users/nfs_j/jkb/ftp/compression/zstd/lib -L/nfs/users/nfs_j/jkb/ftp/compression/zstd/lib -lzstd -Wl,-rpath,/nfs/users/nfs_j/jkb/ftp/compression/zstd/lib

#include <stdio.h>
#include <stdlib.h>
#include "io_lib/cram_block_compression.h"
#include "io_lib/rANS_static.h"

static const char *name(void) {
    return "Qual NIB1";
}

unsigned char *compress_block(int level,
			      cram_slice *s,
			      unsigned char *in,
			      size_t in_size,
			      size_t *out_size) {
    int i_out_size = *out_size;
    unsigned char *comp;
    size_t i,j;

#if 1
    // Reverse complement
    int rec = 0;
    i = 0;
    while (i < in_size) {
	int len = s->crecs[rec].len;
	if (s->crecs[rec].flags & BAM_FREVERSE) {
	    int I, J;
	    unsigned char *cp = in+i;
	    for (I = 0, J = len-1; I < J; I++, J--) {
		unsigned char c;
		c = cp[I];
		cp[I] = cp[J];
		cp[J] = c;
	    }
	}
	rec++;
	i+=len;
    }
#endif

    // Map value to 0..n-1 for n discrete values
    int map[256] = {0}, map_size = 0;
    for (i = 0; i < in_size; i++)
	map[in[i]]++;
    
    for (i = j = 0; i < 256; i++) {
	if (map[i]) {
	    map_size++;
	    map[i]=j++;
	} else {
	    map[i]=-1;
	}
    }

    if (j <= 4) {
	for (i = j = 0; i < in_size; i+=4, j++) {
	    in[j] = (map[in[i+0<in_size?i+0:0]]<<0)
		  | (map[in[i+1<in_size?i+1:0]]<<2)
		  | (map[in[i+2<in_size?i+2:0]]<<4)
		  | (map[in[i+3<in_size?i+3:0]]<<6);
	}
	in_size = j;
    } else if (j <= 16) {
	for (i = j = 0; i < in_size; i+=2, j++) {
	    in[j] = (map[in[i+0<in_size?i+0:0]]<<0)
		  | (map[in[i+1<in_size?i+1:0]]<<4);
	}
	in_size = j;
    }

    // FIXME: write out map

    // FIXME: rewrite this to avoid memcpy
    comp = rans_compress(in, in_size, &i_out_size, 1);
    unsigned char *out_buf = malloc(i_out_size + 1 + map_size);
    unsigned char *cp = out_buf;
    *cp++ = map_size;
    for (i = 0; i < 256; i++) {
	if (map[i] >= 0) *cp++ = i;
    }
    memcpy(cp, comp, i_out_size);
    free(comp);

    *out_size = cp-out_buf + i_out_size;

    return out_buf;;
}

unsigned char *uncompress_block(unsigned char *in,
				size_t in_size,
				size_t *out_size) {
    int i_out_size = *out_size;
    unsigned char *uncomp;
    
    uncomp = rans_uncompress(in, in_size, &i_out_size, 0/*FIXME: unused*/);
    *out_size = i_out_size;

    return uncomp;
}

static cram_compressor c = {
    'n', //FOUR_CC("QNB1"),
    1<<DS_QS, // quality only
    1.0,
    name,
    compress_block,
    uncompress_block,
};

cram_compressor *cram_compressor_init(void) {
    return &c;
}

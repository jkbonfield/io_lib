// FIXME: encoder done as proof of concept to test sizing, but no decoder!

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io_lib/cram_block_compression.h"
#include "io_lib/rANS_static.h"

static const char *name(void) {
    return "Qual NIB1";
}

unsigned char *compress_block(int method,
			      int level,
			      cram_slice *s,
			      unsigned char *in,
			      size_t in_size,
			      size_t *out_size) {
    unsigned int i_out_size = *out_size;
    unsigned char *comp;
    size_t i,j;

#if 0
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

    unsigned char *in2;
    if (j <= 4) {
	in2 = malloc(in_size/4);
	for (i = j = 0; i < in_size; i+=4, j++) {
	    in2[j] = (map[in[i+0<in_size?i+0:0]]<<0)
		  | (map[in[i+1<in_size?i+1:0]]<<2)
		  | (map[in[i+2<in_size?i+2:0]]<<4)
		  | (map[in[i+3<in_size?i+3:0]]<<6);
	}
	in_size = j;
    } else if (j <= 16) {
	in2 = malloc(in_size/2);
	for (i = j = 0; i < in_size; i+=2, j++) {
	    in2[j] = (map[in[i+0<in_size?i+0:0]]<<0)
		  | (map[in[i+1<in_size?i+1:0]]<<4);
	}
	in_size = j;
	in = in2;
    } else {
	return NULL; // Forces use of another codec.
    }

    // FIXME: rewrite this to avoid memcpy
    comp = rans_compress(in2, in_size, &i_out_size, 1);
    unsigned char *out_buf = malloc(i_out_size + 1 + map_size);
    unsigned char *cp = out_buf;
    *cp++ = map_size;
    for (i = 0; i < 256; i++) {
	if (map[i] >= 0) *cp++ = i;
    }
    memcpy(cp, comp, i_out_size);
    free(comp);
    free(in2);

    *out_size = cp-out_buf + i_out_size;

    return out_buf;;
}

unsigned char *uncompress_block(cram_slice *s,
				unsigned char *in,
				size_t in_size,
				size_t *out_size) {
    unsigned int i_out_size = *out_size;
    unsigned char *uncomp, *out;
    
    int map[256] = {0}, map_size = *in++, i, j;
    for (i = 0; i < map_size; i++)
	map[i] = *in++;

    uncomp = rans_uncompress(in, in_size-map_size-1, &i_out_size, 0/*FIXME: unused*/);


    // FIXME: remainder when not multiple of 2 or 4
    if (map_size <= 4) {
	if (!(out = malloc(i_out_size*4))) {
	    free(uncomp);
	    return NULL;
	}
	for (i = j = 0; i < i_out_size; i++) {
	    out[j++] = map[(uncomp[i]>>0)&3];
	    out[j++] = map[(uncomp[i]>>2)&3];
	    out[j++] = map[(uncomp[i]>>4)&3];
	    out[j++] = map[(uncomp[i]>>6)&3];
	}
    } else if (map_size <= 16) {
	if (!(out = malloc(i_out_size*2))) {
	    free(uncomp);
	    return NULL;
	}
	for (i = j = 0; i < i_out_size; i++) {
	    out[j++] = map[(uncomp[i]>>0)&15];
	    out[j++] = map[(uncomp[i]>>4)&15];
	}
    } else {
	free(uncomp);
	return NULL;
    }

    *out_size = j;

#if 0
    // Reverse complement
    int rec = 0;
    i = 0;
    while (i < j) {
	int len = s->crecs[rec].len;
	if (s->crecs[rec].flags & BAM_FREVERSE) {
	    int I, J;
	    unsigned char *cp = out+i;
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

    return out;
}

static cram_compressor c = {
    'N', //FOUR_CC("QNB1"),
    1<<DS_QS, // quality only
    1.0,
    name,
    compress_block,
    uncompress_block,
};

cram_compressor *cram_compressor_init(void) {
    return &c;
}

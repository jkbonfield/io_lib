#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "io_lib/cram_block_compression.h"
#include "io_lib/rANS_static.h"

static const char *name(void) {
    return "Qual rans+dedup";
}

unsigned char *compress_block(int method,
			      int level,
			      cram_slice *s,
			      unsigned char *in,
			      size_t in_size,
			      size_t *out_size) {
    unsigned int i_out_size = *out_size;
    unsigned char *comp;
    size_t i,j, delta;

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

    // Dedup
    unsigned char *in2 = malloc(in_size);
    if (!in2)
	return NULL;
    if (1) {
	int rec = 0;
	int i0 = 0, len0 = 0;
	i = j = 0;
	int last_len = 0;
	int matched = 0;

	while (i < in_size) {
	    int len = s->crecs[rec].len;
	    if (len0 == len && len > 3 && memcmp(in+i0, in+i, len) == 0) {
		matched++;
		if (len == last_len) {
		    in2[j++] = 255;
		} else if (len < 256) {
		    in2[j++] = 254;
		    in2[j++] = len;
		} else if (len < 65536) {
		    in2[j++] = 253;
		    in2[j++] = len>>8;
		    in2[j++] = len;
		} else {
		    in2[j++] = 252;
		    in2[j++] = len>>24;
		    in2[j++] = len>>16;
		    in2[j++] = len>>8;
		    in2[j++] = len;
		}
		last_len = len;
	    } else {
		memcpy(in2+j, in+i, len);
		j += len;
		i0 = i;
		len0 = len;
	    }
	    i += len;
	}

	delta = in_size - j;
	in_size = j;
    }

    // FIXME: rewrite this to avoid memcpy
    comp = rans_compress(in2, in_size, &i_out_size, 1);
    unsigned char *out_buf = malloc(i_out_size + 5);
    unsigned char *cp = out_buf;
    do {
	*cp++ = (delta & 127) + 128*(delta >= 128);
	delta >>= 7;
    } while (delta > 0);
    memcpy(cp, comp, i_out_size);
    free(comp);
    free(in2);

    *out_size = cp-out_buf + i_out_size;

    //printf("%d -> %d\n", (int)in_size, (int)*out_size);

    return out_buf;;
}

unsigned char *uncompress_block(cram_slice *s,
				unsigned char *in,
				size_t in_size,
				size_t *out_size) {
    unsigned int i_out_size = *out_size, i, j;
    unsigned char *uncomp, *out, c;
    int delta = 0;
    
    int bn = 0;
    do {
	c = *in++;
	delta += ((size_t)(c & 0x7f)) << (bn*7);
	bn++;
    } while (c >= 128);

    uncomp = rans_uncompress(in, in_size-bn, &i_out_size, 0/*FIXME: unused*/);


    if (!(out = malloc(i_out_size + delta))) {
	free(uncomp);
	return NULL;
    }

    i = j = 0;
    int len = 0;
    while (j < i_out_size + delta) {
	if (uncomp[i] < 252) {
	    out[j++] = uncomp[i++];
	} else {
	    switch(uncomp[i]) {
	    case 254:
		len = uncomp[i+1];
		i += 2;
		break;
	    case 253:
		len = (uncomp[i+1]<<8) | uncomp[i+2];
		i += 3;
		break;
	    case 252:
		len = (uncomp[i+1]<<24) | (uncomp[i+2]<<16) | (uncomp[i+3]<<8) | uncomp[i+4];
		i += 5;
		break;
	    default:
		i++;
	    }
	    
	    memcpy(&out[j], &out[j-len], len);
	    j += len;
	}
    }

    *out_size = j;
    free(uncomp);

    return out;
}

static cram_compressor c = {
    'D', //FOUR_CC("DDR1"),
    1<<DS_QS, // quality only
    1.0,
    name,
    compress_block,
    uncompress_block,
};

cram_compressor *cram_compressor_init(void) {
    return &c;
}

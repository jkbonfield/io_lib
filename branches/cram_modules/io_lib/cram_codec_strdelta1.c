// cc -g -shared -fPIC -o cram_codec_strdelta1.so cram_codec_strdelta1.c -I/software/badger/include -L/software/badger/lib -lstaden-read

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>

#include "io_lib/hash_table.h"
#include "io_lib/rANS_static.h"
#include "io_lib/cram_block_compression.h"

// FIXME!
int back_buf[1000000] = {};
int prefix[1000000] = {};
int suffix[1000000] = {};

int back_ind = 0;
int pre_ind = 0;

char mid_str[50000000];
int mid_ind = 0;


static const char *name(void) {
    return "String-delta-1";
}

unsigned char *compress_block(int level,
			      unsigned char *data,
			      size_t len,
			      size_t *comp_len) {
    HashTable *h = HashTableCreate(0, HASH_DYNAMIC_SIZE | HASH_NONVOLATILE_KEYS);
    char *s, *last;
    int i, line, last_len;
    unsigned char *comp, *cp, *comp2, *cp2;
    uint32_t comp_len32;

    back_ind = pre_ind = mid_ind = 0;

    line = 0;
    last = NULL;
    for (i = 0; i < len; i++) {
	int j = i, new, name_len;
	HashData hd;
	HashItem *hi;
	unsigned char term_sym;

	while (i < len && data[i] > '\n')
	    i++;
	term_sym = data[i];
	
	name_len = i-j;

	hd.i = line;
	hi = HashTableAdd(h, &data[j], name_len, hd, &new);
	if (new) {
	    back_buf[back_ind++] = 0;
	    if (last) {
		int k = 0, l;
		while (k < name_len && last[k] == data[j+k])
		    k++;

#ifdef FIXED
		k=13;
#endif

		prefix[pre_ind] = k;

		l = name_len-1;
		last_len--;
		while (l > k && last_len >= 0 && last[last_len] == data[j+l])
		    l--, last_len--;
		l++;

#ifdef FIXED
		l=name_len-3;
#endif
		prefix[pre_ind] = k;
		suffix[pre_ind] = name_len-l;
	    } else {
		prefix[pre_ind] = 0;
		suffix[pre_ind] = 0;
	    }

//	    printf("%d %.*s %.*s %.*s %d\n",
//		   prefix[pre_ind], prefix[pre_ind], &data[j],
//		   name_len-suffix[pre_ind] - prefix[pre_ind], &data[j+prefix[pre_ind]],
//		   suffix[pre_ind], &data[j+name_len-suffix[pre_ind]],
//		   suffix[pre_ind]);

	    memcpy(&mid_str[mid_ind], &data[j+prefix[pre_ind]],
		   name_len-suffix[pre_ind] - prefix[pre_ind]);
	    mid_ind += name_len-suffix[pre_ind] - prefix[pre_ind];
	    mid_str[mid_ind++] = term_sym;

	    pre_ind++;
	} else {
	    back_buf[back_ind++] = line - hi->data.i;

	    // Odd  => dup of full name (/2).
	    // Even => prefix/mid/suffix of N-th name back (/2)
	    //back_buf[back_ind++] = (line - hi->data.i)*2+1;

	    //fprintf(stderr, "%d: BACK %d\n", line, (int)(line - hi->data.i));
	    hi->data.i = line;
	}

	line++;
	last = &data[j];
	last_len = name_len;
    }

    fprintf(stderr, "back_ind=%d pre_ind=%d mid_ind=%d\n",
	    back_ind, pre_ind, mid_ind);

    comp = malloc(4 + back_ind*2 + pre_ind*2*2 + mid_ind);
    cp2 = comp2 = malloc(back_ind*2 + pre_ind*2*2 + mid_ind);

    *cp2++ = (line>> 0) & 0xff;
    *cp2++ = (line>> 8) & 0xff;
    *cp2++ = (line>>16) & 0xff;
    *cp2++ = (line>>24) & 0xff;

    //freopen("_0", "w", stdout);
    cp = comp;
    for (i = 0; i < back_ind; i++) {
	int p = back_buf[i];
	do {
	    *cp++ = (p&0x7f) | (p>=128?128:0);
	    //putchar((p&0x7f) | (p>=128?128:0));
	    p >>= 7;
	} while (p);
    }

    cp = rans_compress(comp, cp-comp, &comp_len32, 0);
    *comp_len = comp_len32;
    memcpy(cp2, cp, *comp_len); cp2 += *comp_len;
    free(cp);

    //freopen("_1", "w", stdout);
    cp = comp;
    for (i = 0; i < pre_ind; i++) {
	int p = prefix[i];
	do {
	    *cp++ = (p&0x7f) | (p>=128?128:0);
	    //putchar((p&0x7f) | (p>=128?128:0));
	    p >>= 7;
	} while (p);
    }

    cp = rans_compress(comp, cp-comp, &comp_len32, 0);
    *comp_len = comp_len32;
    memcpy(cp2, cp, *comp_len); cp2 += *comp_len;
    free(cp);

    //freopen("_2", "w", stdout);
    //printf("%.*s", mid_ind, mid_str);
    cp = rans_compress(mid_str, mid_ind, &comp_len32, 1);
    *comp_len = comp_len32;
    memcpy(cp2, cp, *comp_len); cp2 += *comp_len;
    free(cp);

    //freopen("_3", "w", stdout);
    cp = comp;
    for (i = 0; i < pre_ind; i++) {
	int p = suffix[i];
	do {
	    *cp++ = (p&0x7f) | (p>=128?128:0);
	    //putchar((p&0x7f) | (p>=128?128:0));
	    p >>= 7;
	} while (p);
    }

    cp = rans_compress(comp, cp-comp, &comp_len32, 0);
    *comp_len = comp_len32;
    memcpy(cp2, cp, *comp_len); cp2 += *comp_len;
    free(cp);

    free(comp);

    *comp_len = cp2 - comp2;
    return comp2;
}

unsigned char *uncompress_block(unsigned char *data, size_t len,
				size_t *uncomp_len) {
    unsigned char *cp, *udata = NULL;
    int ulen, clen, udata_ind = 0, udata_sz = 0;
    unsigned char *back_buf, *pre_buf, *mid_buf, *suf_buf;
    unsigned char *back_cp, *pre_cp, *mid_cp, *suf_cp;
    int back_sz, pre_sz, mid_sz, suf_sz;
    int i, nlines;
    int *line_ind;

    cp = data;
    nlines = (cp[0]<<0) + (cp[1]<<8) + (cp[2]<<16) + (cp[3]<<24);
    cp += 4;
    fprintf(stderr, "nlines=%d\n", nlines);

    line_ind = calloc(nlines, 4);

    clen = (cp[1]<<0) + (cp[2]<<8) + (cp[3]<<16) + (cp[4]<<24);
    fprintf(stderr, "clen = %d\n", clen);
    back_cp = back_buf = rans_uncompress(cp, clen+9, &back_sz, 0);

    cp += clen+9;
    clen = (cp[1]<<0) + (cp[2]<<8) + (cp[3]<<16) + (cp[4]<<24);
    fprintf(stderr, "clen = %d\n", clen);
    pre_cp = pre_buf = rans_uncompress(cp, clen+9, &pre_sz, 0);

    cp += clen+9;
    clen = (cp[1]<<0) + (cp[2]<<8) + (cp[3]<<16) + (cp[4]<<24);
    fprintf(stderr, "clen = %d\n", clen);
    mid_cp = mid_buf = rans_uncompress(cp, clen+9, &mid_sz, 0);

    cp += clen+9;
    clen = (cp[1]<<0) + (cp[2]<<8) + (cp[3]<<16) + (cp[4]<<24);
    fprintf(stderr, "clen = %d\n", clen);
    suf_cp = suf_buf = rans_uncompress(cp, clen+9, &suf_sz, 0);

    for (i = 0; i < nlines; i++) {
	int t, x;

	line_ind[i] = udata_ind;
	if (udata_ind + 256 >= udata_sz) {
	    udata_sz = udata_sz ? udata_sz*2 : 65536;
	    udata = realloc(udata, udata_sz);
	}

	x = t = 0; do { t |= (*back_cp & 0x7f)<<x; x+=7; } while (*back_cp++ & 128);

	if (t) {
	    //fprintf(stderr, "back_cp=%d\n", t);
	    cp = &udata[line_ind[i-t]];
	    do {
		//putchar(*cp);
		udata[udata_ind++] = *cp;
	    } while (*cp++ > '\n');
	} else {
	    char *last = &udata[line_ind[i-1]];
	    int p, s;

	    x = p = 0; do { p |= (*pre_cp & 0x7f)<<x; x+=7; } while (*pre_cp++ & 128);
	    x = s = 0; do { s |= (*suf_cp & 0x7f)<<x; x+=7; } while (*suf_cp++ & 128);

	    while (p--) {
		//putchar(*last);
		udata[udata_ind++] = *last++;
	    }
	    //putchar(' ');

	    while (*mid_cp > '\n') {
		//putchar(*mid_cp);
		udata[udata_ind++] = *mid_cp++;
	    }
	    //putchar(' ');

	    if (s) {
		while (*last > '\n')
		    last++;

		last -= s;
		while (s--) {
		    //putchar(*last);
		    udata[udata_ind++] = *last++;
		}
	    }
	    //putchar('\n');

	    udata[udata_ind++] = *mid_cp++;
	}
    }
    
    *uncomp_len = udata_ind;
    return udata;
}

static cram_compressor c = {
    'd', //FOUR_CC("STRd"),
    1<<DS_RN, // quality only
    1.0,
    name,
    compress_block,
    uncompress_block,
};

cram_compressor *cram_compressor_init(void) {
    return &c;
}

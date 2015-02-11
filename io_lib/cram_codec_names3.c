#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <io_lib/hash_table.h>
#include <io_lib/string_alloc.h>
#include <io_lib/pooled_alloc.h>
#include <io_lib/cram_block_compression.h>
#include <bzlib.h>


unsigned char *rans_compress(unsigned char *in, unsigned int in_size,
			     unsigned int *out_size, int order);
unsigned char *rans_uncompress(unsigned char *in, unsigned int in_size,
			       unsigned int *out_size);

static const char *name(void) {
    return "Names3 compression";
}

typedef struct trie {
    char c;
    int count;
    struct trie *next[128];
    int sym;
} trie_t;

unsigned char *assign_trie_words(unsigned char trie_prefix[256],
				 unsigned char last_prefix[256],
				 int *trie_counter,
				 unsigned char *cp,
				 trie_t *last, trie_t *t,
				 int depth, int min_val) {
    int j;
    for (j = 0; j < 128; j++) {
	if (!t->next[j])
	    continue;

	trie_prefix[depth] = j;
	trie_prefix[depth+1] = 0;
	if (t->next[j]->count >= min_val) {
	    t->next[j]->sym = 255;
	    cp = assign_trie_words(trie_prefix, last_prefix, trie_counter,
				   cp, t, t->next[j], depth+1, min_val);
	} else {
	    t->sym = 0;
	}
    }

    if (t->count >= min_val) {
	trie_prefix[depth] = 0;
	if (!t->sym && *trie_counter < 254) {
	    t->sym = (*trie_counter)++;
	    if (*trie_prefix) {
		unsigned char *cp2 = trie_prefix;
		unsigned char *cpl = last_prefix;
		int prefix_len = 0;

		while (*cp2 && *cp2 == *cpl) {
		    cp2++;
		    cpl++;
		    prefix_len++;
		}

		*cp++ = prefix_len+1;
		do 
		    *cp++ = *cp2;
		while (*cp2++);

		strncpy((char*)last_prefix, (char*)trie_prefix, 256);
		//fprintf(stderr, ">>%3d %5d %s<<\n", t->sym, t->count, trie_prefix);
	    }
	}
    }

    return cp;
}

void free_trie(trie_t *t) {
    int j;
    for (j = 0; j < 128; j++) {
	if (t->next[j])
	    free_trie(t->next[j]);
    }
    free(t);
}

static int divi = 32;

static void *pool_calloc(pool_alloc_t *p, size_t size) {
    void *v = pool_alloc(p);
    if (v)
	memset(v, 0, size);
    return v;
}


unsigned char *compress_block(int level,
			      unsigned char *data,
			      size_t len,
			      size_t *comp_len) {
    HashTable *h = HashTableCreate(0, HASH_DYNAMIC_SIZE | HASH_NONVOLATILE_KEYS);
    int i;
    unsigned char *comp = NULL, *cp = NULL;
    int nlines = 0;
    pool_alloc_t *t_pool = pool_create(sizeof(trie_t));
    trie_t *t_head, *t;
    //char *last = NULL;
    unsigned char *dict;
    int dict_len;
    unsigned char *back1, *back1p;
    unsigned char *back2, *back2p;
    unsigned char *back3, *back3p;
    unsigned char *back4, *back4p;

    unsigned char trie_prefix[256] = {0};
    unsigned char last_prefix[256] = {0};
    int trie_counter = 128;

    if (!t_pool)
	goto err;

    t_head = (trie_t *)pool_calloc(t_pool, sizeof(*t_head));
    dict = malloc(len+10);

    // Build our trie, also counting input lines
    for (nlines = i = 0; i < len; i++, nlines++) {
	t = t_head;
	t->count++;
	while (i < len && data[i] > '\n') {
	    unsigned char c = data[i++];
	    if (c & 0x80)
		//fprintf(stderr, "8-bit ASCII is unsupported\n");
		goto err;
	    c &= 127;
	    if (!t->next[c]) 
		t->next[c] = (trie_t *)pool_calloc(t_pool, sizeof(*t));
	    t = t->next[c];
	    t->c = c;
	    t->count++;
	}
    }

    comp = malloc(4 + len*2 + nlines*3); // worst case?
    if (!comp) {
	HashTableDestroy(h, 0);
	return NULL;
    }

    // Version
    cp = dict;
    *cp++ = 0;

    // Store line count 
    *cp++ = (nlines>> 0) & 0xff;
    *cp++ = (nlines>> 8) & 0xff;
    *cp++ = (nlines>>16) & 0xff;
    *cp++ = (nlines>>24) & 0xff;

    *cp++ = (len>> 0) & 0xff;
    *cp++ = (len>> 8) & 0xff;
    *cp++ = (len>>16) & 0xff;
    *cp++ = (len>>24) & 0xff;
    //fprintf(stderr, "min_val = %d, nlines = %d\n", t_head->count/divi, nlines);

    // Assign and store tokens
    cp = assign_trie_words(trie_prefix, last_prefix, &trie_counter,
			   cp, NULL, t_head, 0, t_head->count/divi);
    *cp++ = 0; // end of dict.
    dict_len = cp - dict;

    //write(1, comp, cp-comp); cp = comp;

    back1p = back1 = malloc(nlines * sizeof(int));
    back2p = back2 = (nlines >= 0x100 ? malloc(nlines * sizeof(int)) : NULL);
    back3p = back3 = (nlines >= 0x10000 ? malloc(nlines * sizeof(int)) : NULL);
    back4p = back4 = (nlines >= 0x1000000 ? malloc(nlines * sizeof(int)) : NULL);

    // Tokenise and store each line
    cp = comp;
    for (nlines = i = 0; i < len; i++, nlines++) {
	HashData hd;
	HashItem *hi;
	int j = i;
	int token = 0;
	int new;

	// Check if dup line
	while (i < len && data[i] > '\n')
	    i++;
	hd.i = nlines;
	//printf(">>%.*s<<\n", i-j, &data[j]);
	hi = HashTableAdd(h, (char*)&data[j], i-j, hd, &new);
	if (!new) {
	    int p = nlines - hi->data.i;
	    //fprintf(stderr, "%d %lld %d %f\n", nlines, hi->data.i, p, p2/(double)nlines);
	    *back1p++ = (p >>  0) & 0xff;
	    if (back2p) *back2p++ = (p >>  8) & 0xff;
	    if (back3p) *back3p++ = (p >> 16) & 0xff;
	    if (back4p) *back4p++ = (p >> 24) & 0xff;

	    continue;
	} else {
	    *back1p++ = 0;
	    if (back2p) *back2p++ = 0;
	    if (back3p) *back3p++ = 0;
	    if (back4p) *back4p++ = 0;
	}

	t = t_head;

	// Tokenise line
	i = j;
	while (i < len && data[i] > '\n') {
	    unsigned char c = data[i++] & 127;
	    if (t->next[c]) {
		if (t->next[c]->sym == 0) {
		    if (token)
			//printf("<%d>", token);
			//putchar(token);
			*cp++ = token;
		        //*hack2p++ = token;
		    //if (last)
		    //	fprintf(stderr, "%c %.10s %.10s\n", c, &data[i-1], last+i-1-j);
		    *cp++ = c;
		    //if (last)
		    //	*cp++ = data[i-1] ^ last[i-1-j];
		    //else
		    //	*cp++ = c;

		    //putchar(c);
		    token = 0;
		} else {
		    token = t->next[c]->sym;
		}
		t = t->next[c];
	    } else {
		*cp++ = c;
		//putchar(c);
	    }
	}
	if (token)
	    *cp++ = token;
	*cp++ = '\n';
	//putchar('\n');

	//last = &data[j];
    }

    if (1) {
	unsigned char *out = malloc(4 + len*2 + nlines*3 + 6*7);
	unsigned char *x, *outp = out;
	unsigned int len, p;

	x = rans_compress(dict, dict_len, &len, 0);
	p = len; do { *outp++ = (p&0x7f) | (p>=128?128:0); p >>= 7; } while (p);
	memcpy(outp, x, len); outp += len;
	free(x);
	free(dict);
	//fprintf(stderr, "dict\t%d\t%d\n", dict_len, len);

	if (back1) {
	    x = rans_compress(back1, back1p - back1, &len, 0);
	    p = len; do { *outp++ = (p&0x7f) | (p>=128?128:0); p >>= 7; } while (p);
	    memcpy(outp, x, len); outp += len;
	    free(x);
	    free(back1);
	    //fprintf(stderr, "back1\t%d\t%d\n", back1p - back1, len);
	}

	if (back2) {
	    x = rans_compress(back2, back2p - back2, &len, 0);
	    p = len; do { *outp++ = (p&0x7f) | (p>=128?128:0); p >>= 7; } while (p);
	    memcpy(outp, x, len); outp += len;
	    free(x);
	    free(back2);
	    //fprintf(stderr, "back1\t%d\t%d\n", back2p - back2, len);
	}

	if (back3) {
	    x = rans_compress(back3, back3p - back3, &len, 0);
	    p = len; do { *outp++ = (p&0x7f) | (p>=128?128:0); p >>= 7; } while (p);
	    memcpy(outp, x, len); outp += len;
	    free(x);
	    free(back3);
	    //fprintf(stderr, "back1\t%d\t%d\n", back3p - back3, len);
	}

	if (back4) {
	    x = rans_compress(back4, back4p - back4, &len, 0);
	    p = len; do { *outp++ = (p&0x7f) | (p>=128?128:0); p >>= 7; } while (p);
	    memcpy(outp, x, len); outp += len;
	    free(x);
	    free(back4);	
	    //fprintf(stderr, "back1\t%d\t%d\n", back4p - back4, len);
	}

	if (0) {
	    x = rans_compress(comp, cp - comp, &len, 1);
	} else {
	    x = malloc(len = (cp-comp)*1.01 + 600);
	    BZ2_bzBuffToBuffCompress((char *)x, &len, (char *)comp, cp-comp, 9, 0, 30);
	}
	p = len; do { *outp++ = (p&0x7f) | (p>=128?128:0); p >>= 7; } while (p);
	memcpy(outp, x, len); outp += len;
	free(x);
	free(comp);
	//fprintf(stderr, "text\t%d\t%d\n", cp - comp, len);

	cp = out;
	*comp_len = outp - out;
    } else {
	unsigned char *out = malloc(4 + len*2 + nlines*3 + 6*7);
	unsigned char *outp = out;

	memcpy(outp, dict, dict_len); outp += dict_len;
	free(dict);

	if (back1) {
	    memcpy(outp, back1, back1p - back1); outp += back1p - back1;
	    free(back1);
	}
	if (back2) {
	    memcpy(outp, back2, back2p - back2); outp += back2p - back2;
	    free(back2);
	}
	if (back3) {
	    memcpy(outp, back3, back3p - back3); outp += back3p - back3;
	    free(back3);
	}
	if (back4) {
	    memcpy(outp, back4, back4p - back4); outp += back4p - back4;
	    free(back4);
	}

	memcpy(outp, comp, cp - comp); outp += cp - comp;
	free(comp);

	cp = out;
	*comp_len = outp - out;
    }

    //*comp_len = 0; cp = comp;
    //*comp_len = cp-comp; cp = comp;
    //cp = rans_compress(comp, cp-comp, (unsigned int *)comp_len, 1);

 err:
    if (h)
	HashTableDestroy(h, 0);
    if (t_pool)
	pool_destroy(t_pool);

//    if (comp)
//	free(comp);

    return cp;
}

unsigned char *uncompress_block(unsigned char *data,
				size_t len,
				size_t *uncomp_len) {
    int i, x, t;
    unsigned char *word[256] = {0};
    int word_len[256];
    int sym = 128;
    unsigned int nlines, ulen;
    unsigned char **lines;
    unsigned char *uncomp, *cp;
    string_alloc_t *str_pool = string_pool_create(8192);


    unsigned char *dict, *b1, *b2, *b3, *b4, *text;
    unsigned int dict_len, b1_len, b2_len, b3_len, b4_len, text_len;
    unsigned int *back;

    // Uncompress dict
    i = x = t = 0; do { t |= (data[i] & 0x7f)<<x; x+=7; } while (data[i++] & 128);
    dict_len = t;
    dict = rans_uncompress(&data[i], dict_len, &dict_len);
    data += i + t;

    if (dict[0] != 0) {
	fprintf(stderr, "name_dedup: unsupported version number\n");
	string_pool_destroy(str_pool);
	return NULL;
    }


    nlines = (dict[1]<<0) + (dict[2]<<8) + (dict[3]<<16) + (dict[4]<<24);
    ulen   = (dict[5]<<0) + (dict[6]<<8) + (dict[7]<<16) + (dict[8]<<24);
    dict += 9; dict_len -= 9;

    //fprintf(stderr, "nlines = %d, text = %d\n", nlines, ulen);
    lines = calloc(nlines, sizeof(*lines));

    // Decode dictonary
    for (i = 0; i < len; i++) {
	int j, prefix_len;
	if (dict[i] == '\0')
	    break;

	prefix_len = dict[i++]-1;

	j = i;
	while (i < len && dict[i] > '\n')
	    i++;

	word[sym] = (unsigned char *)string_alloc(str_pool, i-j+prefix_len);
	if (prefix_len)
	    memcpy(word[sym], word[sym-1], prefix_len);
	memcpy(word[sym]+prefix_len, &dict[j], i-j);
	word_len[sym] = i-j+prefix_len;
	//printf(stderr, "%3d '%.*s'\n", sym, word_len[sym], word[sym]);
	sym++;
    }
    free(dict-9);


    // Uncompress back refs
    i = x = t = 0; do { t |= (data[i] & 0x7f)<<x; x+=7; } while (data[i++] & 128);
    b1_len = t;
    b1 = rans_uncompress(&data[i], b1_len, &b1_len);
    data += i + t;

    if (nlines >= 0x100) {
	i = x = t = 0; do { t |= (data[i] & 0x7f)<<x; x+=7; } while (data[i++] & 128);
	b2_len = t;
	b2 = rans_uncompress(&data[i], b2_len, &b2_len);
	data += i + t;
    } else {
	b2 = 0;
    }

    if (nlines >= 0x10000) {
	i = x = t = 0; do { t |= (data[i] & 0x7f)<<x; x+=7; } while (data[i++] & 128);
	b3_len = t;
	b3 = rans_uncompress(&data[i], b3_len, &b3_len);
	data += i + t;
    } else {
	b3 = 0;
    }

    if (nlines >= 0x1000000) {
	i = x = t = 0; do { t |= (data[i] & 0x7f)<<x; x+=7; } while (data[i++] & 128);
	b4_len = t;
	b4 = rans_uncompress(&data[i], b4_len, &b4_len);
	data += i + t;
    } else {
	b4 = 0;
    }

    back = malloc(nlines * sizeof(int));
    for (i = 0; i < nlines; i++) {
	back[i] = b1[i] +
	    (b2 ? b2[i]<< 8 : 0) +
	    (b3 ? b3[i]<<16 : 0) +
	    (b4 ? b4[i]<<24 : 0);
    }
    if (b1) free(b1);
    if (b2) free(b2);
    if (b3) free(b3);
    if (b4) free(b4);

    // Uncompress text
    i = x = t = 0; do { t |= (data[i] & 0x7f)<<x; x+=7; } while (data[i++] & 128);
    text_len = t;
    text = malloc(text_len=ulen*1.01+600); // FIXME GUESSWORK
    BZ2_bzBuffToBuffDecompress((char *)text, &text_len, (char *)&data[i], text_len, 0, 0);
    //text = rans_uncompress(&data[i], text_len, &text_len);
    data += i;

    // Decode text
    int j;
    i = j = 0;
    cp = uncomp = (unsigned char *)malloc(ulen);
    //fprintf(stderr, "nlines=%d\n", nlines);
    while (j < nlines) {
	if (back[j]) {
	    unsigned char *cp2 = lines[j-back[j]];
	    //printf("BACK %.20s\n", cp2);
	    do {
		*cp++ = *cp2;
	    } while(*cp2++ > '\n');
	    j++;
	    continue;
	}

	lines[j] = cp;
	do {
	    unsigned char c = text[i];
	    if (c >= 128) {
		//printf(">>%.*s<<", word_len[c], word[c]);
		memcpy(cp, word[c], word_len[c]);
		cp += word_len[c];
	    } else {
		*cp++ = text[i];
	    }
	} while (i < text_len && text[i++] > '\n');

	j++;
    }

    free(lines);
    free(text);
    free(back);
    string_pool_destroy(str_pool);

    *uncomp_len = cp-uncomp;
    return uncomp;
}

static cram_compressor c = {
    'n', //FOUR_CC("\0N3n"),
    1<<DS_RN, // all data series
    1.0,
    name,
    compress_block,
    uncompress_block,
};

cram_compressor *cram_compressor_init(void) {
    return &c;
}

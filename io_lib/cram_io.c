/*
 * TODO: 
 * - Test I/O performance of reading block at a time vs reading
 *   container at a time and decoding the blocks in memory.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "io_lib/cram.h"
#include "io_lib/os.h"
#include "io_lib/deflate_interlaced.h" /* For block_create() */

//Whether CIGAR has just M or uses = and X to indicate match and mismatch
//#define USE_X

/* -----------------------------------------------------------------------------
 * Utility / debugging functions
 */

/* zlib compression code - from Gap5's tg_iface_g.c */
static int tg_zlevel = 6;

static char *zlib_mem_deflate(char *data, size_t size, size_t *cdata_size) {
    z_stream s;
    unsigned char *cdata = NULL; /* Compressed output */
    int cdata_alloc = 0;
    int cdata_pos = 0;
    int err;

    cdata = malloc(cdata_alloc = size*1.05+10);
    cdata_pos = 0;

    /* Initialise zlib stream */
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)data;
    s.avail_in = size;
    s.total_in = 0;
    s.next_out  = cdata;
    s.avail_out = cdata_alloc;
    s.total_out = 0;
    s.data_type = Z_BINARY;

    err = deflateInit2(&s, tg_zlevel == -1 ? Z_DEFAULT_COMPRESSION : tg_zlevel,
		       Z_DEFLATED, 15, 9, Z_DEFAULT_STRATEGY);
    if (err != Z_OK) {
	fprintf(stderr, "zlib deflateInit2 error: %s\n", s.msg);
	return NULL;
    }

    /* Encode to 'cdata' array */
    for (;s.avail_in;) {
	s.next_out = &cdata[cdata_pos];
	s.avail_out = cdata_alloc - cdata_pos;
	if (cdata_alloc - cdata_pos <= 0) {
	    fprintf(stderr, "Deflate produced larger output than expected. Abort\n"); 
	    return NULL;
	}
	err = deflate(&s, Z_NO_FLUSH);
	cdata_pos = cdata_alloc - s.avail_out;
	if (err != Z_OK) {
	    fprintf(stderr, "zlib deflate error: %s\n", s.msg);
	    break;
	}
    }
    if (deflate(&s, Z_FINISH) != Z_STREAM_END) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    *cdata_size = s.total_out;

    if (deflateEnd(&s) != Z_OK) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    return (char *)cdata;
}

static char *zlib_mem_inflate(char *cdata, size_t csize, size_t *size) {
    z_stream s;
    unsigned char *data = NULL; /* Uncompressed output */
    int data_alloc = 0;
    int err;

    /* Starting point at uncompressed size, 4x compressed */
    data = malloc(data_alloc = csize*4+10);

    /* Initialise zlib stream */
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)cdata;
    s.avail_in = csize;
    s.total_in = 0;
    s.next_out  = data;
    s.avail_out = data_alloc;
    s.total_out = 0;

    //err = inflateInit(&s);
    err = inflateInit2(&s, 15 + 32);
    if (err != Z_OK) {
	fprintf(stderr, "zlib inflateInit error: %s\n", s.msg);
	return NULL;
    }

    /* Decode to 'data' array */
    for (;s.avail_in;) {
	s.next_out = &data[s.total_out];
	err = inflate(&s, Z_NO_FLUSH);
	if (err == Z_STREAM_END)
	    break;

	if (err != Z_OK) {
	    fprintf(stderr, "zlib inflate error: %s\n", s.msg);
	    break;
	}

	/* More to come, so realloc */
	data = realloc(data, data_alloc += s.avail_in*4 + 10);
	s.avail_out += s.avail_in*4+10;
    }
    inflateEnd(&s);

    *size = s.total_out;
    return (char *)data;
}

char *cram_block_method2str(enum cram_block_method m) {
    switch(m) {
    case RAW:	return "RAW";
    case GZIP:	return "GZIP";
    case BZIP2:	return "BZIP2";
    }
    return "?";
}

char *cram_content_type2str(enum cram_content_type t) {
    switch (t) {
    case FILE_HEADER:         return "FILE_HEADER";
    case COMPRESSION_HEADER:  return "COMPRESSION_HEADER";
    case MAPPED_SLICE:        return "MAPPED_SLICE";
    case UNMAPPED_SLICE:      return "UNMAPPED_SLICE";
    case EXTERNAL:            return "EXTERNAL";
    case CORE:                return "CORE";
    }
    return "?";
}

/* -----------------------------------------------------------------------------
 * Low level I/O functions for basic encoding and decoding of bits and bytes
 */

/* Static lookup tables */
static int cram_init_done = 0;
static unsigned int bam_flag_swap[512];
static unsigned char L[256];

/*
 * Initialise lookup tables. Used within cram_decode_seq().
 */
static void cram_init(void) {
    int i;

    if (cram_init_done)
	return;

    memset(L, 4, 256);
    L['A'] = 0; L['a'] = 0;
    L['C'] = 1; L['c'] = 1;
    L['G'] = 2; L['g'] = 2;
    L['T'] = 3; L['t'] = 3;

    for (i = 0; i < 512; i++) {
	int f = 0;

	if (i & CRAM_FPAIRED)      f |= BAM_FPAIRED;
	if (i & CRAM_FPROPER_PAIR) f |= BAM_FPROPER_PAIR;
	if (i & CRAM_FUNMAP)       f |= BAM_FUNMAP;
	if (i & CRAM_FREVERSE)     f |= BAM_FREVERSE;
	if (i & CRAM_FREAD1)       f |= BAM_FREAD1;
	if (i & CRAM_FREAD2)       f |= BAM_FREAD2;
	if (i & CRAM_FSECONDARY)   f |= BAM_FSECONDARY;
	if (i & CRAM_FQCFAIL)      f |= BAM_FQCFAIL;
	if (i & CRAM_FDUP)         f |= BAM_FDUP;

	bam_flag_swap[i] = f;
    }

    cram_init_done = 1;
}

/*
 * Reads an integer in ITF-8 encoding from 'cp' and stores it in
 * *val.
 *
 * Returns the number of bytes read on success
 *        -1 on failure
 */
int itf8_decode(cram_fd *fd, int32_t *val_p) {
    static int nbytes[16] = {
	0,0,0,0, 0,0,0,0,                               // 0000xxxx - 0111xxxx
	1,1,1,1,                                        // 1000xxxx - 1011xxxx
	2,2,                                            // 1100xxxx - 1101xxxx
	3,                                              // 1110xxxx
	4,                                              // 1111xxxx
    };

    static int nbits[16] = {
	0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, // 0000xxxx - 0111xxxx
	0x3f, 0x3f, 0x3f, 0x3f,                         // 1000xxxx - 1011xxxx
	0x1f, 0x1f,                                     // 1100xxxx - 1101xxxx
	0x0f,                                           // 1110xxxx
	0x0f,                                           // 1111xxxx
    };

    int sz = 1;
    int32_t val = getc(fd->fp);
    if (val == -1)
	return -1;

    int i = nbytes[val>>4];
    val &= nbits[val>>4];

    while (i) {
	val = (val<<8) | (unsigned char)getc(fd->fp);
	i--;
	sz++;
    }
    *val_p = val;
    
    return sz;
}

/*
 * As above, but decoding from memory
 */
int itf8_get(char *cp, int32_t *val_p) {
    static int nbytes[16] = {
	0,0,0,0, 0,0,0,0,                               // 0000xxxx - 0111xxxx
	1,1,1,1,                                        // 1000xxxx - 1011xxxx
	2,2,                                            // 1100xxxx - 1101xxxx
	3,                                              // 1110xxxx
	4,                                              // 1111xxxx
    };

    static int nbits[16] = {
	0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, // 0000xxxx - 0111xxxx
	0x3f, 0x3f, 0x3f, 0x3f,                         // 1000xxxx - 1011xxxx
	0x1f, 0x1f,                                     // 1100xxxx - 1101xxxx
	0x0f,                                           // 1110xxxx
	0x0f,                                           // 1111xxxx
    };

    int sz = 1;
    int32_t val = (unsigned char)*cp++;

    int i = nbytes[val>>4];
    val &= nbits[val>>4];

    while (i) {
	val = (val<<8) | (unsigned char)*cp++;
	i--;
	sz++;
    }
    *val_p = val;
    
    return sz;
}

/*
 * Stores a value to memory in ITF-8 format.
 *
 * Returns the number of bytes required to store the number.
 * This is a maximum of 5 bytes.
 */
int itf8_put(char *cp, int32_t val) {
    if        (!(val & ~0x00000007f)) { // 1 byte
	*cp = val;
	return 1;
    } else if (!(val & ~0x00003fff)) { // 2 byte
	*cp++ = (val >> 8 ) | 0x80;
	*cp   = val & 0xff;
	return 2;
    } else if (!(val & ~0x01ffffff)) { // 3 byte
	*cp++ = (val >> 16) | 0xc0;
	*cp++ = (val >> 8 ) & 0xff;
	*cp   = val & 0xff;
	return 3;
    } else if (!(val & ~0x0fffffff)) { // 4 byte
	*cp++ = (val >> 24) | 0xe0;
	*cp++ = (val >> 16) & 0xff;
	*cp++ = (val >> 8 ) & 0xff;
	*cp   = val & 0xff;
	return 4;
    } else {                           // 5 byte
	*cp++ = 0xf0;
	*cp++ = (val >> 24) & 0xff;
	*cp++ = (val >> 16) & 0xff;
	*cp++ = (val >> 8 ) & 0xff;
	*cp++ = (val >> 4 ) & 0xff;
	*cp   = val & 0xff;
	return 5;
    }
}


/* Get a single bit, MSB first */
signed int get_bit_MSB(block_t *block) {
    unsigned int val;

    if (block->byte > block->alloc)
	return -1;

    val = block->data[block->byte] >> block->bit;
    if (--block->bit == -1) {
	block->bit = 7;
	block->byte++;
	//printf("(%02X)", block->data[block->byte]);
    }

    //printf("-B%d-", val&1);

    return val & 1;
}

signed int get_bits_MSB(block_t *block, int nbits) {
    unsigned int val = 0;
    int i;

    /* Inefficient implementation! */
    //printf("{");
    for (i = 0; i < nbits; i++)
	//val = (val << 1) | get_bit_MSB(block);
	GET_BIT_MSB(block, val);

    //printf("=0x%x}", val);

    return val;
}

/*
 * Count number of successive 0 and 1 bits
 */
int get_one_bits_MSB(block_t *block) {
    int n = 0, b;
    do {
	b = block->data[block->byte] >> block->bit;
	if (--block->bit == -1) {
	    block->bit = 7;
	    block->byte++;
	}
	n++;
    } while (b&1);

    return n-1;
}

int get_zero_bits_MSB(block_t *block) {
    int n = 0, b;
    do {
	b = block->data[block->byte] >> block->bit;
	if (--block->bit == -1) {
	    block->bit = 7;
	    block->byte++;
	}
	n++;
    } while (!(b&1));

    return n-1;
}


/* -----------------------------------------------------------------------------
 * Mid level I/O functions for manipulating CRAM file structures:
 * Headers, containers, blocks, etc
 */

/*
 * Reads a CRAM file definition structure.
 * Returns file_def ptr on success
 *         NULL on failure
 */
cram_file_def *cram_read_file_def(cram_fd *fd) {
    cram_file_def *def = malloc(sizeof(*def));
    if (!def)
	return NULL;

    if (26 != fread(&def->magic[0], 1, 26, fd->fp)) {
	free(def);
	return NULL;
    }

    if (memcmp(def->magic, "CRAM", 4) != 0) {
	fprintf(stderr, "CRAM magic number failed\n");
	free(def);
	return NULL;
    }

    if (def->major_version != 1 || def->minor_version != 0) {
	fprintf(stderr, "CRAM version number mismatch\n"
		"Expected 1.0, got %d.%d\n",
		def->major_version, def->minor_version);
	free(def);
	return NULL;
    }

    return def;
}

int cram_write_file_def(cram_fd *fd, cram_file_def *def) {
    return (fwrite(&def->magic[0], 1, 26, fd->fp) == 26) ? 0 : -1;
}

void cram_free_file_def(cram_file_def *def) {
    if (def) free(def);
}


/*
 * Reads the SAM header from the first CRAM data block.
 * Also performs minimal parsing to extract read-group
 * and sample information.

 * Returns SAM hdr ptr on success
 *         NULL on failure
 */
cram_SAM_hdr *cram_read_SAM_hdr(cram_fd *fd) {
    cram_SAM_hdr *hdr = calloc(1, sizeof(*hdr));
    if (!hdr)
	return NULL;

    /* Length */
    if (1 != fread(&hdr->length, 4, 1, fd->fp)) {
	free(hdr);
	return NULL;
    }
    hdr->length = le_int4(hdr->length);

    /* Alloc and read */
    if (NULL == (hdr->header = malloc(hdr->length+1))) {
	free(hdr);
	return NULL;
    }

    if (hdr->length != fread(hdr->header, 1, hdr->length, fd->fp)) {
	free(hdr);
	return NULL;
    }

    /* Parse */
    // TODO: Consider reusing bam_parse_header(bam_file_t *b); */
    
    return hdr;
}

/*
 * Creates a CRAM header from a SAM header in string format.
 *
 * FIXME: consider either rejecting this completely and using
 * "char *SAM_hdr" throughout, or instead finishing this off by copying
 * the bam_parse_header() code into here.
 *
 * FIXME 2: check consistency of header. Needs SQ:MD5, HD:SO as POS,
 * RG lines, etc.
 *
 * Returns cram_SAM_hdr* on success
 *         NULL on failure
 */
cram_SAM_hdr *cram_create_SAM_hdr(char *str, size_t len) {
    bam_file_t b;

    cram_SAM_hdr *hdr = calloc(1, sizeof(*hdr));
    if (!hdr)
	return NULL;

    if (NULL == (hdr->header = malloc(len))) {
	free(hdr);
	return NULL;
    }

    memcpy(hdr->header, str, len);
    hdr->length = len;
    hdr->ref_hash = NULL;
    hdr->rg_hash = NULL;

    /* Reuse bam parsing. FIXME: this is ugly and needs merging */
    memset(&b, 0, sizeof(b));
    b.header = str;
    b.header_len = len;
    bam_parse_header(&b);

    hdr->ref_hash = b.ref_hash;
    hdr->rg_hash  = b.rg_hash;
    if (b.rg_id)  free(b.rg_id);
    if (b.rg_len) free(b.rg_len);

    return hdr;
}


/*
 * Writes a CRAM SAM header.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_SAM_hdr(cram_fd *fd, cram_SAM_hdr *hdr) {
    int32_t le_len = le_int4(hdr->length);

    /* Length */
    if (1 != fwrite(&le_len, 4, 1, fd->fp))
	return -1;

    /* Text data */
    if (hdr->length != fwrite(hdr->header, 1, hdr->length, fd->fp))
	return -1;

    return 0;
}


void cram_free_SAM_hdr(cram_SAM_hdr *hdr) {
    if (!hdr)
	return;

    if (hdr->header)
	free(hdr->header);

    if (hdr->ref_hash)
	HashTableDestroy(hdr->ref_hash, 0);

    if (hdr->rg_hash)
	HashTableDestroy(hdr->rg_hash, 1);

    free(hdr);
}


/*
 * Reads a container header.
 * Returns cram_container on success
 *         NULL on failure or no container left (fd->err == 0).
 */
cram_container *cram_read_container(cram_fd *fd) {
    cram_container c2, *c;
    int i, s;
    size_t rd = 0;
    cram_block *b;
    
    fd->err = 0;

    if ((s = itf8_decode(fd, &c2.length))        == -1) return NULL; else rd+=s;
    if ((s = itf8_decode(fd, &c2.ref_seq_id))    == -1) return NULL; else rd+=s;
    if ((s = itf8_decode(fd, &c2.ref_seq_start)) == -1) return NULL; else rd+=s;
    if ((s = itf8_decode(fd, &c2.ref_seq_span))  == -1) return NULL; else rd+=s;
    if ((s = itf8_decode(fd, &c2.num_records))   == -1) return NULL; else rd+=s;
    if ((s = itf8_decode(fd, &c2.num_blocks))    == -1) return NULL; else rd+=s;
    if ((s = itf8_decode(fd, &c2.num_landmarks)) == -1) return NULL; else rd+=s;

    if (!(c = calloc(1, sizeof(*c))))
	return NULL;

    *c = c2;

    if (!(c->landmark = malloc(c->num_landmarks * sizeof(int32_t)))) {
	fd->err = errno;
	cram_free_container(c);
	return NULL;
    }  
    for (i = 0; i < c->num_landmarks; i++) {
	if ((s = itf8_decode(fd, &c->landmark[i])) == -1) {
	    cram_free_container(c);
	    return NULL;
	} else {
	    rd += s;
	}
    }
    c->offset = rd;

    if (!(b = cram_read_block(fd))) {
	cram_free_container(c);
	return NULL;
    }
    assert(b->content_type == COMPRESSION_HEADER);

    /* Keep uncompressed block too as hdr points into it */
    c->comp_hdr = cram_decode_compression_header(b);
    c->comp_hdr_block = b;

    if (!c->comp_hdr) {
	fprintf(stderr, "Unable to decode compression header\n");
	cram_free_container(c);
	return NULL;
    }

    c->slices = NULL;

    return c;
}

/*
 * Creates a new container, specifying the maximum number of slices
 * and records permitted.
 *
 * Returns cram_container ptr on success
 *         NULL on failure
 */
cram_container *cram_new_container(int nrec, int nslice) {
    cram_container *c = calloc(1, sizeof(*c));
    if (!c)
	return NULL;

    c->curr_ref = -2;

    c->max_rec = nrec;
    c->curr_rec = 0;

    c->max_slice = nslice;
    c->curr_slice = 0;

    c->slices = (cram_slice **)calloc(nslice, sizeof(cram_slice *));
    c->slice = NULL;

    c->curr_ctr_rec = 0;

    return c;
}

void cram_free_container(cram_container *c) {
    int i;

    if (!c)
	return;

    if (c->landmark)
	free(c->landmark);

    if (c->comp_hdr)
	cram_free_compression_header(c->comp_hdr);

    if (c->comp_hdr_block)
	cram_free_block(c->comp_hdr_block);

    if (c->slices) {
	for (i = 0; i < c->max_slice; i++)
	    if (c->slices[i])
		cram_free_slice(c->slices[i]);
	free(c->slices);
    }

    free(c);
}

/*
 * Reads a block from a cram file.
 * Returns cram_block pointer on success.
 *         NULL on failure
 */
cram_block *cram_read_block(cram_fd *fd) {
    cram_block *b = malloc(sizeof(*b));
    if (!b)
	return NULL;

    if (-1 == (b->method       = getc(fd->fp))) { free(b); return NULL; }
    if (-1 == (b->content_type = getc(fd->fp))) { free(b); return NULL; }
    if (-1 == itf8_decode(fd, &b->content_id))  { free(b); return NULL; }
    if (-1 == itf8_decode(fd, &b->comp_size))   { free(b); return NULL; }
    if (-1 == itf8_decode(fd, &b->uncomp_size)) { free(b); return NULL; }

    if (!(b->data = malloc(b->comp_size)))      { free(b); return NULL; }
    if (b->comp_size != fread(b->data, 1, b->comp_size, fd->fp)) {
	free(b->data);
	free(b);
	return NULL;
    }

    b->orig_method = b->method;
    b->idx = 0;

    return b;
}

void cram_free_block(cram_block *b) {
    if (!b)
	return;
    if (b->data)
	free(b->data);
    free(b);
}

/*
 * Returns the next 'size' bytes from a block, or NULL if insufficient
 * data left.This is just a pointer into the block data and not an
 * allocated object, so do not free the result.
 */
char *cram_extract_block(cram_block *b, int size) {
    char *cp = b->data + b->idx;
    b->idx += size;
    if (b->idx > b->uncomp_size)
	return NULL;

    return cp;
}

void cram_uncompress_block(cram_block *b) {
    char *uncomp;
    size_t uncomp_size = 0;

    switch (b->method) {
    case RAW:
	return;

    case GZIP:
	uncomp = zlib_mem_inflate(b->data, b->comp_size, &uncomp_size);
	assert((int)uncomp_size == b->uncomp_size);
	free(b->data);
	b->data = uncomp;
	b->method = RAW;
	break;

    case BZIP2:
	fprintf(stderr, "Bzip2 compression not yet implemented\n");
	abort();
	break;
    }
}

/*
 * Decodes a CRAM block compression header.
 * Returns header ptr on success
 *         NULL on failure
 */
cram_block_compression_hdr *cram_decode_compression_header(cram_block *b) {
    char *cp, *cp_copy;
    cram_block_compression_hdr *hdr = calloc(1, sizeof(*hdr));
    int i;
    int32_t map_size, map_count;

    if (!hdr)
	return NULL;

    if (b->method != RAW)
	cram_uncompress_block(b);

    cp = b->data;
    cp += itf8_get(cp, &hdr->ref_seq_id);
    cp += itf8_get(cp, &hdr->ref_seq_start);
    cp += itf8_get(cp, &hdr->ref_seq_span);
    cp += itf8_get(cp, &hdr->num_records);
    cp += itf8_get(cp, &hdr->num_landmarks);
    if (!(hdr->landmark = malloc(hdr->num_landmarks * sizeof(int32_t)))) {
	free(hdr);
	return NULL;
    }
    for (i = 0; i < hdr->num_landmarks; i++) {
	cp += itf8_get(cp, &hdr->landmark[i]);
    }

    hdr->preservation_map = HashTableCreate(4, HASH_NONVOLATILE_KEYS);
    hdr->rec_encoding_map = HashTableCreate(16, HASH_DYNAMIC_SIZE |
					    HASH_NONVOLATILE_KEYS);
    hdr->tag_encoding_map = HashTableCreate(4, HASH_DYNAMIC_SIZE |
					    HASH_NONVOLATILE_KEYS);

    if (!hdr->preservation_map ||
	!hdr->rec_encoding_map ||
	!hdr->tag_encoding_map) {
	cram_free_compression_header(hdr);
	return NULL;
    }

    /* Initialise defaults for preservation map */
    hdr->mapped_qs_included = 0;
    hdr->unmapped_qs_included = 0;
    hdr->unmapped_placed = 0;
    hdr->qs_included = 0;
    hdr->read_names_included = 0;
    memcpy(hdr->substitution_matrix, "CGTNAGTNACTNACGNACGT", 20);

    /* Preservation map */
    cp += itf8_get(cp, &map_size); cp_copy = cp;
    cp += itf8_get(cp, &map_count);
    for (i = 0; i < map_count; i++) {
	HashData hd;

	cp += 2;
	switch(CRAM_KEY(cp[-2],cp[-1])) {
	case CRAM_KEY('M','I'):
	    hd.i = *cp++;
	    HashTableAdd(hdr->preservation_map, "MI", 2, hd, NULL);
	    hdr->mapped_qs_included = hd.i;
	    break;

	case CRAM_KEY('U','I'):
	    hd.i = *cp++;
	    HashTableAdd(hdr->preservation_map, "UI", 2, hd, NULL);
	    hdr->unmapped_qs_included = hd.i;
	    break;

	case CRAM_KEY('P','I'):
	    hd.i = *cp++;
	    HashTableAdd(hdr->preservation_map, "PI", 2, hd, NULL);
	    hdr->unmapped_placed = hd.i;
	    break;

	case CRAM_KEY('R','N'):
	    hd.i = *cp++;
	    HashTableAdd(hdr->preservation_map, "RN", 2, hd, NULL);
	    hdr->read_names_included = hd.i;
	    break;

	case CRAM_KEY('S','M'):
	    hdr->substitution_matrix[0][0] = "CGTN"[(cp[0]>>6)&3];
	    hdr->substitution_matrix[0][1] = "CGTN"[(cp[0]>>4)&3];
	    hdr->substitution_matrix[0][2] = "CGTN"[(cp[0]>>2)&3];
	    hdr->substitution_matrix[0][3] = "CGTN"[(cp[0]>>0)&3];

	    hdr->substitution_matrix[1][0] = "AGTN"[(cp[1]>>6)&3];
	    hdr->substitution_matrix[1][1] = "AGTN"[(cp[1]>>4)&3];
	    hdr->substitution_matrix[1][2] = "AGTN"[(cp[1]>>2)&3];
	    hdr->substitution_matrix[1][3] = "AGTN"[(cp[1]>>0)&3];

	    hdr->substitution_matrix[2][0] = "ACTN"[(cp[2]>>6)&3];
	    hdr->substitution_matrix[2][1] = "ACTN"[(cp[2]>>4)&3];
	    hdr->substitution_matrix[2][2] = "ACTN"[(cp[2]>>2)&3];
	    hdr->substitution_matrix[2][3] = "ACTN"[(cp[2]>>0)&3];

	    hdr->substitution_matrix[3][0] = "ACGN"[(cp[3]>>6)&3];
	    hdr->substitution_matrix[3][1] = "ACGN"[(cp[3]>>4)&3];
	    hdr->substitution_matrix[3][2] = "ACGN"[(cp[3]>>2)&3];
	    hdr->substitution_matrix[3][3] = "ACGN"[(cp[3]>>0)&3];

	    hdr->substitution_matrix[4][0] = "ACGT"[(cp[4]>>6)&3];
	    hdr->substitution_matrix[4][1] = "ACGT"[(cp[4]>>4)&3];
	    hdr->substitution_matrix[4][2] = "ACGT"[(cp[4]>>2)&3];
	    hdr->substitution_matrix[4][3] = "ACGT"[(cp[4]>>0)&3];
	    hd.p = cp;
	    cp += 5;
	    HashTableAdd(hdr->preservation_map, "SM", 2, hd, NULL);
	    break;
	}
    }
    assert(cp - cp_copy == map_size);

    /* Record encoding map */
    cp += itf8_get(cp, &map_size); cp_copy = cp;
    cp += itf8_get(cp, &map_count);
    for (i = 0; i < map_count; i++) {
	HashData hd;
	char *key = cp;
	int32_t encoding;
	int32_t size;
	cram_map *m = malloc(sizeof(*m)); // FIXME: use pooled_alloc

	//hd.p = cp;
	cp += 2;
	cp += itf8_get(cp, &encoding);
	cp += itf8_get(cp, &size);
	//hd.i = ((uint64_t)encoding << 32) | size;

	m->encoding = encoding;
	m->size     = size;
	m->offset   = cp - b->data;
	hd.p = m;

	//printf("%s codes for %.2s\n", cram_encoding2str(encoding), key);

	if (key[0] == 'B' && key[1] == 'F')
	    hdr->BF_codec = cram_decoder_init(encoding, cp, size, E_INT);
	else if (key[0] == 'C' && key[1] == 'F')
	    hdr->CF_codec = cram_decoder_init(encoding, cp, size, E_BYTE);
	else if (key[0] == 'R' && key[1] == 'L')
	    hdr->RL_codec = cram_decoder_init(encoding, cp, size, E_INT);
	else if (key[0] == 'A' && key[1] == 'P')
	    hdr->AP_codec = cram_decoder_init(encoding, cp, size, E_LONG);
	else if (key[0] == 'R' && key[1] == 'G')
	    hdr->RG_codec = cram_decoder_init(encoding, cp, size, E_INT);
	else if (key[0] == 'M' && key[1] == 'F')
	    hdr->MF_codec = cram_decoder_init(encoding, cp, size, E_BYTE);
	else if (key[0] == 'N' && key[1] == 'S')
	    hdr->NS_codec = cram_decoder_init(encoding, cp, size, E_INT);
	else if (key[0] == 'N' && key[1] == 'P')
	    hdr->NP_codec = cram_decoder_init(encoding, cp, size, E_INT);
	else if (key[0] == 'T' && key[1] == 'S')
	    hdr->TS_codec = cram_decoder_init(encoding, cp, size, E_INT);
	else if (key[0] == 'N' && key[1] == 'F')
	    hdr->NF_codec = cram_decoder_init(encoding, cp, size, E_INT);
	else if (key[0] == 'T' && key[1] == 'C')
	    hdr->TC_codec = cram_decoder_init(encoding, cp, size, E_BYTE);
	else if (key[0] == 'T' && key[1] == 'N')
	    hdr->TN_codec = cram_decoder_init(encoding, cp, size, E_BYTE_ARRAY);
	else if (key[0] == 'F' && key[1] == 'N')
	    hdr->FN_codec = cram_decoder_init(encoding, cp, size, E_INT);
	else if (key[0] == 'F' && key[1] == 'C')
	    hdr->FC_codec = cram_decoder_init(encoding, cp, size, E_BYTE);
	else if (key[0] == 'F' && key[1] == 'P')
	    hdr->FP_codec = cram_decoder_init(encoding, cp, size, E_INT);
	else if (key[0] == 'B' && key[1] == 'S')
	    hdr->BS_codec = cram_decoder_init(encoding, cp, size, E_BYTE);
	else if (key[0] == 'I' && key[1] == 'N')
	    hdr->IN_codec = cram_decoder_init(encoding, cp, size, E_BYTE_ARRAY);
	else if (key[0] == 'D' && key[1] == 'L')
	    hdr->DL_codec = cram_decoder_init(encoding, cp, size, E_INT);
	else if (key[0] == 'B' && key[1] == 'A')
	    hdr->BA_codec = cram_decoder_init(encoding, cp, size, E_BYTE);
	else if (key[0] == 'M' && key[1] == 'Q')
	    hdr->MQ_codec = cram_decoder_init(encoding, cp, size, E_INT);
	else if (key[0] == 'R' && key[1] == 'N')
	    hdr->RN_codec = cram_decoder_init(encoding, cp, size, E_BYTE_ARRAY);
	else if (key[0] == 'Q' && key[1] == 'S') {
	    hdr->QS_codec = cram_decoder_init(encoding, cp, size, E_BYTE);
	    hdr->Qs_codec = cram_decoder_init(encoding, cp, size, E_BYTE_ARRAY);
	} else
	    0;//fprintf(stderr, "Unrecognised key: %.2s\n", key);

	cp += size;

	HashTableAdd(hdr->rec_encoding_map, key, 2, hd, NULL);
    }
    assert(cp - cp_copy == map_size);


    /* Tag encoding map */
    cp += itf8_get(cp, &map_size); cp_copy = cp;
    cp += itf8_get(cp, &map_count);
    for (i = 0; i < map_count; i++) {
	HashData hd;
	int32_t encoding;
	int32_t size;
	cram_map *m = malloc(sizeof(*m)); // FIXME: use pooled_alloc
	char *key = cp+1;

	cp += 4; // Strictly ITF8, but this suffices
	cp += itf8_get(cp, &encoding);
	cp += itf8_get(cp, &size);

	m->encoding = encoding;
	m->size     = size;
	m->offset   = cp - b->data;
	m->codec = cram_decoder_init(encoding, cp, size, E_BYTE_ARRAY);
	hd.p = m;
	
	cp += size;

	HashTableAdd(hdr->tag_encoding_map, key, 3, hd, NULL);
    }
    assert(cp - cp_copy == map_size);


    return hdr;
}


void cram_free_compression_header(cram_block_compression_hdr *hdr) {
    if (hdr->landmark)
	free(hdr->landmark);

    if (hdr->preservation_map)
	HashTableDestroy(hdr->preservation_map, 0);

    if (hdr->rec_encoding_map)
	HashTableDestroy(hdr->rec_encoding_map, 1);

    if (hdr->tag_encoding_map) {
	HashItem *hi;
	HashIter *iter = HashTableIterCreate();

	while ((hi = HashTableIterNext(hdr->tag_encoding_map, iter))) {
	    cram_map *m = hi->data.p;
	    m->codec->free(m->codec);
	}

	HashTableDestroy(hdr->tag_encoding_map, 1);

	HashTableIterDestroy(iter);
    }

    if (hdr->BF_codec) hdr->BF_codec->free(hdr->BF_codec);
    if (hdr->CF_codec) hdr->CF_codec->free(hdr->CF_codec);
    if (hdr->RL_codec) hdr->RL_codec->free(hdr->RL_codec);
    if (hdr->AP_codec) hdr->AP_codec->free(hdr->AP_codec);
    if (hdr->RG_codec) hdr->RG_codec->free(hdr->RG_codec);
    if (hdr->MF_codec) hdr->MF_codec->free(hdr->MF_codec);
    if (hdr->NS_codec) hdr->NS_codec->free(hdr->NS_codec);
    if (hdr->NP_codec) hdr->NP_codec->free(hdr->NP_codec);
    if (hdr->TS_codec) hdr->TS_codec->free(hdr->TS_codec);
    if (hdr->NF_codec) hdr->NF_codec->free(hdr->NF_codec);
    if (hdr->TC_codec) hdr->TC_codec->free(hdr->TC_codec);
    if (hdr->TN_codec) hdr->TN_codec->free(hdr->TN_codec);
    if (hdr->FN_codec) hdr->FN_codec->free(hdr->FN_codec);
    if (hdr->FC_codec) hdr->FC_codec->free(hdr->FC_codec);
    if (hdr->FP_codec) hdr->FP_codec->free(hdr->FP_codec);
    if (hdr->BS_codec) hdr->BS_codec->free(hdr->BS_codec);
    if (hdr->IN_codec) hdr->IN_codec->free(hdr->IN_codec);
    if (hdr->DL_codec) hdr->DL_codec->free(hdr->DL_codec);
    if (hdr->BA_codec) hdr->BA_codec->free(hdr->BA_codec);
    if (hdr->MQ_codec) hdr->MQ_codec->free(hdr->MQ_codec);
    if (hdr->RN_codec) hdr->RN_codec->free(hdr->RN_codec);
    if (hdr->QS_codec) hdr->QS_codec->free(hdr->QS_codec);
    if (hdr->Qs_codec) hdr->Qs_codec->free(hdr->Qs_codec);

    free(hdr);
}


/*
 * Decodes a CRAM (un)mapped slice header block.
 * Returns slice header ptr on success
 *         NULL on failure
 */
cram_block_slice_hdr *cram_decode_slice_header(cram_block *b) {
    cram_block_slice_hdr *hdr;
    char *cp = b->data;
    int i;

    if (b->content_type != MAPPED_SLICE &&
	b->content_type != UNMAPPED_SLICE)
	return NULL;

    if (!(hdr  = calloc(1, sizeof(*hdr))))
	return NULL;

    hdr->content_type = b->content_type;

    if (b->content_type == MAPPED_SLICE) {
	cp += itf8_get(cp, &hdr->ref_seq_id);
	cp += itf8_get(cp, &hdr->ref_seq_start);
	cp += itf8_get(cp, &hdr->ref_seq_span);
    }
    cp += itf8_get(cp, &hdr->num_records);
    cp += itf8_get(cp, &hdr->num_blocks);

    hdr->block_content_ids = malloc(hdr->num_blocks * sizeof(int32_t));
    for (i = 0; i < hdr->num_blocks; i++) {
	cp += itf8_get(cp, &hdr->block_content_ids[i]);
    }

    if (b->content_type == MAPPED_SLICE) {
	cp += itf8_get(cp, &hdr->ref_base_id);
    }
    return hdr;
}

void cram_free_slice_header(cram_block_slice_hdr *hdr) {
    if (!hdr)
	return;

    if (hdr->block_content_ids)
	free(hdr->block_content_ids);

    free(hdr);

    return;
}



/*
 * Loads an entire slice.
 * FIXME: In 1.0 the native unit of slices within CRAM is broken
 * as slices contain references to objects in other slices.
 * To work around this while keeping the slice oriented outer loop
 * we read all slices and stitch them together into a fake large
 * slice instead.
 *
 * Returns cram_slice ptr on success
 *         NULL on failure
 */
cram_slice *cram_read_slice(cram_fd *fd) {
    cram_block *b = cram_read_block(fd);
    cram_slice *s = calloc(1, sizeof(*s));
    int i, n;

    if (!b || !s) {
	if (b) cram_free_block(b);
	if (s) free(s);
	return NULL;
    }

    switch (b->content_type) {
    case MAPPED_SLICE:
    case UNMAPPED_SLICE:
	s->hdr = cram_decode_slice_header(b);
	break;

    default:
	fprintf(stderr, "Unexpected block of type %s\n",
		cram_content_type2str(b->content_type));
	cram_free_block(b);
	cram_free_slice(s);
	return NULL;
    }

    s->block = calloc(n = s->hdr->num_blocks, sizeof(*s->block));
    if (!s->block) {
	cram_free_block(b);
	cram_free_slice(s);
	return NULL;
    }

    for (i = 0; i < n; i++) {
	s->block[i] = cram_read_block(fd);
    }

    cram_free_block(b);

    /* Initialise encoding/decoding tables */
    s->cigar = NULL;
    s->cigar_alloc = 0;
    s->ncigar = 0;

    s->name_ds = dstring_create(NULL);
    s->seqs_ds = dstring_create(NULL);
    s->qual_ds = dstring_create(NULL);
    s->aux_ds  = dstring_create(NULL);

    s->crecs = NULL;

    s->last_apos = s->hdr->ref_seq_start;
    
    s->id = fd->slice_num++;

    return s;
}

/*
 * Creates a new empty slice in memory, for subsequent writing to
 * disk.
 *
 * Returns cram_slice ptr on success
 *         NULL on failure
 */
cram_slice *cram_new_slice(enum cram_content_type type, int nrecs) {
    cram_slice *s = calloc(1, sizeof(*s));
    if (!s)
	return NULL;

    s->hdr = (cram_block_slice_hdr *)calloc(1, sizeof(*s->hdr));
    s->hdr->content_type = type;

    s->block = NULL;
    s->last_apos = 0;
    s->id = 0;
    s->crecs = malloc(nrecs * sizeof(cram_record));
    s->cigar = NULL;
    s->cigar_alloc = 0;
    s->ncigar = 0;

    s->name_ds = dstring_create(NULL);
    s->seqs_ds = dstring_create(NULL);
    s->qual_ds = dstring_create(NULL);
    s->aux_ds  = dstring_create(NULL);

    return s;
}

void cram_free_slice(cram_slice *s) {
    if (!s)
	return;

    if (s->block) {
	int i;

	if (s->hdr) {
	    for (i = 0; i < s->hdr->num_blocks; i++) {
		cram_free_block(s->block[i]);
	    }
	}
	free(s->block);
    }

    if (s->hdr)
	cram_free_slice_header(s->hdr);

    if (s->name_ds)
	dstring_destroy(s->name_ds);

    if (s->seqs_ds)
	dstring_destroy(s->seqs_ds);

    if (s->qual_ds)
	dstring_destroy(s->qual_ds);

    if (s->aux_ds)
	dstring_destroy(s->aux_ds);

    if (s->cigar)
	free(s->cigar);

    if (s->crecs)
	free(s->crecs);

    free(s);
}

/*
 * Loads a reference.
 * FIXME: use the sam_comp.cpp get_ref_base() equivalent to allow
 * random access within an indexed reference instead?
 *
 * Returns a ref_seq structure on success
 *         NULL on failure
 */
refs *load_reference(char *fn) {
    struct stat sb;
    char *ref;
    FILE *fp;
    off_t i, j = 0;
    char *name = NULL, *seq = NULL;
    HashData hd;
    HashTable *h;

    refs *r = malloc(sizeof(*r));
    if (!r)
	return NULL;

    h = r->h = HashTableCreate(16, HASH_DYNAMIC_SIZE | HASH_NONVOLATILE_KEYS);

    /* Load reference */
    if (stat(fn, &sb) != 0) {
	perror(fn);
	return NULL;
    }

    if (!(fp = fopen(fn, "r"))) {
	perror(fn);
	return NULL;
    }

    if (!(r->file = ref = malloc(sb.st_size)))
	return NULL;

    if (sb.st_size != fread(ref, 1, sb.st_size, fp)) {
	fprintf(stderr, "Failed to load reference file\n");
	free(ref);
	return NULL;
    }
    fclose(fp);

    /* Parse it */
    for (i = 0; i < sb.st_size; i++) {
	if (ref[i] == '>') {
	    if (name) {
		seq[j] = 0;

		hd.p = seq;
		HashTableAdd(h, name, strlen(name), hd, NULL);
	    }

	    name = &ref[i+1];
	    while (ref[i] != '\n')
		i++;
	    ref[i] = 0;
	    seq = &ref[i+1];
	    j = 0;
	} else {
	    if (ref[i] != '\n')
		seq[j++] = ref[i];
	}
    }
    if (name) {
	seq[j] = 0;
	hd.p = seq;
	HashTableAdd(h, name, strlen(name), hd, NULL);
    }

    r->ref_id = NULL;

    return r;
}

void free_refs(refs *r) {
    if (r->ref_id)
	free(r->ref_id);
    if (r->h)
	HashTableDestroy(r->h, 0);
    if (r->file)
	free(r->file);

    free(r);
}

void refs2id(refs *r, bam_file_t *bfd) {
    int i;
    if (r->ref_id)
	free(r->ref_id);

    r->ref_id = malloc(bfd->nref * sizeof(*r->ref_id));
    for (i = 0; i < bfd->nref; i++) {
	HashItem *hi;
	if ((hi = HashTableSearch(r->h, bfd->ref[i].name, 0))) {
	    r->ref_id[i] = hi->data.p;
	} else {
	    fprintf(stderr, "Unable to find ref name '%s'\n",
		    bfd->ref[i].name);
	}
    }
}

/*
 * Internal part of cram_decode_slice().
 * Generates the sequence, quality and cigar components.
 */
static int cram_decode_seq(cram_container *c, cram_slice *s, block_t *blk,
			   cram_record *cr, bam_file_t *bfd, refs *refs,
			   int bf, int cf, char *seq, char *qual) {
    int prev_pos = 0, f, r = 0, out_sz = 1;
    int seq_pos = 1;
    int cig_len = 0, ref_pos = cr->apos;
    int32_t fn, i32;
    enum cigar_op cig_op = BAM_CMATCH;
    uint32_t *cigar = s->cigar;
    uint32_t ncigar = s->ncigar;
    uint32_t cigar_alloc = s->cigar_alloc;

    r |= c->comp_hdr->FN_codec->decode(s,c->comp_hdr->FN_codec, blk, (char *)&fn, &out_sz);

    ref_pos--; // count from 0
    cr->cigar = ncigar;
    for (f = 0; f < fn; f++) {
	int32_t op, pos;

	if (ncigar+1 >= cigar_alloc) {
	    cigar_alloc = cigar_alloc ? cigar_alloc*2 : 1024;
	    cigar = realloc(cigar, cigar_alloc * sizeof(*cigar));
	}

	r |= c->comp_hdr->FC_codec->decode(s, c->comp_hdr->FC_codec, blk,
					   (char *)&op,  &out_sz);
	r |= c->comp_hdr->FP_codec->decode(s, c->comp_hdr->FP_codec, blk,
					   (char *)&pos, &out_sz);
	pos += prev_pos;

	if (pos > seq_pos) {
	    if (refs) {
		if (ref_pos + pos - seq_pos >= bfd->ref[s->hdr->ref_seq_id].len) {
		    static int whinged = 0;
		    if (!whinged)
			fprintf(stderr, "Ref pos outside of ref sequence boundary\n");
		    whinged = 1;
		} else {
		    memcpy(&seq[seq_pos-1], &refs->ref_id[s->hdr->ref_seq_id][ref_pos],
			   pos - seq_pos);
		}
	    }
#ifdef USE_X
	    if (cig_len && cig_op != BAM_CBASE_MATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    cig_op = BAM_CBASE_MATCH;
#else
	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    cig_op = BAM_CMATCH;
#endif
	    cig_len += pos - seq_pos;
	    ref_pos += pos - seq_pos;
	    seq_pos = pos;
	}

	prev_pos = pos;

	switch(op) {
	case 'S': { // soft clip: IN
	    int32_t out_sz2 = 1;

	    if (cig_len) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    r |= c->comp_hdr->IN_codec
		? c->comp_hdr->IN_codec->decode(s, c->comp_hdr->IN_codec, blk,
						&seq[pos-1], &out_sz2)
		: (seq[pos-1] = 'N', out_sz2 = 1, 0);
	    cigar[ncigar++] = (out_sz2<<4) + BAM_CSOFT_CLIP;
	    cig_op = BAM_CSOFT_CLIP;
	    seq_pos += out_sz2;
	    break;
	}

	case 'X': { // Substitution; BS
#ifdef USE_X
	    if (cig_len && cig_op != BAM_CBASE_MISMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    r |= c->comp_hdr->BS_codec->decode(s, c->comp_hdr->BS_codec, blk,
					       (char *)&i32, &out_sz);
	    seq[pos-1] = 'N'; // FIXME look up BS=i32 value
	    cig_op = BAM_CBASE_MISMATCH;
#else
	    int ref_base;
	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    r |= c->comp_hdr->BS_codec->decode(s, c->comp_hdr->BS_codec, blk,
					       (char *)&i32, &out_sz);
	    if (ref_pos >= bfd->ref[s->hdr->ref_seq_id].len || !refs) {
		seq[pos-1] = 'N';
	    } else {
		ref_base = L[(unsigned char)refs->ref_id[s->hdr->ref_seq_id][ref_pos]];
		seq[pos-1] = c->comp_hdr->substitution_matrix[ref_base][i32];
	    }
	    cig_op = BAM_CMATCH;
#endif
	    cig_len++;
	    seq_pos++;
	    ref_pos++;
	    break;
	}

	case 'D': { // Deletion; DL
	    if (cig_len && cig_op != BAM_CDEL) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    r |= c->comp_hdr->DL_codec->decode(s, c->comp_hdr->DL_codec, blk,
					       (char *)&i32, &out_sz);
	    cig_op = BAM_CDEL;
	    cig_len += i32;
	    ref_pos += i32;
	    //printf("  %d: DL = %d (ret %d)\n", f, i32, r);
	    break;
	}

	case 'I': { // Insertion (several bases); IN
	    int32_t out_sz2 = 1;

	    if (cig_len && cig_op != BAM_CINS) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }

	    r |= c->comp_hdr->IN_codec->decode(s, c->comp_hdr->IN_codec, blk,
					       &seq[pos-1], &out_sz2);
	    cig_op = BAM_CINS;
	    cig_len += out_sz2;
	    seq_pos += out_sz2;
	    //printf("  %d: IN(I) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, dat, r, out_sz2);
	    break;
	}

	case 'i': { // Insertion (single base); BA
	    if (cig_len && cig_op != BAM_CINS) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    r |= c->comp_hdr->BA_codec->decode(s, c->comp_hdr->BA_codec, blk,
					       (char *)&i32, &out_sz);
	    seq[pos-1] = i32;
	    cig_op = BAM_CINS;
	    cig_len++;
	    seq_pos++;
	    //printf("  %d: BA = %c (ret %d)\n", f, i32, r);
	    break;
	}

	case 'B': { // Read base; BA, QS
	    int32_t qc;
#ifdef USE_X
	    if (cig_len && cig_op != BAM_CBASE_MISMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
#else
	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
#endif
	    r  = c->comp_hdr->BA_codec->decode(s, c->comp_hdr->BA_codec, blk,
					       (char *)&i32, &out_sz);
	    r |= c->comp_hdr->QS_codec->decode(s, c->comp_hdr->QS_codec, blk,
					       (char *)&qc, &out_sz);
	    seq[pos-1] = i32;
#ifdef USE_X
	    cig_op = BAM_CBASE_MISMATCH;
#else
	    cig_op = BAM_CMATCH;
#endif
	    cig_len++;
	    seq_pos++;
	    ref_pos++;
	    //printf("  %d: BA/QS(B) = %c/%d (ret %d)\n", f, i32, qc, r);
	    break;
	}

	case 'Q': { // Quality score; QS
	    int32_t qc;
#ifdef USE_X
	    if (cig_len && cig_op != BAM_CBASE_MISMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
#else
	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
#endif
	    r |= c->comp_hdr->QS_codec->decode(s, c->comp_hdr->QS_codec, blk,
					       (char *)&qc, &out_sz);
#ifdef USE_X
	    cig_op = BAM_CBASE_MISMATCH;
#else
	    cig_op = BAM_CMATCH;
#endif
	    cig_len++;
	    seq_pos++;
	    ref_pos++;
	    //printf("  %d: QS = %d (ret %d)\n", f, qc, r);
	}

	default:
	    abort();
	}
    }

    /* An implement match op for any unaccounted for bases */
    if (cr->len >= seq_pos) {
	if (refs) {
	    if (ref_pos + cr->len - seq_pos + 1 > bfd->ref[s->hdr->ref_seq_id].len) {
		static int whinged = 0;
		if (!whinged)
		    fprintf(stderr, "Ref pos outside of ref sequence boundary\n");
		whinged = 1;
	    } else {
		memcpy(&seq[seq_pos-1], &refs->ref_id[s->hdr->ref_seq_id][ref_pos],
		       cr->len - seq_pos + 1);
		ref_pos += cr->len - seq_pos + 1;
	    }
	}

	if (ncigar+1 >= cigar_alloc) {
	    cigar_alloc = cigar_alloc ? cigar_alloc*2 : 1024;
	    cigar = realloc(cigar, cigar_alloc * sizeof(*cigar));
	}
#ifdef USE_X
	if (cig_len && cig_op != BAM_CBASE_MATCH) {
	    cigar[ncigar++] = (cig_len<<4) + cig_op;
	    cig_len = 0;
	}
	cig_op = BAM_CBASE_MATCH;
#else
	if (cig_len && cig_op != BAM_CMATCH) {
	    cigar[ncigar++] = (cig_len<<4) + cig_op;
	    cig_len = 0;
	}
	cig_op = BAM_CMATCH;
#endif
	cig_len += cr->len - seq_pos+1;
    }

    if (cig_len) {
	if (ncigar >= cigar_alloc) {
	    cigar_alloc = cigar_alloc ? cigar_alloc*2 : 1024;
	    cigar = realloc(cigar, cigar_alloc * sizeof(*cigar));
	}

	cigar[ncigar++] = (cig_len<<4) + cig_op;
    }

    cr->ncigar = ncigar - cr->cigar;
    cr->aend = ref_pos;

    //printf("2: %.*s %d .. %d\n", cr->name_len, DSTRING_STR(name_ds) + cr->name, cr->apos, ref_pos);

    r |= c->comp_hdr->MQ_codec->decode(s, c->comp_hdr->MQ_codec, blk,
				       (char *)&cr->mqual, &out_sz);

    if (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES) {
	int32_t out_sz2 = cr->len;

	r |= c->comp_hdr->Qs_codec->decode(s, c->comp_hdr->Qs_codec, blk,
					   qual, &out_sz2);
    } else {
	memset(qual, 30, cr->len);
    }

    s->cigar = cigar;
    s->cigar_alloc = cigar_alloc;
    s->ncigar = ncigar;

    return r;
}

static int cram_decode_aux(cram_container *c, cram_slice *s,
			   block_t *blk, cram_record *cr) {
    int i, r = 0, out_sz = 1;
	    
    r |= c->comp_hdr->TC_codec->decode(s, c->comp_hdr->TC_codec, blk,
				       (char *)&cr->ntags, &out_sz);

    //printf("TC=%d\n", cr->ntags);
    cr->aux_size = 0;
    cr->aux = DSTRING_LEN(s->aux_ds);

    for (i = 0; i < cr->ntags; i++) {
	int32_t id;
	unsigned char tag_data[1024];
	HashItem *hi;

	//printf("Tag %d/%d\n", i+1, cr->ntags);
	r |= c->comp_hdr->TN_codec->decode(s, c->comp_hdr->TN_codec,
					   blk, (char *)&id, &out_sz);
	tag_data[0] = (id>>16) & 0xff;
	tag_data[1] = (id>>8)  & 0xff;
	tag_data[2] = id       & 0xff;

	hi = HashTableSearch(c->comp_hdr->tag_encoding_map,
			     (char *)tag_data, 3);
	if (!hi) {
	    fprintf(stderr, "Unrecognised auxiliary key '%.3s'\n",
		    tag_data);
	    r |= 1;
	} else {
	    cram_map *m = (cram_map *)hi->data.p;
	    int out_sz;

	    r |= m->codec->decode(s, m->codec, blk,
				  (char *)tag_data+3, &out_sz);
	    if (0) {
		int i;
		printf("\tid=%d %.3s\tval = ", id, tag_data);
		for (i=0; i < out_sz; i++) {
		    printf("%02x ", tag_data[i+3]);
		}
		printf("\n");
	    }

	    dstring_nappend(s->aux_ds, (char *)tag_data, out_sz+3);
	    cr->aux_size += out_sz + 3;
	}
    }
    
    return r;
}

int cram_decode_slice(cram_container *c, cram_slice *s, bam_file_t *bfd,
		      refs *refs) {
    cram_block *b = s->block[0];
    int32_t bf, cf, ref_id;
    int out_sz, r = 0;
    block_t *blk = block_create((unsigned char *)b->data, b->uncomp_size);
    int rec;
    char *seq, *qual;
    int unknown_rg = -1, i;

    blk->bit = 7; // MSB first

    /* Look for unknown RG, added as last by Java CRAM? */
    if (bfd->nrg > 0 &&
	!strncmp(bfd->rg_id[bfd->nrg-1], "UNKNOWN", bfd->rg_len[bfd->nrg-1]))
	unknown_rg = bfd->nrg-1;

    cram_init();

    assert(b->content_type == CORE);

    if (s->crecs)
	free(s->crecs);
    s->crecs = malloc(s->hdr->num_records * sizeof(*s->crecs));

    ref_id = s->hdr->ref_seq_id;

    for (rec = 0; rec < s->hdr->num_records; rec++) {
	cram_record *cr = &s->crecs[rec];

	// FIXME: always constant?
	cr->ref_id = ref_id;

	out_sz = 1; /* decode 1 item */
	r |= c->comp_hdr->BF_codec->decode(s, c->comp_hdr->BF_codec, blk,
					   (char *)&bf, &out_sz);
	cr->flags = bam_flag_swap[bf];

	r |= c->comp_hdr->CF_codec->decode(s, c->comp_hdr->CF_codec, blk,
					   (char *)&cf, &out_sz);
	r |= c->comp_hdr->RL_codec->decode(s, c->comp_hdr->RL_codec, blk,
					   (char *)&cr->len, &out_sz);

	r |= c->comp_hdr->AP_codec->decode(s, c->comp_hdr->AP_codec, blk,
					   (char *)&cr->apos, &out_sz);
	cr->apos += s->last_apos;
	s->last_apos=  cr->apos;
		    
	r |= c->comp_hdr->RG_codec->decode(s, c->comp_hdr->RG_codec, blk,
					   (char *)&cr->rg, &out_sz);
	if (cr->rg == unknown_rg)
	    cr->rg = -1;

	cr->name_len = 0;

	if (c->comp_hdr->read_names_included) {
	    int32_t out_sz2 = 1;
	    char *name;

	    // Read directly into dstring
	    dstring_resize(s->name_ds, DSTRING_LEN(s->name_ds) + MAX_NAME_LEN);
	    name = DSTRING_STR(s->name_ds) + DSTRING_LEN(s->name_ds);
	    r |= c->comp_hdr->RN_codec->decode(s, c->comp_hdr->RN_codec, blk,
					       name, &out_sz2);

	    cr->name = DSTRING_LEN(s->name_ds);
	    cr->name_len = out_sz2;

	    /*
	    printf("RN = %d/%.*s\n",
		   cr->name_len, cr->name_len,
		   DSTRING_STR(s->name_ds) + cr->name);
	    */

	    DSTRING_LEN(s->name_ds) += out_sz2;
	}

	cr->mate_line = -1;
	if (cf & CRAM_FLAG_DETACHED) {
	    r |= c->comp_hdr->MF_codec->decode(s, c->comp_hdr->MF_codec, blk,
					       (char *)&cr->mate_flags, &out_sz);

	    if (!c->comp_hdr->read_names_included) {
		int32_t out_sz2 = 1;
		char *name;
	    
		// Read directly into dstring
		dstring_resize(s->name_ds, DSTRING_LEN(s->name_ds) + MAX_NAME_LEN);
		name = DSTRING_STR(s->name_ds) + DSTRING_LEN(s->name_ds);

		r |= c->comp_hdr->RN_codec->decode(s, c->comp_hdr->RN_codec, blk,
						   name, &out_sz2);
		cr->name = DSTRING_LEN(s->name_ds);
		cr->name_len = out_sz2;

		DSTRING_LEN(s->name_ds) += out_sz2;
	    }
		    
	    r |= c->comp_hdr->NS_codec->decode(s, c->comp_hdr->NS_codec, blk,
					       (char *)&cr->mate_ref_id, &out_sz);

	    if (cr->mate_ref_id == -1 && cr->flags & 0x01) {
		/* Paired, but unmapped */
		cr->flags |= 0x08;
	    }

	    r |= c->comp_hdr->NP_codec->decode(s, c->comp_hdr->NP_codec, blk,
					       (char *)&cr->mate_pos, &out_sz);
	    r |= c->comp_hdr->TS_codec->decode(s, c->comp_hdr->TS_codec, blk,
					       (char *)&cr->tlen, &out_sz);
	} else if (cf & CRAM_FLAG_MATE_DOWNSTREAM) {
	    r |= c->comp_hdr->NF_codec->decode(s, c->comp_hdr->NF_codec, blk,
					       (char *)&cr->mate_line, &out_sz);
	    cr->mate_line += rec + 1;

	    //cr->name_len = sprintf(name, "%d", name_id++);
	    //cr->name = DSTRING_LEN(name_ds);
	    //dstring_nappend(name_ds, name, cr->name_len);

	    cr->mate_ref_id = -1;
	    cr->tlen = 0;
	    cr->mate_pos = 0;
	}
	/*
	else if (!name[0]) {
	    //name[0] = '?'; name[1] = 0;
	    //cr->name_len = 1;
	    //cr->name=  DSTRING_LEN(s->name_ds);
	    //dstring_nappend(s->name_ds, "?", 1);

	    cr->mate_ref_id = -1;
	    cr->tlen = 0;
	    cr->mate_pos = 0;
	}
	*/

	/* Auxiliary tags */
	r |= cram_decode_aux(c, s, blk, cr);

	/* Fake up dynamic string growth and appending */
	cr->seq  = DSTRING_LEN(s->seqs_ds);
	cr->qual = DSTRING_LEN(s->qual_ds);
	dstring_resize(s->seqs_ds, DSTRING_LEN(s->seqs_ds) + cr->len);
	dstring_resize(s->qual_ds, DSTRING_LEN(s->qual_ds) + cr->len);
	seq  = DSTRING_STR(s->seqs_ds) + DSTRING_LEN(s->seqs_ds);
	qual = DSTRING_STR(s->qual_ds) + DSTRING_LEN(s->qual_ds);
	DSTRING_LEN(s->seqs_ds) += cr->len;
	DSTRING_LEN(s->qual_ds) += cr->len;

	if (!refs)
	    memset(seq, '=', cr->len);

	if (!(bf & CRAM_FUNMAP)) {
	    /* Decode sequence and generate CIGAR */
	    r |= cram_decode_seq(c, s, blk, cr, bfd, refs, bf, cf, seq, qual);
	} else {
	    int out_sz2 = cr->len;

	    //puts("Unmapped");
	    cr->ncigar = 0;
	    cr->cigar = 0;
	    cr->mqual = 0;

	    r |= c->comp_hdr->BA_codec->decode(s, c->comp_hdr->BA_codec, blk,
					       (char *)seq, &out_sz2);

	    if (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES) {
		out_sz2 = cr->len;
		r |= c->comp_hdr->Qs_codec->decode(s, c->comp_hdr->Qs_codec, blk,
						   qual, &out_sz2);
	    } else {
		memset(qual, 30, cr->len);
	    }
	}
    }

    block_destroy(blk, 1);

    return r;
}

/* -----------------------------------------------------------------------------
 * High level I/O functions for opening, getting or putting sequences and
 * closing CRAM files.
 */

/*
 * Opens a CRAM file for read (mode "rb") or write ("wb").
 * The filename may be "-" to indicate stdin or stdout.
 *
 * Returns file handle on success
 *         NULL on failure.
 */
cram_fd *cram_open(char *filename, char *mode) {
    char *cp;
    cram_fd *fd = calloc(1, sizeof(*fd));
    if (!fd)
	return NULL;

    if (strcmp(filename, "-") == 0) {
	fd->fp = (*mode == 'r') ? stdin : stdout;
    } else {
	fd->fp = fopen(filename, mode);
    }

    if (!fd->fp)
	goto err;
    fd->mode = *mode;

    if (fd->mode == 'r') {
	/* Reader */
	if (!(fd->file_def = cram_read_file_def(fd)))
	    goto err;

	if (!(fd->SAM_hdr = cram_read_SAM_hdr(fd)))
	    goto err;

	//if (0 != parse_SAM_hdr(&fd->SAM_hdr))
	//    goto err;
    } else {
	/* Writer */
	cram_file_def def;
	def.magic[0] = 'C';
	def.magic[1] = 'R';
	def.magic[2] = 'A';
	def.magic[3] = 'M';
	def.major_version = 1;
	def.minor_version = 0;
	memset(def.file_id, 0, 20);
	strncpy(def.file_id, filename, 20);
	if (0 != cram_write_file_def(fd, &def))
	    goto err;

	/* SAM header written later */
    }

    fd->prefix = strdup((cp = strrchr(filename, '/')) ? cp+1 : filename);
    fd->slice_num = 0;
    fd->first_base = fd->last_base = -1;

    fd->ctr = NULL;

    return fd;

 err:
    if (fd->fp)
	fclose(fd->fp);
    if (fd)
	free(fd);

    return NULL;
}


/*
 * Closes a CRAM file.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_close(cram_fd *fd) {
    if (!fd)
	return -1;

    if (fclose(fd->fp) != 0)
	return -1;

    if (fd->file_def)
	cram_free_file_def(fd->file_def);

    if (fd->SAM_hdr)
	cram_free_SAM_hdr(fd->SAM_hdr);

    free(fd->prefix);

    free(fd);
    return 0;
}


/*
 * Sets the read-name prefix for auto-generated sequence names.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_set_prefix(cram_fd *fd, char *prefix) {
    if (!fd)
	return -1;

    if (fd->prefix)
	free(fd->prefix);
    if (!(fd->prefix = strdup(prefix)))
	return -1;

    return 0;
}


/*
 * Converts a cram in-memory record into a bam in-memory record. We
 * pass a pointer to a bam_seq_t pointer along with the a pointer to
 * the allocated size. These can initially be pointers to NULL and zero.
 *
 * This function will reallocate the bam buffer as required and update
 * bam_alloc accordingly, allowing it to be used within a loop
 * efficiently without needing to allocate new bam objects over and
 * over again.
 *
 * Returns the used size of the bam record on success
 *         -1 on failure.
 */
int cram_to_bam(bam_file_t *bfd, cram_fd *fd, cram_slice *s, cram_record *cr,
		int rec, bam_seq_t **bam, size_t *bam_alloc) {
    int bam_idx, bam_len, rg_len, old_idx;
    char *bam_cp;
    char name_a[1024], *name;
    int name_len;

    /* Resolve mate-pair cross-references */
    if (cr->mate_line >= 0) {
	if (cr->mate_line < s->hdr->num_records) {
	    cr->tlen = s->crecs[cr->mate_line].apos >= cr->apos
		?   (s->crecs[cr->mate_line].aend - cr->apos + 1)
		: - (cr->aend - s->crecs[cr->mate_line].apos + 1);
	    cr->mate_pos = s->crecs[cr->mate_line].apos;
	    s->crecs[cr->mate_line].mate_line = rec;
	    cr->mate_ref_id = cr->ref_id;

	    cr->flags |= 0x01; // paired
	    cr->flags &= ~0x08;
	    if (s->crecs[cr->mate_line].flags & 0x10)
		cr->flags |= 0x20;
	    if (cr->flags & 0x10) 
		s->crecs[cr->mate_line].flags |= 0x10;
	} else {
	    fprintf(stderr, "Mate line out of bounds: %d vs [0, %d]\n",
		    cr->mate_line, s->hdr->num_records-1);
	}

	/* FIXME: construct read names here too if needed */
    } else {
	if (cr->mate_flags & 0x10) {
	    cr->flags |= 0x20 | 0x01;
	} else {
	    cr->mate_ref_id = -1;
	}
    }
    
    if (cr->name_len) {
	name = DSTRING_STR(s->name_ds) + cr->name;
	name_len = cr->name_len;
    } else {
	// FIXME: add prefix, container number, slice number, etc
	name = name_a;

	if (cr->mate_line >= 0 && cr->mate_line < rec)
	    name_len = sprintf(name_a, "%s:%"PRId64":%d",
			       fd->prefix, s->id, cr->mate_line);
	else
	    name_len = sprintf(name_a, "%s:%"PRId64":%d",
			       fd->prefix, s->id, rec);
    }

    /* Generate BAM record */
    // FIXME: grow size
    rg_len = (cr->rg != -1) ? bfd->rg_len[cr->rg] + 4 : 0;

    bam_len = cr->name_len + cr->len + (cr->len+1)/2 + 9*36 + cr->ncigar*4
	+ rg_len + cr->aux_size + 1;
    bam_cp = (char *)*bam;
    if (*bam_alloc < bam_len) {
	*bam_alloc = bam_len;
	bam_cp = realloc(bam_cp, *bam_alloc);
	*bam = (bam_seq_t *)bam_cp;
    }
    bam_idx = bam_construct_seq(*bam, *bam_alloc,
				name, name_len,
				cr->flags,
				cr->ref_id,
				cr->apos,
				cr->apos, cr->apos + cr->len, // FIXME
				cr->mqual,
				cr->ncigar, &s->cigar[cr->cigar],
				cr->mate_ref_id,
				cr->mate_pos,
				cr->tlen,
				cr->len,
				DSTRING_STR(s->seqs_ds) + cr->seq,
				DSTRING_STR(s->qual_ds) + cr->qual);
 
    old_idx = bam_idx;

   /* Auxiliary strings */
    if (cr->aux_size != 0) {
	memcpy(&((char *)*bam)[bam_idx], DSTRING_STR(s->aux_ds) + cr->aux,
	       cr->aux_size);
	bam_idx += cr->aux_size;
    }

    /* RG:Z: */
    if (cr->rg != -1) {
	int len = bfd->rg_len[cr->rg];
	((char *)*bam)[bam_idx++] = 'R';
	((char *)*bam)[bam_idx++] = 'G';
	((char *)*bam)[bam_idx++] = 'Z';
	memcpy(&((char *)*bam)[bam_idx], bfd->rg_id[cr->rg], len);
	bam_idx += len;
	((char *)*bam)[bam_idx++] = 0;
    }

    (*bam)->blk_size += bam_idx - old_idx;

    ((char *)*bam)[bam_idx++] = 0;

    return bam_idx;
}


/*
 * Write iterator: put BAM format sequences into a CRAM file.
 * We buffer up a containers worth of data at a time.
 *
 * FIXME: break this into smaller pieces.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_put_bam_seq(cram_fd *fd, bam_seq_t *b) {
    cram_container *c;
    cram_record *cr;
    cram_slice *s;
    int i;
    char *cp, *rg;

    if (!fd->ctr) {
	fd->ctr = cram_new_container(10000, 10);
    }
    c = fd->ctr;

    if (!c->slice || c->curr_rec == c->max_rec ||
	(b->ref != c->curr_ref && c->curr_ref >= -1)) {
	c->curr_ref = b->ref;
	if (c->slice) {
	    s = c->slice;
	    s->hdr->ref_seq_id    = c->curr_ref;
	    s->hdr->ref_seq_start = fd->first_base;
	    s->hdr->ref_seq_span  = fd->last_base - fd->first_base + 1;

	    fprintf(stderr, "End of slice %d, %d/%d..%d\n", c->curr_slice,
		    s->hdr->ref_seq_id, s->hdr->ref_seq_start,
		    s->hdr->ref_seq_start + s->hdr->ref_seq_span-1);

	    if (c->curr_slice == 0) {
		if (c->ref_seq_id != s->hdr->ref_seq_id)
		    c->ref_seq_id  = s->hdr->ref_seq_id;
		c->ref_seq_start = fd->first_base;
	    }

	    c->curr_slice++;
	}

	/* Flush container */
	if (c->curr_slice == c->max_slice) {
	    c->ref_seq_span = fd->last_base - c->ref_seq_start + 1;
	    fprintf(stderr, "Flush container %d/%d..%d\n",
		    c->ref_seq_id, c->ref_seq_start,
		    c->ref_seq_start + c->ref_seq_span -1);

	    /*
	     * Write header
	     * Gather stats
	     * Create container compression header
	     * Write container header
	     * Encode slices
	     * Write slices
	     * Free slices
	     */

	    /* Gather stats */

	    /* Basic strat atm: encoded with NONE */

	    /* Encode slices */
	    for (i = 0; i < c->max_slice; i++) {
		fprintf(stderr, "Encode slice %d\n", i);
	    }
	    

	    // Move to sep func, as we need cram_flush_container for
	    // the closing phase to flush the partial container.
	    for (i = 0; i < c->max_slice; i++) {
		cram_free_slice(c->slices[i]);
		c->slices[i] = NULL;
	    }

	    c->slice = NULL;
	    c->curr_slice = 0;
	}

	fd->first_base = fd->last_base = b->pos;

	/* New slice */
	c->slice = c->slices[c->curr_slice] =
	    cram_new_slice(MAPPED_SLICE, c->max_rec);
	if (!c->slice)
	    return -1;

	c->slice->hdr->ref_seq_start = b->pos;

	c->curr_rec = 0;
    }

    // Create a cram_record
    s = c->slice;
    cr = &s->crecs[c->curr_rec++]; // cache c->slice->crecs as c->rec?

    /* Read group is a mess to implement */
    if ((rg = bam_aux_find(b, "RG"))) {
	HashItem *hi = HashTableSearch(fd->SAM_hdr->rg_hash, rg, 0);
	tag_list_t *tl = hi->data.p;
	cr->rg = hi ? tl->key : -1;
    } else {
	cr->rg = -1;
    }
    //fprintf(stderr, "%s => %d\n", rg ? rg : "\"\"", cr->rg);

    // Fields to resolve later
    //cr->mate_line;    // index to another cram_record
    //cr->mate_flags;   // MF
    //cr->ntags;        // TC

    cr->aux_size = b->blk_size - ((char *)bam_aux(b) - (char *)&b->ref);
    cr->aux = DSTRING_LEN(s->aux_ds);
    dstring_nappend(s->aux_ds, bam_aux(b), cr->aux_size);

    cr->ref_id      = b->ref;
    cr->flags       = bam_flag(b);
    cr->len         = bam_seq_len(b);
    cr->apos        = b->pos;

    cr->name        = DSTRING_LEN(s->name_ds);
    cr->name_len    = bam_name_len(b);
    dstring_nappend(s->name_ds, bam_name(b), bam_name_len(b));

    cr->mate_pos    = b->mate_pos;
    cr->tlen        = b->ins_size;

    cr->seq         = DSTRING_LEN(s->seqs_ds);
    cr->qual        = DSTRING_LEN(s->seqs_ds);
    dstring_resize(s->seqs_ds, DSTRING_LEN(s->seqs_ds) + cr->len);
    dstring_resize(s->qual_ds, DSTRING_LEN(s->seqs_ds) + cr->len);
    cp = DSTRING_STR(s->seqs_ds) + cr->seq;
    for (i = 0; i < cr->len; i++) {
	// do 2 char at a time for efficiency
	cp[i] = bam_nt16_rev_table[bam_seqi(bam_seq(b), i)];
    }

    cp = DSTRING_STR(s->qual_ds) + cr->qual;
    for (i = 0; i < cr->len; i++) {
	cp[i] = bam_qual(b)[i];
    }

    cr->cigar       = s->ncigar;
    cr->ncigar      = bam_cigar_len(b);
    if (cr->cigar + cr->ncigar >= s->cigar_alloc) {
	s->cigar_alloc = s->cigar_alloc ? s->cigar_alloc*2 : 1024;
	s->cigar = realloc(s->cigar, s->cigar_alloc * sizeof(*s->cigar));
    }

    /* Copy and parse */
    {
	int ncigar = cr->ncigar;
	int32_t *cig_to = s->cigar;
	int32_t *cig_from = bam_cigar(b);
	int aend = cr->apos;

	for (i = 0; i < cr->ncigar; i++) {
	    enum cigar_op cig_op = cig_from[i] & BAM_CIGAR_MASK;
	    int cig_len = cig_from[i] > BAM_CIGAR_SHIFT;
	    cig_to[i] = cig_from[i];

	    switch (cig_op) {
	    case BAM_CMATCH:
	    case BAM_CBASE_MATCH:
	    case BAM_CBASE_MISMATCH:
	    case BAM_CDEL:
	    case BAM_CREF_SKIP:
		aend += cig_len;
		break;

	    case BAM_CINS:
	    case BAM_CSOFT_CLIP:
	    case BAM_CHARD_CLIP:
	    case BAM_CPAD:
		break;
	    }
	}
	cr->aend = aend;
    }

    cr->mqual       = bam_map_qual(b);
    cr->mate_ref_id = b->mate_ref;

    c->curr_ctr_rec++;

    // FIXME - use cigar
    if (fd->last_base < cr->aend)
	fd->last_base = cr->aend;

    return 0;
}

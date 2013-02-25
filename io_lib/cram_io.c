/*
 * TODO: 
 *
 * - Add options to control what to store / retrieve.
 *
 * - Investigate generating one block per tag type or encoding tags
 *   via huffman into core?
 *
 * - Replace dstring usage with cram_blocks.
 *
 * - Write an external encoder and use it instead of appending to dstrings.
 *   This means we store the data in the correct place in the code rather than
 *   upfront. But do we need a buffer to copy from anyway?
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
#include <math.h>
#include <ctype.h>

#include "io_lib/cram.h"
#include "io_lib/os.h"
#include "io_lib/deflate_interlaced.h" /* For block_create() */

//Whether CIGAR has just M or uses = and X to indicate match and mismatch
//#define USE_X

/* -----------------------------------------------------------------------------
 * Utility / debugging functions
 */

/* zlib compression code - from Gap5's tg_iface_g.c */
static char *zlib_mem_deflate(char *data, size_t size, size_t *cdata_size,
			      int level, int strat) {
    z_stream s;
    unsigned char *cdata = NULL; /* Compressed output */
    int cdata_alloc = 0;
    int cdata_pos = 0;
    int err;

    cdata = malloc(cdata_alloc = size*1.05+100);
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

    err = deflateInit2(&s, level, Z_DEFLATED, 15|16, 9, strat);
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
static unsigned int bam_flag_swap[0x200], cram_flag_swap[0x800];
static unsigned char L[256];
static unsigned char L2[256];

static char cram_sub_matrix[32][32];

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

    memset(L2, 5, 256);
    L2['A'] = 0; L2['a'] = 0;
    L2['C'] = 1; L2['c'] = 1;
    L2['G'] = 2; L2['g'] = 2;
    L2['T'] = 3; L2['t'] = 3;
    L2['N'] = 4; L2['n'] = 4;

    for (i = 0; i < 0x200; i++) {
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

	bam_flag_swap[i]  = f;
    }

    for (i = 0; i < 0x800; i++) {
	int g = 0;

	if (i & BAM_FPAIRED)	   g |= CRAM_FPAIRED;
	if (i & BAM_FPROPER_PAIR)  g |= CRAM_FPROPER_PAIR;
	if (i & BAM_FUNMAP)        g |= CRAM_FUNMAP;
	if (i & BAM_FREVERSE)      g |= CRAM_FREVERSE;
	if (i & BAM_FREAD1)        g |= CRAM_FREAD1;
	if (i & BAM_FREAD2)        g |= CRAM_FREAD2;
	if (i & BAM_FSECONDARY)    g |= CRAM_FSECONDARY;
	if (i & BAM_FQCFAIL)       g |= CRAM_FQCFAIL;
	if (i & BAM_FDUP)          g |= CRAM_FDUP;

	cram_flag_swap[i] = g;
    }

    memset(cram_sub_matrix, 4, 32*32);
    for (i = 0; i < 32; i++) {
	cram_sub_matrix[i]['A'&0x1f]=0;
	cram_sub_matrix[i]['C'&0x1f]=1;
	cram_sub_matrix[i]['G'&0x1f]=2;
	cram_sub_matrix[i]['T'&0x1f]=3;
	cram_sub_matrix[i]['N'&0x1f]=4;
    }
    for (i = 0; i < 20; i+=4) {
	int j;
	for (j = 0; j < 20; j++) {
	    cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
	    cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
	    cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
	    cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
	}
	cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+0]&0x1f]=0;
	cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+1]&0x1f]=1;
	cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+2]&0x1f]=2;
	cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+3]&0x1f]=3;
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

    int32_t val = getc(fd->fp);
    if (val == -1)
	return -1;

    int i = nbytes[val>>4];
    val &= nbits[val>>4];

    switch(i) {
    case 0:
	*val_p = val;
	return 1;

    case 1:
	val = (val<<8) | (unsigned char)getc(fd->fp);
	*val_p = val;
	return 2;

    case 2:
	val = (val<<8) | (unsigned char)getc(fd->fp);
	val = (val<<8) | (unsigned char)getc(fd->fp);
	*val_p = val;
	return 3;

    case 3:
	val = (val<<8) | (unsigned char)getc(fd->fp);
	val = (val<<8) | (unsigned char)getc(fd->fp);
	val = (val<<8) | (unsigned char)getc(fd->fp);
	*val_p = val;
	return 4;

    case 4: // really 3.5 more, why make it different?
	val = (val<<8) | (unsigned char)getc(fd->fp);
	val = (val<<8) | (unsigned char)getc(fd->fp);
	val = (val<<8) | (unsigned char)getc(fd->fp);
	val = (val<<4) | (((unsigned char)getc(fd->fp)) & 0x0f);
	*val_p = val;
    }

    return 5;
}

/*
 * As above, but decoding from memory
 */
int itf8_get_1(char *cp, int32_t *val_p) {
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

    int32_t val = (unsigned char)*cp++;

    int i = nbytes[val>>4];
    val &= nbits[val>>4];

    switch(i) {
    case 0:
	*val_p = val;
	return 1;

    case 1:
	val = (val<<8) | (unsigned char)*cp++;
	*val_p = val;
	return 2;

    case 2:
	val = (val<<8) | (unsigned char)*cp++;
	val = (val<<8) | (unsigned char)*cp++;
	*val_p = val;
	return 3;

    case 3:
	val = (val<<8) | (unsigned char)*cp++;
	val = (val<<8) | (unsigned char)*cp++;
	val = (val<<8) | (unsigned char)*cp++;
	*val_p = val;
	return 4;

    case 4: // really 4.5, why make it different?
	val = (val<<8) | (unsigned char)*cp++;
	val = (val<<8) | (unsigned char)*cp++;
	val = (val<<8) | (unsigned char)*cp++;
	val = (val<<4) | (((unsigned char)*cp++) & 0x0f);
	*val_p = val;
    }

    return 5;
}

int itf8_get(char *cp, int32_t *val_p) {
    unsigned char *up = (unsigned char *)cp;
    
    if (up[0] < 0x80) {
	*val_p =   up[0];
	return 1;
    } else if (up[0] < 0xc0) {
	*val_p = ((up[0] <<8) |  up[1])                           & 0x3fff;
	return 2;
    } else if (up[0] < 0xe0) {
	*val_p = ((up[0]<<16) | (up[1]<< 8) |  up[2])             & 0x1fffff;
	return 3;
    } else if (up[0] < 0xf0) {
	*val_p = ((up[0]<<24) | (up[1]<<16) | (up[2]<<8) | up[3]) & 0x0fffffff;
	return 4;
    } else {
	*val_p = ((up[0] & 0x0f)<<28) | (up[1]<<20) | (up[2]<<12) | (up[3]<<4) | (up[4] & 0x0f);
	return 5;
    }
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
    } else if (!(val & ~0x01fffff)) { // 3 byte
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
	*cp++ = 0xf0 | ((val>>28) & 0xff);
	*cp++ = (val >> 20) & 0xff;
	*cp++ = (val >> 12) & 0xff;
	*cp++ = (val >> 4 ) & 0xff;
	*cp = val & 0x0f;
	return 5;
    }
}

/*
 * Encodes and writes a single integer in ITF-8 format.
 * Returns 0 on success
 *        -1 on failure
 */
int itf8_encode(cram_fd *fd, int32_t val) {
    char buf[5];
    int len = itf8_put(buf, val);
    return fwrite(buf, 1, len, fd->fp) == len ? 0 : -1;
}


/* -----------------------------------------------------------------------------
 * Mid level I/O functions for manipulating CRAM file structures:
 * Headers, containers, blocks, etc
 */


/*
 * Quick and simple hash lookup for cram_map arrays
 */
static cram_map *map_find(cram_map **map, unsigned char *key, int id) {
    cram_map *m;

    m = map[CRAM_MAP(key[0],key[1])];
    while (m && m->key != id)
	m= m->next;

    assert(m);

    return m;
}

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
    if (1 != fread(&hdr->header_len, 4, 1, fd->fp)) {
	free(hdr);
	return NULL;
    }
    hdr->header_len = le_int4(hdr->header_len);

    /* Alloc and read */
    if (NULL == (hdr->header = malloc(hdr->header_len+100))) {
	free(hdr);
	return NULL;
    }

    if (hdr->header_len != fread(hdr->header, 1, hdr->header_len, fd->fp)) {
	free(hdr);
	return NULL;
    }

    /* Parse */
    bam_parse_header(hdr);
    
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
    cram_SAM_hdr *hdr;
    HashItem *hi;

    if (!(hdr = calloc(1, sizeof(*hdr))))
	return NULL;

    if (NULL == (hdr->header = malloc(len+100))) {
	free(hdr);
	return NULL;
    }

    memcpy(hdr->header, str, len);
    hdr->header_len = len;
    hdr->ref_hash = NULL;
    hdr->rg_hash = NULL;

    bam_parse_header(hdr);

    // If no UNKNOWN read-group, add one.
    if (!(hi = HashTableSearch(hdr->rg_hash, "UNKNOWN", 0))) {
	bam_add_rg(hdr, "UNKNOWN", "UNKNOWN");
    }

    return hdr;
}


/*
 * Writes a CRAM SAM header.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_SAM_hdr(cram_fd *fd, cram_SAM_hdr *hdr) {
    int32_t le_len = le_int4(hdr->header_len);

    /* Length */
    if (1 != fwrite(&le_len, 4, 1, fd->fp))
	return -1;

    /* Text data */
    if (hdr->header_len != fwrite(hdr->header, 1, hdr->header_len, fd->fp))
	return -1;

    fflush(fd->fp);

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
 * Reads a container header plus first block (type COMPRESSION_HEADER).
 * Returns cram_container on success
 *         NULL on failure or no container left (fd->err == 0).
 */
cram_container *cram_read_container(cram_fd *fd) {
    cram_container c2, *c;
    int i, s;
    size_t rd = 0;
    cram_block *b;
    
    fd->err = 0;

    memset(&c2, 0, sizeof(c2));
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
    c->curr_slice = 0;
    c->max_slice = c->num_landmarks;
    c->curr_rec = 0;
    c->max_rec = 0;

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

    c->BF_stats = cram_stats_create();
    c->CF_stats = cram_stats_create();
    c->RN_stats = cram_stats_create();
    c->AP_stats = cram_stats_create();
    c->RG_stats = cram_stats_create();
    c->MQ_stats = cram_stats_create();
    c->NS_stats = cram_stats_create();
    c->NP_stats = cram_stats_create();
    c->TS_stats = cram_stats_create();
    c->MF_stats = cram_stats_create();
    c->NF_stats = cram_stats_create();
    c->RL_stats = cram_stats_create();
    c->FN_stats = cram_stats_create();
    c->FC_stats = cram_stats_create();
    c->FP_stats = cram_stats_create();
    c->DL_stats = cram_stats_create();
    c->BA_stats = cram_stats_create();
    c->QS_stats = cram_stats_create();
    c->BS_stats = cram_stats_create();
    c->TC_stats = cram_stats_create();
    c->TN_stats = cram_stats_create();

    //c->aux_B_stats = cram_stats_create();

    c->tags_used = HashTableCreate(16, HASH_DYNAMIC_SIZE);

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

    if (c->TS_stats) cram_stats_free(c->TS_stats);
    if (c->RG_stats) cram_stats_free(c->RG_stats);
    if (c->FP_stats) cram_stats_free(c->FP_stats);
    if (c->NS_stats) cram_stats_free(c->NS_stats);
    if (c->RN_stats) cram_stats_free(c->RN_stats);
    if (c->CF_stats) cram_stats_free(c->CF_stats);
    if (c->TN_stats) cram_stats_free(c->TN_stats);
    if (c->BA_stats) cram_stats_free(c->BA_stats);
    if (c->TV_stats) cram_stats_free(c->TV_stats);
    if (c->BS_stats) cram_stats_free(c->BS_stats);
    if (c->FC_stats) cram_stats_free(c->FC_stats);
    if (c->BF_stats) cram_stats_free(c->BF_stats);
    if (c->AP_stats) cram_stats_free(c->AP_stats);
    if (c->NF_stats) cram_stats_free(c->NF_stats);
    if (c->MF_stats) cram_stats_free(c->MF_stats);
    if (c->FN_stats) cram_stats_free(c->FN_stats);
    if (c->RL_stats) cram_stats_free(c->RL_stats);
    if (c->DL_stats) cram_stats_free(c->DL_stats);
    if (c->TC_stats) cram_stats_free(c->TC_stats);
    if (c->MQ_stats) cram_stats_free(c->MQ_stats);
    if (c->TM_stats) cram_stats_free(c->TM_stats);
    if (c->IN_stats) cram_stats_free(c->IN_stats);
    if (c->QS_stats) cram_stats_free(c->QS_stats);
    if (c->NP_stats) cram_stats_free(c->NP_stats);

    //if (c->aux_B_stats) cram_stats_free(c->aux_B_stats);
    
    if (c->tags_used) HashTableDestroy(c->tags_used, 0);

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

    //fprintf(stderr, "Block at %d\n", (int)ftell(fd->fp));

    if (-1 == (b->method       = getc(fd->fp))) { free(b); return NULL; }
    if (-1 == (b->content_type = getc(fd->fp))) { free(b); return NULL; }
    if (-1 == itf8_decode(fd, &b->content_id))  { free(b); return NULL; }
    if (-1 == itf8_decode(fd, &b->comp_size))   { free(b); return NULL; }
    if (-1 == itf8_decode(fd, &b->uncomp_size)) { free(b); return NULL; }
    //    fprintf(stderr, "  method %d, ctype %d, cid %d, csize %d, ucsize %d\n",
    //	    b->method, b->content_type, b->content_id, b->comp_size, b->uncomp_size);

    if (b->method == RAW) {
	b->alloc = b->uncomp_size;
	if (!(b->data = malloc(b->uncomp_size))){ free(b); return NULL; }
	if (b->uncomp_size != fread(b->data, 1, b->uncomp_size, fd->fp)) {
	    free(b->data);
	    free(b);
	    return NULL;
	}
    } else {
	b->alloc = b->comp_size;
	if (!(b->data = malloc(b->comp_size)))  { free(b); return NULL; }
	if (b->comp_size != fread(b->data, 1, b->comp_size, fd->fp)) {
	    free(b->data);
	    free(b);
	    return NULL;
	}
    }

    b->orig_method = b->method;
    b->idx = 0;
    b->byte = 0;
    b->bit = 7; // MSB

    return b;
}

/*
 * Writes a CRAM block.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_block(cram_fd *fd, cram_block *b) {
    assert(b->method != RAW || (b->comp_size == b->uncomp_size));

    if (putc(b->method,       fd->fp)   == EOF) return -1;
    if (putc(b->content_type, fd->fp)   == EOF) return -1;
    if (itf8_encode(fd, b->content_id)  ==  -1) return -1;
    if (itf8_encode(fd, b->comp_size)   ==  -1) return -1;
    if (itf8_encode(fd, b->uncomp_size) ==  -1) return -1;

    FILE *tmp = fopen("/dev/null", "wb");
    if (b->method == RAW) {
	fwrite(b->data, 1, b->uncomp_size, tmp); 
	if (b->uncomp_size != fwrite(b->data, 1, b->uncomp_size, fd->fp)) 
	    return -1;
    } else {
	fwrite(b->data, 1, b->comp_size, tmp); 
	if (b->comp_size != fwrite(b->data, 1, b->comp_size, fd->fp)) 
	    return -1;
    }
    fclose(tmp);

    return 0;
}

void cram_free_block(cram_block *b) {
    if (!b)
	return;
    if (b->data)
	free(b->data);
    free(b);
}

cram_block *cram_new_block(enum cram_content_type content_type, int content_id) {
    cram_block *b = malloc(sizeof(*b));
    if (!b)
	return NULL;
    b->method = b->orig_method = RAW;
    b->content_type = content_type;
    b->content_id = content_id;
    b->comp_size = 0;
    b->uncomp_size = 0;
    b->data = NULL;
    b->alloc = 0;
    b->byte = 0;
    b->bit = 7; // MSB

    return b;
}

void cram_uncompress_block(cram_block *b) {
    char *uncomp;
    size_t uncomp_size = 0;

    switch (b->method) {
    case RAW:
	b->uncomp_size = b->comp_size;
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

cram_metrics *cram_new_metrics(void) {
    cram_metrics *m = malloc(sizeof(*m));
    m->m1 = m->m2 = 0;
    m->trial = 2;
    m->next_trial = 100;
    return m;
}

/*
 * Compresses a block using one of two different zlib strategies. If we only
 * want one choice set strat2 to be -1.
 *
 * The logic here is that sometimes Z_RLE does a better job than Z_FILTERED
 * or Z_DEFAULT_STRATEGY on quality data. If so, we'd rather use it as it is
 * significantly faster.
 */
void cram_compress_block(cram_block *b, cram_metrics *metrics,
			 int level,  int strat,
			 int level2, int strat2) {
    char *comp;
    size_t comp_size = 0;

    if (level == 0) {
	b->method = RAW;
	b->comp_size = b->uncomp_size;
	return;
    }

    if (b->method != RAW) {
	fprintf(stderr, "Attempt to compress an already compressed block.\n");
	return;
    }

    if (strat2 >= 0)
	fprintf(stderr, "metrics trial %d, next_trial %d, m1 %d, m2 %d\n",
		metrics->trial, metrics->next_trial, metrics->m1, metrics->m2);

    if (strat2 >= 0 && (metrics->trial || --metrics->next_trial == 0)) {
	char *c1, *c2;
	size_t s1, s2;

	if (metrics->next_trial == 0) {
	    metrics->next_trial = 100;
	    metrics->trial = 2;
	    metrics->m1 = metrics->m2 = 0;
	} else {
	    metrics->trial--;
	}

	c1 = zlib_mem_deflate(b->data, b->uncomp_size, &s1, level, strat);
	c2 = zlib_mem_deflate(b->data, b->uncomp_size, &s2, level2, strat2);
	if (s1 < s2) {
	    fprintf(stderr, "M1 wins\n");
	    comp = c1; comp_size = s1;
	    free(c2);
	    metrics->m1++;
	} else {
	    fprintf(stderr, "M2 wins\n");
	    comp = c2; comp_size = s2;
	    free(c1);
	    metrics->m2++;
	}
    } else if (strat2 >= 0) {
	comp = zlib_mem_deflate(b->data, b->uncomp_size, &comp_size,
				metrics->m1 > metrics->m2 ? level : level2,
				metrics->m1 > metrics->m2 ? strat : strat2);
    } else {
	comp = zlib_mem_deflate(b->data, b->uncomp_size, &comp_size,
				level, strat);
    }

    free(b->data);
    b->data = comp;
    b->method = GZIP;
    b->comp_size = comp_size;

    fprintf(stderr, "Compressed block ID %d from %d to %d\n",
	    b->content_id, b->uncomp_size, b->comp_size);
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

    hdr->preservation_map = HashTableCreate(4, HASH_NONVOLATILE_KEYS |
					    HASH_FUNC_TCL);
    memset(hdr->rec_encoding_map, 0,
	   CRAM_MAP_HASH * sizeof(hdr->rec_encoding_map[0]));
    memset(hdr->tag_encoding_map, 0,
	   CRAM_MAP_HASH * sizeof(hdr->tag_encoding_map[0]));

    if (!hdr->preservation_map) {
	cram_free_compression_header(hdr);
	return NULL;
    }

    /* Initialise defaults for preservation map */
    hdr->mapped_qs_included = 0;
    hdr->unmapped_qs_included = 0;
    hdr->unmapped_placed = 0;
    hdr->qs_included = 0;
    hdr->read_names_included = 0;
    memcpy(hdr->substitution_matrix, "CGTNNAGTNNACTNNACGNNACGTN", 25);

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
	char *key = cp;
	int32_t encoding;
	int32_t size;
	cram_map *m = malloc(sizeof(*m)); // FIXME: use pooled_alloc

	cp += 2;
	cp += itf8_get(cp, &encoding);
	cp += itf8_get(cp, &size);

	// Fill out cram_map purely for cram_dump to dump out.
	m->key = (key[0]<<8)|key[1];
	m->encoding = encoding;
	m->size     = size;
	m->offset   = cp - (char *)b->data;
	m->codec = NULL;

	//printf("%s codes for %.2s\n", cram_encoding2str(encoding), key);

	if (key[0] == 'B' && key[1] == 'F')
	    hdr->BF_codec = cram_decoder_init(encoding, cp, size, E_INT);
	else if (key[0] == 'C' && key[1] == 'F')
	    hdr->CF_codec = cram_decoder_init(encoding, cp, size, E_BYTE);
	else if (key[0] == 'R' && key[1] == 'L')
	    hdr->RL_codec = cram_decoder_init(encoding, cp, size, E_INT);
	else if (key[0] == 'A' && key[1] == 'P')
	    hdr->AP_codec = cram_decoder_init(encoding, cp, size, E_INT);
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
	    hdr->TN_codec = cram_decoder_init(encoding, cp, size, E_INT);
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
	} else if (key[0] == 'T' && key[1] == 'M') {
	} else if (key[0] == 'T' && key[1] == 'V') {
	} else
	    fprintf(stderr, "Unrecognised key: %.2s\n", key);

	cp += size;

	m->next = hdr->rec_encoding_map[CRAM_MAP(key[0], key[1])];
	hdr->rec_encoding_map[CRAM_MAP(key[0], key[1])] = m;
    }
    assert(cp - cp_copy == map_size);


    /* Tag encoding map */
    cp += itf8_get(cp, &map_size); cp_copy = cp;
    cp += itf8_get(cp, &map_count);
    for (i = 0; i < map_count; i++) {
	int32_t encoding;
	int32_t size;
	cram_map *m = malloc(sizeof(*m)); // FIXME: use pooled_alloc
	char *key = cp+1;

	m->key = (key[0]<<16)|(key[1]<<8)|key[2];

	cp += 4; // Strictly ITF8, but this suffices
	cp += itf8_get(cp, &encoding);
	cp += itf8_get(cp, &size);

	m->encoding = encoding;
	m->size     = size;
	m->offset   = cp - (char *)b->data;
	m->codec = cram_decoder_init(encoding, cp, size, E_BYTE_ARRAY);
	
	cp += size;

	m->next = hdr->tag_encoding_map[CRAM_MAP(key[0],key[1])];
	hdr->tag_encoding_map[CRAM_MAP(key[0],key[1])] = m;
    }
    assert(cp - cp_copy == map_size);

    return hdr;
}

static int sub_idx(char *key, char val) {
    int i;

    for (i = 0; *key && *key++ != val; i++);
    return i;
}

/*
 * Encodes a compression header block into a generic cram_block structure.
 *
 * Returns cram_block ptr on success
 *         NULL on failure
 */
cram_block *cram_encode_compression_header(cram_container *c, cram_block_compression_hdr *h) {
    cram_block *cb = cram_new_block(COMPRESSION_HEADER, 0);
    char *buf = malloc(100000); // FIXME, auto-grow this
    char *cp = buf;
    int i, mc;
    //int aux_B_len;
    char map[100000], *mp, cnt_buf[5], aux_B[100000];
    cram_codec *aux_B_codec = NULL;

    if (!cb)
	return NULL;

    //if (c->aux_B_stats->nsamp) {
    //	enum cram_encoding codec = cram_stats_encoding(c->aux_B_stats);
    //	fprintf(stderr, "encoding=%d\n", codec);
    //	aux_B_codec = cram_encoder_init(codec, c->aux_B_stats,
    //					E_BYTE_ARRAY, NULL);
    //	aux_B_len = aux_B_codec->store(aux_B_codec, aux_B, NULL);
    //}

    cp += itf8_put(cp, h->ref_seq_id);
    cp += itf8_put(cp, h->ref_seq_start);
    cp += itf8_put(cp, h->ref_seq_span);
    cp += itf8_put(cp, h->num_records);
    cp += itf8_put(cp, h->num_landmarks);
    for (i = 0; i < h->num_landmarks; i++) {
	cp += itf8_put(cp, h->landmark[i]);
    }

    /* FIXME: should create this when we create the container */
    {
	h->preservation_map = HashTableCreate(4, HASH_NONVOLATILE_KEYS);
	HashData hd;
	hd.i = 0; HashTableAdd(h->preservation_map, "PI", 2, hd, NULL);
	hd.i = 1; HashTableAdd(h->preservation_map, "RN", 2, hd, NULL);
	hd.i = 1; HashTableAdd(h->preservation_map, "UI", 2, hd, NULL);
	hd.i = 1; HashTableAdd(h->preservation_map, "MI", 2, hd, NULL);
    }

    /* Preservation map */
    mp = map; mc = 0;
    if (h->preservation_map) {
        HashItem *hi;
        HashIter *iter = HashTableIterCreate();

        while ((hi = HashTableIterNext(h->preservation_map, iter))) {
            //cram_map *m = hi->data.p;

	    *mp++ = hi->key[0];
	    *mp++ = hi->key[1];
	    switch(CRAM_KEY(hi->key[0], hi->key[1])) {
	    case CRAM_KEY('M','I'):
		*mp++ = hi->data.i;
		break;

	    case CRAM_KEY('U','I'):
		*mp++ = hi->data.i;
		break;

	    case CRAM_KEY('P','I'):
		*mp++ = hi->data.i;
		break;

	    case CRAM_KEY('R','N'):
		*mp++ = hi->data.i;
		break;

	    case CRAM_KEY('S','M'):
		*mp++ =
		    (sub_idx("CGTN", h->substitution_matrix[0][0]) << 6) |
		    (sub_idx("CGTN", h->substitution_matrix[0][1]) << 4) |
		    (sub_idx("CGTN", h->substitution_matrix[0][2]) << 2) |
		    (sub_idx("CGTN", h->substitution_matrix[0][3]) << 0);
		*mp++ =
		    (sub_idx("AGTN", h->substitution_matrix[1][0]) << 6) |
		    (sub_idx("AGTN", h->substitution_matrix[1][1]) << 4) |
		    (sub_idx("AGTN", h->substitution_matrix[1][2]) << 2) |
		    (sub_idx("AGTN", h->substitution_matrix[1][3]) << 0);
		*mp++ =
		    (sub_idx("ACTN", h->substitution_matrix[2][0]) << 6) |
		    (sub_idx("ACTN", h->substitution_matrix[2][1]) << 4) |
		    (sub_idx("ACTN", h->substitution_matrix[2][2]) << 2) |
		    (sub_idx("ACTN", h->substitution_matrix[2][3]) << 0);
		*mp++ =
		    (sub_idx("ACGN", h->substitution_matrix[3][0]) << 6) |
		    (sub_idx("ACGN", h->substitution_matrix[3][1]) << 4) |
		    (sub_idx("ACGN", h->substitution_matrix[3][2]) << 2) |
		    (sub_idx("ACGN", h->substitution_matrix[3][3]) << 0);
		*mp++ =
		    (sub_idx("ACGT", h->substitution_matrix[4][0]) << 6) |
		    (sub_idx("ACGT", h->substitution_matrix[4][1]) << 4) |
		    (sub_idx("ACGT", h->substitution_matrix[4][2]) << 2) |
		    (sub_idx("ACGT", h->substitution_matrix[4][3]) << 0);
		break;

	    default:
		fprintf(stderr, "Unknown preservation key '%.2s'\n", hi->key);
		break;
	    }

	    mc++;
        }

        HashTableIterDestroy(iter);
    }
    cp += itf8_put(cp, mp - map + itf8_put(cnt_buf, mc));
    cp += itf8_put(cp, mc);
    memcpy(cp, map, mp-map);
    cp += mp-map;
    
    /* rec encoding map */
    mp = map;
    mc = 0;
    if (h->BF_codec) mp += h->BF_codec->store(h->BF_codec, mp, "BF"), mc++;
    if (h->CF_codec) mp += h->CF_codec->store(h->CF_codec, mp, "CF"), mc++;
    if (h->RL_codec) mp += h->RL_codec->store(h->RL_codec, mp, "RL"), mc++;
    if (h->AP_codec) mp += h->AP_codec->store(h->AP_codec, mp, "AP"), mc++;
    if (h->RG_codec) mp += h->RG_codec->store(h->RG_codec, mp, "RG"), mc++;
    if (h->MF_codec) mp += h->MF_codec->store(h->MF_codec, mp, "MF"), mc++;
    if (h->NS_codec) mp += h->NS_codec->store(h->NS_codec, mp, "NS"), mc++;
    if (h->NP_codec) mp += h->NP_codec->store(h->NP_codec, mp, "NP"), mc++;
    if (h->TS_codec) mp += h->TS_codec->store(h->TS_codec, mp, "TS"), mc++;
    if (h->NF_codec) mp += h->NF_codec->store(h->NF_codec, mp, "NF"), mc++;
    if (h->TC_codec) mp += h->TC_codec->store(h->TC_codec, mp, "TC"), mc++;
    if (h->TN_codec) mp += h->TN_codec->store(h->TN_codec, mp, "TN"), mc++;
    if (h->FN_codec) mp += h->FN_codec->store(h->FN_codec, mp, "FN"), mc++;
    if (h->FC_codec) mp += h->FC_codec->store(h->FC_codec, mp, "FC"), mc++;
    if (h->FP_codec) mp += h->FP_codec->store(h->FP_codec, mp, "FP"), mc++;
    if (h->BS_codec) mp += h->BS_codec->store(h->BS_codec, mp, "BS"), mc++;
    if (h->IN_codec) mp += h->IN_codec->store(h->IN_codec, mp, "IN"), mc++;
    if (h->DL_codec) mp += h->DL_codec->store(h->DL_codec, mp, "DL"), mc++;
    if (h->BA_codec) mp += h->BA_codec->store(h->BA_codec, mp, "BA"), mc++;
    if (h->MQ_codec) mp += h->MQ_codec->store(h->MQ_codec, mp, "MQ"), mc++;
    if (h->RN_codec) mp += h->RN_codec->store(h->RN_codec, mp, "RN"), mc++;
    if (h->QS_codec) mp += h->QS_codec->store(h->QS_codec, mp, "QS"), mc++;
    if (h->Qs_codec) mp += h->Qs_codec->store(h->Qs_codec, mp, "Qs"), mc++;
    if (h->TM_codec) mp += h->TM_codec->store(h->TM_codec, mp, "TM"), mc++;
    if (h->TV_codec) mp += h->TV_codec->store(h->TV_codec, mp, "TV"), mc++;
    cp += itf8_put(cp, mp - map + itf8_put(cnt_buf, mc));
    cp += itf8_put(cp, mc);
    memcpy(cp, map, mp-map);
    cp += mp-map;
    
    /* tag encoding map */
#if 0
    mp = map; mc = 0;
    if (h->tag_encoding_map) {
        HashItem *hi;
        HashIter *iter = HashTableIterCreate();

        while ((hi = HashTableIterNext(h->tag_encoding_map, iter))) {
            cram_map *m = hi->data.p;
	    mp += itf8_put(mp, (hi->key[0]<<16)|(hi->key[1]<<8)|hi->key[2]);
	    mp += m->codec->store(m->codec, mp, vNULL);
	    mc++;
        }

        HashTableIterDestroy(iter);
    }
#else
    mp = map; mc = 0;
    if (c->tags_used) {
        HashItem *hi;
        HashIter *iter = HashTableIterCreate();

        while ((hi = HashTableIterNext(c->tags_used, iter))) {
	    mc++;
	    mp += itf8_put(mp, (hi->key[0]<<16)|(hi->key[1]<<8)|hi->key[2]);

	    // use block content id 4
	    switch(hi->key[2]) {
	    case 'Z': case 'H':
		// string as byte_array_stop
		*mp++ = 5; // byte_array_stop
		*mp++ = 5;
		*mp++ = '\t';
		*mp++ = 4;
		*mp++ = 0;
		*mp++ = 0;
		*mp++ = 0;
		break;

	    case 'A': case 'c': case 'C':
		// byte array len, 1 byte
		*mp++ = 4;  // byte_array_len
		*mp++ = 9;
		*mp++ = 3; // huffman
		*mp++ = 4;
		*mp++ = 1;
		*mp++ = 1; // 1 byte value
		*mp++ = 1;
		*mp++ = 0;
		*mp++ = 1; // external
		*mp++ = 1;
		*mp++ = 4; // content id
		break;

	    case 's': case 'S':
		// byte array len, 2 byte
		*mp++ = 4;  // byte_array_len
		*mp++ = 9;
		*mp++ = 3; // huffman
		*mp++ = 4;
		*mp++ = 1;
		*mp++ = 2; // 2 byte value
		*mp++ = 1;
		*mp++ = 0;
		*mp++ = 1; // external
		*mp++ = 1;
		*mp++ = 4; // content id
		break;

	    case 'i': case 'I': case 'f':
		// byte array len, 4 byte
		*mp++ = 4; // byte_array_len
		*mp++ = 9;
		*mp++ = 3; // huffman
		*mp++ = 4;
		*mp++ = 1;
		*mp++ = 4; // 4 byte value
		*mp++ = 1;
		*mp++ = 0;
		*mp++ = 1; // external
		*mp++ = 1;
		*mp++ = 4; // content id
		break;

#if 0
	    case 'B':
		// Byte array of variable size so using huffman
		// Not ideal as multiple XX:B:..., YY:B:... ZZ:B:...
		// with different length profiles will all use the same
		// huffman stats. We can deal with optimising this later
		// though.  It'll need one set of statistics per aux type.
		*mp++ = 4; // byte_array_len
		mp += itf8_put(mp, aux_B_len + 3); // len encoding
		memcpy(mp, aux_B, aux_B_len);
		mp += aux_B_len;
		*mp++ = 1; // external
		*mp++ = 1;
		*mp++ = 4; // content id
		break;
#else
	    case 'B':
		// Byte array of variable size, but we generate our tag
		// byte stream at the wrong stage (during reading and not
		// after slice header construction). So we use
		// BYTE_ARRAY_LEN with the length codec being external
		// too.
		*mp++ = 4; // byte_array_len
		*mp++ = 6;
		*mp++ = 1; // len external
		*mp++ = 1;
		*mp++ = 4; // content id;
		*mp++ = 1; // val external
		*mp++ = 1;
		*mp++ = 4; // content id;
		break;
#endif

	    default:
		fprintf(stderr, "Unsupported SAM aux type '%c'\n",
			hi->key[2]);
	    }
	    //mp += m->codec->store(m->codec, mp, NULL);
	}

	HashTableIterDestroy(iter);
    }
#endif
    cp += itf8_put(cp, mp - map + itf8_put(cnt_buf, mc));
    cp += itf8_put(cp, mc);
    memcpy(cp, map, mp-map);
    cp += mp-map;

    fprintf(stderr, "Wrote compression block header in %d bytes\n",
	    (int)(cp-buf));

    cb->data = buf;
    cb->comp_size = cb->uncomp_size = cp - buf;

   //if (aux_B_codec)
   //	aux_B_codec->free(aux_B_codec);

    return cb;
}

/*
 * Creates a new blank container compression header
 *
 * Returns header ptr on success
 *         NULL on failure
 */
cram_block_compression_hdr *cram_new_compression_header(void) {
    cram_block_compression_hdr *hdr = calloc(1, sizeof(*hdr));

    return hdr;
}

int cram_write_container(cram_fd *fd, cram_container *h) {
    char buf[1024], *cp = buf;
    int i;

    cp += itf8_put(cp, h->length);
    cp += itf8_put(cp, h->ref_seq_id);
    cp += itf8_put(cp, h->ref_seq_start);
    cp += itf8_put(cp, h->ref_seq_span);
    cp += itf8_put(cp, h->num_records);
    cp += itf8_put(cp, h->num_blocks);
    cp += itf8_put(cp, h->num_landmarks);
    for (i = 0; i < h->num_landmarks; i++)
	cp += itf8_put(cp, h->landmark[i]);
    if (cp-buf != fwrite(buf, 1, cp-buf, fd->fp))
	return -1;

    return 0;
}


void cram_free_compression_header(cram_block_compression_hdr *hdr) {
    int i;

    if (hdr->landmark)
	free(hdr->landmark);

    if (hdr->preservation_map)
	HashTableDestroy(hdr->preservation_map, 0);

    for (i = 0; i < CRAM_MAP_HASH; i++) {
	cram_map *m, *m2;
	for (m = hdr->rec_encoding_map[i]; m; m = m2) {
	    m2 = m->next;
	    if (m->codec)
		m->codec->free(m->codec);
	    free(m);
	}
    }

    for (i = 0; i < CRAM_MAP_HASH; i++) {
	cram_map *m, *m2;
	for (m = hdr->tag_encoding_map[i]; m; m = m2) {
	    m2 = m->next;
	    if (m->codec)
		m->codec->free(m->codec);
	    free(m);
	}
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

    cp += itf8_get(cp, &hdr->num_content_ids);
    hdr->block_content_ids = malloc(hdr->num_content_ids * sizeof(int32_t));
    for (i = 0; i < hdr->num_content_ids; i++) {
	cp += itf8_get(cp, &hdr->block_content_ids[i]);
    }

    if (b->content_type == MAPPED_SLICE) {
	cp += itf8_get(cp, &hdr->ref_base_id);
    }
    return hdr;
}

/*
 * Encodes a slice compression header. 
 *
 * Returns cram_block on success
 *         NULL on failure
 */
cram_block *cram_encode_slice_header(cram_slice *s) {
    char *buf;
    char *cp;
    cram_block *b = cram_new_block(MAPPED_SLICE, 0);
    int j;

    if (!b)
	return NULL;

    if (NULL == (cp = buf = malloc(5*(6+s->hdr->num_blocks)))) {
	cram_free_block(b);
	return NULL;
    }

    cp += itf8_put(cp, s->hdr->ref_seq_id);
    cp += itf8_put(cp, s->hdr->ref_seq_start);
    cp += itf8_put(cp, s->hdr->ref_seq_span);
    cp += itf8_put(cp, s->hdr->num_records);
    cp += itf8_put(cp, s->hdr->num_blocks);
    cp += itf8_put(cp, s->hdr->num_content_ids);
    for (j = 0; j < s->hdr->num_content_ids; j++) {
	cp += itf8_put(cp, s->hdr->block_content_ids[j]);
    }
    if (s->hdr->content_type == MAPPED_SLICE)
	cp += itf8_put(cp, s->hdr->ref_base_id);
    
    b->data = buf;
    b->comp_size = b->uncomp_size = cp-buf;
    
    return b;
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
    int i, n, max_id;

    if (!b || !s) {
	if (b) cram_free_block(b);
	if (s) free(s);
	return NULL;
    }

    s->hdr_block = b;
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

    for (max_id = i = 0; i < n; i++) {
	s->block[i] = cram_read_block(fd);
	if (s->block[i]->content_type == EXTERNAL &&
	    max_id < s->block[i]->content_id)
	    max_id = s->block[i]->content_id;
    }
    if (max_id < 1024) {
	if (!(s->block_by_id = calloc(max_id+1, sizeof(s->block[0]))))
	    return NULL;

	for (i = 0; i < n; i++) {
	    if (s->block[i]->content_type != EXTERNAL)
		continue;
	    s->block_by_id[s->block[i]->content_id] = s->block[i];
	}
    }

    /* Initialise encoding/decoding tables */
    s->cigar = NULL;
    s->cigar_alloc = 0;
    s->ncigar = 0;

#ifdef DS_SEQ
    s->seqs_ds = dstring_create(NULL);
    s->qual_ds = dstring_create(NULL);
    s->name_ds = dstring_create(NULL);
    s->aux_ds  = dstring_create(NULL);
    s->base_ds = dstring_create(NULL);
#else
    s->seqs_blk = cram_new_block(EXTERNAL, 0);
    s->qual_blk = cram_new_block(EXTERNAL, 1);
    s->name_blk = cram_new_block(EXTERNAL, 2);
    s->aux_blk  = cram_new_block(EXTERNAL, 4);
    s->base_blk = cram_new_block(EXTERNAL, 0);
#endif


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

    s->hdr_block = NULL;
    s->block = NULL;
    s->block_by_id = NULL;
    s->last_apos = 0;
    s->id = 0;
    s->crecs = malloc(nrecs * sizeof(cram_record));
    s->cigar = NULL;
    s->cigar_alloc = 0;
    s->ncigar = 0;

#ifdef DS_SEQ
    s->seqs_ds = dstring_create(NULL);
    s->qual_ds = dstring_create(NULL);
    s->name_ds = dstring_create(NULL);
    s->aux_ds  = dstring_create(NULL);
    s->base_ds = dstring_create(NULL);
#else
    s->seqs_blk = cram_new_block(EXTERNAL, 0);
    s->qual_blk = cram_new_block(EXTERNAL, 1);
    s->name_blk = cram_new_block(EXTERNAL, 2);
    s->aux_blk  = cram_new_block(EXTERNAL, 4);
    s->base_blk = cram_new_block(EXTERNAL, 0);
#endif

    s->features = NULL;
    s->nfeatures = s->afeatures = 0;

#ifndef TN_external
    s->TN = NULL;
    s->nTN = s->aTN = 0;
#endif

    // Volatile keys as we do realloc in dstring
    s->pair = HashTableCreate(10000, HASH_DYNAMIC_SIZE);

#ifdef BA_external
    s->BA_len = 0;
#endif

    return s;
}

void cram_free_slice(cram_slice *s) {
    if (!s)
	return;

    if (s->hdr_block)
	cram_free_block(s->hdr_block);

    if (s->block) {
	int i;

	if (s->hdr) {
	    for (i = 0; i < s->hdr->num_blocks; i++) {
		cram_free_block(s->block[i]);
	    }
	}
	free(s->block);
    }

    if (s->block_by_id)
	free(s->block_by_id);

    if (s->hdr)
	cram_free_slice_header(s->hdr);

#ifdef DS_SEQ
    if (s->seqs_ds)
	dstring_destroy(s->seqs_ds);

    if (s->qual_ds)
	dstring_destroy(s->qual_ds);

    if (s->name_ds)
	dstring_destroy(s->name_ds);

    if (s->aux_ds)
	dstring_destroy(s->aux_ds);

    if (s->base_ds)
	dstring_destroy(s->base_ds);
#else
    if (s->seqs_blk)
	cram_free_block(s->seqs_blk);

    if (s->qual_blk)
	cram_free_block(s->qual_blk);

    if (s->name_blk)
	cram_free_block(s->name_blk);

    if (s->aux_blk)
	cram_free_block(s->aux_blk);

    if (s->base_blk)
	cram_free_block(s->base_blk);
#endif

    if (s->cigar)
	free(s->cigar);

    if (s->crecs)
	free(s->crecs);

    if (s->features)
	free(s->features);

#ifndef TN_external
    if (s->TN)
	free(s->TN);
#endif

    if (s->pair)
	HashTableDestroy(s->pair, 0);

    free(s);
}


cram_stats *cram_stats_create(void) {
    cram_stats *st = calloc(1, sizeof(cram_stats));
    return st;
}

void cram_stats_add(cram_stats *st, int32_t val) {
    st->nsamp++;

    //assert(val >= 0);

    if (val < MAX_STAT_VAL && val >= 0) {
	st->freqs[val]++;
    } else {
	HashItem *hi;

	if (!st->h) {
	    st->h = HashTableCreate(2048, HASH_DYNAMIC_SIZE|HASH_NONVOLATILE_KEYS|HASH_INT_KEYS);
	}

	if ((hi = HashTableSearch(st->h, (char *)val, 4))) {
	    hi->data.i++;
	} else {
	    HashData hd;
	    hd.i = 1;
	    HashTableAdd(st->h, (char *)val, 4, hd, NULL);
	}
    }
}

void cram_stats_del(cram_stats *st, int32_t val) {
    st->nsamp--;

    //assert(val >= 0);

    if (val < MAX_STAT_VAL && val >= 0) {
	st->freqs[val]--;
    } else if (st->h) {
	HashItem *hi;

	if ((hi = HashTableSearch(st->h, (char *)val, 4))) {
	    if (--hi->data.i == 0)
		HashTableDel(st->h, hi, 0);
	} else {
	    fprintf(stderr, "Failed to remove val %d from cram_stats\n", val);
	    st->nsamp++;
	}
    } else {
	fprintf(stderr, "Failed to remove val %d from cram_stats\n", val);
	st->nsamp++;
    }
}

void cram_stats_dump(cram_stats *st) {
    int i;
    fprintf(stderr, "cram_stats:\n");
    for (i = 0; i < MAX_STAT_VAL; i++) {
	if (!st->freqs[i])
	    continue;
	fprintf(stderr, "\t%d\t%d\n", i, st->freqs[i]);
    }
    if (st->h) {
	HashIter *iter=  HashTableIterCreate();
	HashItem *hi;

	while ((hi = HashTableIterNext(st->h, iter))) {
	    fprintf(stderr, "\t%d\t%d\n", (int)hi->key, (int)hi->data.i);
	}
	HashTableIterDestroy(iter);
    }
}

/* Returns the number of bits set in val; it the highest bit used */
static int nbits(int v) {
    static const int MultiplyDeBruijnBitPosition[32] = {
	1, 10, 2, 11, 14, 22, 3, 30, 12, 15, 17, 19, 23, 26, 4, 31,
	9, 13, 21, 29, 16, 18, 25, 8, 20, 28, 24, 7, 27, 6, 5, 32
    };

    v |= v >> 1; // first up to set all bits 1 after the first 1 */
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;

    // DeBruijn magic to find top bit
    return MultiplyDeBruijnBitPosition[(uint32_t)(v * 0x07C4ACDDU) >> 27];
}

#if 0
static int sort_freqs(const void *vp1, const void *vp2) {
    const int i1 = *(const int *)vp1;
    const int i2 = *(const int *)vp2;
    return i1-i2;
}
#endif

/*
 * Computes entropy from integer frequencies for various encoding methods and
 * picks the best encoding.
 *
 * FIXME: we could reuse some of the code here for the actual encoding parameters
 * too. Eg the best 'k' for SUBEXP or the code lengths for huffman.
 *
 * Returns the best codec to use.
 */
enum cram_encoding cram_stats_encoding(cram_stats *st) {
    enum cram_encoding best_encoding = E_NULL;
    int best_size = INT_MAX, bits;
    int nvals, i, ntot = 0, max_val = 0, min_val = INT_MAX, k;
    int *vals = NULL, *freqs = NULL, vals_alloc = 0, *codes;

    //cram_stats_dump(st);

    /* Count number of unique symbols */
    for (nvals = i = 0; i < MAX_STAT_VAL; i++) {
	if (!st->freqs[i])
	    continue;
	if (nvals >= vals_alloc) {
	    vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
	    vals  = realloc(vals,  vals_alloc * sizeof(int));
	    freqs = realloc(freqs, vals_alloc * sizeof(int));
	}
	vals[nvals] = i;
	freqs[nvals] = st->freqs[i];
	ntot += freqs[nvals];
	if (max_val < i) max_val = i;
	if (min_val > i) min_val = i;
	nvals++;
    }
    if (st->h) {
	HashIter *iter=  HashTableIterCreate();
	HashItem *hi;
	int i;

	while ((hi = HashTableIterNext(st->h, iter))) {
	    if (nvals >= vals_alloc) {
		vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
		vals  = realloc(vals,  vals_alloc * sizeof(int));
		freqs = realloc(freqs, vals_alloc * sizeof(int));
	    }
	    i = (int)hi->key;
	    vals[nvals]=i;
	    freqs[nvals] = hi->data.i;
	    ntot += freqs[nvals];
	    if (max_val < i) max_val = i;
	    if (min_val > i) min_val = i;
	    nvals++;
	}
	HashTableIterDestroy(iter);
    }

    st->nvals = nvals;
    assert(ntot == st->nsamp);

    if (nvals <= 1) {
	free(vals);
	free(freqs);
	return E_HUFFMAN;
    }

    /* We only support huffman now anyway... */
    free(vals); free(freqs); return E_HUFFMAN;

    fprintf(stderr, "Range = %d..%d, nvals=%d, ntot=%d\n",
	    min_val, max_val, nvals, ntot);

    /* Theoretical entropy */
    {
	double dbits = 0;
	for (i = 0; i < nvals; i++) {
	    dbits += freqs[i] * log((double)freqs[i]/ntot);
	}
	dbits /= -log(2);
	fprintf(stderr, "Entropy = %f\n", dbits);
    }

#if 0
    /* Unary */
    if (min_val >= 0) {
	for (bits = i = 0; i < nvals; i++)
	    bits += freqs[i]*(vals[i]+1);
	fprintf(stderr, "UNARY   = %d\n", bits);
	if (best_size > bits)
	    best_size = bits, best_encoding = E_NULL; //E_UNARY;
    }

    /* Beta */
    bits = nbits(max_val - min_val) * ntot;
    fprintf(stderr, "BETA    = %d\n", bits);
    if (best_size > bits)
	best_size = bits, best_encoding = E_BETA;

    /* Gamma */
    for (bits = i = 0; i < nvals; i++)
	bits += ((nbits(vals[i]-min_val+1)-1) + nbits(vals[i]-min_val+1)) * freqs[i];
    fprintf(stderr, "GAMMA   = %d\n", bits);
    if (best_size > bits)
	best_size = bits, best_encoding = E_GAMMA;

    /* Subexponential */
    for (k = 0; k < 10; k++) {
	for (bits = i = 0; i < nvals; i++) {
	    if (vals[i]-min_val < (1<<k))
		bits += (1 + k)*freqs[i];
	    else
		bits += (nbits(vals[i]-min_val)*2-k)*freqs[i];
	}

	fprintf(stderr, "SUBEXP%d = %d\n", k, bits);
	if (best_size > bits)
	    best_size = bits, best_encoding = E_SUBEXP;
    }
#endif

    /* byte array len */

    /* byte array stop */

    /* External? Guesswork! */

    /* Huffman */
//    qsort(freqs, nvals, sizeof(freqs[0]), sort_freqs);
//    for (i = 0; i < nvals; i++) {
//	fprintf(stderr, "%d = %d\n", i, freqs[i]);
//	vals[i] = 0;
//    }

    /* Grow freqs to 2*freqs, to store sums */
    /* Vals holds link data */
    freqs = realloc(freqs, 2*nvals*sizeof(*freqs));
    codes = calloc(2*nvals, sizeof(*codes));

    /* Inefficient, use pointers to form chain so we can insert and maintain
     * a sorted list? This is currently O(nvals^2) complexity.
     */
    for (;;) {
	int low1 = INT_MAX, low2 = INT_MAX;
	int ind1 = 0, ind2 = 0;
	for (i = 0; i < nvals; i++) {
	    if (freqs[i] < 0)
		continue;
	    if (low1 > freqs[i]) 
		low2 = low1, ind2 = ind1, low1 = freqs[i], ind1 = i;
	    else if (low2 > freqs[i])
		low2 = freqs[i], ind2 = i;
	}
	if (low2 == INT_MAX)
	    break;

	//fprintf(stderr, "Merge ind %d (%d), %d (%d) = %d+%d, => %d=%d\n",
	//	ind1, vals[ind1], ind2, vals[ind2], low1, low2,
	//	nvals, low1+low2);

	freqs[nvals] = low1 + low2;
	codes[ind1] = nvals;
	codes[ind2] = nvals;
	freqs[ind1] *= -1;
	freqs[ind2] *= -1;
	nvals++;
    }
    nvals = nvals/2+1;

    for (i = 0; i < nvals; i++) {
	int code_len = 0;
	for (k = codes[i]; k; k = codes[k])
	    code_len++;
	codes[i] = code_len;
	freqs[i] *= -1;
	//fprintf(stderr, "%d / %d => %d\n", vals[i], freqs[i], codes[i]);
    }

    for (bits = i = 0; i < nvals; i++) {
	bits += freqs[i] * codes[i];
    }
    fprintf(stderr, "HUFFMAN = %d\n", bits);
    if (best_size >= bits)
	best_size = bits, best_encoding = E_HUFFMAN;
    free(codes);

    free(vals);
    free(freqs);

    return best_encoding;
}

void cram_stats_free(cram_stats *st) {
    if (st->h)
	HashTableDestroy(st->h, 0);
    free(st);
}

/*
 * Encodes a single slice from a container
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int cram_encode_slice(cram_fd *fd, cram_container *c,
			     cram_block_compression_hdr *h, cram_slice *s) {
    int rec, r = 0, last_pos;
    cram_block *core;

    /*
     * Slice external blocks:
     * ID 0 => base calls (insertions, soft-clip)
     * ID 1 => qualities
     * ID 2 => names
     * ID 3 => TS (insert size), NP (next frag)
     * ID 4 => tags
     */

    /* Create cram slice header, num_blocks etc */
    s->hdr->ref_base_id = 0;
#ifdef BA_external
    s->hdr->num_content_ids = 6;
#else
    s->hdr->num_content_ids = 5;
#endif
    s->hdr->num_blocks = s->hdr->num_content_ids+1;
    s->block = calloc(s->hdr->num_blocks, sizeof(s->block[0]));
    s->hdr->block_content_ids = malloc(s->hdr->num_content_ids *
				       sizeof(int32_t));
    s->hdr->block_content_ids[0] = 0;
    s->hdr->block_content_ids[1] = 1;
    s->hdr->block_content_ids[2] = 2;
    s->hdr->block_content_ids[3] = 3;
    s->hdr->block_content_ids[4] = 4;
#ifdef BA_external
    s->hdr->block_content_ids[5] = 5;
#endif

    s->block[0] = cram_new_block(CORE, 0);     // Core 
    s->block[1] = cram_new_block(EXTERNAL, 0); // IN Bases
    s->block[2] = cram_new_block(EXTERNAL, 1); // Qual
    s->block[3] = cram_new_block(EXTERNAL, 2); // Names
    s->block[4] = cram_new_block(EXTERNAL, 3); // TS/NP
    s->block[5] = cram_new_block(EXTERNAL, 4); // Tags
#ifdef BA_external
    s->block[6] = cram_new_block(EXTERNAL, 5); // BA bases
#endif

    core = s->block[0];
		 
    /* Create a formal method for stealing from dstrings! */
    //s->block[4]->data = calloc(5, s->hdr->num_records);
    s->block[4]->data = calloc(10, s->hdr->num_records); // NP TS
    s->block[4]->comp_size = s->block[4]->uncomp_size = 0;

#ifdef BA_external
    s->block[6]->data = calloc(1, s->BA_len);
    s->block[6]->comp_size = s->block[6]->uncomp_size = 0;
#endif

    /* Generate core block */
    s->hdr_block = cram_encode_slice_header(s);

    last_pos = s->hdr->ref_seq_start;
    for (rec = 0; rec < s->hdr->num_records; rec++) {
	cram_record *cr = &s->crecs[rec];
	int32_t i32;
	unsigned char uc;

	//fprintf(stderr, "Encode seq %d, %d/%d FN=%d, %s\n", rec, core->byte, core->bit, cr->nfeature, s->name_ds->str + cr->name);

	//printf("BF=0x%x\n", cr->flags);
	//	    bf = cram_flag_swap[cr->flags];
	i32 = cram_flag_swap[cr->flags & 0x7ff];
	r |= h->BF_codec->encode(s, h->BF_codec, core, (char *)&i32, 1);

	uc = cr->cram_flags;
	r |= h->CF_codec->encode(s, h->CF_codec, core,
				 (char *)&uc, 1);

	r |= h->RL_codec->encode(s, h->RL_codec, core,
				 (char *)&cr->len, 1);

	i32 = cr->apos - last_pos;
	r |= h->AP_codec->encode(s, h->AP_codec, core, (char *)&i32, 1);
	last_pos = cr->apos;

	r |= h->RG_codec->encode(s, h->RG_codec, core,
				 (char *)&cr->rg, 1);

	if (c->comp_hdr->read_names_included) {
	    // FIXME: RN codec
	    // Already stored in block[3].
	}

	if (cr->cram_flags & CRAM_FLAG_DETACHED) {
	    char mf = cr->mate_flags;
	    r |= h->MF_codec->encode(s, h->MF_codec, core, &mf, 1);

	    if (!c->comp_hdr->read_names_included) {
		// FIXME: RN codec
		// Already stored in block[3].
	    }

#ifndef NS_external
	    r |= h->NS_codec->encode(s, h->NS_codec, core,
				     (char *)&cr->mate_ref_id, 1);
#else
	    s->block[4]->uncomp_size +=
		itf8_put(&s->block[4]->data[s->block[4]->uncomp_size],
			 cr->mate_ref_id);
#endif

#ifndef TS_external
	    r |= h->NP_codec->encode(s, h->NP_codec, core,
				     (char *)&cr->mate_pos, 1);

	    r |= h->TS_codec->encode(s, h->TS_codec, core,
				     (char *)&cr->tlen, 1);
#else
	    s->block[4]->uncomp_size +=
		itf8_put(&s->block[4]->data[s->block[4]->uncomp_size],
			 cr->mate_pos);
	    s->block[4]->uncomp_size +=
		itf8_put(&s->block[4]->data[s->block[4]->uncomp_size],
			 cr->tlen);
#endif
	} else if (cr->cram_flags & CRAM_FLAG_MATE_DOWNSTREAM) {
	    r |= h->NF_codec->encode(s, h->NF_codec, core,
				     (char *)&cr->mate_line, 1);
	}

	/* Aux tags */
#if 1
	uc = cr->ntags;
	r |= h->TC_codec->encode(s, h->TC_codec, core, (char *)&uc, 1);
#ifndef TN_external
	{
	    int j;
	    for (j = 0; j < cr->ntags; j++) {
		uint32_t i32 = s->TN[cr->TN_idx + j]; // id
		r |= h->TN_codec->encode(s, h->TN_codec, core,
					 (char *)&i32, 1);
	    }
	}
#endif
#else
	uc = 0; // no tags
	r |= h->TC_codec->encode(s, h->TC_codec, core, (char *)&uc, 1);
#endif

	// qual
	// FIXME: QS codec
	// Already stored in block[2].

	// features (diffs)
	if (!(cr->flags & BAM_FUNMAP)) {
	    int prev_pos = 0, j;

	    r |= h->FN_codec->encode(s, h->FN_codec, core,
				     (char *)&cr->nfeature, 1);
	    for (j = 0; j < cr->nfeature; j++) {
		cram_feature *f = &s->features[cr->feature + j];

		uc = f->X.code;
		r |= h->FC_codec->encode(s, h->FC_codec, core,
					 (char *)&uc, 1);
		i32 = f->X.pos - prev_pos;
		r |= h->FP_codec->encode(s, h->FP_codec, core,
					 (char *)&i32, 1);
		prev_pos = f->X.pos;

		switch(f->X.code) {
		    //char *seq;

		case 'X':
		    //fprintf(stderr, "    FC=%c FP=%d base=%d\n", f->X.code, i32, f->X.base);
		
		    uc = f->X.base;
		    r |= h->BS_codec->encode(s, h->BS_codec, core,
					     (char *)&uc, 1);
		    break;
		case 'S':
		    //seq = DSTRING_STR(s->seqs_ds) + f->S.seq_idx;
		    //r |= h->IN_codec->encode(s, h->IN_codec, core,
		    //			     seq, f->S.len);
		    break;
		case 'I':
		    //seq = DSTRING_STR(s->seqs_ds) + f->S.seq_idx;
		    //r |= h->IN_codec->encode(s, h->IN_codec, core,
		    //			     seq, f->S.len);
		    break;
		case 'i':
		    uc = f->i.base;
#ifdef BA_external
		    s->block[6]->data[s->block[6]->uncomp_size++] = uc;
#else
		    r |= h->BA_codec->encode(s, h->BA_codec, core,
					     (char *)&uc, 1);
#endif
		    //seq = DSTRING_STR(s->seqs_ds) + f->S.seq_idx;
		    //r |= h->IN_codec->encode(s, h->IN_codec, core,
		    //			     seq, 1);
		    break;
		case 'D':
		    i32 = f->D.len;
		    r |= h->DL_codec->encode(s, h->DL_codec, core,
					     (char *)&i32, 1);
		    break;

		case 'B':
//		    // Used when we try to store a non ACGTN base or an N
//		    // that aligns against a non ACGTN reference

		    uc  = f->B.base;
#ifdef BA_external
		    s->block[6]->data[s->block[6]->uncomp_size++] = uc;
#else
		    r |= h->BA_codec->encode(s, h->BA_codec, core,
					     (char *)&uc, 1);
#endif

//                  Already added
//		    uc  = f->B.qual;
//		    r |= h->QS_codec->encode(s, h->QS_codec, core,
//					     (char *)&uc, 1);
		    break;

		default:
		    fprintf(stderr, "unhandled feature code %c\n",
			    f->X.code);
		    return -1;
		}
	    }

	    r |= h->MQ_codec->encode(s, h->MQ_codec, core,
				     (char *)&cr->mqual, 1);
	} else {
#ifdef DS_SEQ
	    char *seq = DSTRING_STR(s->seqs_ds) + cr->seq;
#else
	    char *seq = BLOCK_DATA(s->seqs_blk) + cr->seq;
#endif
#ifdef BA_external
	    memcpy(&s->block[6]->data[s->block[6]->uncomp_size],
		   seq, cr->len);
	    s->block[6]->uncomp_size += cr->len;
#else
	    r |= h->BA_codec->encode(s, h->BA_codec, core, seq, cr->len);
#endif
	}
    }
    s->block[0]->uncomp_size = s->block[0]->byte + (s->block[0]->bit < 7);
    s->block[0]->comp_size = s->block[0]->uncomp_size;

#ifdef DS_SEQ
    s->block[1]->data = s->base_ds->str; s->base_ds->str = NULL;
    s->block[1]->comp_size = s->block[1]->uncomp_size = s->base_ds->length;

    s->block[2]->data = s->qual_ds->str; s->qual_ds->str = NULL;
    s->block[2]->comp_size = s->block[2]->uncomp_size = s->qual_ds->length;

    s->block[3]->data = s->name_ds->str; s->name_ds->str = NULL;
    s->block[3]->comp_size = s->block[3]->uncomp_size = s->name_ds->length;

    s->block[5]->data = s->aux_ds->str; s->aux_ds->str = NULL;
    s->block[5]->comp_size = s->block[5]->uncomp_size = s->aux_ds->length;
#else
    cram_free_block(s->block[1]);
    cram_free_block(s->block[2]);
    cram_free_block(s->block[3]);
    cram_free_block(s->block[5]);
    BLOCK_UPLEN(s->base_blk); s->block[1] = s->base_blk; s->base_blk = NULL;
    BLOCK_UPLEN(s->qual_blk); s->block[2] = s->qual_blk; s->qual_blk = NULL;
    BLOCK_UPLEN(s->name_blk); s->block[3] = s->name_blk; s->name_blk = NULL;
    BLOCK_UPLEN(s->aux_blk);  s->block[5] = s->aux_blk;  s->aux_blk  = NULL;
#endif

    s->block[4]->comp_size = s->block[4]->uncomp_size;
    
#ifdef BA_external
    s->block[6]->comp_size = s->block[6]->uncomp_size;
#endif

    /* Compress the other blocks */
    cram_compress_block(s->block[1], fd->m[0], fd->level,Z_FILTERED, -1, -1);
    cram_compress_block(s->block[2], fd->m[1], fd->level,Z_FILTERED,  1,Z_RLE);
    cram_compress_block(s->block[3], fd->m[2], fd->level,Z_FILTERED, -1, -1);
    cram_compress_block(s->block[4], fd->m[3], fd->level,Z_FILTERED, -1, -1);
    cram_compress_block(s->block[5], fd->m[4], fd->level,Z_FILTERED, -1, -1);
#ifdef BA_external
    cram_compress_block(s->block[6], fd->m[5], fd->level,Z_FILTERED, -1, -1);
#endif

    return r ? -1 : 0;
}

/*
 * Encodes all slices in a container into blocks.
 * Returns 0 on success
 *        -1 on failure
 *
 * FIXME: separate into encode_container and write_container. Ideally
 * we should be able to do read_container / write_container or
 * decode_container / encode_container.
 */
int cram_encode_container(cram_fd *fd, cram_container *c) {
    int i, j, slice_offset;
    cram_block_compression_hdr *h = cram_new_compression_header();
    cram_block *c_hdr;

    /* Complete partial slice header */
    if (c->slice) {
	cram_slice *s = c->slice;
	//s->hdr->ref_seq_id    = c->curr_ref;
	assert(s->hdr->ref_seq_start == fd->first_base);
	s->hdr->ref_seq_start = fd->first_base;
	s->hdr->ref_seq_span  = fd->last_base - fd->first_base + 1;
	s->hdr->num_records   = c->curr_rec;
    }

    c->num_records = 0;
    c->num_blocks = 0;
    c->length = 0;

//    /* Gather stats on slices, for huffman codes etc */
//    for (i = 0; i < c->curr_slice; i++) {
//	cram_slice *s = c->slices[i];
//
//	last_pos = s->hdr->ref_seq_start;
//	for (j = 0; j < s->hdr->num_records; j++) { 
//	    cram_record *cr = &s->crecs[j];
//	}
//    }

    //fprintf(stderr, "=== BF ===\n");
    h->BF_codec = cram_encoder_init(cram_stats_encoding(c->BF_stats),
				    c->BF_stats, E_INT, NULL);

    //fprintf(stderr, "=== CF ===\n");
    h->CF_codec = cram_encoder_init(cram_stats_encoding(c->CF_stats),
				    c->CF_stats, E_BYTE, NULL);

//    fprintf(stderr, "=== RN ===\n");
//    h->RN_codec = cram_encoder_init(cram_stats_encoding(c->RN_stats),
//				    c->RN_stats, E_BYTE_ARRAY, NULL);

    //fprintf(stderr, "=== AP ===\n");
    h->AP_codec = cram_encoder_init(cram_stats_encoding(c->AP_stats),
				    c->AP_stats, E_INT, NULL);

    //fprintf(stderr, "=== RG ===\n");
    h->RG_codec = cram_encoder_init(cram_stats_encoding(c->RG_stats),
				    c->RG_stats, E_INT, NULL);

    //fprintf(stderr, "=== MQ ===\n");
    h->MQ_codec = cram_encoder_init(cram_stats_encoding(c->MQ_stats),
				    c->MQ_stats, E_INT, NULL);

    //fprintf(stderr, "=== NS ===\n");
#ifdef NS_external
    h->NS_codec = cram_encoder_init(E_EXTERNAL, NULL, E_INT, (void *)3);
#else
    h->NS_codec = cram_encoder_init(cram_stats_encoding(c->NS_stats),
				    c->NS_stats, E_INT, NULL);
#endif

    //fprintf(stderr, "=== MF ===\n");
    h->MF_codec = cram_encoder_init(cram_stats_encoding(c->MF_stats),
				    c->MF_stats, E_BYTE, NULL);

#ifdef TS_external
    h->TS_codec = cram_encoder_init(E_EXTERNAL, NULL, E_INT, (void *)3);
    h->NP_codec = cram_encoder_init(E_EXTERNAL, NULL, E_INT, (void *)3);
#else
    //fprintf(stderr, "=== TS ===\n");
    h->TS_codec = cram_encoder_init(cram_stats_encoding(c->TS_stats),
				    c->TS_stats, E_INT, NULL);
    //fprintf(stderr, "=== NP ===\n");
    h->NP_codec = cram_encoder_init(cram_stats_encoding(c->NP_stats),
				    c->NP_stats, E_INT, NULL);
#endif

    //fprintf(stderr, "=== NF ===\n");
    h->NF_codec = cram_encoder_init(cram_stats_encoding(c->NF_stats),
				    c->NF_stats, E_INT, NULL);

    //fprintf(stderr, "=== RL ===\n");
    h->RL_codec = cram_encoder_init(cram_stats_encoding(c->RL_stats),
				    c->RL_stats, E_INT, NULL);

    //fprintf(stderr, "=== FN ===\n");
    h->FN_codec = cram_encoder_init(cram_stats_encoding(c->FN_stats),
				    c->FN_stats, E_INT, NULL);

    //fprintf(stderr, "=== FC ===\n");
    h->FC_codec = cram_encoder_init(cram_stats_encoding(c->FC_stats),
				    c->FC_stats, E_BYTE, NULL);

    //fprintf(stderr, "=== FP ===\n");
    h->FP_codec = cram_encoder_init(cram_stats_encoding(c->FP_stats),
				    c->FP_stats, E_INT, NULL);

    //fprintf(stderr, "=== DL ===\n");
    h->DL_codec = cram_encoder_init(cram_stats_encoding(c->DL_stats),
				    c->DL_stats, E_INT, NULL);

#ifdef BA_external
    h->BA_codec = cram_encoder_init(E_EXTERNAL, NULL, E_BYTE, (void *)5);
#else
    //fprintf(stderr, "=== BA ===\n");
    h->BA_codec = cram_encoder_init(cram_stats_encoding(c->BA_stats),
				    c->BA_stats, E_BYTE, NULL);
#endif

    //fprintf(stderr, "=== BS ===\n");
    h->BS_codec = cram_encoder_init(cram_stats_encoding(c->BS_stats),
				    c->BS_stats, E_BYTE, NULL);

    //fprintf(stderr, "=== TC ===\n");
    h->TC_codec = cram_encoder_init(cram_stats_encoding(c->TC_stats),
				    c->TC_stats, E_BYTE, NULL);

    //fprintf(stderr, "=== TN ===\n");
    if (0) {
	//h->TN_codec = cram_encoder_init(cram_stats_encoding(c->TN_stats),
	//				    c->TN_stats, E_INT, NULL);
	cram_byte_array_len_encoder e;
	e.len_len = 6;
	e.len_dat = (uc *)"\003\004\001\003\001\000";
	e.val_len = 3;
	e.val_dat = (uc *)"\001\001\004\001\003\001\000";
	h->TN_codec = cram_encoder_init(E_BYTE_ARRAY_LEN, NULL,
					E_INT, (void *)&e);
    } else {
#ifdef TN_external
	h->TN_codec = cram_encoder_init(E_EXTERNAL, NULL, E_INT, (void *)4);
#else
	h->TN_codec = cram_encoder_init(cram_stats_encoding(c->TN_stats),
					c->TN_stats, E_INT, NULL);
#endif
    }
    
    if (1) {
	// HACK
	int i2[2] = {0, 0};
	h->IN_codec = cram_encoder_init(E_BYTE_ARRAY_STOP, NULL,
					E_BYTE_ARRAY, (void *)i2);
    } else {
	// Do we know the length in soft-clips? Maybe this is impossible.
	h->IN_codec = cram_encoder_init(E_EXTERNAL, NULL,
					E_BYTE_ARRAY, (void *)0);
    }
    {
	// HACK
	//int i2[2] = {0, 1};
	//h->QS_codec = cram_encoder_init(E_BYTE_ARRAY_STOP, NULL, (void *)i2);
	h->QS_codec = cram_encoder_init(E_EXTERNAL, NULL, E_BYTE, (void *)1);
    }
    {
	// HACK
	int i2[2] = {0, 2};
	h->RN_codec = cram_encoder_init(E_BYTE_ARRAY_STOP, NULL,
					E_BYTE_ARRAY, (void *)i2);
    }


    /* Encode slices */
    for (i = 0; i < c->curr_slice; i++) {
	fprintf(stderr, "Encode slice %d\n", i);
	if (cram_encode_slice(fd, c, h, c->slices[i]) != 0)
	    return -1;
    }

    /* Create compression header */
    {
	h->ref_seq_id    = c->ref_seq_id;
	h->ref_seq_start = c->ref_seq_start;
	h->ref_seq_span  = c->ref_seq_span;
	h->num_records   = c->num_records;
	
	h->mapped_qs_included = 0;   // fixme
	h->unmapped_qs_included = 0; // fixme
	// h->...  fixme
	memcpy(h->substitution_matrix, CRAM_SUBST_MATRIX, 20);

	c_hdr = cram_encode_compression_header(c, h);
    }

    /* Compute landmarks */
    /* Fill out slice landmarks */
    c->num_landmarks = c->curr_slice;
    c->landmark = malloc(c->num_landmarks * sizeof(*c->landmark));
    if (!c->landmark)
	return -1;

    /*
     * Slice offset starts after the first block, so we need to simulate
     * writing it to work out the correct offset
     */
    {
	char tmp[5];
	slice_offset = c_hdr->method == RAW
	    ? c_hdr->uncomp_size
	    : c_hdr->comp_size;
	slice_offset += 2 +
	    itf8_put(tmp, c_hdr->content_id) +
	    itf8_put(tmp, c_hdr->comp_size) +
	    itf8_put(tmp, c_hdr->uncomp_size);
    }

    c->ref_seq_id = c->slices[0]->hdr->ref_seq_id;
    c->ref_seq_start = c->slices[0]->hdr->ref_seq_start;
    c->ref_seq_span = c->slices[0]->hdr->ref_seq_span;
    for (i = 0; i < c->curr_slice; i++) {
	cram_slice *s = c->slices[i];
	
	c->num_records += s->hdr->num_records;
	c->num_blocks += s->hdr->num_blocks + 2;
	c->landmark[i] = slice_offset;

	if (s->hdr->ref_seq_start + s->hdr->ref_seq_span >
	    c->ref_seq_start + c->ref_seq_span) {
	    c->ref_seq_span = s->hdr->ref_seq_start + s->hdr->ref_seq_span
		- c->ref_seq_start;
	}
	
	slice_offset += s->hdr_block->method == RAW
	    ? s->hdr_block->uncomp_size
	    : s->hdr_block->comp_size;

	for (j = 0; j < s->hdr->num_blocks; j++) {
	    slice_offset += s->block[j]->method == RAW
		? s->block[j]->uncomp_size
		: s->block[j]->comp_size;
	}

	c->length += slice_offset;
    }

    /* Write the container header */
    if (0 != cram_write_container(fd, c))
	return -1;


    /*
     * Write the blocks
     */
    cram_write_block(fd, c_hdr);
    cram_free_block(c_hdr);

    for (i = 0; i < c->curr_slice; i++) {
	cram_slice *s = c->slices[i];

	/* No compression for now... */

	cram_write_block(fd, s->hdr_block);

	for (j = 0; j < s->hdr->num_blocks; j++) {
	    if (-1 == cram_write_block(fd, s->block[j]))
		return -1;
	}
    }

    cram_free_compression_header(h);

    return fflush(fd->fp) == 0 ? 0 : -1;
}

/*
 * Adds a feature code to a read within a slice. For purposes of minimising
 * memory allocations and fragmentation we have one array of features for all
 * reads within the slice. We return the index into this array for this new
 * feature.
 *
 * Returns feature index on success
 *         -1 on failure.
 */
static int counter =0;
int cram_add_feature(cram_container *c, cram_slice *s,
		     cram_record *r, cram_feature *f) {
    if (s->nfeatures >= s->afeatures) {
	s->afeatures = s->afeatures ? s->afeatures*2 : 1024;
	s->features = realloc(s->features, s->afeatures * sizeof(*s->features));
	if (!s->features)
	    return -1;
    }

    if (!r->nfeature++) {
	r->feature = s->nfeatures;
	cram_stats_add(c->FP_stats, f->X.pos);
    } else {
	cram_stats_add(c->FP_stats,
		       f->X.pos - s->features[r->feature + r->nfeature-2].X.pos);
    }
    cram_stats_add(c->FC_stats, f->X.code);
    counter++;

    s->features[s->nfeatures++] = *f;

    return 0;
}

static int cram_add_substitution(cram_container *c, cram_slice *s, cram_record *r,
				 int pos, char base, char qual, char ref) {
    cram_feature f;

    // seq=ACGTN vs ref=ACGT or seq=ACGT vs ref=ACGTN
    if (L2[base]<4 || (L2[base]<5 && L2[ref]<4)) {
	f.X.pos = pos+1;
	f.X.code = 'X';
	f.X.base = cram_sub_matrix[ref&0x1f][base&0x1f];
	cram_stats_add(c->BS_stats, f.X.base);
    } else {
	f.B.pos = pos+1;
	f.B.code = 'B';
	f.B.base = base;
	f.B.qual = qual;
	fprintf(stderr, "Add %c %c\n", base, ref);
	cram_stats_add(c->BA_stats, f.B.base);
	cram_stats_add(c->QS_stats, f.B.qual);
#ifdef DS_SEQ
	dstring_append_char(s->qual_ds, qual);
#else
	BLOCK_APPEND_CHAR(s->qual_blk, qual);
#endif
    }
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_deletion(cram_container *c, cram_slice *s, cram_record *r,
			     int pos, int len, char *base) {
    cram_feature f;
    f.D.pos = pos+1;
    f.D.code = 'D';
    f.D.len = len;
    cram_stats_add(c->DL_stats, len);
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_softclip(cram_container *c, cram_slice *s, cram_record *r,
			     int pos, int len, char *base) {
    cram_feature f;
    f.S.pos = pos+1;
    f.S.code = 'S';
    f.S.len = len;
#ifdef DS_SEQ
    f.S.seq_idx = DSTRING_LEN(s->base_ds);
    dstring_nappend(s->base_ds, base, len);
    dstring_append_char(s->base_ds, '\0');
#else
    f.S.seq_idx = BLOCK_SIZE(s->base_blk);
    BLOCK_APPEND(s->base_blk, base, len);
    BLOCK_APPEND_CHAR(s->base_blk, '\0');
#endif
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_insertion(cram_container *c, cram_slice *s, cram_record *r,
			      int pos, int len, char *base) {
    cram_feature f;
    f.I.pos = pos+1;
    if (len == 1) {
	f.i.code = 'i';
	f.i.base = *base;
#ifdef BA_external
	s->BA_len++;
#else
	cram_stats_add(c->BA_stats, *base);
#endif
    } else {
	f.I.code = 'I';
	f.I.len = len;
#ifdef DS_SEQ
	f.I.seq_idx = DSTRING_LEN(s->base_ds);
	dstring_nappend(s->base_ds, base, len);
	dstring_append_char(s->base_ds, '\0');
#else
	f.S.seq_idx = BLOCK_SIZE(s->base_blk);
	BLOCK_APPEND(s->base_blk, base, len);
	BLOCK_APPEND_CHAR(s->base_blk, '\0');
#endif
    }
    return cram_add_feature(c, s, r, &f);
}

/*
 * Loads a reference.
 * FIXME: use the sam_comp.cpp get_ref_base() equivalent to allow
 * random access within an indexed reference instead?
 *
 * Returns a ref_seq structure on success
 *         NULL on failure
 */
static refs *load_reference(char *fn) {
    struct stat sb;
    char *ref;
    FILE *fp;
    off_t i, j = 0;
    char *name = NULL, *seq = NULL;
    HashData hd;
    char fai_fn[PATH_MAX];
    char line[1024];

    refs *r = malloc(sizeof(*r));
    if (!r)
	return NULL;

    /* Open reference, for later use */
    if (stat(fn, &sb) != 0) {
	perror(fn);
	return NULL;
    }

    if (!(r->fp = fopen(fn, "r"))) {
	perror(fn);
	return NULL;
    }

    r->ref_id = NULL;
    r->h_meta = HashTableCreate(16, HASH_DYNAMIC_SIZE | HASH_NONVOLATILE_KEYS);
    //r->h_seq  = HashTableCreate(16, HASH_DYNAMIC_SIZE | HASH_NONVOLATILE_KEYS);

    /* Parse .fai file and load meta-data */
    sprintf(fai_fn, "%.*s.fai", PATH_MAX-5, fn);
    if (stat(fai_fn, &sb) != 0) {
	perror(fai_fn);
	return NULL;
    }
    if (!(fp = fopen(fai_fn, "r"))) {
	perror(fai_fn);
	return NULL;
    }
    while (fgets(line, 1024, fp) != NULL) {
	ref_entry *e = malloc(sizeof(*e));
	char *cp;

	if (!r)
	    return NULL;

	// id
	for (cp = line; *cp && !isspace(*cp); cp++)
	    ;
	*cp++ = 0;
	strncpy(e->name, line, 255); e->name[255] = 0;
	
	// length
	while (*cp && isspace(*cp))
	    cp++;
	e->length = strtoll(cp, &cp, 10);

	// offset
	while (*cp && isspace(*cp))
	    cp++;
	e->offset = strtoll(cp, &cp, 10);

	// bases per line
	while (*cp && isspace(*cp))
	    cp++;
	e->bases_per_line = strtol(cp, &cp, 10);

	// line length
	while (*cp && isspace(*cp))
	    cp++;
	e->line_length = strtol(cp, &cp, 10);

	hd.p = e;
	HashTableAdd(r->h_meta, e->name, strlen(e->name), hd, NULL);
    }

    return r;
}

void free_refs(refs *r) {
    if (r->ref_id)
	free(r->ref_id);
    //if (r->h_seq)
    //    HashTableDestroy(r->h_seq, 0);
    if (r->h_meta)
	HashTableDestroy(r->h_meta, 1);

    if (r->fp)
	fclose(r->fp);

    free(r);
}

/*
 * Indexes references by the order they appear in a BAM file. This may not
 * necessarily be the same order they appear in the fasta reference file.
 */
void refs2id(refs *r, bam_file_t *bfd) {
    int i;
    if (r->ref_id)
	free(r->ref_id);

    r->ref_id = malloc(bfd->nref * sizeof(*r->ref_id));
    for (i = 0; i < bfd->nref; i++) {
	HashItem *hi;
	if ((hi = HashTableSearch(r->h_meta, bfd->ref[i].name, 0))) {
	    r->ref_id[i] = hi->data.p;
	} else {
	    fprintf(stderr, "Unable to find ref name '%s'\n",
		    bfd->ref[i].name);
	}
    }
}

/*
 * Returns a portion of a reference sequence from start to end inclusive.
 * The returned pointer is owned by the cram_file fd and should not be freed
 * by the caller. It is valid only until the next cram_get_ref is called
 * with the same fd parameter (so is thread-safe if given multiple files).
 *
 * To return the entire reference sequence, specify start as 1 and end
 * as 0.
 *
 * Returns reference on success
 *         NULL on failure
 */
char *cram_get_ref(cram_fd *fd, int id, int start, int end) {
    ref_entry *r;
    off_t offset, len;
    char *cp_from, *cp_to;

    //fd->first_base = start;
    //fd->last_base = end;

    if (id < 0) {
	if (fd->ref)
	    free(fd->ref);
	fd->ref = NULL;
	return NULL;
    }

    if (!fd->refs->ref_id[id])
	return NULL;

    if (!(r = fd->refs->ref_id[id])) {
	fprintf(stderr, "No reference found for id %d\n", id);
	return NULL;
    }
    if (end < 1)
	end = r->length;
    if (end >= r->length)
	end  = r->length; 
    assert(start >= 1);
    
    /*
     * Check if asking for same reference. A common issue in sam_to_cram
     */
    if (id == fd->ref_id &&
	start == fd->ref_start &&
	end == fd->ref_end) {
	return fd->ref;
    }

    // Compute location in file
    offset = r->offset + (start-1)/r->bases_per_line * r->line_length + 
	(start-1)%r->bases_per_line;
    if (0 != fseeko(fd->refs->fp, offset, SEEK_SET)) {
	perror("fseeko() on reference file");
	return NULL;
    }

    len = r->offset + (end-1)/r->bases_per_line * r->line_length + 
	(end-1)%r->bases_per_line - offset + 1;

    // Load the data enmasse and then strip whitespace as we go 
    fd->ref = realloc(fd->ref, len);     // FIXME: grow only?

    if (len != fread(fd->ref, 1, len, fd->refs->fp)) {
	perror("fread() on reference file");
	return NULL;
    }

    for (cp_from = cp_to = fd->ref; len; len--, cp_from++) {
	if (!isspace(*cp_from))
	    *cp_to++ = toupper(*cp_from);
    }
    if (cp_to - fd->ref != end-start+1) {
	fprintf(stderr, "Malformed reference file?\n");
	return NULL;
    }

    fd->ref_id    = id;
    fd->ref_start = start;
    fd->ref_end   = end;

    return fd->ref;
}

/*
 * Internal part of cram_decode_slice().
 * Generates the sequence, quality and cigar components.
 */
static int cram_decode_seq(cram_fd *fd, cram_container *c, cram_slice *s,
			   cram_block *blk, cram_record *cr, bam_file_t *bfd,
			   int bf, int cf, char *seq, char *qual) {
    int prev_pos = 0, f, r = 0, out_sz = 1;
    int seq_pos = 1;
    int cig_len = 0, ref_pos = cr->apos;
    int32_t fn, i32;
    enum cigar_op cig_op = BAM_CMATCH;
    uint32_t *cigar = s->cigar;
    uint32_t ncigar = s->ncigar;
    uint32_t cigar_alloc = s->cigar_alloc;
    uint32_t nm = 0, md_dist = 0;
    int orig_aux;
    int decode_md = fd->decode_md;
    char buf[20];

#ifdef DS_SEQ
    if (decode_md) {
	orig_aux = DSTRING_LEN(s->aux_ds);
	dstring_nappend(s->aux_ds, "MDZ", 3);
    }
#else
    if (decode_md) {
	orig_aux = BLOCK_SIZE(s->aux_blk);
	BLOCK_APPEND(s->aux_blk, "MDZ", 3);
    }
#endif

    r |= c->comp_hdr->FN_codec->decode(s,c->comp_hdr->FN_codec, blk, (char *)&fn, &out_sz);

    ref_pos--; // count from 0
    cr->cigar = ncigar;
    for (f = 0; f < fn; f++) {
	int32_t pos;
	char op;

	if (ncigar+1 >= cigar_alloc) {
	    cigar_alloc = cigar_alloc ? cigar_alloc*2 : 1024;
	    cigar = realloc(cigar, cigar_alloc * sizeof(*cigar));
	}

	r |= c->comp_hdr->FC_codec->decode(s, c->comp_hdr->FC_codec, blk,
					   &op,  &out_sz);
	r |= c->comp_hdr->FP_codec->decode(s, c->comp_hdr->FP_codec, blk,
					   (char *)&pos, &out_sz);
	pos += prev_pos;

	if (pos > seq_pos) {
	    if (s->ref && s->hdr->ref_seq_id >= 0) {
		if (ref_pos + pos - seq_pos > bfd->ref[s->hdr->ref_seq_id].len) {
		    static int whinged = 0;
		    if (!whinged)
			fprintf(stderr, "Ref pos outside of ref sequence boundary\n");
		    whinged = 1;
		} else {
		    memcpy(&seq[seq_pos-1], &s->ref[ref_pos - s->ref_start +1],
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
	    md_dist += pos - seq_pos;
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
	    unsigned char base;
#ifdef USE_X
	    if (cig_len && cig_op != BAM_CBASE_MISMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    r |= c->comp_hdr->BS_codec->decode(s, c->comp_hdr->BS_codec, blk,
					       (char *)&base, &out_sz);
	    seq[pos-1] = 'N'; // FIXME look up BS=base value
	    cig_op = BAM_CBASE_MISMATCH;
#else
	    int ref_base;
	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    r |= c->comp_hdr->BS_codec->decode(s, c->comp_hdr->BS_codec, blk,
					       (char *)&base, &out_sz);
	    if (ref_pos >= bfd->ref[s->hdr->ref_seq_id].len || !s->ref) {
		seq[pos-1] = 'N';
	    } else {
		ref_base = L[(unsigned char)s->ref[ref_pos - s->ref_start +1]];
		seq[pos-1] = c->comp_hdr->substitution_matrix[ref_base][base];
		if (decode_md) {
#ifdef DS_SEQ
		    dstring_appendf(s->aux_ds, "%d%c", md_dist,
				    s->ref[ref_pos - s->ref_start +1]);
#else
		    BLOCK_APPENDF_2(s->aux_blk, buf, "%d%c",
				    md_dist, s->ref[ref_pos-s->ref_start +1]);
#endif
		    md_dist = 0;
		}
	    }
	    cig_op = BAM_CMATCH;
#endif
	    nm++;
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
	    if (decode_md) {
#ifdef DS_SEQ
		dstring_appendf(s->aux_ds, "%d^%.*s", md_dist, i32,
				&s->ref[ref_pos - s->ref_start +1]);
#else
		BLOCK_APPENDF_1(s->aux_blk, buf, "%d^", md_dist);
		BLOCK_APPEND(s->aux_blk, &s->ref[ref_pos - s->ref_start +1],
			     i32);
#endif
		md_dist = 0;
	    }
	    cig_op = BAM_CDEL;
	    cig_len += i32;
	    ref_pos += i32;
	    nm      += i32;
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
	    nm      += out_sz2;
	    //printf("  %d: IN(I) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, dat, r, out_sz2);
	    break;
	}

	case 'i': { // Insertion (single base); BA
	    if (cig_len && cig_op != BAM_CINS) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    r |= c->comp_hdr->BA_codec->decode(s, c->comp_hdr->BA_codec, blk,
					       (char *)&seq[pos-1], &out_sz);
	    cig_op = BAM_CINS;
	    cig_len++;
	    seq_pos++;
	    nm++;
	    //printf("  %d: BA = %c (ret %d)\n", f, seq[pos-1], r);
	    break;
	}

	case 'B': { // Read base; BA, QS
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
					       (char *)&seq[pos-1], &out_sz);
	    r |= c->comp_hdr->QS_codec->decode(s, c->comp_hdr->QS_codec, blk,
					       (char *)&qual[pos-1], &out_sz);
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
					       (char *)&qual[pos-1], &out_sz);
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
	if (s->ref) {
	    if (ref_pos + cr->len - seq_pos + 1 > bfd->ref[s->hdr->ref_seq_id].len) {
		static int whinged = 0;
		if (!whinged)
		    fprintf(stderr, "Ref pos outside of ref sequence boundary\n");
		whinged = 1;
	    } else {
		memcpy(&seq[seq_pos-1], &s->ref[ref_pos - s->ref_start +1],
		       cr->len - seq_pos + 1);
		ref_pos += cr->len - seq_pos + 1;
		md_dist += cr->len - seq_pos + 1;
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
    if (decode_md) {
#ifdef DS_SEQ
	dstring_appendf(s->aux_ds, "%d", md_dist);
#else
	BLOCK_APPENDF_1(s->aux_blk, buf, "%d", md_dist);
#endif
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

    if (decode_md) {
	char buf[7];
#ifdef DS_SEQ
	dstring_append_char(s->aux_ds, '\0'); // null terminate MD:Z:
	cr->aux_size += DSTRING_LEN(s->aux_ds) - orig_aux;
	buf[0] = 'N'; buf[1] = 'M'; buf[2] = 'I';
	buf[3] = (nm>> 0) & 0xff;
	buf[4] = (nm>> 8) & 0xff;
	buf[5] = (nm>>16) & 0xff;
	buf[6] = (nm>>24) & 0xff;
	dstring_nappend(s->aux_ds, buf, 7);
	cr->aux_size += 7;
#else
	BLOCK_APPEND_CHAR(s->aux_blk, '\0'); // null terminate MD:Z:
	cr->aux_size += BLOCK_SIZE(s->aux_blk) - orig_aux;
	buf[0] = 'N'; buf[1] = 'M'; buf[2] = 'I';
	buf[3] = (nm>> 0) & 0xff;
	buf[4] = (nm>> 8) & 0xff;
	buf[5] = (nm>>16) & 0xff;
	buf[6] = (nm>>24) & 0xff;
	BLOCK_APPEND(s->aux_blk, buf, 7);
	cr->aux_size += 7;
#endif
    }

    return r;
}

static int cram_decode_aux(cram_container *c, cram_slice *s,
			   cram_block *blk, cram_record *cr) {
    int i, r = 0, out_sz = 1;
    unsigned char ntags;
	    
    r |= c->comp_hdr->TC_codec->decode(s, c->comp_hdr->TC_codec, blk,
				       (char *)&ntags, &out_sz);
    cr->ntags = ntags;

    //printf("TC=%d\n", cr->ntags);
    cr->aux_size = 0;
#ifdef DS_SEQ
    cr->aux = DSTRING_LEN(s->aux_ds);
#else
    cr->aux = BLOCK_SIZE(s->aux_blk);
#endif

    for (i = 0; i < cr->ntags; i++) {
	int32_t id, out_sz = 1;
	unsigned char tag_data[1024], key[2];;
	HashItem *hi;
	cram_map *m;

	//printf("Tag %d/%d\n", i+1, cr->ntags);
	r |= c->comp_hdr->TN_codec->decode(s, c->comp_hdr->TN_codec,
					   blk, (char *)&id, &out_sz);
	if (out_sz == 3) {
	    tag_data[0] = ((char *)&id)[0];
	    tag_data[1] = ((char *)&id)[1];
	    tag_data[2] = ((char *)&id)[2];
	} else {
	    tag_data[0] = (id>>16) & 0xff;
	    tag_data[1] = (id>>8)  & 0xff;
	    tag_data[2] = id       & 0xff;
	} 

	m = map_find(c->comp_hdr->tag_encoding_map, tag_data, id);
	assert(m);

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

#ifdef DS_SEQ
	dstring_nappend(s->aux_ds, (char *)tag_data, out_sz+3);
#else
	BLOCK_APPEND(s->aux_blk, (char *)tag_data, out_sz+3);
#endif
	cr->aux_size += out_sz + 3;
    }
    
    return r;
}

int cram_decode_slice(cram_fd *fd, cram_container *c, cram_slice *s,
		      bam_file_t *bfd) {
    cram_block *blk = s->block[0];
    int32_t bf, ref_id;
    unsigned char cf;
    int out_sz, r = 0;
    int rec;
    char *seq, *qual;
    int unknown_rg = -1;

    blk->bit = 7; // MSB first

    /* Look for unknown RG, added as last by Java CRAM? */
    if (bfd->nrg > 0 &&
	!strncmp(bfd->rg_id[bfd->nrg-1], "UNKNOWN", bfd->rg_len[bfd->nrg-1]))
	unknown_rg = bfd->nrg-1;

    assert(blk->content_type == CORE);

    if (s->crecs)
	free(s->crecs);
    s->crecs = malloc(s->hdr->num_records * sizeof(*s->crecs));

    ref_id = s->hdr->ref_seq_id;

#if 0
    s->ref = cram_get_ref(fd, s->hdr->ref_seq_id, s->hdr->ref_seq_start,
			  s->hdr->ref_seq_start + s->hdr->ref_seq_span -1);
    s->ref_start = s->hdr->ref_seq_start;
#else
    // Avoid Java cramtools bug by loading entire reference seq
    s->ref = cram_get_ref(fd, s->hdr->ref_seq_id, 1, 0);
    s->ref_start = 1;
#endif

    for (rec = 0; rec < s->hdr->num_records; rec++) {
	cram_record *cr = &s->crecs[rec];

	//fprintf(stderr, "Decode seq %d, %d/%d\n", rec, blk->byte, blk->bit);

	// FIXME: always constant?
	cr->ref_id = ref_id;
	cr->s = s;

	out_sz = 1; /* decode 1 item */
	r |= c->comp_hdr->BF_codec->decode(s, c->comp_hdr->BF_codec, blk,
					   (char *)&bf, &out_sz);
	cr->flags = bam_flag_swap[bf];

	r |= c->comp_hdr->CF_codec->decode(s, c->comp_hdr->CF_codec, blk,
					   (char *)&cf, &out_sz);
	cr->cram_flags = cf;

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

#ifdef DS_SEQ
	    // Read directly into dstring
	    DSTRING_RESIZE(s->name_ds, DSTRING_LEN(s->name_ds) + MAX_NAME_LEN);
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
#else
	    // Read directly into name cram_block
	    BLOCK_GROW(s->name_blk, MAX_NAME_LEN);
	    name = BLOCK_END(s->name_blk);
	    r |= c->comp_hdr->RN_codec->decode(s, c->comp_hdr->RN_codec, blk,
					       name, &out_sz2);

	    cr->name = BLOCK_SIZE(s->name_blk);
	    BLOCK_SIZE(s->name_blk) += out_sz2;
	    cr->name_len = out_sz2;
#endif
	}

	cr->mate_line = -1;
	cr->mate_ref_id = -1;
	if (cf & CRAM_FLAG_DETACHED) {
	    unsigned char mf;
	    r |= c->comp_hdr->MF_codec->decode(s, c->comp_hdr->MF_codec, blk,
					       (char *)&mf, &out_sz);
	    cr->mate_flags = mf;

	    if (!c->comp_hdr->read_names_included) {
		int32_t out_sz2 = 1;
		char *name;
	    
#ifdef DS_SEQ
		// Read directly into dstring
		DSTRING_RESIZE(s->name_ds, DSTRING_LEN(s->name_ds) + MAX_NAME_LEN);
		name = DSTRING_STR(s->name_ds) + DSTRING_LEN(s->name_ds);

		r |= c->comp_hdr->RN_codec->decode(s, c->comp_hdr->RN_codec, blk,
						   name, &out_sz2);
		cr->name = DSTRING_LEN(s->name_ds);
		cr->name_len = out_sz2;

		DSTRING_LEN(s->name_ds) += out_sz2;
#else
		// Read directly into name cram_block
		BLOCK_GROW(s->name_blk, MAX_NAME_LEN);
		name = BLOCK_END(s->name_blk);
		r |= c->comp_hdr->RN_codec->decode(s, c->comp_hdr->RN_codec, blk,
						   name, &out_sz2);

		cr->name = BLOCK_SIZE(s->name_blk);
		BLOCK_SIZE(s->name_blk) += out_sz2;
		cr->name_len = out_sz2;
#endif
	    }
		    
	    r |= c->comp_hdr->NS_codec->decode(s, c->comp_hdr->NS_codec, blk,
					       (char *)&cr->mate_ref_id, &out_sz);

// Skip as mate_ref of "*" is legit. It doesn't mean unmapped, just unknown.
//	    if (cr->mate_ref_id == -1 && cr->flags & 0x01) {
//		/* Paired, but unmapped */
//		cr->flags |= BAM_FMUNMAP;
//	    }

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
	    cr->tlen = -1;
	    cr->mate_pos = 0;
	} else {
	    cr->tlen = -1;
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
#ifdef DS_SEQ
	cr->qual = DSTRING_LEN(s->qual_ds);
	DSTRING_RESIZE(s->qual_ds, DSTRING_LEN(s->qual_ds) + cr->len);
	qual = DSTRING_STR(s->qual_ds) + DSTRING_LEN(s->qual_ds);
	DSTRING_LEN(s->qual_ds) += cr->len;

	cr->seq  = DSTRING_LEN(s->seqs_ds);
	DSTRING_RESIZE(s->seqs_ds, DSTRING_LEN(s->seqs_ds) + cr->len);
	seq  = DSTRING_STR(s->seqs_ds) + DSTRING_LEN(s->seqs_ds);
	DSTRING_LEN(s->seqs_ds) += cr->len;
#else
	cr->seq = BLOCK_SIZE(s->seqs_blk);
	BLOCK_GROW(s->seqs_blk, cr->len);
	seq = BLOCK_END(s->seqs_blk);
	BLOCK_SIZE(s->seqs_blk) += cr->len;

	cr->qual = BLOCK_SIZE(s->qual_blk);
	BLOCK_GROW(s->qual_blk, cr->len);
	qual = BLOCK_END(s->qual_blk);
	BLOCK_SIZE(s->qual_blk) += cr->len;
#endif

	if (!s->ref)
	    memset(seq, '=', cr->len);

	if (!(bf & CRAM_FUNMAP)) {
	    /* Decode sequence and generate CIGAR */
	    r |= cram_decode_seq(fd, c, s, blk, cr, bfd, bf, cf, seq, qual);
	} else {
	    int out_sz2 = cr->len;

	    //puts("Unmapped");
	    cr->cigar = 0;
	    cr->ncigar = 0;
	    cr->aend = -1;
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
    int i;
    char *cp;
    cram_fd *fd = calloc(1, sizeof(*fd));
    if (!fd)
	return NULL;

    fd->level = 5;
    if (strlen(mode) > 2 && mode[2] >= '0' && mode[2] <= '9')
	fd->level = mode[2] - '0';

    cram_init();

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
    fd->refs = NULL; // refs meta-data structure
    fd->ref  = NULL; // current ref as char*

    fd->decode_md = 0;

    for (i = 0; i < 6; i++)
	fd->m[i] = cram_new_metrics();

    return fd;

 err:
    if (fd->fp)
	fclose(fd->fp);
    if (fd)
	free(fd);

    return NULL;
}

/* 
 * Sets options on the cram_fd. See CRAM_OPT_* definitions in cram_structs.h.
 * Use this immediately after opening.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_set_option(cram_fd *fd, enum cram_option opt, cram_opt *val) {
    switch (opt) {
    case CRAM_OPT_DECODE_MD:
	fd->decode_md = val->i;
	break;

    case CRAM_OPT_PREFIX:
	if (fd->prefix)
	    free(fd->prefix);
	if (!(fd->prefix = strdup(val->s)))
	    return -1;
	break;
    }

    return 0;
}


void cram_load_reference(cram_fd *fd, char *fn) {
    fd->refs = load_reference(fn);
    refs2id(fd->refs, fd->SAM_hdr);
}

/*
 * Closes a CRAM file.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_close(cram_fd *fd) {
    int i;

    if (!fd)
	return -1;

    if (fd->mode == 'w' && fd->ctr) {
	if(fd->ctr->slice)
	    fd->ctr->curr_slice++;
	if (-1 == cram_encode_container(fd, fd->ctr))
	    return -1;
    }

    if (fclose(fd->fp) != 0)
	return -1;

    if (fd->file_def)
	cram_free_file_def(fd->file_def);

    if (fd->SAM_hdr)
	cram_free_SAM_hdr(fd->SAM_hdr);

    free(fd->prefix);

    if (fd->ctr)
	cram_free_container(fd->ctr);

    if (fd->refs)
	free_refs(fd->refs);
    if (fd->ref)
	free(fd->ref);

    for (i = 0; i < 6; i++)
	if (fd->m[i])
	    free(fd->m[i]);

    free(fd);
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
	    /*
	     * On the first read, loop through computing lengths.
	     * It's not perfect as we have one slice per reference so we
	     * cannot detect when TLEN should be zero due to seqs that
	     * map to multiple references.
	     *
	     * We also cannot set tlen correct when it spans a slice for
	     * other reasons. This may make tlen too small. Should we
	     * fix this by forcing TLEN to be stored verbatim in such cases?
	     *
	     * Or do we just admit defeat and output 0 for tlen? It's the
	     * safe option...
	     */
	    if (cr->tlen == -1) {
		int id1 = rec, id2 = rec;
		int apos = cr->apos, aend;
		int tlen;

		do {
		    aend = s->crecs[id2].aend;
		    if (s->crecs[id2].mate_line == -1) {
			s->crecs[id2].mate_line = rec;
			break;
		    }
		    id2 = s->crecs[id2].mate_line;
		} while (id2 != id1);

		tlen = aend - apos + 1;
		id1 = id2 = rec;

		// leftmost is +ve, rightmost -ve, all others undefined
		s->crecs[id2].tlen = tlen;
		tlen *= -1;
		id2 = s->crecs[id2].mate_line;
		while (id2 != id1) {
		    s->crecs[id2].tlen = tlen;
		    id2 = s->crecs[id2].mate_line;
		}
	    }

	    cr->mate_pos = s->crecs[cr->mate_line].apos;
	    //if (s->crecs[cr->mate_line].mate_line == -1)
	    //	s->crecs[cr->mate_line].mate_line = rec;
	    cr->mate_ref_id = cr->ref_id;

	    // paired
	    cr->flags |= BAM_FPAIRED;
	    s->crecs[cr->mate_line].flags |= BAM_FPAIRED;

	    // set mate unmapped if needed
	    if (s->crecs[cr->mate_line].flags & BAM_FUNMAP) {
		cr->flags |= BAM_FMUNMAP;
		cr->tlen = s->crecs[cr->mate_line].tlen = 0;
	    }
	    if (cr->flags & BAM_FUNMAP) {
		s->crecs[cr->mate_line].flags |= BAM_FMUNMAP;
		cr->tlen = s->crecs[cr->mate_line].tlen = 0;
	    }

	    // set mate reversed if needed
	    if (s->crecs[cr->mate_line].flags & BAM_FREVERSE)
		cr->flags |= BAM_FMREVERSE;
	    if (cr->flags & BAM_FREVERSE) 
		s->crecs[cr->mate_line].flags |= BAM_FMREVERSE;
	} else {
	    fprintf(stderr, "Mate line out of bounds: %d vs [0, %d]\n",
		    cr->mate_line, s->hdr->num_records-1);
	}

	/* FIXME: construct read names here too if needed */
    } else {
	if (cr->mate_flags & CRAM_M_REVERSE) {
	    cr->flags |= BAM_FPAIRED | BAM_FMREVERSE;
	} else if (cr->mate_flags & CRAM_M_UNMAP) {
	    cr->flags |= BAM_FMUNMAP;
	    //cr->mate_ref_id = -1;
	}
    }
    
    if (cr->name_len) {
#ifdef DS_SEQ
	name = DSTRING_STR(s->name_ds) + cr->name;
#else
	name = BLOCK_DATA(s->name_blk) + cr->name;
#endif
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
#ifdef DS_SEQ
				DSTRING_STR(s->seqs_ds) + cr->seq,
				DSTRING_STR(s->qual_ds) + cr->qual);
#else
                                BLOCK_DATA(s->seqs_blk) + cr->seq,
				BLOCK_DATA(s->qual_blk) + cr->qual);
#endif
 
    old_idx = bam_idx;

   /* Auxiliary strings */
    if (cr->aux_size != 0) {
#ifdef DS_SEQ
	memcpy(&((char *)*bam)[bam_idx], DSTRING_STR(s->aux_ds) + cr->aux,
	       cr->aux_size);
#else
	memcpy(&((char *)*bam)[bam_idx], BLOCK_DATA(s->aux_blk) + cr->aux,
	       cr->aux_size);
#endif
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
 * Encodes auxiliary data.
 * Returns the read-group parsed out of the BAM aux fields on success
 *         NULL on failure or no rg present (FIXME)
 */
static char *cram_encode_aux(cram_fd *fd, bam_seq_t *b, cram_container *c,
			     cram_slice *s, cram_record *cr) {
    char *aux, *tmp, *rg = NULL;
    int aux_size = b->blk_size - ((char *)bam_aux(b) - (char *)&b->ref);
	
    /* Worst case is 1 nul char on every ??:Z: string, so +33% */
#ifdef DS_SEQ
    DSTRING_RESIZE(s->aux_ds, DSTRING_LEN(s->aux_ds) + aux_size*1.34 + 1);
    tmp = DSTRING_STR(s->aux_ds) + DSTRING_LEN(s->aux_ds);
#else
    BLOCK_GROW(s->aux_blk, aux_size*1.34+1);
    tmp = BLOCK_END(s->aux_blk);
#endif

    aux = bam_aux(b);
#ifndef TN_external
    cr->TN_idx = s->nTN;
#endif
    while (aux[0] != 0) {
	HashData hd; hd.i = 0;
	int32_t i32;

	if (aux[0] == 'R' && aux[1] == 'G' && aux[2] == 'Z') {
	    rg = &aux[3];
	    while (*aux++);
	    continue;
	}
	if (aux[0] == 'M' && aux[1] == 'D' && aux[2] == 'Z') {
	    while (*aux++);
	    continue;
	}
	if (aux[0] == 'N' && aux[1] == 'M') {
	    switch(aux[2]) {
	    case 'A': case 'C': case 'c': aux+=4; break;
	    case 'I': case 'i': case 'f': aux+=7; break;
	    default:
		fprintf(stderr, "Unhandled type code for NM tag\n");
		return NULL;
	    }
	    continue;
	}

	cr->ntags++;
	// replace with fast hash too
	HashTableAdd(c->tags_used, aux, 3, hd, NULL);

	i32 = (aux[0]<<16) | (aux[1]<<8) | aux[2];
#ifndef TN_external
	if (s->nTN >= s->aTN) {
	    s->aTN = s->aTN ? s->aTN*2 : 1024;
	    s->TN = realloc(s->TN, s->aTN * sizeof(*s->TN));
	}
	s->TN[s->nTN++] = i32;
	cram_stats_add(c->TN_stats, i32);
#else
	tmp += itf8_put(tmp, i32);
#endif

	switch(aux[2]) {
	case 'A': case 'C': case 'c':
	    aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    *tmp++=*aux++;
	    break;

	case 'S': case 's':
	    aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    *tmp++=*aux++; *tmp++=*aux++;
	    break;

	case 'I': case 'i': case 'f':
	    aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    break;

	case 'd':
	    aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    break;

	case 'Z': case 'H':
	    aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    while ((*tmp++=*aux++));
	    *tmp++ = '\t'; // stop byte
	    break;

	case 'B': {
	    int type = aux[3], blen;
	    uint32_t count = (uint32_t)((((unsigned char *)aux)[4]<< 0) +
					(((unsigned char *)aux)[5]<< 8) +
					(((unsigned char *)aux)[6]<<16) +
					(((unsigned char *)aux)[7]<<24));
	    // skip TN field
	    aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;

	    // We use BYTE_ARRAY_LEN with external length, so store that first
	    switch (type) {
	    case 'c': case 'C':
		blen = count;
		break;
	    case 's': case 'S':
		blen = 2*count;
		break;
	    case 'i': case 'I': case 'f':
		blen = 4*count;
		break;
	    default:
		fprintf(stderr, "Unknown sub-type '%c' for aux type 'B'\n",
			type);
		return NULL;
		    
	    }

	    tmp += itf8_put(tmp, blen+5);

	    *tmp++=*aux++; // sub-type & length
	    *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;

	    // The tag data itself
	    memcpy(tmp, aux, blen); tmp += blen; aux += blen;

	    //cram_stats_add(c->aux_B_stats, blen);
	    break;
	}
	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", aux[2]);
	    return NULL;
	}
    }
    cram_stats_add(c->TC_stats, cr->ntags);

#ifdef DS_SEQ
    cr->aux = DSTRING_LEN(s->aux_ds);
    cr->aux_size = tmp - (DSTRING_STR(s->aux_ds) + cr->aux);
    DSTRING_LEN(s->aux_ds) = tmp - DSTRING_STR(s->aux_ds);
    assert(s->aux_ds->length <= s->aux_ds->allocated);
#else
    cr->aux = BLOCK_SIZE(s->aux_blk);
    cr->aux_size = (uc *)tmp - (BLOCK_DATA(s->aux_blk) + cr->aux);
    BLOCK_SIZE(s->aux_blk) = (uc *)tmp - BLOCK_DATA(s->aux_blk);
    assert(s->aux_blk->byte <= s->aux_blk->alloc);
#endif

    return rg;
}

/*
 * Handles creation of a new container, flushing any existing slices as
 * appropriate.
 *
 * Really this is next slice, which may or may not lead to a new container.
 *
 * Returns cram_container pointer on success
 *         NULL on failure.
 */
static cram_container *cram_next_container(cram_fd *fd, bam_seq_t *b) {
    cram_container *c = fd->ctr;
    cram_slice *s;
    int i;

    /* First occurence */
    if (c->curr_ref == -2)
	c->curr_ref = b->ref;

    if (c->slice) {
	s = c->slice;
	s->hdr->ref_seq_id    = c->curr_ref;
	s->hdr->ref_seq_start = fd->first_base;
	s->hdr->ref_seq_span  = fd->last_base - fd->first_base + 1;
	s->hdr->num_records   = c->curr_rec;

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

	/* Encode slices */
	if (-1 == cram_encode_container(fd, c))
	    return NULL;

	// Move to sep func, as we need cram_flush_container for
	// the closing phase to flush the partial container.
	for (i = 0; i < c->max_slice; i++) {
	    cram_free_slice(c->slices[i]);
	    c->slices[i] = NULL;
	}

	c->slice = NULL;
	c->curr_slice = 0;

	/* Easy approach for purposes of freeing stats */
	cram_free_container(c);
	c = fd->ctr = cram_new_container(SEQS_PER_SLICE, SLICE_PER_CNT);
	c->curr_ref = b->ref;
    }

    c->last_pos = fd->first_base = fd->last_base = b->pos+1;

    /* New slice */
    c->slice = c->slices[c->curr_slice] =
	cram_new_slice(MAPPED_SLICE, c->max_rec);
    if (!c->slice)
	return NULL;

    c->slice->hdr->ref_seq_id = b->ref;
    c->slice->hdr->ref_seq_start = b->pos+1;
    c->slice->last_apos = b->pos+1;

    c->curr_rec = 0;

    return c;
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
    char *ref, *seq, *qual;

    if (!fd->ctr)
	fd->ctr = cram_new_container(SEQS_PER_SLICE, SLICE_PER_CNT);

    c = fd->ctr;

    if (!c->slice || c->curr_rec == c->max_rec ||
	(b->ref != c->curr_ref && c->curr_ref >= -1)) {

	if (NULL == (c = cram_next_container(fd, b)))
	    return -1;
	
	cram_get_ref(fd, b->ref, 1, 0);
    }

    ref = fd->ref;

    // Create a cram_record
    s = c->slice;
    cr = &s->crecs[c->curr_rec++]; // cache c->slice->crecs as c->rec?


    //fprintf(stderr, "%s => %d\n", rg ? rg : "\"\"", cr->rg);

    // Fields to resolve later
    //cr->mate_line;    // index to another cram_record
    //cr->mate_flags;   // MF
    //cr->ntags;        // TC
    cr->ntags      = 0; //cram_stats_add(c->TC_stats, cr->ntags);
    rg = cram_encode_aux(fd, b, c, s, cr);

    //cr->aux_size = b->blk_size - ((char *)bam_aux(b) - (char *)&b->ref);
    //cr->aux = DSTRING_LEN(s->aux_ds);
    //dstring_nappend(s->aux_ds, bam_aux(b), cr->aux_size);

    /* Read group, identified earlier */
    if (rg) {
	HashItem *hi = HashTableSearch(fd->SAM_hdr->rg_hash, rg, 0);
	tag_list_t *tl = hi->data.p;
	cr->rg = hi ? tl->key : -1;
    } else {
	HashItem *hi = HashTableSearch(fd->SAM_hdr->rg_hash, "UNKNOWN", 0);
	assert(hi);
	cr->rg = ((tag_list_t *)hi->data.p)->key;
    }
    cram_stats_add(c->RG_stats, cr->rg);

    
    cr->ref_id      = b->ref;
    cr->flags       = bam_flag(b);
    if (bam_cigar_len(b) == 0)
	cr->flags |= BAM_FUNMAP;
    cram_stats_add(c->BF_stats, cram_flag_swap[cr->flags & 0x7ff]);

    cr->cram_flags  = CRAM_FLAG_PRESERVE_QUAL_SCORES; // FIXME
    //cram_stats_add(c->CF_stats, cr->cram_flags);

    cr->len         = bam_seq_len(b);  cram_stats_add(c->RL_stats, cr->len);
    cr->apos        = b->pos+1;
    cram_stats_add(c->AP_stats, cr->apos - s->last_apos);
    s->last_apos = cr->apos;

#ifdef DS_SEQ
    cr->name        = DSTRING_LEN(s->name_ds);
    cr->name_len    = bam_name_len(b); cram_stats_add(c->RN_stats, cr->name_len);
    dstring_nappend(s->name_ds, bam_name(b), bam_name_len(b));
#else
    cr->name        = BLOCK_SIZE(s->name_blk);
    cr->name_len    = bam_name_len(b); cram_stats_add(c->RN_stats, cr->name_len);
    BLOCK_APPEND(s->name_blk, bam_name(b), bam_name_len(b));
#endif


    /*
     * This seqs_ds is largely pointless and it could reuse the same memory
     * over and over.
     * s->base_ds is what we need for encoding.
     */
#ifdef DS_SEQ
    cr->seq         = DSTRING_LEN(s->seqs_ds);
    cr->qual        = DSTRING_LEN(s->qual_ds);
    DSTRING_RESIZE(s->seqs_ds, DSTRING_LEN(s->seqs_ds) + cr->len);
    DSTRING_RESIZE(s->qual_ds, DSTRING_LEN(s->qual_ds) + cr->len);
    seq = cp = DSTRING_STR(s->seqs_ds) + cr->seq;

    for (i = 0; i < cr->len; i++) {
	// FIXME: do 2 char at a time for efficiency
	cp[i] = bam_nt16_rev_table[bam_seqi(bam_seq(b), i)];
    }
    DSTRING_LEN(s->seqs_ds) += cr->len;

    qual = cp = bam_qual(b);
#else
    cr->seq         = BLOCK_SIZE(s->seqs_blk);
    cr->qual        = BLOCK_SIZE(s->qual_blk);
    BLOCK_GROW(s->seqs_blk, cr->len);
    BLOCK_GROW(s->qual_blk, cr->len);
    seq = cp = BLOCK_END(s->seqs_blk);

    for (i = 0; i < cr->len; i++) {
	// FIXME: do 2 char at a time for efficiency
	cp[i] = bam_nt16_rev_table[bam_seqi(bam_seq(b), i)];
    }
    BLOCK_SIZE(s->seqs_blk) += cr->len;

    qual = cp = bam_qual(b);
#endif

    /* Copy and parse */
    if (!(cr->flags & BAM_FUNMAP)) {
	int32_t *cig_to, *cig_from;
	int apos = cr->apos-1, spos = 0;

	cr->cigar       = s->ncigar;
	cr->ncigar      = bam_cigar_len(b);
	if (cr->cigar + cr->ncigar >= s->cigar_alloc) {
	    s->cigar_alloc = s->cigar_alloc ? s->cigar_alloc*2 : 1024;
	    s->cigar = realloc(s->cigar, s->cigar_alloc * sizeof(*s->cigar));
	}

	cig_to = (int32_t *)s->cigar;
	cig_from = (int32_t *)bam_cigar(b);

	cr->feature = 0;
	cr->nfeature = 0;
	for (i = 0; i < cr->ncigar; i++) {
	    enum cigar_op cig_op = cig_from[i] & BAM_CIGAR_MASK;
	    int cig_len = cig_from[i] >> BAM_CIGAR_SHIFT;
	    cig_to[i] = cig_from[i];

	    /* Can also generate events from here for CRAM diffs */

	    switch (cig_op) {
		int l;

	    // Don't trust = and X ops to be correct.
	    case BAM_CMATCH:
	    case BAM_CBASE_MATCH:
	    case BAM_CBASE_MISMATCH:
		//fprintf(stderr, "\nBAM_CMATCH\nR: %.*s\nS: %.*s\n",
		//	cig_len, &ref[apos], cig_len, &seq[spos]);
		for (l = 0; l < cig_len; l++, apos++, spos++) {
		    if (ref[apos] != seq[spos]) {
			//fprintf(stderr, "Subst: %d; %c vs %c\n",
			//	spos, ref[apos], seq[spos]);
			cram_add_substitution(c, s, cr, spos,
					      seq[spos], qual[spos],
					      ref[apos]);
		    }
		}
		break;
#if 0
	    case BAM_CBASE_MATCH:
		//fprintf(stderr, "\nBAM_CBASE_MATCH\nR: %.*s\nS: %.*s\n",
		//	cig_len, &ref[apos], cig_len, &seq[spos]);
		apos += cig_len;
		spos += cig_len;
		break;

	    case BAM_CBASE_MISMATCH:
		//fprintf(stderr, "\nBAM_CBASE_MISMATCH\nR: %.*s\nS: %.*s\n",
		//	cig_len, &ref[apos], cig_len, &seq[spos]);
		for (l = 0; l < cig_len; l++, apos++, spos++) {
		    cram_add_substitution(c, s, cr, spos,
					  seq[spos], qual[spos], ref[apos]);
		}
		break;
#endif
		
	    case BAM_CDEL:
		cram_add_deletion(c, s, cr, spos, cig_len, &seq[spos]);
		apos += cig_len;
		break;

	    case BAM_CREF_SKIP:
		fprintf(stderr, "BAM_CREF_SKIP unimplemented\n");
		apos += cig_len;
		break;

	    case BAM_CINS:
		cram_add_insertion(c, s, cr, spos, cig_len, &seq[spos]);
		spos += cig_len;
		break;

	    case BAM_CSOFT_CLIP:
		cram_add_softclip(c, s, cr, spos, cig_len, &seq[spos]);
		spos += cig_len;
		break;

	    case BAM_CHARD_CLIP:
		fprintf(stderr, "BAM_HARD_CLIP unimplemented\n");
		break;
	
	    case BAM_CPAD:
		fprintf(stderr, "BAM_HARD_CLIP unimplemented\n");
		break;
	    }
	}
	cr->aend = apos;
	cram_stats_add(c->FN_stats, cr->nfeature);
    } else {
	// Unmapped
	cr->cigar  = 0;
	cr->ncigar = 0;
	cr->nfeature = 0;
	cr->aend = cr->apos;
#ifdef BA_external
	s->BA_len += cr->len;
#else
	for (i = 0; i < cr->len; i++)
	    cram_stats_add(c->BA_stats, seq[i]);
#endif
    }

    /*
     * Append to the qual block now. We do this here as
     * cram_add_substitution() can generate BA/QS events which need to 
     * be in the qual block before we append the rest of the data.
     */
#ifdef DS_SEQ
    DSTRING_RESIZE(s->qual_ds, DSTRING_LEN(s->qual_ds) + cr->len);
    qual = cp = DSTRING_STR(s->qual_ds) + cr->qual;
    for (i = 0; i < cr->len; i++) {
	cp[i] = bam_qual(b)[i];
    }
    DSTRING_LEN(s->qual_ds) += cr->len;
#else
    BLOCK_GROW(s->qual_blk, cr->len);
    qual = cp = BLOCK_END(s->qual_blk);
    for (i = 0; i < cr->len; i++) {
	cp[i] = bam_qual(b)[i];
    }
    BLOCK_SIZE(s->qual_blk) += cr->len;
#endif    

    /* Now we know apos and aend both, update mate-pair information */
    {
	int new;
	HashData hd;
	HashItem *hi;

	hd.i = c->curr_rec-1;
	//fprintf(stderr, "Checking %"PRId64"/%.*s\t", hd.i,
	//	cr->name_len, DSTRING_STR(s->name_ds)+cr->name);
	if (cr->flags & BAM_FPAIRED) {
#ifdef DS_SEQ
	    hi = HashTableAdd(s->pair, DSTRING_STR(s->name_ds)+cr->name,
			      cr->name_len, hd, &new);
#else
	    hi = HashTableAdd(s->pair, BLOCK_DATA(s->name_blk)+cr->name,
			      cr->name_len, hd, &new);
	} else {
	    new = 1;
	}
#endif
	if (!new) {
	    cram_record *p = &s->crecs[hi->data.i];
	    
	    //fprintf(stderr, "paired %"PRId64"\n", hi->data.i);

	    // copy from p to cr
	    cr->mate_pos = p->apos;
	    cram_stats_add(c->NP_stats, cr->mate_pos);

	    cr->tlen = cr->aend - p->apos;
	    cram_stats_add(c->TS_stats, cr->tlen);

	    cr->mate_flags =
		((p->flags & BAM_FMUNMAP)   == BAM_FMUNMAP)   * CRAM_M_UNMAP +
		((p->flags & BAM_FMREVERSE) == BAM_FMREVERSE) * CRAM_M_REVERSE;
	    cram_stats_add(c->MF_stats, cr->mate_flags);

	    // copy from cr to p
	    cram_stats_del(c->NP_stats, p->mate_pos);
	    p->mate_pos = cr->apos;
	    cram_stats_add(c->NP_stats, p->mate_pos);

	    cram_stats_del(c->MF_stats, p->mate_flags);
	    p->mate_flags =
		((cr->flags & BAM_FMUNMAP)   == BAM_FMUNMAP)  * CRAM_M_UNMAP +
		((cr->flags & BAM_FMREVERSE) == BAM_FMREVERSE)* CRAM_M_REVERSE;
	    cram_stats_add(c->MF_stats, p->mate_flags);

	    cram_stats_del(c->TS_stats, p->tlen);
	    p->tlen = p->apos - cr->aend;
	    cram_stats_add(c->TS_stats, p->tlen);

	    // Clear detached from cr flags
	    //cram_stats_del(c->CF_stats, cr->cram_flags);
	    cr->cram_flags &= ~CRAM_FLAG_DETACHED;
	    cram_stats_add(c->CF_stats, cr->cram_flags);

	    // Clear detached from p flags and set downstream
	    cram_stats_del(c->CF_stats, p->cram_flags);
	    p->cram_flags  &= ~CRAM_FLAG_DETACHED;
	    p->cram_flags  |=  CRAM_FLAG_MATE_DOWNSTREAM;
	    cram_stats_add(c->CF_stats, p->cram_flags);

	    p->mate_line = hd.i - (hi->data.i + 1);
	    cram_stats_add(c->NF_stats, p->mate_line);

	    hi->data.i = c->curr_rec-1;
	    //HashTableDel(s->pair, hi, 0);
	} else {
	    //fprintf(stderr, "unpaired\n");

	    /* Derive mate flags from this flag */
	    cr->mate_flags = 0;
	    if (bam_flag(b) & BAM_FMUNMAP)
		cr->mate_flags |= CRAM_M_UNMAP;
	    if (bam_flag(b) & BAM_FMREVERSE)
		cr->mate_flags |= CRAM_M_REVERSE;

	    cram_stats_add(c->MF_stats, cr->mate_flags);

	    cr->mate_pos    = MAX(b->mate_pos+1, 0);
	    cram_stats_add(c->NP_stats, cr->mate_pos);

	    cr->tlen        = b->ins_size;
	    cram_stats_add(c->TS_stats, cr->tlen);

	    cr->cram_flags |= CRAM_FLAG_DETACHED;
	    cram_stats_add(c->CF_stats, cr->cram_flags);
	}
    }


    cr->mqual       = bam_map_qual(b); cram_stats_add(c->MQ_stats, cr->mqual);
    cr->mate_ref_id = b->mate_ref;     cram_stats_add(c->NS_stats, b->mate_ref);

    c->curr_ctr_rec++;

    // FIXME - use cigar
    if (fd->last_base < cr->aend)
	fd->last_base = cr->aend;

    return 0;
}


/*
 * Read the next cram record and return it.
 * Note that to decode cram_record the caller will need to look up some data
 * in the current slice, pointed to by fd->ctr->slice. This is valid until
 * the next call to cram_get_seq (which may invalidate it).
 *
 * Returns record pointer on success (do not free)
 *        NULL on failure
 */
cram_record *cram_get_seq(cram_fd *fd) {
    cram_container *c;
    cram_slice *s;

    if (!(c = fd->ctr)) {
	// Load first container. Do as part of cram_open?
	if (!(c = fd->ctr = cram_read_container(fd)))
	    return NULL;
    }

    s = c->slice;
    if (!s || c->curr_rec == c->max_rec) {
	int id;

	// new slice
	if (s)
	    cram_free_slice(s);

	if (c->curr_slice == c->max_slice) {
	    // new container
	    cram_free_container(c);

	    if (!(c = fd->ctr = cram_read_container(fd)))
		return NULL;
	}

	s = c->slice = cram_read_slice(fd);
	c->curr_slice++;
	c->curr_rec = 0;
	c->max_rec = s->hdr->num_records;

	// FIXME: this should be correct, but we have a bug in the
	// Java CRAM implementation?
	s->last_apos = s->hdr->ref_seq_start;
	    
	for (id = 0; id < s->hdr->num_blocks; id++)
	    cram_uncompress_block(s->block[id]);

	/* Test decoding of 1st seq */
	if (cram_decode_slice(fd, c, s, fd->SAM_hdr) != 0) {
	    fprintf(stderr, "Failure to decode slice\n");
	    return NULL;
	}
    }

    return &s->crecs[c->curr_rec++];
}

/*
 * Read the next cram record and convert it to a bam_seq_t struct.
 *
 * Returns 0 on success
 *        -1 on EOF or failure (check fd->err)
 */
int cram_get_bam_seq(cram_fd *fd, bam_seq_t **bam, size_t *bam_alloc) {
    cram_record *cr;
    cram_container *c;
    cram_slice *s;

    if (!(cr = cram_get_seq(fd)))
	return -1;

    c = fd->ctr;
    s = c->slice;
    cram_to_bam(fd->SAM_hdr, fd, s, cr, c->curr_rec-1, bam, bam_alloc);

    return 0;
}


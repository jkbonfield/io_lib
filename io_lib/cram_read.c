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
static unsigned int bam_flag_swap[0x200];
static unsigned char L[256];

/*
 * Initialise lookup tables. Used within cram_decode_seq().
 */
void cram_reader_init(void) {
    int i;

    if (cram_init_done)
	return;

    memset(L, 4, 256);
    L['A'] = 0; L['a'] = 0;
    L['C'] = 1; L['c'] = 1;
    L['G'] = 2; L['g'] = 2;
    L['T'] = 3; L['t'] = 3;

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

    cram_init_done = 1;
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

    cp = (char *)b->data;
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
    char *cp = (char *)b->data;
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

    s->seqs_blk = cram_new_block(EXTERNAL, 0);
    s->qual_blk = cram_new_block(EXTERNAL, 1);
    s->name_blk = cram_new_block(EXTERNAL, 2);
    s->aux_blk  = cram_new_block(EXTERNAL, 4);
    s->base_blk = cram_new_block(EXTERNAL, 0);
#ifdef TN_external
    s->tn_blk   = cram_new_block(EXTERNAL, 6);
#endif


    s->crecs = NULL;

    s->last_apos = s->hdr->ref_seq_start;
    
    s->id = fd->slice_num++;

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
#ifdef TN_external
    if (s->tn_blk)
	cram_free_block(s->tn_blk);
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


#if 0
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
#endif

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
enum cram_encoding cram_stats_encoding(cram_fd *fd, cram_stats *st) {
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

    if (fd->verbose > 1)
	fprintf(stderr, "Range = %d..%d, nvals=%d, ntot=%d\n",
		min_val, max_val, nvals, ntot);

    /* Theoretical entropy */
    {
	double dbits = 0;
	for (i = 0; i < nvals; i++) {
	    dbits += freqs[i] * log((double)freqs[i]/ntot);
	}
	dbits /= -log(2);
	if (fd->verbose > 1)
	    fprintf(stderr, "Entropy = %f\n", dbits);
    }

#if 0
    /* Unary */
    if (min_val >= 0) {
	for (bits = i = 0; i < nvals; i++)
	    bits += freqs[i]*(vals[i]+1);
	if (fd->verbose > 1)
	    fprintf(stderr, "UNARY   = %d\n", bits);
	if (best_size > bits)
	    best_size = bits, best_encoding = E_NULL; //E_UNARY;
    }

    /* Beta */
    bits = nbits(max_val - min_val) * ntot;
    if (fd->verbose > 1)
	fprintf(stderr, "BETA    = %d\n", bits);
    if (best_size > bits)
	best_size = bits, best_encoding = E_BETA;

    /* Gamma */
    for (bits = i = 0; i < nvals; i++)
	bits += ((nbits(vals[i]-min_val+1)-1) + nbits(vals[i]-min_val+1)) * freqs[i];
    if (fd->verbose > 1)
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

	if (fd->verbose > 1)
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
    if (fd->verbose > 1)
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
 * Loads a reference.
 * FIXME: use the sam_comp.cpp get_ref_base() equivalent to allow
 * random access within an indexed reference instead?
 *
 * Returns a ref_seq structure on success
 *         NULL on failure
 */
static refs *load_reference(char *fn) {
    struct stat sb;
    FILE *fp;
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

    if (decode_md) {
	orig_aux = BLOCK_SIZE(s->aux_blk);
	BLOCK_APPEND(s->aux_blk, "MDZ", 3);
    }

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
		    BLOCK_APPENDF_2(s->aux_blk, buf, "%d%c",
				    md_dist, s->ref[ref_pos-s->ref_start +1]);
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
		BLOCK_APPENDF_1(s->aux_blk, buf, "%d^", md_dist);
		BLOCK_APPEND(s->aux_blk, &s->ref[ref_pos - s->ref_start +1],
			     i32);
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
	BLOCK_APPENDF_1(s->aux_blk, buf, "%d", md_dist);
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
	BLOCK_APPEND_CHAR(s->aux_blk, '\0'); // null terminate MD:Z:
	cr->aux_size += BLOCK_SIZE(s->aux_blk) - orig_aux;
	buf[0] = 'N'; buf[1] = 'M'; buf[2] = 'I';
	buf[3] = (nm>> 0) & 0xff;
	buf[4] = (nm>> 8) & 0xff;
	buf[5] = (nm>>16) & 0xff;
	buf[6] = (nm>>24) & 0xff;
	BLOCK_APPEND(s->aux_blk, buf, 7);
	cr->aux_size += 7;
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
    cr->aux = BLOCK_SIZE(s->aux_blk);

    for (i = 0; i < cr->ntags; i++) {
	int32_t id, out_sz = 1;
	unsigned char tag_data[1024];
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

	BLOCK_APPEND(s->aux_blk, (char *)tag_data, out_sz+3);
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

	    // Read directly into name cram_block
	    BLOCK_GROW(s->name_blk, MAX_NAME_LEN);
	    name = (char *)BLOCK_END(s->name_blk);
	    r |= c->comp_hdr->RN_codec->decode(s, c->comp_hdr->RN_codec, blk,
					       name, &out_sz2);

	    cr->name = BLOCK_SIZE(s->name_blk);
	    BLOCK_SIZE(s->name_blk) += out_sz2;
	    cr->name_len = out_sz2;
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
	    
		// Read directly into name cram_block
		BLOCK_GROW(s->name_blk, MAX_NAME_LEN);
		name = (char *)BLOCK_END(s->name_blk);
		r |= c->comp_hdr->RN_codec->decode(s, c->comp_hdr->RN_codec, blk,
						   name, &out_sz2);

		cr->name = BLOCK_SIZE(s->name_blk);
		BLOCK_SIZE(s->name_blk) += out_sz2;
		cr->name_len = out_sz2;
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
	cr->seq = BLOCK_SIZE(s->seqs_blk);
	BLOCK_GROW(s->seqs_blk, cr->len);
	seq = (char *)BLOCK_END(s->seqs_blk);
	BLOCK_SIZE(s->seqs_blk) += cr->len;

	cr->qual = BLOCK_SIZE(s->qual_blk);
	BLOCK_GROW(s->qual_blk, cr->len);
	qual = (char *)BLOCK_END(s->qual_blk);
	BLOCK_SIZE(s->qual_blk) += cr->len;

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
	cram_reader_init();

	if (!(fd->file_def = cram_read_file_def(fd)))
	    goto err;

	if (!(fd->SAM_hdr = cram_read_SAM_hdr(fd)))
	    goto err;

	//if (0 != parse_SAM_hdr(&fd->SAM_hdr))
	//    goto err;
    } else {
	/* Writer */
	cram_file_def def;

	cram_writer_init();

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
    fd->verbose = 0;

    for (i = 0; i < 7; i++)
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

    case CRAM_OPT_VERBOSITY:
	fd->verbose = val->i;
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

    for (i = 0; i < 7; i++)
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
	name = (char *)BLOCK_DATA(s->name_blk) + cr->name;
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
                                (char *)BLOCK_DATA(s->seqs_blk) + cr->seq,
				(char *)BLOCK_DATA(s->qual_blk) + cr->qual);
 
    old_idx = bam_idx;

   /* Auxiliary strings */
    if (cr->aux_size != 0) {
	memcpy(&((char *)*bam)[bam_idx], BLOCK_DATA(s->aux_blk) + cr->aux,
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


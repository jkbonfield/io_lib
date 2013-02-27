/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 *
 * - In-memory decoding of CRAM data structures.
 * - Iterator for reading CRAM record by record.
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


/* ----------------------------------------------------------------------
 * CRAM compression headers
 */

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

/* ----------------------------------------------------------------------
 * CRAM slices
 */

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

/* ----------------------------------------------------------------------
 * Primary CRAM sequence decoder
 */

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
		ref_base = fd->L1[(uc)s->ref[ref_pos - s->ref_start +1]];
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

/*
 * Decode an entire slice from container blocks. Fills out s->crecs[] array.
 * Returns 0 on success
 *        -1 on failure
 */
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
	cr->flags = fd->bam_flag_swap[bf];

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

/* ----------------------------------------------------------------------
 * CRAM sequence iterators.
 */

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
static int cram_to_bam(bam_file_t *bfd, cram_fd *fd, cram_slice *s,
		       cram_record *cr, int rec, bam_seq_t **bam,
		       size_t *bam_alloc) {
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


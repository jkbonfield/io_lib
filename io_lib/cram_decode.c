/*
 * Copyright (c) 2013, 2014, 2015 Genome Research Ltd.
 * Author(s): James Bonfield, Rob Davies
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 * 
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *    Institute nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific
 *    prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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
#include "io_lib/md5.h"
#include "io_lib/crc32.h"

//Whether CIGAR has just M or uses = and X to indicate match and mismatch
//#define USE_X


/* ----------------------------------------------------------------------
 * CRAM compression headers
 */

/*
 * Decodes the Tag Dictionary record in the preservation map
 * Updates the cram compression header.
 * 
 * Returns number of bytes decoded on success
 *        -1 on failure
 */
int cram_decode_TD(cram_fd *fd, char *cp, const char *endp,
		   cram_block_compression_hdr *h) {
    char *op = cp;
    unsigned char *dat;
    cram_block *b;
    int32_t blk_size = 0;
    int nTL, i, sz, err = 0;

    if (!(b = cram_new_block(0, 0)))
	return -1;

    /* Decode */
    blk_size = fd->vv.varint_get32(&cp, endp, &err);
    if (!blk_size) {
	h->nTL = 0;
	h->TL = NULL;
	cram_free_block(b);
	return cp - op;
    }

    if (err || blk_size < 0 || endp - cp < blk_size) {
        cram_free_block(b);
	return -1;
    }

    BLOCK_APPEND(b, cp, blk_size);
    cp += blk_size;
    sz = cp - op;

    // Force nul termination if missing
    if (BLOCK_DATA(b)[BLOCK_SIZE(b)-1])
	BLOCK_APPEND_CHAR(b, '\0');

    /* Set up TL lookup table */
    dat = BLOCK_DATA(b);

    // Count
    for (nTL = i = 0; i < BLOCK_SIZE(b); i++) {
	nTL++;
	while (dat[i])
	    i++;
    }

    // Copy
    h->nTL = nTL;
    if (!(h->TL = calloc(h->nTL, sizeof(unsigned char *)))) {
        cram_free_block(b);
        return -1;
    }
    for (nTL = i = 0; i < BLOCK_SIZE(b); i++) {
	h->TL[nTL++] = &dat[i];
	while (dat[i])
	    i++;
    }
    h->TD_blk = b;
    
    return sz;
}

/*
 * Decodes a CRAM block compression header.
 * Returns header ptr on success
 *         NULL on failure
 */
cram_block_compression_hdr *cram_decode_compression_header(cram_fd *fd,
							   cram_block *b) {
    char *cp, *endp, *cp_copy;
    cram_block_compression_hdr *hdr = calloc(1, sizeof(*hdr));
    int i, err = 0;
    int32_t map_size = 0, map_count = 0;

    if (!hdr)
	return NULL;

    if (b->method != RAW) {
	if (cram_uncompress_block(b)) {
	    free(hdr);
	    return NULL;
	}
    }

    cp = (char *)b->data;
    endp = cp + b->uncomp_size;

    if (IS_CRAM_1_VERS(fd)) {
	hdr->ref_seq_id = fd->vv.varint_get32(&cp, endp, &err);
	if (CRAM_MAJOR_VERS(fd->version) >= 4) {
	    hdr->ref_seq_start = fd->vv.varint_get64(&cp, endp, &err);
	    hdr->ref_seq_span  = fd->vv.varint_get64(&cp, endp, &err);
	} else {
	    hdr->ref_seq_start = fd->vv.varint_get32(&cp, endp, &err);
	    hdr->ref_seq_span  = fd->vv.varint_get32(&cp, endp, &err);
	}
	hdr->num_records   = fd->vv.varint_get32(&cp, endp, &err);
	hdr->num_landmarks = fd->vv.varint_get32(&cp, endp, &err);
	if (!(hdr->landmark = malloc(hdr->num_landmarks * sizeof(int32_t)))) {
	    free(hdr);
	    return NULL;
	}
	for (i = 0; i < hdr->num_landmarks; i++)
	    hdr->landmark[i] = fd->vv.varint_get32(&cp, endp, &err);;
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
    hdr->read_names_included = 0;
    hdr->AP_delta = 1;
    hdr->qs_seq_orient = 1;
    memcpy(hdr->substitution_matrix, "CGTNAGTNACTNACGNACGT", 20);

    /* Preservation map */
    map_size  = fd->vv.varint_get32(&cp, endp, &err); cp_copy = cp;
    map_count = fd->vv.varint_get32(&cp, endp, &err);
    for (i = 0; i < map_count; i++) {
	HashData hd;
	if (endp - cp < 2) {
	    cram_free_compression_header(hdr);
	    return NULL;
	}
	cp += 2;
	switch(CRAM_KEY(cp[-2],cp[-1])) {
	case CRAM_KEY('M','I'): // was mapped QS included in V1.0
	case CRAM_KEY('U','I'): // was unmapped QS included in V1.0
	case CRAM_KEY('P','I'): // was unmapped placed in V1.0
	    hd.i = *cp++;
	    break;

	case CRAM_KEY('R','N'):
	    hd.i = *cp++;
	    if (!HashTableAdd(hdr->preservation_map, "RN", 2, hd, NULL)) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	    hdr->read_names_included = hd.i;
	    break;

	case CRAM_KEY('A','P'):
	    hd.i = *cp++;
	    if (!HashTableAdd(hdr->preservation_map, "AP", 2, hd, NULL)) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	    hdr->AP_delta = hd.i;
	    break;

	case CRAM_KEY('R','R'):
	    hd.i = *cp++;
	    if (!HashTableAdd(hdr->preservation_map, "RR", 2, hd, NULL)) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	    hdr->no_ref = !hd.i;
	    break;

	case CRAM_KEY('Q','O'):
	    hd.i = *cp++;
	    if (!HashTableAdd(hdr->preservation_map, "QO", 2, hd, NULL)) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	    hdr->qs_seq_orient = hd.i;
	    break;

	case CRAM_KEY('S','M'):
	    if (endp - cp < 5) {
	        cram_free_compression_header(hdr);
		return NULL;
	    }
	    hdr->substitution_matrix[0][(cp[0]>>6)&3] = 'C';
	    hdr->substitution_matrix[0][(cp[0]>>4)&3] = 'G';
	    hdr->substitution_matrix[0][(cp[0]>>2)&3] = 'T';
	    hdr->substitution_matrix[0][(cp[0]>>0)&3] = 'N';

	    hdr->substitution_matrix[1][(cp[1]>>6)&3] = 'A';
	    hdr->substitution_matrix[1][(cp[1]>>4)&3] = 'G';
	    hdr->substitution_matrix[1][(cp[1]>>2)&3] = 'T';
	    hdr->substitution_matrix[1][(cp[1]>>0)&3] = 'N';

	    hdr->substitution_matrix[2][(cp[2]>>6)&3] = 'A';
	    hdr->substitution_matrix[2][(cp[2]>>4)&3] = 'C';
	    hdr->substitution_matrix[2][(cp[2]>>2)&3] = 'T';
	    hdr->substitution_matrix[2][(cp[2]>>0)&3] = 'N';

	    hdr->substitution_matrix[3][(cp[3]>>6)&3] = 'A';
	    hdr->substitution_matrix[3][(cp[3]>>4)&3] = 'C';
	    hdr->substitution_matrix[3][(cp[3]>>2)&3] = 'G';
	    hdr->substitution_matrix[3][(cp[3]>>0)&3] = 'N';

	    hdr->substitution_matrix[4][(cp[4]>>6)&3] = 'A';
	    hdr->substitution_matrix[4][(cp[4]>>4)&3] = 'C';
	    hdr->substitution_matrix[4][(cp[4]>>2)&3] = 'G';
	    hdr->substitution_matrix[4][(cp[4]>>0)&3] = 'T';
	    hd.p = cp;
	    cp += 5;
	    if (!HashTableAdd(hdr->preservation_map, "SM", 2, hd, NULL)) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	    break;

	case CRAM_KEY('T','D'): {
	    int sz = cram_decode_TD(fd, cp, endp, hdr); // tag dictionary
	    if (sz < 0) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	    hd.p = cp;
	    cp += sz;
	    if (!HashTableAdd(hdr->preservation_map, "TD", 2, hd, NULL)) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	    break;
	}

	default:
	    fprintf(stderr, "Unrecognised preservation map key %c%c\n",
		    cp[-2], cp[-1]);
	    // guess byte;
	    cp++;
	    break;
	}
    }
    if (cp - cp_copy != map_size) {
	cram_free_compression_header(hdr);
	return NULL;
    }

    /* Record encoding map */
    map_size  = fd->vv.varint_get32(&cp, endp, &err); cp_copy = cp;
    map_count = fd->vv.varint_get32(&cp, endp, &err);
    int is_v4 = CRAM_MAJOR_VERS(fd->version) >= 4 ? 1 : 0;
    for (i = 0; i < map_count; i++) {
	char *key = cp;
	int32_t encoding = E_NULL;
	int32_t size = 0;
	cram_map *m = malloc(sizeof(*m)); // FIXME: use pooled_alloc

	if (!m || endp - cp < 4) {
	    free(m);
	    cram_free_compression_header(hdr);
	    return NULL;
	}

	cp += 2;
	encoding = fd->vv.varint_get32(&cp, endp, &err);
	size     = fd->vv.varint_get32(&cp, endp, &err);

	// Fill out cram_map purely for cram_dump to dump out.
	m->key = (key[0]<<8)|key[1];
	m->encoding = encoding;
	m->size     = size;
	m->offset   = cp - (char *)b->data;
	m->codec = NULL;

	if (m->encoding == E_NULL)
	    continue;

	if (size < 0 || endp - cp < size) {
	    free(m);
	    cram_free_compression_header(hdr);
	    return NULL;
	}

	//printf("%s codes for %.2s\n", cram_encoding2str(encoding), key);

	/*
	 * For CRAM1.0 CF and BF are Byte and not Int.
	 * Practically speaking it makes no difference unless we have a
	 * 1.0 format file that stores these in EXTERNAL as only then
	 * does Byte vs Int matter.
	 *
	 * Neither this C code nor Java reference implementations did this,
	 * so we gloss over it and treat them as int.
	 */

	if (key[0] == 'B' && key[1] == 'F') {
	    if (!(hdr->codecs[DS_BF] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'C' && key[1] == 'F') {
	    if (!(hdr->codecs[DS_CF] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'I') {
	    if (!(hdr->codecs[DS_RI] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'L') {
	    if (!(hdr->codecs[DS_RL] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'A' && key[1] == 'P') {
	    if (!(hdr->codecs[DS_AP] = cram_decoder_init(hdr, encoding, cp, size,
							 is_v4 ? E_LONG : E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'G') {
	    if (!(hdr->codecs[DS_RG] = cram_decoder_init(hdr, encoding, cp, size,
							 E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'M' && key[1] == 'F') {
	    if (!(hdr->codecs[DS_MF] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'N' && key[1] == 'S') {
	    if (!(hdr->codecs[DS_NS] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'N' && key[1] == 'P') {
	    if (!(hdr->codecs[DS_NP] = cram_decoder_init(hdr, encoding, cp, size,
							 is_v4 ? E_LONG : E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'S') {
	    if (!(hdr->codecs[DS_TS] = cram_decoder_init(hdr, encoding, cp, size,
							 is_v4 ? E_LONG : E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'N' && key[1] == 'F') {
	    if (!(hdr->codecs[DS_NF] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'C') {
	    if (!(hdr->codecs[DS_TC] = cram_decoder_init(hdr, encoding, cp, size, E_BYTE,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'N') {
	    if (!(hdr->codecs[DS_TN] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'F' && key[1] == 'N') {
	    if (!(hdr->codecs[DS_FN] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'F' && key[1] == 'C') {
	    if (!(hdr->codecs[DS_FC] = cram_decoder_init(hdr, encoding, cp, size, E_BYTE,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'F' && key[1] == 'P') {
	    if (!(hdr->codecs[DS_FP] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'B' && key[1] == 'S') {
	    if (!(hdr->codecs[DS_BS] = cram_decoder_init(hdr, encoding, cp, size, E_BYTE,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'I' && key[1] == 'N') {
	    if (!(hdr->codecs[DS_IN] = cram_decoder_init(hdr, encoding, cp, size,
							 E_BYTE_ARRAY,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'S' && key[1] == 'C') {
	    if (!(hdr->codecs[DS_SC] = cram_decoder_init(hdr, encoding, cp, size,
							 E_BYTE_ARRAY,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'D' && key[1] == 'L') {
	    if (!(hdr->codecs[DS_DL] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'B' && key[1] == 'A') {
	    if (!(hdr->codecs[DS_BA] = cram_decoder_init(hdr, encoding, cp, size, E_BYTE,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'B' && key[1] == 'B') {
	    if (!(hdr->codecs[DS_BB] = cram_decoder_init(hdr, encoding, cp, size, E_BYTE_ARRAY,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'S') {
	    if (!(hdr->codecs[DS_RS] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'P' && key[1] == 'D') {
	    if (!(hdr->codecs[DS_PD] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'H' && key[1] == 'C') {
	    if (!(hdr->codecs[DS_HC] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'M' && key[1] == 'Q') {
	    if (!(hdr->codecs[DS_MQ] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'N') {
	    if (!(hdr->codecs[DS_RN] = cram_decoder_init(hdr, encoding, cp, size,
							 E_BYTE_ARRAY_BLOCK,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'Q' && key[1] == 'S') {
	    if (!(hdr->codecs[DS_QS] = cram_decoder_init(hdr, encoding, cp, size, E_BYTE,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'Q' && key[1] == 'Q') {
	    if (!(hdr->codecs[DS_QQ] = cram_decoder_init(hdr, encoding, cp, size, E_BYTE_ARRAY,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'L') {
	    if (!(hdr->codecs[DS_TL] = cram_decoder_init(hdr, encoding, cp, size, E_INT,
							 fd->version, &fd->vv))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'M') {
	} else if (key[0] == 'T' && key[1] == 'V') {
	} else
	    fprintf(stderr, "Unrecognised key: %.2s\n", key);

	cp += size;

	m->next = hdr->rec_encoding_map[CRAM_MAP(key[0], key[1])];
	hdr->rec_encoding_map[CRAM_MAP(key[0], key[1])] = m;
    }
    if (cp - cp_copy != map_size) {
	cram_free_compression_header(hdr);
	return NULL;
    }

    /* Tag encoding map */
    map_size  = fd->vv.varint_get32(&cp, endp, &err); cp_copy = cp;
    map_count = fd->vv.varint_get32(&cp, endp, &err);
    for (i = 0; i < map_count; i++) {
	int32_t encoding = E_NULL;
	int32_t size = 0;
	cram_map *m = malloc(sizeof(*m)); // FIXME: use pooled_alloc
	uint8_t key[3];

	if (!m || endp - cp < 6) {
	    free(m);
	    cram_free_compression_header(hdr);
	    return NULL;
	}

	m->key = fd->vv.varint_get32(&cp, endp, &err);
	key[0] = m->key>>16;
	key[1] = m->key>>8;
	key[2] = m->key;
	encoding = fd->vv.varint_get32(&cp, endp, &err);
	size     = fd->vv.varint_get32(&cp, endp, &err);

	m->encoding = encoding;
	m->size     = size;
	m->offset   = cp - (char *)b->data;
	if (size < 0 || endp - cp < size ||
	    !(m->codec = cram_decoder_init(hdr, encoding, cp, size,
					   E_BYTE_ARRAY_BLOCK, fd->version,
					   &fd->vv))) {
	    cram_free_compression_header(hdr);
	    free(m);
	    return NULL;
	}
	
	cp += size;

	m->next = hdr->tag_encoding_map[CRAM_MAP(key[0],key[1])];
	hdr->tag_encoding_map[CRAM_MAP(key[0],key[1])] = m;
    }
    if (err || cp - cp_copy != map_size) {
	cram_free_compression_header(hdr);
	return NULL;
    }

    return hdr;
}


/*
 * Note we also need to scan through the record encoding map to
 * see which data series share the same block, either external or
 * CORE. For example if we need the BF data series but MQ and CF
 * are also encoded in the same block then we need to add those in
 * as a dependency in order to correctly decode BF.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_dependent_data_series(cram_fd *fd,
			       cram_block_compression_hdr *hdr,
			       cram_slice *s) {
    int *block_used;
    int core_used = 0;
    int i;
    static int i_to_id[] = {
	DS_BF, DS_AP, DS_FP, DS_RL, DS_DL, DS_NF, DS_BA, DS_QS,
	DS_FC, DS_FN, DS_BS, DS_IN, DS_RG, DS_MQ, DS_TL, DS_RN,
	DS_NS, DS_NP, DS_TS, DS_MF, DS_CF, DS_RI, DS_RS, DS_PD,
	DS_HC, DS_SC, DS_BB, DS_QQ,
    };
    uint32_t orig_ds;

    /*
     * Set the data_series bit field based on fd->required_fields
     * contents.
     */
    if (fd->required_fields && fd->required_fields != INT_MAX) {
	s->data_series = 0;

	if (fd->required_fields & SAM_QNAME)
	    s->data_series |= CRAM_RN;

	if (fd->required_fields & SAM_FLAG)
	    s->data_series |= CRAM_BF;

	if (fd->required_fields & SAM_RNAME)
	    s->data_series |= CRAM_RI | CRAM_BF;

	if (fd->required_fields & SAM_POS)
	    s->data_series |= CRAM_AP | CRAM_BF;

	if (fd->required_fields & SAM_MAPQ)
	    s->data_series |= CRAM_MQ;

	if (fd->required_fields & SAM_CIGAR)
	    s->data_series |= CRAM_CIGAR;

	if (fd->required_fields & SAM_RNEXT)
	    s->data_series |= CRAM_CF | CRAM_NF | CRAM_RI | CRAM_NS |CRAM_BF;

	if (fd->required_fields & SAM_PNEXT)
	    s->data_series |= CRAM_CF | CRAM_NF | CRAM_AP | CRAM_NP | CRAM_BF;

	if (fd->required_fields & SAM_TLEN)
	    s->data_series |= CRAM_CF | CRAM_NF | CRAM_AP | CRAM_TS |
		CRAM_BF | CRAM_MF | CRAM_RI | CRAM_CIGAR;

	if (fd->required_fields & SAM_SEQ)
	    s->data_series |= CRAM_SEQ;

	if (!(fd->required_fields & SAM_AUX))
	    // No easy way to get MD/NM without other tags at present
	    s->decode_md = 0;

	if (fd->required_fields & SAM_QUAL) {
	    s->data_series |= CRAM_QUAL;
	    if (CRAM_MAJOR_VERS(fd->version) >= 4)
		s->data_series |= CRAM_BF;
	}

	if (fd->required_fields & SAM_AUX)
	    s->data_series |= CRAM_RG | CRAM_TL | CRAM_aux;

	if (fd->required_fields & SAM_RGAUX)
	    s->data_series |= CRAM_RG | CRAM_BF;

	// Always uncompress CORE block
	if (cram_uncompress_block(s->block[0]))
	    return -1;
    } else {
	s->data_series = CRAM_ALL;

	for (i = 0; i < s->hdr->num_blocks; i++) {
	    if (cram_uncompress_block(s->block[i]))
		return -1;
	}

	return 0;
    }

    block_used = calloc(s->hdr->num_blocks+1, sizeof(int));
    if (!block_used)
	return -1;

    do {
	/*
	 * Also set data_series based on code prerequisites. Eg if we need
	 * CRAM_QS then we also need to know CRAM_RL so we know how long it
	 * is, or if we need FC/FP then we also need FN (number of features).
	 *
	 * It's not reciprocal though. We may be needing to decode FN
	 * but have no need to decode FC, FP and cigar ops.
	 */
	if (s->data_series & CRAM_RS)    s->data_series |= CRAM_FC|CRAM_FP;
	if (s->data_series & CRAM_PD)    s->data_series |= CRAM_FC|CRAM_FP;
	if (s->data_series & CRAM_HC)    s->data_series |= CRAM_FC|CRAM_FP;
	if (s->data_series & CRAM_QS)    s->data_series |= CRAM_FC|CRAM_FP;
	if (s->data_series & CRAM_IN)    s->data_series |= CRAM_FC|CRAM_FP;
	if (s->data_series & CRAM_SC)    s->data_series |= CRAM_FC|CRAM_FP;
	if (s->data_series & CRAM_BS)    s->data_series |= CRAM_FC|CRAM_FP;
	if (s->data_series & CRAM_DL)    s->data_series |= CRAM_FC|CRAM_FP;
	if (s->data_series & CRAM_BA)    s->data_series |= CRAM_FC|CRAM_FP;
	if (s->data_series & CRAM_BB)    s->data_series |= CRAM_FC|CRAM_FP;
	if (s->data_series & CRAM_QQ)    s->data_series |= CRAM_FC|CRAM_FP;

	// cram_decode_seq() needs seq[] array
	if (s->data_series & (CRAM_SEQ|CRAM_CIGAR)) s->data_series |= CRAM_RL;

	if (s->data_series & CRAM_FP)    s->data_series |= CRAM_FC;
	if (s->data_series & CRAM_FC)    s->data_series |= CRAM_FN;
	if (s->data_series & CRAM_aux)   s->data_series |= CRAM_TL;
	if (s->data_series & CRAM_MF)    s->data_series |= CRAM_CF;
	if (s->data_series & CRAM_MQ)    s->data_series |= CRAM_BF;
	if (s->data_series & CRAM_BS)    s->data_series |= CRAM_RI;
	if (s->data_series & (CRAM_MF |CRAM_NS |CRAM_NP |CRAM_TS |CRAM_NF))
	    s->data_series |= CRAM_CF;
	if (!hdr->read_names_included && s->data_series & CRAM_RN)
	    s->data_series |= CRAM_CF | CRAM_NF;
	if (s->data_series & (CRAM_BA | CRAM_QS | CRAM_BB | CRAM_QQ))
	    s->data_series |= CRAM_BF | CRAM_CF | CRAM_RL;

	orig_ds = s->data_series;

	// Find which blocks are in use.
	for (i = 0; i < sizeof(i_to_id)/sizeof(*i_to_id); i++) {
	    int bnum1, bnum2, j;
	    cram_codec *c = hdr->codecs[i_to_id[i]];

	    if (!(s->data_series & (1<<i)))
		continue;

	    if (!c)
		continue;

	    bnum1 = cram_codec_to_id(c, &bnum2);

	    for (;;) {
		switch (bnum1) {
		case -2:
		    break;

		case -1:
		    core_used = 1;
		    break;

		default:
		    for (j = 0; j < s->hdr->num_blocks; j++) {
			if (s->block[j]->content_type == EXTERNAL &&
			    s->block[j]->content_id == bnum1) {
			    block_used[j] = 1;
			    if (cram_uncompress_block(s->block[j])) {
				free(block_used);
				return -1;
			    }
			}
		    }
		    break;
		}

		if (bnum2 == -2 || bnum1 == bnum2)
		    break;

		bnum1 = bnum2; // 2nd pass
	    }
	}

	// Tags too
	if ((fd->required_fields & SAM_AUX) ||
	    (s->data_series & CRAM_aux)) {
	    for (i = 0; i < CRAM_MAP_HASH; i++) {
		int bnum1, bnum2, j;
		cram_map *m = hdr->tag_encoding_map[i];

		while (m) {
		    cram_codec *c = m->codec;
		    if (!c) {
			m = m->next;
			continue;
		    }

		    bnum1 = cram_codec_to_id(c, &bnum2);

		    for (;;) {
			switch (bnum1) {
			case -2:
			    break;

			case -1:
			    core_used = 1;
			    break;

			default:
			    for (j = 0; j < s->hdr->num_blocks; j++) {
				if (s->block[j]->content_type == EXTERNAL &&
				    s->block[j]->content_id == bnum1) {
				    block_used[j] = 1;
				    if (cram_uncompress_block(s->block[j])) {
					free(block_used);
					return -1;
				    }
				}
			    }
			    break;
			}

			if (bnum2 == -2 || bnum1 == bnum2)
			    break;

			bnum1 = bnum2; // 2nd pass
		    }

		    m = m->next;
		}
	    }
	}

	// We now know which blocks are in used, so repeat and find
	// which other data series need to be added.
	for (i = 0; i < sizeof(i_to_id)/sizeof(*i_to_id); i++) {
	    int bnum1, bnum2, j;
	    cram_codec *c = hdr->codecs[i_to_id[i]];

	    if (!c)
		continue;

	    bnum1 = cram_codec_to_id(c, &bnum2);

	    for (;;) {
		switch (bnum1) {
		case -2:
		    break;

		case -1:
		    if (core_used) {
			//printf(" + data series %08x:\n", 1<<i);
			s->data_series |= 1<<i;
		    }
		    break;

		default:
		    for (j = 0; j < s->hdr->num_blocks; j++) {
			if (s->block[j]->content_type == EXTERNAL &&
			    s->block[j]->content_id == bnum1) {
			    if (block_used[j]) {
				//printf(" + data series %08x:\n", 1<<i);
				s->data_series |= 1<<i;
			    }
			}
		    }
		    break;
		}

		if (bnum2 == -2 || bnum1 == bnum2)
		    break;

		bnum1 = bnum2; // 2nd pass
	    }
	}

	// Tags too
	for (i = 0; i < CRAM_MAP_HASH; i++) {
	    int bnum1, bnum2, j;
	    cram_map *m = hdr->tag_encoding_map[i];

	    while (m) {
		cram_codec *c = m->codec;
		if (!c) {
		    m = m->next;
		    continue;
		}

		bnum1 = cram_codec_to_id(c, &bnum2);
		
		for (;;) {
		    switch (bnum1) {
		    case -2:
			break;

		    case -1:
			//printf(" + data series %08x:\n", CRAM_aux);
			s->data_series |= CRAM_aux;
			break;

		    default:
			for (j = 0; j < s->hdr->num_blocks; j++) {
			    if (s->block[j]->content_type == EXTERNAL &&
				s->block[j]->content_id == bnum1) {
				if (block_used[j]) {
				    //printf(" + data series %08x:\n",
				    //       CRAM_aux);
				    s->data_series |= CRAM_aux;
				}
			    }
			}
			break;
		    }

		    if (bnum2 == -2 || bnum1 == bnum2)
			break;

		    bnum1 = bnum2; // 2nd pass
		}

		m = m->next;
	    }
	}
    } while (orig_ds != s->data_series);

    free(block_used);
    return 0;
}

/*
 * Checks whether an external block is used solely by a single data series.
 * Returns the codec type if so (EXTERNAL, BYTE_ARRAY_LEN, BYTE_ARRAY_STOP)
 *         or 0 if not (E_NULL).
 */
static int cram_ds_unique(cram_block_compression_hdr *hdr, cram_codec *c,
			  int id) {
    int i, n_id = 0;
    enum cram_encoding e_type = 0;

    for (i = 0; i < DS_END; i++) {
	cram_codec *c;
	int bnum1, bnum2, old_n_id;

	if (!(c = hdr->codecs[i]))
	    continue;

	bnum1 = cram_codec_to_id(c, &bnum2);

	old_n_id = n_id;
	if (bnum1 == id) {
	    n_id++;
	    e_type = c->codec;
	}
	if (bnum2 == id) {
	    n_id++;
	    e_type = c->codec;
	}

	if (n_id == old_n_id+2)
	    n_id--; // len/val in same place counts once only.
    }

    return n_id == 1 ? e_type : 0;
}

/*
 * Attempts to estimate the size of some blocks so we can preallocate them
 * before decoding.  Although decoding will automatically grow the blocks,
 * it is typically more efficient to preallocate.
 */
void cram_decode_estimate_sizes(cram_block_compression_hdr *hdr, cram_slice *s,
				int *qual_size, int *name_size,
				int *q_id) {
    int bnum1, bnum2;
    cram_codec *cd;

    /* Qual */
    if (!(cd = hdr->codecs[DS_QS]))
	return;

    bnum1 = cram_codec_to_id(cd, &bnum2);
    if (bnum1 < 0 && bnum2 >= 0) bnum1 = bnum2;
    if (cram_ds_unique(hdr, cd, bnum1)) {
	cram_block *b = cram_get_block_by_id(s, bnum1);
	if (b && *qual_size < b->uncomp_size) *qual_size = b->uncomp_size;
	if (q_id && cd->codec == E_EXTERNAL)
	    *q_id = bnum1;
    }

    /* Name */
    cd = hdr->codecs[DS_RN];
    bnum1 = cram_codec_to_id(cd, &bnum2);
    if (bnum1 < 0 && bnum2 >= 0) bnum1 = bnum2;
    if (cram_ds_unique(hdr, cd, bnum1)) {
	cram_block *b = cram_get_block_by_id(s, bnum1);
	if (b && *name_size < b->uncomp_size) *name_size = b->uncomp_size;
    }
}


/* ----------------------------------------------------------------------
 * CRAM slices
 */

/*
 * Decodes a CRAM (un)mapped slice header block.
 * Returns slice header ptr on success
 *         NULL on failure
 */
cram_block_slice_hdr *cram_decode_slice_header(cram_fd *fd, cram_block *b) {
    cram_block_slice_hdr *hdr;
    unsigned char *cp;
    unsigned char *cp_end;
    int i, err = 0;

    if (b->method != RAW) {
        /* Spec. says slice header should be RAW, but we can future-proof
	   by trying to decode it if it isn't. */
        if (cram_uncompress_block(b) < 0)
            return NULL;
    }

    cp = (unsigned char *)BLOCK_DATA(b);
    cp_end = cp + b->uncomp_size;

    if (b->content_type != MAPPED_SLICE &&
	b->content_type != UNMAPPED_SLICE)
	return NULL;

    if (!(hdr  = calloc(1, sizeof(*hdr))))
	return NULL;

    hdr->content_type = b->content_type;

    if (b->content_type == MAPPED_SLICE) {
        hdr->ref_seq_id = fd->vv.varint_get32s((char **)&cp, (char *)cp_end, &err);
	if (CRAM_MAJOR_VERS(fd->version) >= 4) {
	    hdr->ref_seq_start = fd->vv.varint_get64((char **)&cp, (char *)cp_end, &err);
	    hdr->ref_seq_span  = fd->vv.varint_get64((char **)&cp, (char *)cp_end, &err);
	} else {
	    hdr->ref_seq_start = fd->vv.varint_get32((char **)&cp, (char *)cp_end, &err);
	    hdr->ref_seq_span  = fd->vv.varint_get32((char **)&cp, (char *)cp_end, &err);
	}
    }
    hdr->num_records = fd->vv.varint_get32((char **)&cp, (char *) cp_end, &err);
    hdr->record_counter = 0;
    if (CRAM_MAJOR_VERS(fd->version) == 2) {
	hdr->record_counter = fd->vv.varint_get32((char **)&cp, (char *)cp_end, &err);
    } else if (CRAM_MAJOR_VERS(fd->version) >= 3) {
	hdr->record_counter = fd->vv.varint_get64((char **)&cp, (char *)cp_end, &err);
    }
    hdr->num_blocks      = fd->vv.varint_get32((char **)&cp, (char *)cp_end, &err);
    hdr->num_content_ids = fd->vv.varint_get32((char **)&cp, (char *)cp_end, &err);
    if (hdr->num_content_ids < 1 ||
	hdr->num_content_ids >= SIZE_MAX / sizeof(int32_t)) {
        /* Slice must have at least one data block,
	   and malloc'd size shouldn't wrap. */
        free(hdr);
        return NULL;
    }
    hdr->block_content_ids = malloc(hdr->num_content_ids * sizeof(int32_t));
    if (!hdr->block_content_ids) {
	free(hdr);
	return NULL;
    }

    for (i = 0; i < hdr->num_content_ids; i++) {
	hdr->block_content_ids[i] = fd->vv.varint_get32((char **)&cp,
							 (char *)cp_end,
							 &err);
    }
    if (err) {
	free(hdr->block_content_ids);
	free(hdr);
	return NULL;
    }

    if (b->content_type == MAPPED_SLICE) {
        hdr->ref_base_id = fd->vv.varint_get32((char **)&cp, (char *) cp_end, &err);
    }

    if (!IS_CRAM_1_VERS(fd)) {
        if (cp_end - cp < 16) {
            free(hdr->block_content_ids);
            free(hdr);
            return NULL;
        }
	memcpy(hdr->md5, cp, 16);
	cp += 16;
    } else {
	memset(hdr->md5, 0, 16);
    }

    // Decode any optional tag:type:value fields
    if (cp == cp_end)
	return hdr;

    hdr->tags = HashTableCreate(4, HASH_FUNC_TCL);
    while (cp <= cp_end) {
	unsigned int sub_len;
	unsigned char id[3];
	HashData hd;

	if (cp_end - cp < 4)
	    return hdr;

	id[0] = cp[0];
	id[1] = cp[1];
	id[2] = '\0';

	switch (cp[2]) {
	case 'c':
	    id[2] = 'i';
	    hd.i = (int8_t)(cp[3]);
	    cp += 4;
	    break;

	case 'C':
	    id[2] = 'i';
	    hd.i = (uint8_t)(cp[3]);
	    cp += 4;
	    break;

	case 's':
	    if (cp_end - cp < 5) break;
	    id[2] = 'i';
	    hd.i = (int16_t)(cp[3] + (cp[4]<<8));
	    cp += 5;
	    break;

	case 'S':
	    if (cp_end - cp < 5) break;
	    id[2] = 'i';
	    hd.i = (uint16_t)(cp[3] + (cp[4]<<8));
	    cp += 5;
	    break;

	case 'i':
	    if (cp_end - cp < 7) break;
	    id[2] = 'i';
	    hd.i = (int32_t)(cp[3] + (cp[4]<<8) + (cp[5]<<16) + (cp[6]<<24));
	    cp += 7;
	    break;

	case 'I':
	    if (cp_end - cp < 7) break;
	    id[2] = 'i';
	    hd.i = (uint32_t)(cp[3] + (cp[4]<<8) + (cp[5]<<16) + (cp[6]<<24));
	    cp += 7;
	    break;

	case 'f':
	    if (cp_end - cp < 7) break;
	    id[2] = cp[2];
	    hd.f = bam_aux_f(cp+2);
	    cp += 7;
	    break;

	case 'A':
	    id[2] = cp[2];
	    hd.i = cp[3];
	    cp += 4;
	    break;

	case 'Z': case 'H':
	    id[2] = cp[2];
	    hd.p = &cp[3];
	    cp += 3;
	    while (cp < cp_end && *cp != '\0') cp++;
	    if (cp < cp_end) {
	        /* Skip NUL */
	        cp++;
	    } else {
	        /* Add missing NUL termination */
	        assert(cp == BLOCK_DATA(b) + b->uncomp_size);
	        BLOCK_RESIZE(b, b->uncomp_size + 1);
		cp = cp_end = BLOCK_DATA(b) + b->uncomp_size;
		*cp = '\0';
	    }
	    break;

	case 'B':
	    // <type>B<code><4-len><...>
	    if (cp_end - cp < 8) break;
	    hd.p = &cp[3];
	    sub_len = cp[4] + (cp[5]<<8) + (cp[6]<<16) + (cp[7]<<24);
	    switch (cp[3]) {
	    case 'c': case 'C':
	        if (cp_end - cp < 8 + sub_len) break;
	        id[2] = cp[2];
	        cp += 8 + sub_len;
	        break;
	    case 's': case 'S':
	        if (cp_end - cp < 8 + 2 * sub_len) break;
	        id[2] = cp[2];
	        cp += 8 + 2*sub_len;
	        break;
	    case 'i': case 'I': case 'f':
	        if (cp_end - cp < 8 + 4 * sub_len) break;
	        id[2] = cp[2];
	        cp += 8 + 4*sub_len;
	        break;
	    default:
	        fprintf(stderr, "Unknown aux type 'B' sub-code.\n");
	        cp = cp_end;
	    }
	    break;

        default:
	    fprintf(stderr, "Unknown aux type.\n");
	    cp = cp_end;
	}

	if (id[2] != '\0') {
	    HashTableAdd(hdr->tags, (char *)id, 3, hd, NULL);
	} else {
	    cp = cp_end;
	    break;
	}

	if (id[0] == 'B' && id[1] == 'D' && id[2] == 'B') {
	    unsigned char *p = hd.p;
	    hdr->BD_crc =  p[5] | (p[6]<<8) | (p[7]<<16) | (p[8]<<24);
	}

	if (id[0] == 'S' && id[1] == 'D' && id[2] == 'B') {
	    unsigned char *p = hd.p;
	    hdr->SD_crc =  p[5] | (p[6]<<8) | (p[7]<<16) | (p[8]<<24);
	}
    }

    if (!err)
	return hdr;

    free(hdr);
    return NULL;
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
			   cram_block *blk, cram_record *cr, SAM_hdr *bfd,
			   int cf, char *seq, char *qual,
			   int has_MD, int has_NM) {
    int prev_pos = 0, f, r = 0, out_sz = 1;
    int seq_pos = 1;
    int cig_len = 0;
    int64_t ref_pos = cr->apos;
    int32_t fn, i32;
    enum cigar_op cig_op = BAM_CMATCH;
    uint32_t *cigar = s->cigar;
    uint32_t ncigar = s->ncigar;
    uint32_t cigar_alloc = s->cigar_alloc;
    uint32_t nm = 0;
    int32_t md_dist = 0;
    int orig_aux = 0;
    int decode_md = s->ref && ((s->decode_md && !has_MD) || has_MD < 0);
    int decode_nm = s->ref && ((s->decode_md && !has_NM) || has_NM < 0);
    uint32_t ds = s->data_series;

    if ((ds & CRAM_QS) && !(cf & CRAM_FLAG_PRESERVE_QUAL_SCORES)) {
	memset(qual, 255, cr->len);
    }

    if (cr->cram_flags & CRAM_FLAG_NO_SEQ)
        decode_md = decode_nm = 0;

    if (decode_md) {
	orig_aux = BLOCK_SIZE(s->aux_blk); // end of existing aux block
	if (has_MD == 0)
	    BLOCK_APPEND(s->aux_blk, "MDZ", 3);
    }
    
    if (ds & CRAM_FN) {
	if (!c->comp_hdr->codecs[DS_FN]) return -1;
	r |= c->comp_hdr->codecs[DS_FN]->decode(s,c->comp_hdr->codecs[DS_FN],
						blk, (char *)&fn, &out_sz);
        if (r) return r;
    } else {
	fn = 0;
    }

    ref_pos--; // count from 0
    cr->cigar = ncigar;

    if (!(ds & (CRAM_FC | CRAM_FP)))
	goto skip_cigar;

    for (f = 0; f < fn; f++) {
	int32_t pos = 0;
	char op;

	if (ncigar+2 >= cigar_alloc) {
	    cigar_alloc = cigar_alloc ? cigar_alloc*2 : 1024;
	    s->cigar = cigar;
	    if (!(cigar = realloc(cigar, cigar_alloc * sizeof(*cigar))))
		return -1;
	}

	if (ds & CRAM_FC) {
	    if (!c->comp_hdr->codecs[DS_FC]) return -1;
	    r |= c->comp_hdr->codecs[DS_FC]->decode(s,
						    c->comp_hdr->codecs[DS_FC],
						    blk,
						    &op,  &out_sz);
	    if (r) return r;
	}

	if (!(ds & CRAM_FP))
	    continue;

	if (!c->comp_hdr->codecs[DS_FP]) return -1;
	r |= c->comp_hdr->codecs[DS_FP]->decode(s,
						c->comp_hdr->codecs[DS_FP],
						blk,
						(char *)&pos, &out_sz);
	if (r) return r;
	pos += prev_pos;

	if (pos <= 0) {
	    fprintf(stderr, "Error: feature position %d before start of read.\n",
		    pos);
	    return -1;
	}

	if (pos > seq_pos) {
	    if (pos > cr->len+1)
		return -1;

	    if (s->ref && cr->ref_id >= 0) {
		if (ref_pos + pos - seq_pos > bfd->ref[cr->ref_id].len) {
		    static int whinged = 0;
		    int rlen;
		    if (!whinged)
			fprintf(stderr, "Ref pos outside of ref "
				"sequence boundary\n");
		    whinged = 1;
		    rlen = bfd->ref[cr->ref_id].len - ref_pos;
		    if (rlen > 0) {
			memcpy(&seq[seq_pos-1],
			       &s->ref[ref_pos - s->ref_start +1], rlen);
			if ((pos - seq_pos) - rlen > 0)
			    memset(&seq[seq_pos-1+rlen], 'N',
				   (pos - seq_pos) - rlen);
		    } else {
		        memset(&seq[seq_pos-1], 'N', cr->len - seq_pos + 1);
		    }
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
	    if (md_dist >= 0)
		md_dist += pos - seq_pos;
	    seq_pos = pos;
	}

	prev_pos = pos;

	if (!(ds & CRAM_FC))
	    goto skip_cigar;

	switch(op) {
	case 'S': { // soft clip: IN
	    int32_t out_sz2 = 1;
	    int have_sc = 0;

	    if (cig_len) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    switch (CRAM_MAJOR_VERS(fd->version)) {
	    case 1:
	        if (ds & CRAM_IN) {
		    r |= c->comp_hdr->codecs[DS_IN]
			? c->comp_hdr->codecs[DS_IN]
			             ->decode(s, c->comp_hdr->codecs[DS_IN],
					      blk,
					      cr->len ? &seq[pos-1] : NULL,
					      &out_sz2)
			: (seq[pos-1] = 'N', out_sz2 = 1, 0);
		    have_sc = 1;
		}
		break;
	    case 2:
	    default:
	        if (ds & CRAM_SC) {
		    r |= c->comp_hdr->codecs[DS_SC]
			? c->comp_hdr->codecs[DS_SC]
			             ->decode(s, c->comp_hdr->codecs[DS_SC],
					      blk,
					      cr->len ? &seq[pos-1] : NULL,
					      &out_sz2)
			: (seq[pos-1] = 'N', out_sz2 = 1, 0);
		    have_sc = 1;
		}
		break;

//		default:
//		    r |= c->comp_hdr->codecs[DS_BB]
//			? c->comp_hdr->codecs[DS_BB]
//			             ->decode(s, c->comp_hdr->codecs[DS_BB],
//					      blk, &seq[pos-1], &out_sz2)
//			: (seq[pos-1] = 'N', out_sz2 = 1, 0);
	    }
	    if (have_sc) {
		if (r) return r;
		cigar[ncigar++] = (out_sz2<<4) + BAM_CSOFT_CLIP;
		cig_op = BAM_CSOFT_CLIP;
		seq_pos += out_sz2;
	    }
	    break;
	}

	case 'X': { // Substitution; BS
	    unsigned char base;
#ifdef USE_X
	    if (cig_len && cig_op != BAM_CBASE_MISMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_BS) {
		if (!c->comp_hdr->codecs[DS_BS]) return -1;
		r |= c->comp_hdr->codecs[DS_BS]
		                ->decode(s, c->comp_hdr->codecs[DS_BS], blk,
					 (char *)&base, &out_sz);
		if (pos-1 < cr->len)
		    seq[pos-1] = 'N'; // FIXME look up BS=base value
	    }
	    cig_op = BAM_CBASE_MISMATCH;
#else
	    int ref_base;
	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_BS) {
		if (!c->comp_hdr->codecs[DS_BS]) return -1;
		r |= c->comp_hdr->codecs[DS_BS]
		                ->decode(s, c->comp_hdr->codecs[DS_BS], blk,
					 (char *)&base, &out_sz);
		if (r) return -1;
		if (ref_pos >= bfd->ref[cr->ref_id].len || !s->ref) {
		    if (pos-1 < cr->len)
			seq[pos-1] = c->comp_hdr->
			    substitution_matrix[fd->L1['N']][base];
		    if (decode_md || decode_nm) {
			if (md_dist >= 0 && decode_md)
			    BLOCK_APPEND_UINT(s->aux_blk, md_dist);
			md_dist = -1;
			nm--;
		    }
		} else {
                    unsigned char ref_call = ref_pos <= s->ref_end
                        ? (uc)s->ref[ref_pos - s->ref_start +1]
                        : 'N';
                    ref_base = fd->L1[ref_call];
		    if (pos-1 < cr->len)
			seq[pos-1] = c->comp_hdr->
			    substitution_matrix[ref_base][base];
		    if (decode_md) {
                        BLOCK_APPEND_UINT(s->aux_blk, md_dist);
                        BLOCK_APPEND_CHAR(s->aux_blk, ref_call);
			md_dist = 0;
		    }
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
	    if (ds & CRAM_DL) {
		if (!c->comp_hdr->codecs[DS_DL]) return -1;
		r |= c->comp_hdr->codecs[DS_DL]
		                ->decode(s, c->comp_hdr->codecs[DS_DL], blk,
					 (char *)&i32, &out_sz);
		if (r) return r;
		if (decode_md || decode_nm) {
		    if (md_dist >= 0 && decode_md)
			BLOCK_APPEND_UINT(s->aux_blk, md_dist);
		    if (ref_pos + i32 <= bfd->ref[cr->ref_id].len) {
			if (decode_md) {
			    BLOCK_APPEND_CHAR(s->aux_blk, '^');
			    BLOCK_APPEND(s->aux_blk,
					 &s->ref[ref_pos - s->ref_start +1],
					 i32);
			    md_dist = 0;
			}
			nm += i32;
		    } else {
			uint32_t dlen;
			if (bfd->ref[cr->ref_id].len >= ref_pos) {
			    if (decode_md) {
				BLOCK_APPEND_CHAR(s->aux_blk, '^');
				BLOCK_APPEND(s->aux_blk,
					     &s->ref[ref_pos - s->ref_start+1],
					     bfd->ref[cr->ref_id].len-ref_pos);
				BLOCK_APPEND_UINT(s->aux_blk, 0);
			    }
			    dlen = i32 - (bfd->ref[cr->ref_id].len - ref_pos);
			    nm += i32 - dlen;
			} else {
			    dlen = i32;
			}

			md_dist = -1;
		    }
		}
		cig_op = BAM_CDEL;
		cig_len += i32;
		ref_pos += i32;
		//printf("  %d: DL = %d (ret %d)\n", f, i32, r);
	    }
	    break;
	}

	case 'I': { // Insertion (several bases); IN
	    int32_t out_sz2 = 1;

	    if (cig_len && cig_op != BAM_CINS) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }

	    if (ds & CRAM_IN) {
		if (!c->comp_hdr->codecs[DS_IN]) return -1;
		r |= c->comp_hdr->codecs[DS_IN]
		                ->decode(s, c->comp_hdr->codecs[DS_IN], blk,
					 cr->len ? &seq[pos-1] : NULL,
					 &out_sz2);
		if (r) return r;
		cig_op = BAM_CINS;
		cig_len += out_sz2;
		seq_pos += out_sz2;
		nm      += out_sz2;
		//printf("  %d: IN(I) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, dat, r, out_sz2);
	    }
	    break;
	}

	case 'i': { // Insertion (single base); BA
	    if (cig_len && cig_op != BAM_CINS) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_BA) {
		if (!c->comp_hdr->codecs[DS_BA]) return -1;
		r |= c->comp_hdr->codecs[DS_BA]
		                ->decode(s, c->comp_hdr->codecs[DS_BA], blk,
					 cr->len ? &seq[pos-1] : NULL,
					 &out_sz);
		if (r) return r;
		//printf("  %d: BA = %c (ret %d)\n", f, seq[pos-1], r);
	    }
	    cig_op = BAM_CINS;
	    cig_len++;
	    seq_pos++;
	    nm++;
	    break;
	}

	case 'b': { // Several bases
	    int32_t len = 1;

	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }

	    if (ds & CRAM_BB) {
		if (!c->comp_hdr->codecs[DS_BB]) return -1;
		r |= c->comp_hdr->codecs[DS_BB]
		    ->decode(s, c->comp_hdr->codecs[DS_BB], blk,
			     cr->len ? &seq[pos-1] : NULL,
			     &len);
		if (r) return r;

		if (decode_md || decode_nm) {
		    int x;
		    if (md_dist >= 0 && decode_md)
			BLOCK_APPEND_UINT(s->aux_blk, md_dist);

		    for (x = 0; x < len; x++) {
			if (x && decode_md)
			    BLOCK_APPEND_UINT(s->aux_blk, 0);
			if (ref_pos+x >= bfd->ref[cr->ref_id].len || !s->ref) {
			    md_dist = -1;
			    break;
			} else {
			    if (decode_md) {
				char r = s->ref[ref_pos+x-s->ref_start +1];
				BLOCK_APPEND_CHAR(s->aux_blk, r);
			    }
			}
		    }

		    nm += x;
		    md_dist = 0;
		}
	    }

	    cig_op = BAM_CMATCH;

	    cig_len+=len;
	    seq_pos+=len;
	    ref_pos+=len;
	    //prev_pos+=len;
	    break;
	}

	case 'q': { // Several quality values
	    int32_t len = 1;

	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }

	    if (ds & CRAM_QQ) {
                if (!c->comp_hdr->codecs[DS_QQ]) return -1;
                if ((ds & CRAM_QS) && !(cf & CRAM_FLAG_PRESERVE_QUAL_SCORES)
		    && (unsigned char)*qual == 255)
                    memset(qual, 30, cr->len); // ASCII '?' = qual 30
		if (!c->comp_hdr->codecs[DS_QQ]) return -1;
		r |= c->comp_hdr->codecs[DS_QQ]
		    ->decode(s, c->comp_hdr->codecs[DS_QQ], blk,
			     (char *)&qual[pos-1], &len);
		if (r) return r;
	    }

	    cig_op = BAM_CMATCH;

	    //prev_pos+=len;
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
	    if (ds & CRAM_BA) {
		if (!c->comp_hdr->codecs[DS_BA]) return -1;
		r |= c->comp_hdr->codecs[DS_BA]
		                ->decode(s, c->comp_hdr->codecs[DS_BA], blk,
					 cr->len ? &seq[pos-1] : NULL,
					 &out_sz);

		if (decode_md || decode_nm) {
		    if (md_dist >= 0 && decode_md)
			BLOCK_APPEND_UINT(s->aux_blk, md_dist);
		    if (ref_pos >= bfd->ref[cr->ref_id].len || !s->ref) {
			md_dist = -1;
		    } else {
			if (decode_md)
			    BLOCK_APPEND_CHAR(s->aux_blk,
					      s->ref[ref_pos-s->ref_start +1]);
			nm++;
			md_dist = 0;
		    }
		}
	    }
	    if (ds & CRAM_QS) {
		if (!c->comp_hdr->codecs[DS_QS]) return -1;
                if (!(cf & CRAM_FLAG_PRESERVE_QUAL_SCORES)
		    && (unsigned char)*qual == 255)
                    memset(qual, 30, cr->len); // ASCII '?' = qual 30
		r |= c->comp_hdr->codecs[DS_QS]
		                ->decode(s, c->comp_hdr->codecs[DS_QS], blk,
					 (char *)&qual[pos-1], &out_sz);
	    }
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
	    if (ds & CRAM_QS) {
		if (!c->comp_hdr->codecs[DS_QS]) return -1;
                if (!(cf & CRAM_FLAG_PRESERVE_QUAL_SCORES)
		    && (unsigned char)*qual == 255)
                    memset(qual, 30, cr->len); // ASCII '?' = qual 30
		r |= c->comp_hdr->codecs[DS_QS]
		                ->decode(s, c->comp_hdr->codecs[DS_QS], blk,
					 (char *)&qual[pos-1], &out_sz);
		//printf("  %d: QS = %d (ret %d)\n", f, qc, r);
	    }
	    break;
	}

	case 'H': { // hard clip; HC
	    if (cig_len && cig_op != BAM_CHARD_CLIP) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_HC) {
		if (!c->comp_hdr->codecs[DS_HC]) return -1;
		r |= c->comp_hdr->codecs[DS_HC]
		                ->decode(s, c->comp_hdr->codecs[DS_HC], blk,
					 (char *)&i32, &out_sz);
		if (r) return r;
		cig_op = BAM_CHARD_CLIP;
		cig_len += i32;
	    }
	    break;
	}

	case 'P': { // padding; PD
	    if (cig_len && cig_op != BAM_CPAD) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_PD) {
		if (!c->comp_hdr->codecs[DS_PD]) return -1;
		r |= c->comp_hdr->codecs[DS_PD]
		                ->decode(s, c->comp_hdr->codecs[DS_PD], blk,
					 (char *)&i32, &out_sz);
		if (r) return r;
		cig_op = BAM_CPAD;
		cig_len += i32;
	    }
	    break;
	}

	case 'N': { // Ref skip; RS
	    if (cig_len && cig_op != BAM_CREF_SKIP) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_RS) {
		if (!c->comp_hdr->codecs[DS_RS]) return -1;
		r |= c->comp_hdr->codecs[DS_RS]
		                ->decode(s, c->comp_hdr->codecs[DS_RS], blk,
					 (char *)&i32, &out_sz);
		if (r) return r;
		cig_op = BAM_CREF_SKIP;
		cig_len += i32;
		ref_pos += i32;
	    }
	    break;
	}

	default:
            fprintf(stderr, "Error: Unknown feature code '%c'\n", op);
	    return -1;
	}
    }

    if (!(ds & CRAM_FC))
	goto skip_cigar;

    /* An implicit match op for any unaccounted for bases */
    if ((ds & CRAM_FN) && cr->len >= seq_pos) {
	if (s->ref) {
	    if (ref_pos + cr->len - seq_pos + 1 > bfd->ref[cr->ref_id].len) {
		static int whinged = 0;
		int rlen;
		if (!whinged)
		    fprintf(stderr, "Ref pos outside of ref sequence boundary\n");
		whinged = 1;
		rlen = bfd->ref[cr->ref_id].len - ref_pos;
		if (rlen > 0) {
		    if (seq_pos-1 + rlen < cr->len)
			memcpy(&seq[seq_pos-1],
			       &s->ref[ref_pos - s->ref_start +1], rlen);
		    if ((cr->len - seq_pos + 1) - rlen > 0)
		        memset(&seq[seq_pos-1+rlen], 'N',
                               (cr->len - seq_pos + 1) - rlen);
		} else {
		    if (cr->len - seq_pos + 1 > 0)
			memset(&seq[seq_pos-1], 'N', cr->len - seq_pos + 1);
		}
	    } else {
		if (cr->len - seq_pos + 1 > 0)
		    memcpy(&seq[seq_pos-1], &s->ref[ref_pos - s->ref_start +1],
			   cr->len - seq_pos + 1);
		ref_pos += cr->len - seq_pos + 1;
		if (md_dist >= 0)
		    md_dist += cr->len - seq_pos + 1;
	    }
	} else {
	    // So alignment end can be computed even when not decoding sequence
	    ref_pos += cr->len - seq_pos + 1;
	}

	if (ncigar+1 >= cigar_alloc) {
	    cigar_alloc = cigar_alloc ? cigar_alloc*2 : 1024;
	    s->cigar = cigar;
	    if (!(cigar = realloc(cigar, cigar_alloc * sizeof(*cigar))))
		return -1;
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

 skip_cigar:

    if ((ds & CRAM_FN) && decode_md) {
	if (md_dist >= 0)
	    BLOCK_APPEND_UINT(s->aux_blk, md_dist);
    }

    if (cig_len) {
	if (ncigar >= cigar_alloc) {
	    cigar_alloc = cigar_alloc ? cigar_alloc*2 : 1024;
	    s->cigar = cigar;
	    if (!(cigar = realloc(cigar, cigar_alloc * sizeof(*cigar))))
		return -1;
	}

	cigar[ncigar++] = (cig_len<<4) + cig_op;
    }

    cr->ncigar = ncigar - cr->cigar;
    cr->aend = ref_pos;

    //printf("2: %.*s %d .. %d %d\n", cr->name_len, (char *)BLOCK_DATA(s->name_blk) + cr->name, cr->apos, ref_pos, seq_pos);

    if (ds & CRAM_MQ) {
	if (!c->comp_hdr->codecs[DS_MQ]) return -1;
	r |= c->comp_hdr->codecs[DS_MQ]
	                ->decode(s, c->comp_hdr->codecs[DS_MQ], blk,
				 (char *)&cr->mqual, &out_sz);
    } else {
	cr->mqual = 40;
    }

    if ((ds & CRAM_QS) && (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES)) {
	int32_t out_sz2 = cr->len;

	if (ds & CRAM_QS) {
	    if (!c->comp_hdr->codecs[DS_QS]) return -1;
	    r |= c->comp_hdr->codecs[DS_QS]
		            ->decode(s, c->comp_hdr->codecs[DS_QS], blk,
				     qual, &out_sz2);
	}
    }

    s->cigar = cigar;
    s->cigar_alloc = cigar_alloc;
    s->ncigar = ncigar;

    if (cr->cram_flags & CRAM_FLAG_NO_SEQ)
	cr->len = 0;

    if (decode_md) {
	BLOCK_APPEND_CHAR(s->aux_blk, '\0'); // null terminate MD:Z:
	size_t sz = BLOCK_SIZE(s->aux_blk) - orig_aux;
	if (has_MD < 0) {
	    // has_MD < 0; already have MDZ allocated in aux at -has_MD,
	    // but wrote MD to end of aux (at orig_aux).
	    // We need some memmoves to shuffle it around.
	    char tmp_MD_[1024], *tmp_MD = tmp_MD_;
	    unsigned char *orig_aux_p = BLOCK_DATA(s->aux_blk) + orig_aux;
	    if (sz > 1024) {
		tmp_MD = malloc(sz);
		if (!tmp_MD)
		    return -1;
	    }
	    memcpy(tmp_MD, orig_aux_p, sz);
	    memmove(&BLOCK_DATA(s->aux_blk)[-has_MD] + sz,
		    &BLOCK_DATA(s->aux_blk)[-has_MD],
		    orig_aux_p - &BLOCK_DATA(s->aux_blk)[-has_MD]);
	    memcpy(&BLOCK_DATA(s->aux_blk)[-has_MD], tmp_MD, sz);
	    if (tmp_MD != tmp_MD_)
		free(tmp_MD);

	    if (-has_NM > -has_MD)
		// we inserted before NM, so move it up a bit
		has_NM -= sz;
	}
	// else has_MD == 0 and we've already appended MD to the end.

	cr->aux_size += sz;
    }

    if (decode_nm) {
	if (has_NM == 0) {
	    char buf[7];
	    buf[0] = 'N'; buf[1] = 'M'; buf[2] = 'I';
	    buf[3] = (nm>> 0) & 0xff;
	    buf[4] = (nm>> 8) & 0xff;
	    buf[5] = (nm>>16) & 0xff;
	    buf[6] = (nm>>24) & 0xff;
	    BLOCK_APPEND(s->aux_blk, buf, 7);
	    cr->aux_size += 7;
	} else {
	    // Preallocated space for NM at -has_NM into aux block
	    unsigned char *buf = BLOCK_DATA(s->aux_blk) + -has_NM;
	    buf[0] = (nm>> 0) & 0xff;
	    buf[1] = (nm>> 8) & 0xff;
	    buf[2] = (nm>>16) & 0xff;
	    buf[3] = (nm>>24) & 0xff;
	}
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

    return m;
}

//#define map_find(M,K,I) M[CRAM_MAP(K[0],K[1])];while (m && m->key != I);m= m->next


static int cram_decode_aux_1_0(cram_container *c, cram_slice *s,
			       cram_block *blk, cram_record *cr) {
    int i, r = 0, out_sz = 1;
    unsigned char ntags;
	    
    if (!c->comp_hdr->codecs[DS_TC]) return -1;
    r |= c->comp_hdr->codecs[DS_TC]->decode(s, c->comp_hdr->codecs[DS_TC], blk,
				       (char *)&ntags, &out_sz);
    cr->ntags = ntags;

    //printf("TC=%d\n", cr->ntags);
    cr->aux_size = 0;
    cr->aux = BLOCK_SIZE(s->aux_blk);

    for (i = 0; i < cr->ntags; i++) {
	int32_t id, out_sz = 1;
	unsigned char tag_data[3];
	cram_map *m;

	//printf("Tag %d/%d\n", i+1, cr->ntags);
	if (!c->comp_hdr->codecs[DS_TN]) return -1;
	r |= c->comp_hdr->codecs[DS_TN]->decode(s, c->comp_hdr->codecs[DS_TN],
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
	if (!m)
	    return -1;
	BLOCK_APPEND(s->aux_blk, (char *)tag_data, 3);

	if (!m->codec) return -1;
	r |= m->codec->decode(s, m->codec, blk, (char *)s->aux_blk, &out_sz);

	cr->aux_size += out_sz + 3;
    }
    
    return r;
}

// has_MD and has_NM are filled out with 0 for none present,
// 1 for present and verbatim, and -pos for present as placeholder
// (MD*, NM*) to be generated and filled out at offset +pos.
static int cram_decode_aux(cram_fd *fd, cram_container *c, cram_slice *s,
			   cram_block *blk, cram_record *cr,
			   int *has_MD, int *has_NM) {
    int i, r = 0, out_sz = 1;
    int32_t TL = 0;
    unsigned char *TN;
    uint32_t ds = s->data_series;
	    
    if (!(ds & (CRAM_TL|CRAM_aux))) {
	cr->aux = 0;
	cr->aux_size = 0;
	return 0;
    }

    if (!c->comp_hdr->codecs[DS_TL]) return -1;
    r |= c->comp_hdr->codecs[DS_TL]->decode(s, c->comp_hdr->codecs[DS_TL], blk,
					    (char *)&TL, &out_sz);
    if (r || TL < 0 || TL >= c->comp_hdr->nTL)
	return -1;

    TN = c->comp_hdr->TL[TL];
    cr->ntags = strlen((char *)TN)/3; // optimise to remove strlen

    //printf("TC=%d\n", cr->ntags);
    cr->aux_size = 0;
    cr->aux = BLOCK_SIZE(s->aux_blk);

    if (!(ds & CRAM_aux))
	return 0;

    for (i = 0; i < cr->ntags; i++) {
	int32_t id, out_sz = 1;
	unsigned char tag_data[7];
	cram_map *m;

	if (TN[0] == 'M' && TN[1] == 'D' && has_MD)
	    *has_MD = (BLOCK_SIZE(s->aux_blk)+3) * (TN[2] == '*' ? -1 : 1);
	if (TN[0] == 'N' && TN[1] == 'M' && has_NM)
	    *has_NM = (BLOCK_SIZE(s->aux_blk)+3) * (TN[2] == '*' ? -1 : 1);;

	//printf("Tag %d/%d\n", i+1, cr->ntags);
	tag_data[0] = TN[0];
	tag_data[1] = TN[1];
	tag_data[2] = TN[2];
	if (IS_CRAM_4_VERS(fd) && TN[2] == '*') {
	    // Place holder, fill out contents later.
	    int tag_data_size;
	    if (TN[0] == 'N' && TN[1] == 'M') {
		// Use a fixed size, so we can allocate room for it now.
		memcpy(&tag_data[2], "I\0\0\0\0", 5);
		tag_data_size = 7;
	    } else if (TN[0] == 'R' && TN[1] == 'G') {
		// RG is variable size, but known already.  Insert now
		TN += 3;
		if (cr->rg < 0 || cr->rg >= fd->header->nrg)
		    continue;
		tag_data[2] = 'Z';
		BLOCK_APPEND(s->aux_blk, (char *)tag_data, 3);
		BLOCK_APPEND(s->aux_blk, fd->header->rg[cr->rg].name,
			     fd->header->rg[cr->rg].name_len);
		BLOCK_APPEND_CHAR(s->aux_blk, '\0');
		cr->aux_size += 3 + fd->header->rg[cr->rg].name_len + 1;
		cr->rg = -1; // prevents auto-add later
		continue;
	    } else {
		// Unknown size.  We'll insert MD into stream later.
		tag_data[2] = 'Z';
		tag_data_size = 3;
	    }
	    BLOCK_APPEND(s->aux_blk, (char *)tag_data, tag_data_size);
	    cr->aux_size += tag_data_size;
	    TN += 3;
	} else {
	    TN += 3;
	    // Actual tag, so copy into BAM tag data stream
	    id = (tag_data[0]<<16) | (tag_data[1]<<8) | tag_data[2];

	    m = map_find(c->comp_hdr->tag_encoding_map, tag_data, id);
	    if (!m)
		return -1;
	    BLOCK_APPEND(s->aux_blk, (char *)tag_data, 3);

	    if (!m->codec) return -1;
	    r |= m->codec->decode(s, m->codec, blk, (char *)s->aux_blk, &out_sz);
	    if (r) break;
	    cr->aux_size += out_sz + 3;
	}
    }
    
    return r;
}

/* Resolve mate pair cross-references between recs within this slice */
static int cram_decode_slice_xref(cram_slice *s, int required_fields) {
    int rec;

    if (!(required_fields & (SAM_RNEXT | SAM_PNEXT | SAM_TLEN))) {
	for (rec = 0; rec < s->hdr->num_records; rec++) {
	    cram_record *cr = &s->crecs[rec];

	    cr->tlen = 0;
	    cr->mate_pos = 0;
	    cr->mate_ref_id = -1;
	}

	return 0;
    }

    for (rec = 0; rec < s->hdr->num_records; rec++) {
	cram_record *cr = &s->crecs[rec];

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
		if (cr->tlen == INT64_MIN) {
		    int id1 = rec, id2 = rec;
		    int64_t aleft = cr->apos, aright = cr->aend;
		    int64_t tlen;
		    int ref = cr->ref_id;

		    // number of segments starting at the same point.
		    int left_cnt = 0;

		    do {
			if (aleft > s->crecs[id2].apos)
			    aleft = s->crecs[id2].apos, left_cnt = 1;
			else if (aleft == s->crecs[id2].apos)
			    left_cnt++;
			if (aright < s->crecs[id2].aend)
			    aright = s->crecs[id2].aend;
			if (s->crecs[id2].mate_line == -1) {
			    s->crecs[id2].mate_line = rec;
			    break;
			}
			if (s->crecs[id2].mate_line <= id2 ||
                            s->crecs[id2].mate_line >= s->hdr->num_records)
			    return -1;
			id2 = s->crecs[id2].mate_line;

			if (s->crecs[id2].ref_id != ref)
			    ref = -1;
		    } while (id2 != id1);

		    if (ref != -1) {
			tlen = aright - aleft + 1;
			id1 = id2 = rec;

			/*
			 * When we have two seqs with identical start and
			 * end coordinates, set +/- tlen based on 1st/last
			 * bit flags instead, as a tie breaker.
			 */
			if (s->crecs[id2].apos == aleft) {
			    if (left_cnt == 1 || 
				(s->crecs[id2].flags & BAM_FREAD1))
				s->crecs[id2].tlen = tlen;
			    else
				s->crecs[id2].tlen = -tlen;
			} else {
			    s->crecs[id2].tlen = -tlen;
			}

			id2 = s->crecs[id2].mate_line;
			while (id2 != id1) {
			    if (s->crecs[id2].apos == aleft) {
				if (left_cnt == 1 || 
				    (s->crecs[id2].flags & BAM_FREAD1))
				    s->crecs[id2].tlen = tlen;
				else
				    s->crecs[id2].tlen = -tlen;
			    } else {
				s->crecs[id2].tlen = -tlen;
			    }
			    id2 = s->crecs[id2].mate_line;
			}
		    } else {
			id1 = id2 = rec;

			s->crecs[id2].tlen = 0;
			id2 = s->crecs[id2].mate_line;
			while (id2 != id1) {
			    s->crecs[id2].tlen = 0;
			    id2 = s->crecs[id2].mate_line;
			}
		    }
		}

		cr->mate_pos = s->crecs[cr->mate_line].apos;
		cr->mate_ref_id = s->crecs[cr->mate_line].ref_id;

		// paired
		cr->flags |= BAM_FPAIRED;

		// set mate unmapped if needed
		if (s->crecs[cr->mate_line].flags & BAM_FUNMAP) {
		    cr->flags |= BAM_FMUNMAP;
		    cr->tlen = 0;
		}
		if (cr->flags & BAM_FUNMAP) {
		    cr->tlen = 0;
		}

		// set mate reversed if needed
		if (s->crecs[cr->mate_line].flags & BAM_FREVERSE)
		    cr->flags |= BAM_FMREVERSE;
	    } else {
		fprintf(stderr, "Mate line out of bounds: %d vs [0, %d]\n",
			cr->mate_line, s->hdr->num_records-1);
	    }

	    /* FIXME: construct read names here too if needed */
	} else {
	    if (cr->mate_flags & CRAM_M_REVERSE) {
		cr->flags |= BAM_FPAIRED | BAM_FMREVERSE;
	    }
	    if (cr->mate_flags & CRAM_M_UNMAP) {
		cr->flags |= BAM_FMUNMAP;
		//cr->mate_ref_id = -1;
	    }
	    if (!(cr->flags & BAM_FPAIRED))
		cr->mate_ref_id = -1;
	}

	if (cr->tlen == INT64_MIN)
	    cr->tlen = 0; // Just incase
    }

    for (rec = 0; rec < s->hdr->num_records; rec++) {
	cram_record *cr = &s->crecs[rec];
	if (cr->explicit_tlen != INT64_MIN)
	    cr->tlen = cr->explicit_tlen;
    }

    return 0;
}

static char *md5_print(unsigned char *md5, char *out) {
    int i;
    for (i = 0; i < 16; i++) {
	out[i*2+0] = "0123456789abcdef"[md5[i]>>4];
	out[i*2+1] = "0123456789abcdef"[md5[i]&15];
    }
    out[32] = 0;

    return out;
}

static int cram_to_bam(SAM_hdr *bfd, cram_fd *fd, cram_slice *s,
		       cram_record *cr, int rec, bam_seq_t **bam);

/*
 * Bulk conversion of an entire cram slice to an array of bam objects.
 * (Assumption that fd->required_fields will not change from one
 * cram_get_bam_seq() call to the  next.)
 *
 * Returns 0 on success
 *        -1 on failure
 */

static int bam_size(SAM_hdr *bfd, cram_fd *fd, cram_record *cr) {
    int name_len, rg_len;

    if (fd->required_fields & SAM_QNAME) {
	if (cr->name_len)
	    name_len = cr->name_len;
	else
	    name_len = strlen(fd->prefix) + 20; // overestimate
    } else {
	name_len = 1;
    }

    rg_len = (cr->rg != -1) ? bfd->rg[cr->rg].name_len + 4 : 0;

    return sizeof(bam_seq_t)
	+ name_len + 1
	+ round4(name_len+1)
	+ 4 * cr->ncigar
	+ (cr->len+1)/2
	+ cr->len
	+ cr->aux_size + rg_len + 1;
}

/*
 * Returns 0 on success,
 *        -1 on failure
 */
static int bulk_cram_to_bam(SAM_hdr *bfd, cram_fd *fd, cram_slice *s) {
    int i;
    int r = 0;
    int sizes[10000];

    size_t len = 0;
    for (i = 0; i < s->hdr->num_records; i++) {
	int sz = bam_size(bfd, fd, &s->crecs[i]);
	if (i < 10000)
	    sizes[i] = sz;
	len += round8(sz);
    }

    s->bl = (bam_seq_t **)malloc(s->hdr->num_records * sizeof(*s->bl) + len + 8);
    if (!s->bl)
	return -1;

    // Round up to next multiple of 8, to ensure bam structs are 8-byte aligned.
    char *x = ((char *)s->bl) + round8(s->hdr->num_records * sizeof(*s->bl));
    for (i = 0; i < s->hdr->num_records; i++) {
	bam_seq_t *b = (bam_seq_t *)x, *o = b;
	int bsize = i < 10000 ? sizes[i] : bam_size(bfd, fd, &s->crecs[i]);
	b->alloc = bsize;
	r |= (cram_to_bam(fd->header, fd, s, &s->crecs[i], i, &b) < 0);
	// if we allocated enough, the above won't have resized b
	assert(o == b && o->alloc == bsize);
	x += round8(bsize);
	s->bl[i] = b;
    }

    return r?-1:0;
}

/*
 * Utility function to decode tlen (ISIZE), as it's called
 * in multiple places.
 *
 * Returns codec return value (0 on success).
 */
static int cram_decode_tlen(cram_fd *fd, cram_container *c, cram_slice *s,
			    cram_block *blk, int64_t *tlen) {
    int out_sz = 1, r = 0;

    if (!c->comp_hdr->codecs[DS_TS]) return -1;
    if (CRAM_MAJOR_VERS(fd->version) < 4) {
	int32_t i32;
	r |= c->comp_hdr->codecs[DS_TS]
	    ->decode(s, c->comp_hdr->codecs[DS_TS], blk,
		     (char *)&i32, &out_sz);
	*tlen = i32;
    } else {
	r |= c->comp_hdr->codecs[DS_TS]
	    ->decode(s, c->comp_hdr->codecs[DS_TS], blk,
		     (char *)tlen, &out_sz);
    }
    return r;
}

/*
 * Decode an entire slice from container blocks. Fills out s->crecs[] array.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_decode_slice(cram_fd *fd, cram_container *c, cram_slice *s,
		      SAM_hdr *bfd) {
    cram_block *blk = s->block[0];
    int32_t bf, ref_id;
    unsigned char cf;
    int out_sz, r = 0;
    int rec;
    char *seq = NULL, *qual = NULL;
    int unknown_rg = -1;
    int embed_ref;
    char **refs = NULL;
    uint32_t ds;

    if (cram_dependent_data_series(fd, c->comp_hdr, s) != 0)
	return -1;

    ds = s->data_series;
    //printf("%08x\n", ds);

    blk->bit = 7; // MSB first

    // Study the blocks and estimate approx sizes to preallocate.
    // This looks to speed up decoding by around 8-9%.
    // We can always shrink back down at the end if we overestimated.
    // However it's likely that this also saves memory as own growth
    // factor (*=1.5) is never applied.
    {
	int qsize, nsize, q_id;
	qsize = c->num_bases / c->max_slice;
	nsize = c->num_records / c->max_slice * 50;
	cram_decode_estimate_sizes(c->comp_hdr, s, &qsize, &nsize, &q_id);

	if (qsize && (ds & CRAM_RL)) BLOCK_RESIZE_EXACT(s->seqs_blk, qsize+1);
	if (qsize && (ds & CRAM_RL)) BLOCK_RESIZE_EXACT(s->qual_blk, qsize+1);
	if (nsize && (ds & CRAM_NS)) BLOCK_RESIZE_EXACT(s->name_blk, nsize+1);

	// To do - consider using q_id here to usurp the quality block and
	// avoid a memcpy during decode.
	// Specifically when quality is an external block uniquely used by
	// DS_QS only, then we can set s->qual_blk directly to this
	// block and save the codec->decode() calls. (Approx 3% cpu saving)
    }

    /* Look for unknown RG, added as last by Java CRAM? */
    if (bfd->nrg > 0 &&
	bfd->rg[bfd->nrg-1].name != NULL &&
	!strcmp(bfd->rg[bfd->nrg-1].name, "UNKNOWN"))
	unknown_rg = bfd->nrg-1;

    if (blk->content_type != CORE)
	return -1;

    if (s->crecs)
	free(s->crecs);
    if (!(s->crecs = malloc(s->hdr->num_records * sizeof(*s->crecs))))
	return -1;

    ref_id = s->hdr->ref_seq_id;
    if (CRAM_MAJOR_VERS(fd->version) < 4)
	embed_ref = s->hdr->ref_base_id >= 0 ? 1 : 0;
    else
	embed_ref = s->hdr->ref_base_id > 0 ? 1 : 0;

    // Embedded references, or more specifically embedded consensuses, mean
    // that we may be reading, modifying (filter, subsample, etc), writing
    // and the written file may have a different embedded reference (cons)
    // to the read file.  In order to know whether or not this changes the
    // MD tag, we must first know what the MD tag is.
    //
    // Until we have a way of knowing whether the MD tag was originally in
    // the file (this is a CRAM flaw), we have to assume it was and was
    // silently dropped.  The only work around is to forcibly enable the
    // generation of MD tags in this scenario.  Ugly, and something to
    // fix for CRAM v4.0.
    if (embed_ref)
	s->decode_md = 1;

    if (ref_id >= 0) {
	if (embed_ref) {
	    cram_block *b;
	    if (s->hdr->ref_base_id < 0) {
		fprintf(stderr, "No reference specified and "
			"no embedded reference is available.\n");
		return -1;
	    }
	    b = cram_get_block_by_id(s, s->hdr->ref_base_id);
	    if (!b)
		return -1;
            if (cram_uncompress_block(b) != 0)
                return -1;
	    s->ref = (char *)BLOCK_DATA(b);
	    s->ref_start = s->hdr->ref_seq_start;
	    s->ref_end   = s->hdr->ref_seq_start + s->hdr->ref_seq_span-1;
	    if (s->ref_end - s->ref_start > b->uncomp_size) {
		fprintf(stderr, "Embedded reference is too small.\n");
		return -1;
	    }
	} else if (!c->comp_hdr->no_ref) {
	    //// Avoid Java cramtools bug by loading entire reference seq 
	    //s->ref = cram_get_ref(fd, s->hdr->ref_seq_id, 1, 0);
	    //s->ref_start = 1;

	    if (fd->required_fields & SAM_SEQ)
		s->ref =
		cram_get_ref(fd, s->hdr->ref_seq_id,
			     s->hdr->ref_seq_start,
			     s->hdr->ref_seq_start + s->hdr->ref_seq_span -1);
	    s->ref_start = s->hdr->ref_seq_start;
	    s->ref_end   = s->hdr->ref_seq_start + s->hdr->ref_seq_span-1;

	    /* Sanity check */
	    if (s->ref_start < 0) {
		fprintf(stderr, "Slice starts before base 1.\n");
		s->ref_start = 0;
	    }
	    if (fd->ref_lock) pthread_mutex_lock(fd->ref_lock);
	    pthread_mutex_lock(&fd->refs->lock);
	    if ((fd->required_fields & SAM_SEQ) &&
		s->ref_end > fd->refs->ref_id[ref_id]->length) {
		s->ref_end = fd->refs->ref_id[ref_id]->length;
	    }
	    pthread_mutex_unlock(&fd->refs->lock);
	    if (fd->ref_lock) pthread_mutex_unlock(fd->ref_lock);
	}
    }

    if ((fd->required_fields & SAM_SEQ) &&
	s->ref == NULL && s->hdr->ref_seq_id >= 0 && !c->comp_hdr->no_ref) {
	fprintf(stderr, "Unable to fetch reference #%d %"PRId64"..%"PRId64"\n",
		s->hdr->ref_seq_id, s->hdr->ref_seq_start,
		s->hdr->ref_seq_start + s->hdr->ref_seq_span-1);
	return -1;
    }

    if (!IS_CRAM_1_VERS(fd)
	&& (fd->required_fields & SAM_SEQ)
	&& s->hdr->ref_seq_id >= 0
	&& !fd->ignore_md5
	&& memcmp(s->hdr->md5, "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", 16)) {
	MD5_CTX md5;
	unsigned char digest[16];

	if (s->ref && s->hdr->ref_seq_id >= 0) {
	    int start, len;

	    if (s->hdr->ref_seq_start >= s->ref_start) {
		start = s->hdr->ref_seq_start - s->ref_start;
	    } else {
		fprintf(stderr, "Slice starts before base 1.\n");
		start = 0;
	    }

	    if (s->hdr->ref_seq_span <= s->ref_end - s->ref_start + 1) {
		len = s->hdr->ref_seq_span;
	    } else {
		fprintf(stderr, "Slice ends beyond reference end.\n");
		len = s->ref_end - s->ref_start + 1;
	    }

	    MD5_Init(&md5);
	    if (start + len > s->ref_end - s->ref_start + 1)
		len = s->ref_end - s->ref_start + 1 - start;
	    if (len >= 0)
		MD5_Update(&md5, s->ref + start, len);
	    MD5_Final(digest, &md5);
	} else if (!s->ref && s->hdr->ref_base_id >= 0) {
	    cram_block *b = cram_get_block_by_id(s, s->hdr->ref_base_id);
	    if (b) {
		MD5_Init(&md5);
		MD5_Update(&md5, b->data, b->uncomp_size);
		MD5_Final(digest, &md5);
	    }
	}

	if ((!s->ref && s->hdr->ref_base_id < 0)
	    || memcmp(digest, s->hdr->md5, 16) != 0) {
	    char M[33];
	    fprintf(stderr, "ERROR: md5sum reference mismatch for ref "
		    "%d pos %d..%d\n", ref_id, s->ref_start, s->ref_end);
	    fprintf(stderr, "CRAM: %s\n", md5_print(s->hdr->md5, M));
	    fprintf(stderr, "Ref : %s\n", md5_print(digest, M));
	    return -1;
	}
    }

    if (ref_id == -2) {
	if (fd->ref_lock) pthread_mutex_lock(fd->ref_lock);
	pthread_mutex_lock(&fd->refs->lock);
	refs = calloc(fd->refs->nref, sizeof(char *));
	pthread_mutex_unlock(&fd->refs->lock);
	if (fd->ref_lock) pthread_mutex_unlock(fd->ref_lock);
	if (!refs)
	    return -1;
    }

#define RETURN return printf("Fail to decode CRAM rec %d at %s:%d\n", rec, __FILE__, __LINE__),
    int last_ref_id = -9; // Arbitrary -ve marker for not-yet-set
    for (rec = 0; rec < s->hdr->num_records; rec++) {
	cram_record *cr = &s->crecs[rec];
	int has_MD, has_NM;

	//fprintf(stderr, "Decode seq %d, %d/%d\n", rec, blk->byte, blk->bit);

	cr->s = s;

	out_sz = 1; /* decode 1 item */
	if (ds & CRAM_BF) {
	    if (!c->comp_hdr->codecs[DS_BF]) RETURN -1;
	    r |= c->comp_hdr->codecs[DS_BF]
		            ->decode(s, c->comp_hdr->codecs[DS_BF], blk,
				     (char *)&bf, &out_sz);
	    if (r || bf < 0 ||
		bf >= sizeof(fd->bam_flag_swap)/sizeof(*fd->bam_flag_swap))
		RETURN -1;
	    bf = fd->bam_flag_swap[bf];
	    cr->flags = bf;
	} else {
	    cr->flags = bf = 0x4; // unmapped
	}

	if (ds & CRAM_CF) {
	    if (IS_CRAM_1_VERS(fd)) {
		/* CF is byte in 1.0, int32 in 2.0 */
		if (!c->comp_hdr->codecs[DS_CF]) RETURN -1;
		r |= c->comp_hdr->codecs[DS_CF]
		                ->decode(s, c->comp_hdr->codecs[DS_CF], blk,
				 	 (char *)&cf, &out_sz);
		if (r) RETURN -1;
		cr->cram_flags = cf;
	    } else {
		if (!c->comp_hdr->codecs[DS_CF]) RETURN -1;
		r |= c->comp_hdr->codecs[DS_CF]
		                ->decode(s, c->comp_hdr->codecs[DS_CF], blk,
					 (char *)&cr->cram_flags, &out_sz);
		if (r) RETURN -1;
		cf = cr->cram_flags;
	    }
	} else {
	    cf = cr->cram_flags = 0;
	}

	if (!IS_CRAM_1_VERS(fd) && ref_id == -2) {
	    if (ds & CRAM_RI) {
		if (!c->comp_hdr->codecs[DS_RI]) RETURN -1;
		r |= c->comp_hdr->codecs[DS_RI]
		                ->decode(s, c->comp_hdr->codecs[DS_RI], blk,
					 (char *)&cr->ref_id, &out_sz);
		if (r) RETURN -1;
		if ((fd->required_fields & (SAM_SEQ|SAM_TLEN))
		    && cr->ref_id >= 0
		    && cr->ref_id != last_ref_id) {
		    if (!c->comp_hdr->no_ref) {
			if (!refs[cr->ref_id])
			    refs[cr->ref_id] = cram_get_ref(fd, cr->ref_id,
							    1, 0);
			s->ref = refs[cr->ref_id];

			if (!fd->unsorted && last_ref_id >= 0 && refs[last_ref_id]) {
			    cram_ref_decr(fd->refs, last_ref_id);
			    refs[last_ref_id] = NULL;
			}
		    }
		    s->ref_start = 1;
		    if (fd->ref_lock) pthread_mutex_lock(fd->ref_lock);
		    pthread_mutex_lock(&fd->refs->lock);
		    s->ref_end = fd->refs->ref_id[cr->ref_id]->length;
		    pthread_mutex_unlock(&fd->refs->lock);
		    if (fd->ref_lock) pthread_mutex_unlock(fd->ref_lock);

		    last_ref_id = cr->ref_id;
		}
	    } else {
		cr->ref_id = -1;
	    }
	} else {
	    cr->ref_id = ref_id; // Forced constant in CRAM 1.0
	}
	if (cr->ref_id >= bfd->nref) {
	    fprintf(stderr, "Requested unknown reference ID %d\n", cr->ref_id);
            RETURN -1;
	}

	if (ds & CRAM_RL) {
	    if (!c->comp_hdr->codecs[DS_RL]) RETURN -1;
	    r |= c->comp_hdr->codecs[DS_RL]
		            ->decode(s, c->comp_hdr->codecs[DS_RL], blk,
				     (char *)&cr->len, &out_sz);
	    if (r) RETURN r;
	    if (cr->len < 0) {
	        fprintf(stderr, "Read has negative length\n");
		RETURN -1;
	    }
	}

	if (ds & CRAM_AP) {
	    if (!c->comp_hdr->codecs[DS_AP]) RETURN -1;
	    if (CRAM_MAJOR_VERS(fd->version) < 4) {
		int32_t i32;
		r |= c->comp_hdr->codecs[DS_AP]
		                ->decode(s, c->comp_hdr->codecs[DS_AP], blk,
					 (char *)&i32, &out_sz);
		cr->apos = i32;
	    } else {
		r |= c->comp_hdr->codecs[DS_AP]
		                ->decode(s, c->comp_hdr->codecs[DS_AP], blk,
					 (char *)&cr->apos, &out_sz);
	    }
	    if (r) RETURN r;
	    if (c->comp_hdr->AP_delta)
		cr->apos += s->last_apos;
	    s->last_apos=  cr->apos;
	} else {
	    cr->apos = c->ref_seq_start;
	}
		    
	if (ds & CRAM_RG) {
	    if (!c->comp_hdr->codecs[DS_RG]) RETURN -1;
	    r |= c->comp_hdr->codecs[DS_RG]
		           ->decode(s, c->comp_hdr->codecs[DS_RG], blk,
				    (char *)&cr->rg, &out_sz);
	    if (r) RETURN r;
	    if (cr->rg == unknown_rg)
		cr->rg = -1;
	} else {
	    cr->rg = -1;
	}

	cr->name_len = 0;

	if (c->comp_hdr->read_names_included) {
	    int32_t out_sz2 = 1;

	    // Read directly into name cram_block
	    cr->name = BLOCK_SIZE(s->name_blk);
	    if (ds & CRAM_RN) {
		if (!c->comp_hdr->codecs[DS_RN]) RETURN -1;
		r |= c->comp_hdr->codecs[DS_RN]
		                ->decode(s, c->comp_hdr->codecs[DS_RN], blk,
					 (char *)s->name_blk, &out_sz2);
		if (r) RETURN r;
		cr->name_len = out_sz2;
	    }
	}

	cr->mate_pos = 0;
	cr->mate_line = -1;
	cr->mate_ref_id = -1;
	cr->explicit_tlen = INT64_MIN;
	if ((ds & CRAM_CF) && (cf & CRAM_FLAG_DETACHED)) {
	    if (ds & CRAM_MF) {
		if (IS_CRAM_1_VERS(fd)) {
		    /* MF is byte in 1.0, int32 in 2.0 */
		    unsigned char mf;
		    if (!c->comp_hdr->codecs[DS_MF]) RETURN -1;
		    r |= c->comp_hdr->codecs[DS_MF]
			            ->decode(s, c->comp_hdr->codecs[DS_MF],
					     blk, (char *)&mf, &out_sz);
		    if (r) RETURN r;
		    cr->mate_flags = mf;
		} else {
		    if (!c->comp_hdr->codecs[DS_MF]) RETURN -1;
		    r |= c->comp_hdr->codecs[DS_MF]
			            ->decode(s, c->comp_hdr->codecs[DS_MF],
					     blk,
					     (char *)&cr->mate_flags,
					     &out_sz);
		    if (r) RETURN r;
		}
	    } else {
		cr->mate_flags = 0;
	    }

	    if (!c->comp_hdr->read_names_included) {
		int32_t out_sz2 = 1;
	    
		// Read directly into name cram_block
		cr->name = BLOCK_SIZE(s->name_blk);
		if (ds & CRAM_RN) {
		    if (!c->comp_hdr->codecs[DS_RN]) RETURN -1;
		    r |= c->comp_hdr->codecs[DS_RN]
			            ->decode(s, c->comp_hdr->codecs[DS_RN],
					     blk, (char *)s->name_blk,
					     &out_sz2);
		    if (r) RETURN r;
		    cr->name_len = out_sz2;
		}
	    }
		    
	    if (ds & CRAM_NS) {
		if (!c->comp_hdr->codecs[DS_NS]) RETURN -1;
		r |= c->comp_hdr->codecs[DS_NS]
		                ->decode(s, c->comp_hdr->codecs[DS_NS], blk,
					 (char *)&cr->mate_ref_id, &out_sz);
		if (r) RETURN r;
	    }

// Skip as mate_ref of "*" is legit. It doesn't mean unmapped, just unknown.
//	    if (cr->mate_ref_id == -1 && cr->flags & 0x01) {
//		/* Paired, but unmapped */
//		cr->flags |= BAM_FMUNMAP;
//	    }

	    if (ds & CRAM_NP) {
		if (!c->comp_hdr->codecs[DS_NP]) RETURN -1;
		if (CRAM_MAJOR_VERS(fd->version) < 4) {
		    int32_t i32;
		    r |= c->comp_hdr->codecs[DS_NP]
			            ->decode(s, c->comp_hdr->codecs[DS_NP], blk,
					     (char *)&i32, &out_sz);
		    cr->mate_pos = i32;
		} else {
		    r |= c->comp_hdr->codecs[DS_NP]
			->decode(s, c->comp_hdr->codecs[DS_NP], blk,
				 (char *)&cr->mate_pos, &out_sz);
		}
		if (r) RETURN r;
	    }

	    if (ds & CRAM_TS) {
		r = cram_decode_tlen(fd, c, s, blk, &cr->tlen);
		if (r) RETURN r;
	    } else {
		cr->tlen = INT64_MIN;
	    }

	} else if ((ds & CRAM_CF) && (cf & CRAM_FLAG_MATE_DOWNSTREAM)) {
	    // else not detached
	    if (ds & CRAM_NF) {
		if (!c->comp_hdr->codecs[DS_NF]) RETURN -1;
		r |= c->comp_hdr->codecs[DS_NF]
		                ->decode(s, c->comp_hdr->codecs[DS_NF], blk,
					 (char *)&cr->mate_line, &out_sz);
		if (r) RETURN r;
		cr->mate_line += rec + 1;

		//cr->name_len = sprintf(name, "%d", name_id++);
		//cr->name = DSTRING_LEN(name_ds);
		//dstring_nappend(name_ds, name, cr->name_len);

		cr->mate_ref_id = -1;
		cr->tlen = INT64_MIN;
		cr->mate_pos = 0;
	    } else  {
		cr->mate_flags = 0;
		cr->tlen = INT64_MIN;
	    }
	    if ((ds & CRAM_CF) && (cf & CRAM_FLAG_EXPLICIT_TLEN)) {
		if (ds & CRAM_TS) {
		    r = cram_decode_tlen(fd, c, s, blk, &cr->explicit_tlen);
		    if (r) RETURN r;
		} else {
		    cr->mate_flags = 0;
		    cr->tlen = INT64_MIN;
		}
	    }
	} else if ((ds & CRAM_CF) && (cf & CRAM_FLAG_EXPLICIT_TLEN)) {
	    if (ds & CRAM_TS) {
		r = cram_decode_tlen(fd, c, s, blk, &cr->explicit_tlen);
		if (r) RETURN r;
	    } else {
		cr->mate_flags = 0;
		cr->tlen = INT64_MIN;
	    }
	} else {
	    cr->mate_flags = 0;
	    cr->tlen = INT64_MIN;
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
	has_MD = has_NM = 0;
	if (IS_CRAM_1_VERS(fd))
	    r |= cram_decode_aux_1_0(c, s, blk, cr);
	else
	    r |= cram_decode_aux(fd, c, s, blk, cr, &has_MD, &has_NM);
	if (r) RETURN r;

	/* Fake up dynamic string growth and appending */
	if (ds & CRAM_RL) {
	    cr->seq = BLOCK_SIZE(s->seqs_blk);
	    BLOCK_GROW(s->seqs_blk, cr->len);
	    seq = (char *)BLOCK_END(s->seqs_blk);
	    BLOCK_SIZE(s->seqs_blk) += cr->len;

	    if (!seq)
		RETURN -1;
	
	    cr->qual = BLOCK_SIZE(s->qual_blk);
	    BLOCK_GROW(s->qual_blk, cr->len);
	    qual = (char *)BLOCK_END(s->qual_blk);
	    BLOCK_SIZE(s->qual_blk) += cr->len;

	    if (!s->ref)
		memset(seq, '=', cr->len);
	}

	if (!(bf & BAM_FUNMAP)) {
            if ((ds & CRAM_AP) && cr->apos <= 0) {
                fprintf(stderr,
			"Read has alignment position %"PRId64
			" but no unmapped flag\n",
			cr->apos);
		RETURN -1;
	    }
	    /* Decode sequence and generate CIGAR */
	    if (ds & (CRAM_SEQ | CRAM_MQ)) {
		r |= cram_decode_seq(fd, c, s, blk, cr, bfd, cf, seq, qual,
				     has_MD, has_NM);
		if (r) RETURN r;
	    } else {
		cr->cigar = 0;
		cr->ncigar = 0;
		cr->aend = cr->apos;
		cr->mqual = 0;
	    }
	} else {
	    int out_sz2 = cr->len;

	    //puts("Unmapped");
	    cr->cigar = 0;
	    cr->ncigar = 0;
	    cr->aend = cr->apos;
	    cr->mqual = 0;

	    if (ds & CRAM_BA && cr->len) {
		if (!c->comp_hdr->codecs[DS_BA]) RETURN -1;
		r |= c->comp_hdr->codecs[DS_BA]
		                ->decode(s, c->comp_hdr->codecs[DS_BA], blk,
					 (char *)seq, &out_sz2);
		if (r) RETURN r;
	    }

	    if ((ds & CRAM_CF) && (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES)) {
		out_sz2 = cr->len;
		if (ds & CRAM_QS && cr->len >= 0) {
		    if (!c->comp_hdr->codecs[DS_QS]) RETURN -1;
		    r |= c->comp_hdr->codecs[DS_QS]
			            ->decode(s, c->comp_hdr->codecs[DS_QS],
					     blk, qual, &out_sz2);
		    if (r) RETURN r;
		}
	    } else {
		if (ds & CRAM_RL)
		    memset(qual, 255, cr->len);
	    }
	}

	if (!c->comp_hdr->qs_seq_orient && (ds & CRAM_QS) && (cr->flags & BAM_FREVERSE)) {
	    int i, j;
	    for (i = 0, j = cr->len-1; i < j; i++, j--) {
		unsigned char c;
		c = qual[i];
		qual[i] = qual[j];
		qual[j] = c;
	    }
	}

	if (!fd->ignore_chksum) {
	    if (s->hdr->BD_crc && ds & CRAM_BA && s->ref)
		s->BD_crc += iolib_crc32(0L, (Bytef *) seq, cr->len);
	    
	    if (s->hdr->SD_crc &&
		(ds & CRAM_QS) &&
		(cf & CRAM_FLAG_PRESERVE_QUAL_SCORES)) {
		s->SD_crc += iolib_crc32(0L, (Bytef *) qual, cr->len);
	    }
	}
    }

    if (fd->ref_lock) pthread_mutex_lock(fd->ref_lock);
    if (refs) {
	int i;
	for (i = 0; i < fd->refs->nref; i++) {
	    if (refs[i])
		cram_ref_decr(fd->refs, i);
	}
	free(refs);
    } else if (ref_id >= 0 && s->ref != fd->ref_free && !embed_ref) {
	cram_ref_decr(fd->refs, ref_id);
    }
    if (fd->ref_lock) pthread_mutex_unlock(fd->ref_lock);

    /* Resolve mate pair cross-references between recs within this slice */
    r |= cram_decode_slice_xref(s, fd->required_fields);

    // Free the original blocks as we no longer need these.
    {
	int i;
	for (i = 0; i < s->hdr->num_blocks; i++) {
	    cram_block *b = s->block[i];
	    cram_free_block(b);
	    s->block[i] = NULL;
	}
    }

    // Also see initial BLOCK_RESIZE_EXACT at top of function.
    // As we grow blocks we overallocate by up to 50%. So shrink
    // back to their final sizes here.
    //
//    fprintf(stderr, "%d %d // %d %d // %d %d // %d %d\n",
//    	    (int)s->seqs_blk->byte, (int)s->seqs_blk->alloc, 
//    	    (int)s->qual_blk->byte, (int)s->qual_blk->alloc, 
//    	    (int)s->name_blk->byte, (int)s->name_blk->alloc, 
//    	    (int)s->aux_blk->byte,  (int)s->aux_blk->alloc);
    BLOCK_RESIZE_EXACT(s->seqs_blk, BLOCK_SIZE(s->seqs_blk)+1);
    BLOCK_RESIZE_EXACT(s->qual_blk, BLOCK_SIZE(s->qual_blk)+1);
    BLOCK_RESIZE_EXACT(s->name_blk, BLOCK_SIZE(s->name_blk)+1);
    BLOCK_RESIZE_EXACT(s->aux_blk,  BLOCK_SIZE(s->aux_blk)+1);

    /* Checksum */
    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
	if (ds & CRAM_BA && s->BD_crc && s->hdr->BD_crc) {
	    if (s->BD_crc != s->hdr->BD_crc) {
		fprintf(stderr, "BD checksum failure: %08x vs %08x\n",
			(uint32_t)s->BD_crc, (uint32_t)s->hdr->BD_crc);
		r |= -1;
	    }
	}

	if (ds & CRAM_QS && s->SD_crc && s->hdr->SD_crc) {
	    if (s->SD_crc != s->hdr->SD_crc) {
		fprintf(stderr, "SD checksum failure: %08x vs %08x\n",
			(uint32_t)s->SD_crc, (uint32_t)s->hdr->SD_crc);
		r |= -1;
	    }
	}
    }

    // If we're wanting BAM records, convert these up-front too.
    // This is useful when we're streaming lots of data in a
    // multi-threaded environment as the cram to bam conversion is
    // then threaded too.
    //
    // Possible future optimisation - check range query and don't
    // convert all reads to BAM.

    if (fd->pool)
	r |= bulk_cram_to_bam(bfd, fd, s);

    return r;
}

typedef struct {
    cram_fd *fd;
    cram_container *c;
    cram_slice *s;
    SAM_hdr *h;
    int exit_code;
} cram_decode_job;

void *cram_decode_slice_thread(void *arg) {
    cram_decode_job *j = (cram_decode_job *)arg;

    j->exit_code = cram_decode_slice(j->fd, j->c, j->s, j->h);

    return j;
}

/*
 * Spawn a multi-threaded version of cram_decode_slice().
 */
int cram_decode_slice_mt(cram_fd *fd, cram_container *c, cram_slice *s,
			 SAM_hdr *bfd) {
    cram_decode_job *j;
    int nonblock;

    if (!fd->pool)
	return cram_decode_slice(fd, c, s, bfd);

    if (!(j = malloc(sizeof(*j))))
	return -1;

    j->fd = fd;
    j->c  = c;
    j->s  = s;
    j->h  = bfd;
    
    nonblock = t_pool_results_queue_sz(fd->rqueue) ? 1 : 0;

    if (-1 == t_pool_dispatch2(fd->pool, fd->rqueue, cram_decode_slice_thread,
			       j, nonblock)) {
	/* Would block */
	fd->job_pending = j;
    } else {
	fd->job_pending = NULL;
    }

    // flush too
    return 0;
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
 * (*bam)->alloc accordingly, allowing it to be used within a loop
 * efficiently without needing to allocate new bam objects over and
 * over again.
 *
 * Returns the used size of the bam record on success
 *         -1 on failure.
 */
static int cram_to_bam(SAM_hdr *bfd, cram_fd *fd, cram_slice *s,
		       cram_record *cr, int rec, bam_seq_t **bam) {
    int bam_idx, rg_len;
    char name_a[1024], *name;
    int name_len;
    char *aux, *aux_orig;
    char *seq, *qual;

    /* Assign names if not explicitly set */
    if (fd->required_fields & SAM_QNAME) {
	if (cr->name_len) {
	    name = (char *)BLOCK_DATA(s->name_blk) + cr->name;
	    name_len = cr->name_len;
	} else {
	    name = name_a;
	    if (cr->mate_line >= 0 && cr->mate_line < s->max_rec &&
		s->crecs[cr->mate_line].name_len > 0) {
		// Copy our mate if non-zero.
		memcpy(name_a, BLOCK_DATA(s->name_blk)+s->crecs[cr->mate_line].name,
		       s->crecs[cr->mate_line].name_len);
		name = name_a + s->crecs[cr->mate_line].name_len;
	    } else {
		// Otherwise generate a name based on prefix
		name_len = strlen(fd->prefix);
		memcpy(name, fd->prefix, name_len);
		name += name_len;
		*name++ = ':';
		if (cr->mate_line >= 0 && cr->mate_line < rec) {
		    name = (char *)append_uint64((unsigned char *)name,
						 s->hdr->record_counter +
						 cr->mate_line + 1);
		} else {
		    name = (char *)append_uint64((unsigned char *)name,
						 s->hdr->record_counter +
						 rec + 1);
		}
	    }
	    name_len = name - name_a;
	    name = name_a;
	}
    } else {
	name = "?";
	name_len = 1;
    }

    /* Generate BAM record */
    if (cr->rg < -1 || cr->rg >= bfd->nrg)
	return -1;
    rg_len = (cr->rg != -1) ? bfd->rg[cr->rg].name_len + 4 : 0;

    if (fd->required_fields & (SAM_SEQ | SAM_QUAL)) {
	if (!BLOCK_DATA(s->seqs_blk))
	    return -1;
	seq = (char *)BLOCK_DATA(s->seqs_blk) + cr->seq;
    } else {
	seq = "*";
	cr->len = 0;
    }

    if (fd->required_fields & SAM_QUAL) {
	if (!BLOCK_DATA(s->qual_blk))
	    return -1;
	qual = (char *)BLOCK_DATA(s->qual_blk) + cr->qual;
    } else {
	qual = NULL;
    }

    bam_idx = bam_construct_seq(bam, cr->aux_size + rg_len,
				name, name_len,
				cr->flags,
				cr->ref_id,
				cr->apos,
				cr->aend,
				cr->mqual,
				cr->ncigar, &s->cigar[cr->cigar],
				cr->mate_ref_id,
				cr->mate_pos,
				cr->tlen,
				cr->len,
				seq,
				qual);
    if (bam_idx == -1)
	return -1;

    aux = aux_orig = (char *)bam_aux(*bam);

    /* Auxiliary strings */
    if (cr->aux_size != 0) {
	memcpy(aux, BLOCK_DATA(s->aux_blk) + cr->aux, cr->aux_size);
	aux += cr->aux_size;
    }

    /* RG:Z: */
    if (cr->rg != -1) {
	int len = bfd->rg[cr->rg].name_len;
	*aux++ = 'R'; *aux++ = 'G'; *aux++ = 'Z';
	memcpy(aux, bfd->rg[cr->rg].name, len);
	aux += len;
	*aux++ = 0;
    }
    
#ifndef SAMTOOLS
    bam_set_blk_size(*bam, bam_blk_size(*bam) + (aux - aux_orig));
#endif

    *aux++ = 0;

    return bam_idx + (aux - aux_orig);
}

/*
 * Here be dragons! The multi-threading code in this is crufty beyond belief.
 */
/*
 * Load first container.
 * Called when fd->ctr is NULL>
 *
 * Returns container on success
 *        NULL on failure.
 */
static cram_container *cram_first_slice(cram_fd *fd) {
    cram_container *c;

    do {
	if (!(c = fd->ctr = cram_read_container(fd)))
	    return NULL;
	c->curr_slice_mt = c->curr_slice;
    } while (c->length == 0);

    /*
     * The first container may be a result of a sub-range query.
     * In which case it may still not be the optimal starting point
     * due to skipped containers/slices in the index. 
     */
    // No need for locks here as we're in the main thread.
    if (fd->range.refid != -2) {
	while (c->ref_seq_id != -2 &&
	       (c->ref_seq_id < fd->range.refid ||
		(fd->range.refid >= 0 && c->ref_seq_id == fd->range.refid
		 && c->ref_seq_start + c->ref_seq_span-1 < fd->range.start))) {
	    if (0 != cram_seek(fd, c->length, SEEK_CUR))
		return NULL;
	    cram_free_container(fd->ctr);
	    do {
		if (!(c = fd->ctr = cram_read_container(fd)))
		    return NULL;
	    } while (c->length == 0);
	}

	if (c->ref_seq_id != -2 && c->ref_seq_id != fd->range.refid) {
	    fd->eof = 1;
	    return NULL;
	}
    }

    if (!(c->comp_hdr_block = cram_read_block(fd)))
	return NULL;
    if (c->comp_hdr_block->content_type != COMPRESSION_HEADER)
	return NULL;

    c->comp_hdr = cram_decode_compression_header(fd, c->comp_hdr_block);
    if (!c->comp_hdr)
	return NULL;
    if (!c->comp_hdr->AP_delta &&
	sam_hdr_sort_order(fd->header) != ORDER_COORD) {
	if (fd->ref_lock) pthread_mutex_lock(fd->ref_lock);
	fd->unsorted = 1;
	if (fd->ref_lock) pthread_mutex_unlock(fd->ref_lock);
    }

    return c;
}

static cram_slice *cram_next_slice(cram_fd *fd, cram_container **cp) {
    cram_container *c_curr;  // container being consumed via cram_get_seq()
    cram_slice *s_curr = NULL;

    // Populate the first container if unknown.
    if (!(c_curr = fd->ctr)) {
	if (!(c_curr = cram_first_slice(fd)))
	    return NULL;
    }

    // Discard previous slice
    if ((s_curr = c_curr->slice)) {
	c_curr->slice = NULL;
	cram_free_slice(s_curr);
	s_curr = NULL;
    }

    // If we've consumed all slices in this container, also discard
    // the container too.
    if (c_curr->curr_slice == c_curr->max_slice) {
	if (fd->ctr == c_curr)
	    fd->ctr = NULL;
	if (fd->ctr_mt == c_curr)
	    fd->ctr_mt = NULL;
	cram_free_container(c_curr);
	c_curr = NULL;
    }

    if (!fd->ctr_mt)
	fd->ctr_mt = c_curr;

    // Fetch the next slice (and the container if necessary).
    //
    // If single threaded this loop bails out as soon as it finds
    // a slice in range.  In this case c_next and c_curr end up being
    // the same thing.
    //
    // If multi-threaded, we loop until we have filled out
    // thread pool input queue.  Here c_next and c_curr *may* differ, as
    // can fd->ctr and fd->ctr_mt.
    for (;;) {
	cram_container *c_next = fd->ctr_mt;
	cram_slice *s_next = NULL;

	// Next slice; either from the last job we failed to push
	// to the input queue or via more I/O.
	if (fd->job_pending) {
	    cram_decode_job *j = (cram_decode_job *)fd->job_pending;
	    c_next = j->c;
	    s_next = j->s;
	    free(fd->job_pending);
	    fd->job_pending = NULL;
	} else if (!fd->ooc) {
	empty_container:
	    if (!c_next || c_next->curr_slice_mt == c_next->max_slice) {
		// new container
		for(;;) {
		    if (!(c_next = cram_read_container(fd))) {
			if (fd->pool) {
			    fd->ooc = 1;
			    break;
			}

			return NULL;
		    }
		    c_next->curr_slice_mt = c_next->curr_slice;

		    if (c_next->length != 0)
			break;

		    cram_free_container(c_next);
		}
		if (fd->ooc)
		    break;

		/* Skip containers not yet spanning our range */
		if (fd->range.refid != -2 && c_next->ref_seq_id != -2) {
		    // ref_id beyond end of range; bail out
		    if (c_next->ref_seq_id != fd->range.refid) {
			cram_free_container(c_next);
			fd->ctr_mt = NULL;
			fd->ooc = 1;
			break;
		    }

		    // position beyond end of range; bail out
		    if (c_next->ref_seq_start > fd->range.end) {
			cram_free_container(c_next);
			fd->ctr_mt = NULL;
			fd->ooc = 1;
			break;
		    }

		    // before start of range; skip to next container
		    if (c_next->ref_seq_start + c_next->ref_seq_span-1 <
			fd->range.start) {
			c_next->curr_slice_mt = c_next->max_slice;
			cram_seek(fd, c_next->length, SEEK_CUR);
			cram_free_container(c_next);
			c_next = NULL;
			continue;
		    }
		}

		// Container is valid range, so remember it for restarting
		// this function.
		fd->ctr_mt = c_next;

		if (!(c_next->comp_hdr_block = cram_read_block(fd)))
		    return NULL;
		if (c_next->comp_hdr_block->content_type != COMPRESSION_HEADER)
		    return NULL;

		c_next->comp_hdr =
		    cram_decode_compression_header(fd, c_next->comp_hdr_block);
		if (!c_next->comp_hdr)
		    return NULL;

		if (!c_next->comp_hdr->AP_delta &&
		    sam_hdr_sort_order(fd->header) != ORDER_COORD) {
		    if (fd->ref_lock) pthread_mutex_lock(fd->ref_lock);
		    fd->unsorted = 1;
		    if (fd->ref_lock) pthread_mutex_unlock(fd->ref_lock);
		}
	    }

	    if (c_next->num_records == 0) {
		cram_free_container(c_next);
		if (fd->ctr_mt == c_next)
		    fd->ctr_mt = NULL;
		c_next = NULL;
		goto empty_container;
	    }

	    if (!(s_next = c_next->slice = cram_read_slice(fd)))
		return NULL;

	    s_next->slice_num = ++c_next->curr_slice_mt;
	    s_next->curr_rec = 0;
	    s_next->max_rec = s_next->hdr->num_records;

	    s_next->last_apos = s_next->hdr->ref_seq_start;
	    
	    // We know the container overlaps our range, but with multi-slice
	    // containers we may have slices that do not.  Skip these also.
	    if (fd->range.refid != -2 && s_next->hdr->ref_seq_id != -2) {
		// ref_id beyond end of range; bail out
		if (s_next->hdr->ref_seq_id != fd->range.refid) {
		    fd->ooc = 1;
		    cram_free_slice(s_next);
		    c_next->slice = s_next = NULL;
		    break;
		}

		// position beyond end of range; bail out
		if (s_next->hdr->ref_seq_start > fd->range.end) {
		    fd->ooc = 1;
		    cram_free_slice(s_next);
		    c_next->slice = s_next = NULL;
		    break;
		}

		// before start of range; skip to next slice
		if (s_next->hdr->ref_seq_start + s_next->hdr->ref_seq_span-1 <
		    fd->range.start) {
		    cram_free_slice(s_next);
		    c_next->slice = s_next = NULL;
		    continue;
		}
	    }
	} // end: if (!fd->ooc)

	if (!c_next || !s_next)
	    break;

	// Decode the slice, either right now (non-threaded) or by pushing
	// it to the a decode queue (threaded).
	if (cram_decode_slice_mt(fd, c_next, s_next, fd->header) != 0) {
	    fprintf(stderr, "Failure to decode slice\n");
	    cram_free_slice(s_next);
	    c_next->slice = NULL;
	    return NULL;
	}

	// No thread pool, so don't loop again
	if (!fd->pool) {
	    c_curr = c_next;
	    s_curr = s_next;
	    break;
	}

	// With thread pool, but we have a job pending so our decode queue
	// is full.
	if (fd->job_pending)
	    break;

	// Otherwise we're threaded with room in the decode input queue, so
	// keep reading slices for decode.
	// Push it a bit far, to qsize in queue rather than pending arrival,
	// as cram tends to be a bit bursty in decode timings.
	if (t_pool_results_queue_len(fd->rqueue) > fd->pool->qsize)
	    break;
    } // end of for(;;)


    // When not threaded we've already have c_curr and s_curr.
    // Otherwise we need get them by pulling off the decode output queue.
    if (fd->pool) {
	t_pool_result *res;
	cram_decode_job *j;
	
	if (fd->ooc && t_pool_results_queue_empty(fd->rqueue)) {
	    fd->eof = 1;
	    return NULL;
	}

	res = t_pool_next_result_wait(fd->rqueue);

	if (!res || !res->data) {
	    fprintf(stderr, "Call to t_pool_next_result failed\n");
	    return NULL;
	}

	j = (cram_decode_job *)res->data;
	c_curr = j->c;
	s_curr = j->s;

	if (j->exit_code != 0) {
	    fprintf(stderr, "Slice decode failure\n");
	    fd->eof = 0;
	    t_pool_delete_result(res, 1);
	    return NULL;
	}

	t_pool_delete_result(res, 1);
    }

    *cp = c_curr;

    // Update current slice being processed (as opposed to current
    // slice in the multi-threaded reahead.
    fd->ctr = c_curr;
    if (c_curr) {
	c_curr->slice = s_curr;
	if (s_curr)
	    c_curr->curr_slice = s_curr->slice_num;
    }
    if (s_curr)
	s_curr->curr_rec = 0;
    else
	fd->eof = 1;

    return s_curr;
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

    for (;;) {
	c = fd->ctr;
	if (c && c->slice && c->slice->curr_rec < c->slice->max_rec) {
	    s = c->slice;
	} else {
	    if (!(s = cram_next_slice(fd, &c)))
		return NULL;
	    continue; /* In case slice contains no records */
	}

	if (fd->range.refid != -2) {
	    if (fd->range.refid == -1 && s->crecs[s->curr_rec].ref_id != -1) {
		// Special case when looking for unmapped blocks at end.
		// If these are mixed in with mapped data (c->ref_id == -2)
		// then we need skip until we find the unmapped data, if at all
		s->curr_rec++;
		continue;
	    }
	    if (s->crecs[s->curr_rec].ref_id < fd->range.refid &&
		s->crecs[s->curr_rec].ref_id != -1) {
		// Looking for a mapped read, but not there yet.  Special case
		// as -1 (unmapped) shouldn't be considered < refid.
		s->curr_rec++;
		continue;
	    }

	    if (s->crecs[s->curr_rec].ref_id != fd->range.refid) {
		fd->eof = 1;
		cram_free_slice(s);
		c->slice = NULL;
		return NULL;
	    }

	    if (fd->range.refid != -1 && s->crecs[s->curr_rec].apos > fd->range.end) {
		fd->eof = 1;
		cram_free_slice(s);
		c->slice = NULL;
		return NULL;
	    }

	    if (fd->range.refid != -1 && s->crecs[s->curr_rec].aend < fd->range.start) {
		s->curr_rec++;
		continue;
	    }
	}

	break;
    }

    fd->ctr = c;
    c->slice = s;
    return &s->crecs[s->curr_rec++];
}

/*
 * Read the next cram record and convert it to a bam_seq_t struct.
 *
 * Returns 0 on success
 *        -1 on EOF or failure (check fd->err)
 */
int cram_get_bam_seq(cram_fd *fd, bam_seq_t **bam) {
    cram_record *cr;
    cram_container *c;
    cram_slice *s;

    if (!(cr = cram_get_seq(fd))) {
	//*bam=0;
	return -1;
    }

    c = fd->ctr;
    s = c->slice;

    if (s->bl) {
	//*bam = s->bl[s->curr_rec-1]; return 0;

	// Ideally we'd just do: *bam = s->bl[s->curr_rec-1];
	// That works, but it changes the API as the bam object is
	// no longer a malloced block of memory and cannot be
	// freed by the caller.  (Possibly we can do *bam=0
	// in the case where cram_get_seq hits EOF, but this is
	// also assuming that the *bam object was ours and not a
	// result of, say, a merge with a bam file.)
	//
	// Hence instead we laboriously manage the memory and do a
	// memcpy each time.  (This is around an extra 40% time taken
	// in main to decode a CRAM file, harming parallel execution.)
	int sz = s->bl[s->curr_rec-1]->alloc;
	if (!*bam) {
	    if (!(*bam = malloc(sz)))
		return -1;
	    (*bam)->alloc = sz;
	} else if ((*bam)->alloc < sz) {
	    if (!(*bam = realloc(*bam, sz)))
		return -1;
	    (*bam)->alloc = sz;
	}
	memcpy(*bam, s->bl[s->curr_rec-1], sz);
	return 0;
    }

    return cram_to_bam(fd->header, fd, s, cr, s->curr_rec-1, bam) >= 0 ? 0 : -1;
}

/*
 * Copyright (c) 2013-2015 Genome Research Ltd.
 * Author(s): James Bonfield
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
 * A cut down version of cram_dump.c to accumulate only size
 * information per data series.  This is much faster than cram_dump as
 * it does not require uncompression of data blocks.
 */

#include "io_lib_config.h"

#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <getopt.h>

#include <io_lib/cram.h>

// Lifted out of cram_io.c.
// Maybe make it public as cram_write_full_container.
static int cram_flush_container2(cram_fd *fd, cram_container *c) {
    int i, j;

    //fprintf(stderr, "Writing container %d, sum %u\n", c->record_counter, sum);

    /* Write the container struct itself */
    if (0 != cram_write_container(fd, c))
	return -1;

    /* And the compression header */
    if (0 != cram_write_block(fd, c->comp_hdr_block))
	return -1;

    /* Followed by the slice blocks */
    for (i = 0; i < c->curr_slice; i++) {
	cram_slice *s = c->slices[i];

	if (0 != cram_write_block(fd, s->hdr_block))
	    return -1;

	for (j = 0; j < s->hdr->num_blocks; j++) {
	    if (0 != cram_write_block(fd, s->block[j]))
		return -1;
	}
    }
    
    return CRAM_IO_FLUSH(fd) == 0 ? 0 : -1;
}


// Lifted from htslib/cram/cram_external.c
/*
 * Converts a cram_block_compression_hdr struct used for decoding to
 * one used for encoding.  Maybe this should be a transparent
 * operation applied on-demand.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_block_compression_hdr_decoder2encoder(cram_fd *fd, cram_container *c,
					       cram_block_compression_hdr *ch) {
    int i;

    if (!ch || !ch->codecs)
	return -1;

    for (i = 0; i < DS_END; i++) {
	cram_codec *co = ch->codecs[i];
	if (!co)
	    continue;

	if (-1 == cram_codec_decoder2encoder(fd, co))
	    return -1;
    }

    // Fix tag encoding map.
    if (ch->tag_encoding_map) {
	if (!(c->tags_used = HashTableCreate(16, HASH_DYNAMIC_SIZE)))
	    return -1;

	for (i = 0; i < CRAM_MAP_HASH; i++) {
	    cram_map *m;
	    HashData hd;
	    for (m = ch->tag_encoding_map[i]; m; m = m->next) {
		cram_tag_map *tm = calloc(1, sizeof(*tm));
		unsigned char key[3];
		tm->codec = m->codec;
		if (-1 == cram_codec_decoder2encoder(fd, tm->codec))
		    return -1;
		hd.p = tm;
		key[0] = (m->key>>16)&0xff;
		key[1] = (m->key>> 8)&0xff;
		key[2] = (m->key>> 0)&0xff;
		
		HashTableAdd(c->tags_used, key, 3, hd, NULL);
	    }
	}
    }

    // Fix various container bits which get generated as compression header bits
    c->pos_sorted = ch->AP_delta;

    return 0;
}

typedef struct {
    int id1, id2;
} ds_keys;

// Extracts content ids for all of the data series mentioned
// in the ds_h hash, storing the id in the top and bottom
// 32-bit values of hd->data.i (i32[2]).
//
// Also updates the content ID hash, ci_h, indicating which
// content IDs we to keep and to remove.
int ds_to_id(cram_map **ma, char *data, HashTable *ds_h, HashTable *ci_h) {
    int i;
    uintptr_t k;

    for (i = 0; i < CRAM_MAP_HASH; i++) {
	cram_map *m;
	for (m = ma[i]; m; m = m->next) {
	    HashItem *hi;

	    // Crude, only works with single byte ITF8 values
	    if (m->encoding == E_EXTERNAL ||
		m->encoding == E_BYTE_ARRAY_STOP ||
		m->encoding == E_BYTE_ARRAY_LEN) {
		HashData hd;
		// content id, iff cram_decoder_init fails below
		hd.i = (unsigned char)data[m->offset + m->size-1];

		// 2 byte keys are data series.
		// 3 byte keys are aux tags, but we hack the last byte
		// to be all permutations of valid data types.
		k = m->key;
		if (k>>16)
		    k &= ~0xff;

		cram_codec *c = cram_decoder_init(m->encoding,
						  data + m->offset,
						  m->size, E_BYTE_ARRAY, 0);
		int id1, id2;
		if (c) {
		    id1 = cram_codec_to_id(c, &id2);
		    c->free(c);
		    if (id1 >= 0) {
			int drop = 0;
			if ((hi = HashTableSearch(ds_h, (char *)k, sizeof(k))))
			    drop = 1, hi->data.i32[0] = id1;
			hd.i = 0;
			uintptr_t k2 = id1;
			hi = HashTableAdd(ci_h, (char *)k2, sizeof(k2), hd, NULL);
			hi->data.i32[drop]++;
		    }
		    if (id2 >= 0) {
			int drop = 0;
			if ((hi = HashTableSearch(ds_h, (char *)k, sizeof(k))))
			    drop = 1, hi->data.i32[1] = id2;
			hd.i = 0;
			uintptr_t k2 = id2;
			hi = HashTableAdd(ci_h, (char *)k2, sizeof(k2), hd, NULL);
			hi->data.i32[drop]++;
		    }
		} else {
		    int drop = 0;
		    if ((hi = HashTableSearch(ds_h, (char *)k, 4)))
			drop = 1, hi->data.i32[0] = id1;
		    hd.i = 0;
		    uintptr_t k2 = id1;
		    hi = HashTableAdd(ci_h, (char *)k2, sizeof(k2), hd, NULL);
		    hi->data.i32[drop]++;
		}
	    }
	}
    }

    return 0;
}


/*
 * Returns 0 on success;
 *        -1 on failure
 */
int filter_blocks(cram_fd *fd_in, cram_fd *fd_out, HashTable *ds_h, int drop_qs, char *keep_aux) {
    cram_container *c;
    off_t pos, pos2, hpos;
    char tag_to_del[128][128] = {0};
    char tag_to_keep[128][128] = {0};
    //int TC_id;

    if (keep_aux) {
	while (*keep_aux) {
	    tag_to_keep[keep_aux[0]&0x7f][keep_aux[1]&0x7f] = 1;
	    while (*keep_aux && *keep_aux != ',')
		keep_aux++;
	    if (*keep_aux)
		keep_aux++;
	}
    }

    // 1. Load container struct
    // 2. Load compression header into c->comp_hdr
    // 3. Load all slices for this container into c->slices[i]
    // 4.    Filter slice blocks and edit slice header.
    // 5. Edit compression header
    // 6. Edit container num_blocks and size.
    // 7. Write container
    // 8. Write compression header
    // 9. Write slices.

    // Load container struct
    pos = CRAM_IO_TELLO(fd_in);
    while ((c = cram_read_container(fd_in))) {
	int i, j;

	if (fd_in->empty_container) {
	    cram_free_container(c);
	    if (0 != cram_write_eof_block(fd_out))
		return -1;
	    break;
	}

	// Maps content IDs to keep/remove counters.
	// We use this to spot content IDs that clash (we asked to remove
	// them but something else is in this block that we must keep).
	HashItem *hi;
	HashTable *ci_h = HashTableCreate(128, HASH_DYNAMIC_SIZE|
					  HASH_NONVOLATILE_KEYS |
					  HASH_INT_KEYS);

	// Reset counters for ds_h items.
	HashIter *iter;
	iter = HashTableIterCreate();
	while ((hi = HashTableIterNext(ds_h, iter))) {
	    hi->data.i32[0] = UINT_MAX;
	    hi->data.i32[1] = UINT_MAX;
	}
	HashTableIterDestroy(iter);

	if (fd_in->err) {
	    perror("Cram container read");
	    return 1;
	}

	hpos = CRAM_IO_TELLO(fd_in);

	if (!c->length) {
	    pos = CRAM_IO_TELLO(fd_in);
	    continue;
	}

	// Load compression header and parse the content IDs.
	if (!(c->comp_hdr_block = cram_read_block(fd_in)))
	    return 1;
	assert(c->comp_hdr_block->content_type == COMPRESSION_HEADER);

	c->comp_hdr = cram_decode_compression_header(fd_in, c->comp_hdr_block);
	if (!c->comp_hdr)
	    return 1;

	cram_block_compression_hdr_decoder2encoder(fd_in, c, c->comp_hdr);

	if (keep_aux) {
	    iter = HashTableIterCreate();
	    while ((hi = HashTableIterNext(c->tags_used, iter))) {
		// If we have any tags listed in keep_aux, then explicitly
		// consider all others as candidate for removal.
		if (!tag_to_keep[hi->key[0]&0x7f][hi->key[1]&0x7f]) {
		    uintptr_t k = (uintptr_t)((hi->key[0]<<16)|
					      (hi->key[1]<<8));
		    HashData hd;
		    hd.i32[0] = UINT_MAX; hd.i32[1] = UINT_MAX;
		    HashTableAdd(ds_h, (char *)k, sizeof(k), hd, NULL);
		}
	    }
	    HashTableIterDestroy(iter);
	}

	ds_to_id(c->comp_hdr->rec_encoding_map,
		 (char *)c->comp_hdr_block->data, ds_h, ci_h);
	ds_to_id(c->comp_hdr->tag_encoding_map,
		 (char *)c->comp_hdr_block->data, ds_h, ci_h);

	// Determine which IDs matter (the ones we wish to remove), and check
	// that they are not shared by another data-series.
	int cid1 = -1, cid2 = -1;

	// Work out which tags we will be removing.
	// This is based on the ones we requested in ds_h and the
	// ones we can according to ci_h.
	iter = HashTableIterCreate();
	while ((hi = HashTableIterNext(c->tags_used, iter))) {
	    uintptr_t k = (uintptr_t)((hi->key[0]<<16)|
				      (hi->key[1]<<8));

	    HashItem *ds_hi;
	    if (!(ds_hi = HashTableSearch(ds_h, (char *)k, sizeof(k)))) {
		//printf("tag_to_del[%c][%c]=0\n", hi->key[0], hi->key[1]);
		tag_to_del[hi->key[0]][hi->key[1]] = 0;
		continue;
	    }

	    int keep = 1;
	    if (ds_hi->data.i32[0] != UINT_MAX) {
		HashItem *ci_hi;
		uintptr_t k2 = ds_hi->data.i32[0];
		if ((ci_hi = HashTableSearch(ci_h, (char *)k2, sizeof(k2)))) {
		    //if (TC_id < 0) ci_hi->data.i32[0] = 1;
		    if (ci_hi->data.i32[0] == 0 && ci_hi->data.i32[1] > 0)
			keep = 0;
		}
	    }
	    if (ds_hi->data.i32[1] != UINT_MAX) {
		HashItem *ci_hi;
		uintptr_t k2 = ds_hi->data.i32[1];
		if ((ci_hi = HashTableSearch(ci_h, (char *)k2, sizeof(k2)))) {
		    //if (TC_id < 0) ci_hi->data.i32[0] = 1;
		    if (ci_hi->data.i32[0] == 0 && ci_hi->data.i32[1] > 0)
			keep = 0;
		}
	    }

	    //printf("tag_to_del[%c][%c]=%d(*)\n", hi->key[0], hi->key[1], 1-keep);
	    tag_to_del[hi->key[0]][hi->key[1]] = 1-keep;
	}
	HashTableIterDestroy(iter);


	// Fix TD map to remove any tags.
	//int *nTL_dec = calloc(c->comp_hdr->nTL, sizeof(int));
	//int n_tag_dec = 0;
	char *cp1, *cp2;
	cp1 = cp2 = c->comp_hdr->TL[0];
	i = 0;
	while (i < c->comp_hdr->nTL) {
	    //char *tmp=cp1;printf(">  %s\n", tmp);
	    while (*cp1) {
		if (tag_to_del[cp1[0]][cp1[1]]) {
		    //printf("Del %.3s\n", cp1);
		    //nTL_dec[i]++;
		    //n_tag_dec++;
		} else {
		    cp2[0] = cp1[0];
		    cp2[1] = cp1[1];
		    cp2[2] = cp1[2];
		    cp2 += 3;
		}
		cp1 += 3;
	    }
	    *cp2++ = 0;
	    //printf("<  %s\n", tmp);
	    cp1++;
	    i++;
	}
	BLOCK_SIZE(c->comp_hdr->TD_blk) = cp2 - (char *)c->comp_hdr->TL[0];

	c->curr_slice = 0;
	c->slices = calloc(c->num_landmarks, sizeof(*c->slices));
	uintptr_t k;
	// assume slices in container have same no. blocks
	int *blocks_removed = NULL;
	int n_removed = 0;
	for (j = 0; j < c->num_landmarks; j++) {
	    cram_slice *s;
	    int id, id2;
	    
	    s = cram_read_slice(fd_in);
	    c->slices[c->curr_slice++] = s;
	    
	    if (!blocks_removed)
		blocks_removed = malloc(s->hdr->num_blocks * sizeof(int));

	    // Filter slice
	    for (id = id2 = 0; id < s->hdr->num_blocks; id++) {
		k = s->block[id]->content_id;
		
		hi = HashTableSearch(ci_h, (char *)k, sizeof(k));
		if (hi && hi->data.i32[0] == 0 && hi->data.i32[1] > 0) {
		    cram_free_block(s->block[id]);
		    blocks_removed[n_removed++] = k;
		} else {
		    s->block[id2] = s->block[id];
		    if (id > 0)
			s->hdr->block_content_ids[id2-1] = s->hdr->block_content_ids[id-1];
		    id2++;
		}
	    }
	    s->hdr->num_blocks = id2;

	    cram_free_block(s->hdr_block);
	    s->hdr_block = cram_encode_slice_header(fd_out, s);
	}

	// Compute new compression header
	cram_block *c_hdr;
	if (drop_qs) {
	    // Edit QS codec to be HUFFMAN constant value 255 (no qual).
	    cram_stats *stats = cram_stats_create();
	    cram_stats_add(stats, 255);
	    cram_stats_encoding(fd_out, stats);
	    if (c->comp_hdr->codecs[DS_QS])
		c->comp_hdr->codecs[DS_QS]->free(c->comp_hdr->codecs[DS_QS]);
	    c->comp_hdr->codecs[DS_QS] = cram_encoder_init(E_HUFFMAN, stats, E_BYTE,
							   NULL, fd_out->version);
	    cram_stats_free(stats);
	}

	free(blocks_removed);

	c_hdr = cram_encode_compression_header(fd_out, c, c->comp_hdr);
	if (c->comp_hdr_block)
	    cram_free_block(c->comp_hdr_block);
	c->comp_hdr_block = c_hdr;


	// New slice offsets.
	// This entire section should perhaps be a function. See end of the
	// cram_encode_container() function, where this came from.
	uint32_t slice_offset;
	{
	    slice_offset = c_hdr->method == RAW
		? c_hdr->uncomp_size
		: c_hdr->comp_size;
	    slice_offset += 2 + 4*IS_CRAM_3_VERS(fd_out) +
		itf8_size(c_hdr->content_id) +
		itf8_size(c_hdr->comp_size) +
		itf8_size(c_hdr->uncomp_size);
	}
	c->num_blocks = 1; // compression header
	c->length = 0;
	for (i = 0; i < c->curr_slice; i++) {
	    int j;
	    cram_slice *s = c->slices[i];
	
	    c->num_blocks += s->hdr->num_blocks + 1; // slice header
	    c->landmark[i] = slice_offset;

	    if (s->hdr->ref_seq_start + s->hdr->ref_seq_span >
		c->ref_seq_start + c->ref_seq_span) {
		c->ref_seq_span = s->hdr->ref_seq_start + s->hdr->ref_seq_span
		    - c->ref_seq_start;
	    }
	
	    slice_offset += s->hdr_block->method == RAW
		? s->hdr_block->uncomp_size
		: s->hdr_block->comp_size;

	    slice_offset += 2 + 4*IS_CRAM_3_VERS(fd_out) + 
		itf8_size(s->hdr_block->content_id) +
		itf8_size(s->hdr_block->comp_size) +
		itf8_size(s->hdr_block->uncomp_size);

	    for (j = 0; j < s->hdr->num_blocks; j++) {
		slice_offset += 2 + 4*IS_CRAM_3_VERS(fd_out) + 
		    itf8_size(s->block[j]->content_id) +
		    itf8_size(s->block[j]->comp_size) +
		    itf8_size(s->block[j]->uncomp_size);

		slice_offset += s->block[j]->method == RAW
		    ? s->block[j]->uncomp_size
		    : s->block[j]->comp_size;
	    }
	}
	c->length += slice_offset; // just past the final slice

	cram_flush_container2(fd_out, c);
	    
	HashTableDestroy(c->tags_used, 1);
	c->tags_used = NULL; // Avoids freeing codecs twice.
	cram_free_container(c);

	HashTableDestroy(ci_h, 0);
    }

    return 0;
}

void usage(int err) {
    fprintf(err ? stderr : stdout,
	    "Usage: cram_filter [options] in.cram out.cram\n"
	    "Valid options:\n"
	    "    -q            Drop quality strings (CRAM QS).\n"
	    "    -t tag-list   Discard comma separated list of tag types.\n"
	    "    -T tag-list   Keep only aux. tag types in the specified list.\n");
    exit(err);
}

int main(int argc, char **argv) {
    cram_fd *fd_in, *fd_out;
    int drop_qs = 0;
    char *keep_aux = NULL;
    int c;

    // Map of data series 2 or 3 byte code to content_id(s).
    HashTable *ds_h = HashTableCreate(128, HASH_DYNAMIC_SIZE|
				      HASH_NONVOLATILE_KEYS |
				      HASH_INT_KEYS);

    while ((c = getopt(argc, argv, "hqt:T:")) != -1) {
	switch (c) {
	case 't': {
	    while (*optarg) {
		HashData hd;
		uintptr_t k = (optarg[0]<<16) | (optarg[1]<<8);
		hd.i32[0] = UINT_MAX; hd.i32[1] = UINT_MAX;
		HashTableAdd(ds_h, (char *)k, sizeof(k), hd, NULL);
		while (*optarg && *optarg != ',')
		    optarg++;
		if (*optarg == ',')
		    optarg++;
	    }
	    break;
	}
	case 'T': keep_aux = optarg; break;
	case 'q': drop_qs = 1; break;
	case 'h': usage(0);
	default:  usage(1);
	}
    }

    if (argc - optind != 2) {
	fprintf(stderr, "Usage: cram_size [options] in.cram out.cram\n");
	return 1;
    }

    if (NULL == (fd_in = cram_open(argv[optind], "rb"))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", argv[optind]);
	return 1;
    }

    if (NULL == (fd_out = cram_open(argv[optind+1], "wb"))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", argv[optind+1]);
	return 1;
    }

    fd_out->header = fd_in->header;
    sam_hdr_incr_ref(fd_in->header);
    cram_write_SAM_hdr(fd_out, fd_out->header);


    // Mark the things we wish to remove.
    if (drop_qs) {
	HashData hd;
	uintptr_t k = ('Q'<<8) | 'S';
	hd.i32[0] = UINT_MAX;
	hd.i32[1] = UINT_MAX;
	HashTableAdd(ds_h, (char *)k, sizeof(k), hd, NULL);
    }

    if (0 != filter_blocks(fd_in, fd_out, ds_h, drop_qs, keep_aux))
	return 1;
    cram_close(fd_in);
    cram_close(fd_out);

    HashTableDestroy(ds_h, 0);

    return 0;
}

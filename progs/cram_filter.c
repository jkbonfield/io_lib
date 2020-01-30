/*
 * Copyright (c) 2016, 2019 Genome Research Ltd.
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
 * A tool to slice-n-dice cram files at the container / block level,
 * for efficient production of a subset without needing to uncompress
 * and recompress.
 */

#include "io_lib_config.h"

#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <getopt.h>

#include <io_lib/cram.h>

// Variable sized integer function pointers.
varint_vec vv;

// Lifted out of cram_io.c.
// Maybe make it public as cram_write_full_container.
static int cram_flush_container2(cram_fd *fd, cram_container *c) {
    int i, j;

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


// Lifted from htslib/cram/cram_external.c with modifications.
/*
 * Converts a cram_block_compression_hdr struct used for decoding to
 * one used for encoding.  Maybe this should be a transparent
 * operation applied on-demand.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int
cram_block_compression_hdr_decoder2encoder(cram_fd *fd_in, cram_fd *fd_out, cram_container *c,
					   cram_block_compression_hdr *ch) {
    int i;

    if (!ch)
	return -1;

    for (i = 0; i < DS_END; i++) {
	cram_codec *co = ch->codecs[i];
	if (!co)
	    continue;

	if (-1 == cram_codec_decoder2encoder(fd_in, co))
	    return -1;
    }

    // Fix tag encoding map.
    if (!(c->tags_used = HashTableCreate(16, HASH_DYNAMIC_SIZE)))
	return -1;

    for (i = 0; i < CRAM_MAP_HASH; i++) {
	cram_map *m;
	HashData hd;
	for (m = ch->tag_encoding_map[i]; m; m = m->next) {
	    cram_tag_map *tm = calloc(1, sizeof(*tm));
	    if (!tm) return -1;
	    unsigned char key[3];
	    tm->codec = m->codec;
	    if (-1 == cram_codec_decoder2encoder(fd_in, tm->codec))
		return -1;
	    hd.p = tm;
	    key[0] = (m->key>>16)&0xff;
	    key[1] = (m->key>> 8)&0xff;
	    key[2] = (m->key>> 0)&0xff;
		
	    HashTableAdd(c->tags_used, (char *)key, 3, hd, NULL);
	}
    }

    // Migrate misc. container header bits into the container itself.
    c->pos_sorted = ch->AP_delta;
    if (ch->read_names_included == 0)
	fd_out->lossy_read_names = 1;

    return 0;
}


// Extracts content ids for all of the data series mentioned
// in the ds_h hash, storing the id in the top and bottom
// 32-bit values of hd->data.i (i32[2]).
//
// Also updates the content ID hash, ci_h, indicating which
// content IDs we to keep and to remove.
int ds_to_id(cram_block_compression_hdr *hdr,
	     cram_map **ma, char *data, HashTable *ds_h, HashTable *ci_h) {
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

		cram_codec *c = cram_decoder_init(hdr, m->encoding,
						  data + m->offset,
						  m->size, E_BYTE_ARRAY, 0, &vv);
		int id1 = 0, id2;
		if (c) {
		    id1 = cram_codec_to_id(c, &id2);
		    c->free(c);
		    if (id1 >= 0) {
			int drop = 0;
			if ((hi = HashTableSearch(ds_h, (char *)k, sizeof(k))))
			    drop = 1, hi->data.i32[0] = id1;
			hd.i = 0;
			uintptr_t k2 = id1;
			hi = HashTableAdd(ci_h, (char *)k2, sizeof(k2),
					  hd, NULL);
			hi->data.i32[drop]++;
		    }
		    if (id2 >= 0) {
			int drop = 0;
			if ((hi = HashTableSearch(ds_h, (char *)k, sizeof(k))))
			    drop = 1, hi->data.i32[1] = id2;
			hd.i = 0;
			uintptr_t k2 = id2;
			hi = HashTableAdd(ci_h, (char *)k2, sizeof(k2),
					  hd, NULL);
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
 * Reads the compression header, turns it into an encoder suitable
 * structure, and updates ds_h (the hash of data series to discard)
 * based on keep_aux if set.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
static int process_comp_hdr(cram_fd *fd_in, cram_fd *fd_out, HashTable *ds_h,
			    cram_container *c, char *keep_aux,
			    char (*tag_to_keep)[128]) {
    if (!(c->comp_hdr_block = cram_read_block(fd_in)))
	return -1;
    assert(c->comp_hdr_block->content_type == COMPRESSION_HEADER);

    c->comp_hdr = cram_decode_compression_header(fd_in, c->comp_hdr_block);
    if (!c->comp_hdr)
	return -1;

    if (cram_block_compression_hdr_decoder2encoder(fd_in, fd_out, c, c->comp_hdr))
	return -1;

    // If we have any tags listed in keep_aux, then explicitly
    // consider all others as candidate for removal.
    if (keep_aux) {
	HashItem *hi;
	HashIter *iter = HashTableIterCreate();
	if (!iter)
	    return -1;
	while ((hi = HashTableIterNext(c->tags_used, iter))) {
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
    
    return 0;
}

/*
 * Given a list of data series / tags we desire to remove (ds_h)
 * we fill out the ci_h hash, indexed on content-id, allowing us to
 * work out which blocks we can actually delete and which are shared
 * with another data series that we are not also removing.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
static int find_tags_to_del(cram_fd *fd_in,
			    HashTable *ds_h, HashTable *ci_h,
			    cram_container *c,
			    char (*tag_to_keep)[128],
			    char (*tag_to_del)[128]) {
    HashItem *hi;
    HashIter *iter;

    if (ds_to_id(c->comp_hdr, c->comp_hdr->rec_encoding_map,
		 (char *)c->comp_hdr_block->data, ds_h, ci_h))
	return -1;
    if (ds_to_id(c->comp_hdr, c->comp_hdr->tag_encoding_map,
		 (char *)c->comp_hdr_block->data, ds_h, ci_h))
	return -1;

    // Work out which tags we will be removing.
    // This is based on the ones we requested in ds_h and the
    // ones we can according to ci_h.
    if (!(iter = HashTableIterCreate()))
	return -1;

    while ((hi = HashTableIterNext(c->tags_used, iter))) {
	uintptr_t k = (uintptr_t)((hi->key[0]<<16)|
				  (hi->key[1]<<8));

	HashItem *ds_hi;
	if (!(ds_hi = HashTableSearch(ds_h, (char *)k, sizeof(k)))) {
	    //printf("tag_to_del[%c][%c]=0\n", hi->key[0], hi->key[1]);
	    tag_to_del[hi->key[0]&0x7f][hi->key[1]&0x7f] = 0;
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
	tag_to_del[hi->key[0]&0x7f][hi->key[1]&0x7f] = 1-keep;
    }
    HashTableIterDestroy(iter);

    return 0;
}


// Fix TD map to remove any tags.
static void fix_TD_map(cram_container *c, char (*tag_to_del)[128]) {
    char *cp1, *cp2;
    cp1 = cp2 = (char *)c->comp_hdr->TL[0];
    int i;

    i = 0;
    while (i < c->comp_hdr->nTL) {
	while (*cp1) {
	    if (!tag_to_del[cp1[0]&0x7f][cp1[1]&0x7f]) {
		cp2[0] = cp1[0];
		cp2[1] = cp1[1];
		cp2[2] = cp1[2];
		cp2 += 3;
	    }
	    cp1 += 3;
	}
	*cp2++ = 0;
	cp1++;
	i++;
    }
    BLOCK_SIZE(c->comp_hdr->TD_blk) = cp2 - (char *)c->comp_hdr->TL[0];
    c->comp_hdr->TD_blk->crc32 = 0; // force recompute
}

/*
 * Filters all slices in the current container based on the list
 * of blocks to remove present in ci_h.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
static int filter_container(cram_fd *fd_in, cram_fd *fd_out,
			    HashTable *ci_h, cram_container *c,
			    int *eor) {
    int j;

    c->curr_slice = 0;
    c->slices = calloc(c->num_landmarks, sizeof(*c->slices));
    if (!c->slices)
	return -1;
    // assume slices in container have same no. blocks
    for (j = 0; j < c->num_landmarks; j++) {
	cram_slice *s;
	int id, id2;
	    
	s = cram_read_slice(fd_in);
	c->slices[c->curr_slice++] = s;

	// Check if first slice is beyond range.
	// This is complicated.  For fixed RI we've already checked for HUFFMAN
	// headers, but E_EXTERNAL blocks we need to actually uncompressed and
	// pull out the first value from the first slice in this container.
	//
	// We uncompress a duplicate of the appropriate block, so we don't
	// need to recompress it afterwards.
	*eor = 0;
	if (j == 0 && fd_in->range.refid != -2 && c->ref_seq_id == -2) {
	    if (c->comp_hdr->codecs[DS_RI] &&
		c->comp_hdr->codecs[DS_RI]->codec == E_EXTERNAL) {
		int cid = c->comp_hdr->codecs[DS_RI]->e_external.content_id;
		cram_block *b = cram_get_block_by_id(s, cid);
		cram_block *dup = malloc(sizeof(*dup));
		*dup = *b;
		dup->data = malloc(b->comp_size);
		memcpy(dup->data, b->data, b->comp_size);
		
		cram_uncompress_block(dup);
		int32_t rid;
		char *cp = (char *)BLOCK_DATA(dup);
		rid = fd_in->vv.varint_get32(&cp, NULL, NULL);
		cram_free_block(dup);
		if (rid > fd_in->range.refid) {
		    *eor = 1;
		    return 0;
		}
	    }
	}
	    
	// Filter slice
	for (id = id2 = 0; id < s->hdr->num_blocks; id++) {
	    uintptr_t k = s->block[id]->content_id;
		
	    HashItem *hi = HashTableSearch(ci_h, (char *)k, sizeof(k));
	    if (hi && hi->data.i32[0] == 0 && hi->data.i32[1] > 0) {
		cram_free_block(s->block[id]);
	    } else {
		s->block[id2] = s->block[id];
		if (id > 0)
		    s->hdr->block_content_ids[id2-1] =
			s->hdr->block_content_ids[id-1];
		id2++;
	    }
	}
	s->hdr->num_blocks = id2;

	cram_free_block(s->hdr_block);
	s->hdr_block = cram_encode_slice_header(fd_out, s);
    }

    return 0;
}


/*
 * Rebuilds the container compression header (& block).
 * This rewrites the QS data series if we're removing this block so
 * that all quality values are 255.  Runs of qual 255 is how BAM
 * handles quality "*" in SAM.
 */
void correct_compression_header(cram_fd *fd_out,
				cram_container *c,
				int drop_qs) {
    cram_block *c_hdr;

    if (drop_qs) {
	// Edit QS codec to be HUFFMAN constant value 255 (no qual).
	cram_stats *stats = cram_stats_create();
	cram_stats_add(stats, 255);
	cram_stats_encoding(fd_out, stats);
	if (c->comp_hdr->codecs[DS_QS])
	    c->comp_hdr->codecs[DS_QS]->free(c->comp_hdr->codecs[DS_QS]);
	c->comp_hdr->codecs[DS_QS] = cram_encoder_init(E_HUFFMAN, stats,
						       E_BYTE, NULL,
						       fd_out->version, &vv);
	cram_stats_free(stats);
    }

    c_hdr = cram_encode_compression_header(fd_out, c, c->comp_hdr);
    if (c->comp_hdr_block)
	cram_free_block(c->comp_hdr_block);
    c->comp_hdr_block = c_hdr;
}


/*
 * Recomputes the number of blocks and slice "landmarks".
 */
void update_slice_offsets(cram_fd *fd_out, cram_container *c) {
    cram_block *c_hdr = c->comp_hdr_block;
    uint32_t slice_offset;
    int i;

    slice_offset = c_hdr->method == RAW
	? c_hdr->uncomp_size
	: c_hdr->comp_size;
    slice_offset += 2 + 4*IS_CRAM_3_VERS(fd_out) +
	fd_out->vv.varint_size(c_hdr->content_id) +
	fd_out->vv.varint_size(c_hdr->comp_size) +
	fd_out->vv.varint_size(c_hdr->uncomp_size);
    
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
	    fd_out->vv.varint_size(s->hdr_block->content_id) +
	    fd_out->vv.varint_size(s->hdr_block->comp_size) +
	    fd_out->vv.varint_size(s->hdr_block->uncomp_size);

	for (j = 0; j < s->hdr->num_blocks; j++) {
	    slice_offset += 2 + 4*IS_CRAM_3_VERS(fd_out) + 
		fd_out->vv.varint_size(s->block[j]->content_id) +
		fd_out->vv.varint_size(s->block[j]->comp_size) +
		fd_out->vv.varint_size(s->block[j]->uncomp_size);

	    slice_offset += s->block[j]->method == RAW
		? s->block[j]->uncomp_size
		: s->block[j]->comp_size;
	}
    }
    c->length += slice_offset; // just past the final slice
}

/*
 * The heart of the CRAM block filtering algorithm.
 *
 * 1. Load container struct
 * 2. Load compression header into c->comp_hdr
 * 3. Load all slices for this container into c->slices[i]
 * 4.    Filter slice blocks and edit slice header.
 * 5. Edit compression header
 * 6. Edit container num_blocks and size.
 * 7. Write container
 * 8. Write compression header
 * 9. Write slices.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
int filter_blocks(cram_fd *fd_in, cram_fd *fd_out, HashTable *ds_h,
		  int drop_qs, char *keep_aux, int n_containers) {
    cram_container *c;
    char tag_to_del[128][128] = {{0}};
    char tag_to_keep[128][128] = {{0}};
    HashTable *ci_h = NULL;

    if (keep_aux) {
	while (*keep_aux) {
	    tag_to_keep[keep_aux[0]&0x7f][keep_aux[1]&0x7f] = 1;
	    while (*keep_aux && *keep_aux != ',')
		keep_aux++;
	    if (*keep_aux)
		keep_aux++;
	}
    }

    // Load container struct
    while ((c = cram_read_container(fd_in))) {
	if (fd_in->empty_container) {
	    cram_free_container(c);
	    continue;
	}

	if (fd_in->range.refid != -2) {
	    // It's possible we may have multiple references in one container.
	    // In theory we seeked to the first one we could, so just go with
	    // it. (NB: our index may have missing entries, so not optimal.)

	    // Beyond chr.
	    if (c->ref_seq_id != -2 && c->ref_seq_id > fd_in->range.refid)
		goto tidy;

	    // Beyond seq in chr
	    if (c->ref_seq_id != -2 && c->ref_seq_id == fd_in->range.refid &&
		c->ref_seq_start > fd_in->range.end)
		goto tidy;

	    // We may also have mixed length sequences causing range X-Y be
	    // covering disjoint containers on disk. Eg seqs in blocks
	    // labeled with - = . and ;.
	    //
	    //                       |<--range--->
	    // ----------------------------------------
	    //  -----     ====     ....     ;;;;     
	    //    -----     ====     ....     ;;;;     
	    //      -----     ====     ....     ;;;;     
	    //
	    // The "-" "." and ";" blocks cover range, but not "=" block.
	    // => Skip the "=" type blocks.
	    if (c->ref_seq_id != -2 &&
		!(c->ref_seq_id == fd_in->range.refid &&
		  c->ref_seq_start <= fd_in->range.end &&
		  c->ref_seq_start + c->ref_seq_span-1 >= fd_in->range.start))
		continue;
	} else {
	    // How do we deal with multi-ref containers?  We'll have seeked to
	    // the ideal starting point, but detecting the end point is tricky.
	    // We need to fetch the first RI value, if able, and check this.
	    //
	    // See later in this function for the implementation.
	}

	// Maps content IDs to keep/remove counters.
	// We use this to spot content IDs that clash (we asked to remove
	// them but something else is in this block that we must keep).
	HashItem *hi;
	if (!(ci_h = HashTableCreate(128, HASH_DYNAMIC_SIZE|
				     HASH_NONVOLATILE_KEYS |
				     HASH_INT_KEYS)))
	    return -1;

	// Reset counters for ds_h items.
	HashIter *iter;
	if (!(iter = HashTableIterCreate()))
	    return -1;
	while ((hi = HashTableIterNext(ds_h, iter))) {
	    hi->data.i32[0] = UINT_MAX;
	    hi->data.i32[1] = UINT_MAX;
	}
	HashTableIterDestroy(iter);

	if (fd_in->err) {
	    perror("Cram container read");
	    return -1;
	}

	if (!c->length)
	    continue;

	// Load compression header and parse the content IDs.
	if (process_comp_hdr(fd_in, fd_out, ds_h, c, keep_aux, tag_to_keep))
	    return -1;

	// Check RI tag if we're doing a range query and in multi-ref mode.
	if (fd_in->range.refid != -2 && c->ref_seq_id == -2) {
	    int id = -1;
	    if (c->comp_hdr->codecs[DS_RI]) {
		if (c->comp_hdr->codecs[DS_RI]->codec == E_HUFFMAN &&
		    // 1 item huffman is fine (zero bits used).
		    c->comp_hdr->codecs[DS_RI]->e_huffman.nvals == 1) {
		    id = c->comp_hdr->codecs[DS_RI]->e_huffman.codes[0].symbol;
		}
	    }
	    if (id > fd_in->range.refid)
		goto tidy;
	}

	// Check which blocks we are permitted to delete (eg shared).
	if (find_tags_to_del(fd_in, ds_h, ci_h, c, tag_to_keep, tag_to_del))
	    return -1;
	fix_TD_map(c, tag_to_del);

	// Filter all slices in this container.
	int eor;
	if (filter_container(fd_in, fd_out, ci_h, c, &eor))
	    return -1;
	if (eor == 1)
	    goto tidy;

	// Compute new compression header
	correct_compression_header(fd_out, c, drop_qs);

	// New slice offsets.
	update_slice_offsets(fd_out, c);

	// Write out the container and all slices.
	if (cram_flush_container2(fd_out, c) != 0)
	    return -1;
	    
	HashTableDestroy(c->tags_used, 1);
	c->tags_used = NULL; // Avoids freeing codecs twice.
	cram_free_container(c);

	HashTableDestroy(ci_h, 0); ci_h = NULL;

	if (n_containers && --n_containers <= 0)
	    break;
    }

    return cram_write_eof_block(fd_out);

 tidy:
    HashTableDestroy(c->tags_used, 1);
    c->tags_used = NULL; // Avoids freeing codecs twice.
    cram_free_container(c);
    if (ci_h)
	HashTableDestroy(ci_h, 0);

    return cram_write_eof_block(fd_out);
}


/*
 * -----------------------------------------------------------------------------
 */

/*
 * Loads the cram index and seeks to the Nth container.
 * Returns 0 on success
 *        -1 on failure
 */
int index_start(cram_fd *fd, char *fn, int container) {
    char fn_idx[PATH_MAX];
    FILE *fp;
    size_t len, buf_alloc = 0, buf_sz = 0;
    char *buf = NULL;

    snprintf(fn_idx, PATH_MAX, "%s.crai", fn);
    
    if (!(fp = fopen(fn_idx, "r"))) {
	perror(fn_idx);
	return -1;
    }
    
    // Load the entire index into memory
    buf = malloc((buf_alloc = 65536));
    if (!buf)
	return -1;
    while ((len = fread(buf + buf_sz, 1, 65536, fp)) > 0) {
	buf_sz += len;
	if (buf_alloc < buf_sz + 65536) {
	    buf_alloc *= 2;
	    buf = realloc(buf, buf_alloc);
	    if (!buf)
		return -1;
	}
    }
    
    // Uncompress if required
    if (buf_sz >= 2 && buf[0] == 31 && (unsigned char)buf[1] == 139) {
	char *u = zlib_mem_inflate(buf, buf_sz, &buf_sz);
	free(buf);
	if (!u)
	    return -1;
	buf = u;
    }

    // Skip <container> lines
    char *cp = buf;
    while (container--) {
	while (*cp && *cp != '\n')
	    cp++;
	if (!*cp)
	    break;
	cp++;
    }

    if (container != -1) {
	free(buf);
	return -1;
    }

    long dummy = 0;
    dummy |= strtol(cp, &cp, 10);
    dummy |= strtol(cp, &cp, 10);
    dummy |= strtol(cp, &cp, 10);
    if (dummy == LONG_MIN || dummy == LONG_MAX) {
	free(buf);
	return -1;
    }
    if (cram_seek(fd, strtoll(cp, &cp, 10), SEEK_SET) != 0) {
	free(buf);
	return -1;
    }

    free(buf);
    return 0;
}

void usage(int err) {
    fprintf(err ? stderr : stdout,
	"Usage: cram_filter [options] in.cram out.cram\n\n"
	"Valid options:\n"
	"    -n start[-end]    Only emit containers 'start' to 'end' inclusive.\n"
	"                      '-n 100' is equivalent to '-n 100-100'.\n"
	"    -r range          Limit output to containers overlapping 'range'.\n"
        "                      '-r chr1' matches all of chr1.\n"
	"                      '-r chr1:1000' is equivalent to '-r chr1:1000-1000'.\n"
	"    -q                Drop quality strings (CRAM QS).\n"
	"    -t tag-list       Discard comma separated list of tag types.\n"
	"    -T tag-list       Keep only aux. tag types in the specified list.\n"
	"    -!                Disable all checking of checksums.\n"
	"    -h                Show this help.\n"
	);
    exit(err);
}


int main(int argc, char **argv) {
    cram_fd *fd_in, *fd_out;
    int drop_qs = 0, ignore_md5 = 0;
    char *keep_aux = NULL, *range = NULL;
    int c, c_start = 0, c_end = -1, require_index = 0;

    // Map of data series 2 or 3 byte code to content_id(s).
    HashTable *ds_h = HashTableCreate(128, HASH_DYNAMIC_SIZE|
				      HASH_NONVOLATILE_KEYS |
				      HASH_INT_KEYS);
    if (!ds_h)
	return 1;

    // Parse arguments
    while ((c = getopt(argc, argv, "hqt:T:!n:r:")) != -1) {
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
	case '!': ignore_md5 = 1; break;

	case 'n':
	    c_start = strtol(optarg, &optarg, 0);
	    if (optarg && *optarg++)
		c_end = strtol(optarg, &optarg, 0);
	    else
		c_end = c_start;

	    require_index = 1;
	    break;

	case 'r':
	    range = optarg;
	    require_index = 2;
	    break;

	case 'h': usage(0);
	default:  usage(1);
	}
    }

    if (argc - optind != 2)
	usage(1);

    if (NULL == (fd_in = cram_open(argv[optind], "rb"))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", argv[optind]);
	return 1;
    }

    cram_init_varint(&vv, fd_in->file_def->major_version);
    
    // Parse index if required by -n and -r options.
    if (require_index == 1) {
	// For -n we parse the index manually as we need to track
	// the order of containers.
	if (index_start(fd_in, argv[optind], c_start) != 0) {
	    fprintf(stderr, "Failed to seek to container #%d\n", c_start);
	    return 1;
	}
    } if (require_index == 2) {
	// Need this for -r only.
	if (cram_index_load(fd_in, argv[optind])) {
	    fprintf(stderr, "Unable to load .crai index\n");
	    return 1;
	}

	int refid, start, end;
	cram_range r;
	char *cp = strchr(range, ':');
	if (cp) {
	    *cp = 0;
	    switch(sscanf(cp+1, "%d-%d", &start, &end)) {
	    case 1:
		end = start;
		break;
	    case 2:
		break;
	    default:
		fprintf(stderr, "Malformed range format\n");
		return 1;
	    }
	} else {
	    start = INT_MIN;
	    end   = INT_MAX;
	}

	if ((refid = sam_hdr_name2ref(fd_in->header, range)) == -1
	    && *range != '*') {
	    fprintf(stderr, "Unknown reference name '%s'\n", range);
	    return 1;
	}
	r.refid = refid;
	r.start = start;
	r.end   = end;
	if (cram_set_option(fd_in, CRAM_OPT_RANGE, &r) != 0) {
	    fprintf(stderr, "Failed to seek to range.\n");
	    return 1;
	}
    }

    if (NULL == (fd_out = cram_open(argv[optind+1], "wb"))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", argv[optind+1]);
	return 1;
    }

    if (ignore_md5) {
	if (cram_set_option(fd_in, CRAM_OPT_IGNORE_MD5, ignore_md5))
	    return 1;
	if (cram_set_option(fd_in, CRAM_OPT_IGNORE_CHKSUM, ignore_md5))
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

    if (0 != filter_blocks(fd_in, fd_out, ds_h, drop_qs, keep_aux,
			   c_end - c_start+1)) {
	fprintf(stderr, "Filter blocks failed\n");
	return 1;
    }
    cram_close(fd_in);
    cram_close(fd_out);

    HashTableDestroy(ds_h, 0);

    return 0;
}

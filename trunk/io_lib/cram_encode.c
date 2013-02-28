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

/*
 * Returns index of val into key.
 * Basically strchr(key, val)-key;
 */
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
cram_block *cram_encode_compression_header(cram_fd *fd, cram_container *c,
					   cram_block_compression_hdr *h) {
    cram_block *cb = cram_new_block(COMPRESSION_HEADER, 0);
    char *buf = malloc(100000); // FIXME, auto-grow this
    char *cp = buf;
    int i, mc;
    //int aux_B_len;
    char map[100000], *mp, cnt_buf[5];
    //char aux_B[100000];
    //cram_codec *aux_B_codec = NULL;

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

    if (fd->verbose)
	fprintf(stderr, "Wrote compression block header in %d bytes\n",
		(int)(cp-buf));

    cb->data = (unsigned char *)buf;
    cb->comp_size = cb->uncomp_size = cp - buf;

   //if (aux_B_codec)
   //	aux_B_codec->free(aux_B_codec);

    return cb;
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
    
    b->data = (unsigned char *)buf;
    b->comp_size = b->uncomp_size = cp-buf;
    
    return b;
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
    int nblk;

    /*
     * Slice external blocks:
     * ID 0 => base calls (insertions, soft-clip)
     * ID 1 => qualities
     * ID 2 => names
     * ID 3 => TS (insert size), NP (next frag)
     * ID 4 => tag values
     * ID 5 => BA, ifdef BA_external
     * ID 6 => tag IDs (TN), ifdef TN_external
     */

    /* Create cram slice header, num_blocks etc */
    s->hdr->ref_base_id = 0;
    nblk = 5;
#ifdef BA_external
    nblk++;
#endif
#ifdef TN_external
    nblk++;
#endif
    s->hdr->num_content_ids = nblk;
    s->hdr->num_blocks = s->hdr->num_content_ids+1;
    s->block = calloc(s->hdr->num_blocks, sizeof(s->block[0]));
    s->hdr->block_content_ids = malloc(s->hdr->num_content_ids *
				       sizeof(int32_t));
    s->hdr->block_content_ids[0] = 0;
    s->hdr->block_content_ids[1] = 1;
    s->hdr->block_content_ids[2] = 2;
    s->hdr->block_content_ids[3] = 3;
    s->hdr->block_content_ids[4] = 4;
    nblk = 5;
#ifdef BA_external
    s->hdr->block_content_ids[(s->ba_id = ++nblk)-1] = 5;
#endif
#ifdef TN_external
    s->hdr->block_content_ids[(s->tn_id = ++nblk)-1] = 6;
#endif

    s->block[0] = cram_new_block(CORE, 0);     // Core 
    s->block[1] = cram_new_block(EXTERNAL, 0); // IN Bases
    s->block[2] = cram_new_block(EXTERNAL, 1); // Qual
    s->block[3] = cram_new_block(EXTERNAL, 2); // Names
    s->block[4] = cram_new_block(EXTERNAL, 3); // TS/NP
    s->block[5] = cram_new_block(EXTERNAL, 4); // Tags
#ifdef BA_external
    s->block[s->ba_id] = cram_new_block(EXTERNAL, 5); // BA bases
#endif
#ifdef TN_external
    s->block[s->tn_id] = cram_new_block(EXTERNAL, 6); // TN ids
#endif

    core = s->block[0];
		 
    /* Create a formal method for stealing from dstrings! */
    //s->block[4]->data = calloc(5, s->hdr->num_records);
    s->block[4]->data = calloc(10, s->hdr->num_records); // NP TS
    s->block[4]->comp_size = s->block[4]->uncomp_size = 0;

#ifdef BA_external
    s->block[s->ba_id]->data = calloc(1, s->BA_len);
    s->block[s->ba_id]->comp_size = s->block[s->ba_id]->uncomp_size = 0;
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
	i32 = fd->cram_flag_swap[cr->flags & 0x7ff];
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
	    // RN codec: Already stored in block[3].
	}

	if (cr->cram_flags & CRAM_FLAG_DETACHED) {
	    char mf = cr->mate_flags;
	    r |= h->MF_codec->encode(s, h->MF_codec, core, &mf, 1);

	    if (!c->comp_hdr->read_names_included) {
		// RN codec: Already stored in block[3].
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
		itf8_put((char *)&s->block[4]->data[s->block[4]->uncomp_size],
			 cr->mate_pos);
	    s->block[4]->uncomp_size +=
		itf8_put((char *)&s->block[4]->data[s->block[4]->uncomp_size],
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
	// QS codec : Already stored in block[2].

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
		    s->block[s->ba_id]->data[s->block[s->ba_id]->uncomp_size++] = uc;
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
		    s->block[s->ba_id]->data[s->block[s->ba_id]->uncomp_size++] = uc;
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
	    char *seq = (char *)BLOCK_DATA(s->seqs_blk) + cr->seq;
#ifdef BA_external
	    memcpy(&s->block[s->ba_id]->data[s->block[s->ba_id]->uncomp_size],
		   seq, cr->len);
	    s->block[s->ba_id]->uncomp_size += cr->len;
#else
	    r |= h->BA_codec->encode(s, h->BA_codec, core, seq, cr->len);
#endif
	}
    }
    s->block[0]->uncomp_size = s->block[0]->byte + (s->block[0]->bit < 7);
    s->block[0]->comp_size = s->block[0]->uncomp_size;

    cram_free_block(s->block[1]);
    cram_free_block(s->block[2]);
    cram_free_block(s->block[3]);
    cram_free_block(s->block[5]);
    BLOCK_UPLEN(s->base_blk); s->block[1] = s->base_blk; s->base_blk = NULL;
    BLOCK_UPLEN(s->qual_blk); s->block[2] = s->qual_blk; s->qual_blk = NULL;
    BLOCK_UPLEN(s->name_blk); s->block[3] = s->name_blk; s->name_blk = NULL;
    BLOCK_UPLEN(s->aux_blk);  s->block[5] = s->aux_blk;  s->aux_blk  = NULL;

#ifdef TN_external
    cram_free_block(s->block[s->tn_id]);
    BLOCK_UPLEN(s->tn_blk); s->block[s->tn_id] = s->tn_blk; s->tn_blk = NULL;
#endif

    s->block[4]->comp_size = s->block[4]->uncomp_size;
    
#ifdef BA_external
    s->block[s->ba_id]->comp_size = s->block[s->ba_id]->uncomp_size;
#endif

    /* Compress the CORE Block too, with minimal zlib level */
    if (fd->level != 0)
	cram_compress_block(fd, s->block[0], NULL, 1, Z_FILTERED, -1, -1);

    /* Compress the other blocks */
    cram_compress_block(fd, s->block[1], fd->m[0], fd->level, Z_FILTERED,
			-1, -1);			      
    cram_compress_block(fd, s->block[2], fd->m[1], fd->level, Z_FILTERED, 
			1,Z_RLE);			      
    cram_compress_block(fd, s->block[3], fd->m[2], fd->level, Z_FILTERED,
			-1, -1);			      
    cram_compress_block(fd, s->block[4], fd->m[3], fd->level, Z_FILTERED,
			-1, -1);			      
    cram_compress_block(fd, s->block[5], fd->m[4], fd->level, Z_FILTERED,
			-1, -1);
#ifdef BA_external
    cram_compress_block(fd, s->block[s->ba_id], fd->m[5],
			fd->level, Z_FILTERED, -1, -1);
#endif
#ifdef TN_external
    cram_compress_block(fd, s->block[s->tn_id], fd->m[6],
			fd->level, Z_DEFAULT_STRATEGY, -1, -1);
#endif

    return r ? -1 : 0;
}

/*
 * Encodes all slices in a container into blocks.
 * Returns 0 on success
 *        -1 on failure
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
    h->BF_codec = cram_encoder_init(cram_stats_encoding(fd, c->BF_stats),
				    c->BF_stats, E_INT, NULL);

    //fprintf(stderr, "=== CF ===\n");
    h->CF_codec = cram_encoder_init(cram_stats_encoding(fd, c->CF_stats),
				    c->CF_stats, E_BYTE, NULL);

//    fprintf(stderr, "=== RN ===\n");
//    h->RN_codec = cram_encoder_init(cram_stats_encoding(fd, c->RN_stats),
//				    c->RN_stats, E_BYTE_ARRAY, NULL);

    //fprintf(stderr, "=== AP ===\n");
    h->AP_codec = cram_encoder_init(cram_stats_encoding(fd, c->AP_stats),
				    c->AP_stats, E_INT, NULL);

    //fprintf(stderr, "=== RG ===\n");
    h->RG_codec = cram_encoder_init(cram_stats_encoding(fd, c->RG_stats),
				    c->RG_stats, E_INT, NULL);

    //fprintf(stderr, "=== MQ ===\n");
    h->MQ_codec = cram_encoder_init(cram_stats_encoding(fd, c->MQ_stats),
				    c->MQ_stats, E_INT, NULL);

    //fprintf(stderr, "=== NS ===\n");
#ifdef NS_external
    h->NS_codec = cram_encoder_init(E_EXTERNAL, NULL, E_INT, (void *)3);
#else
    h->NS_codec = cram_encoder_init(cram_stats_encoding(fd, c->NS_stats),
				    c->NS_stats, E_INT, NULL);
#endif

    //fprintf(stderr, "=== MF ===\n");
    h->MF_codec = cram_encoder_init(cram_stats_encoding(fd, c->MF_stats),
				    c->MF_stats, E_BYTE, NULL);

#ifdef TS_external
    h->TS_codec = cram_encoder_init(E_EXTERNAL, NULL, E_INT, (void *)3);
    h->NP_codec = cram_encoder_init(E_EXTERNAL, NULL, E_INT, (void *)3);
#else
    //fprintf(stderr, "=== TS ===\n");
    h->TS_codec = cram_encoder_init(cram_stats_encoding(fd, c->TS_stats),
				    c->TS_stats, E_INT, NULL);
    //fprintf(stderr, "=== NP ===\n");
    h->NP_codec = cram_encoder_init(cram_stats_encoding(fd, c->NP_stats),
				    c->NP_stats, E_INT, NULL);
#endif

    //fprintf(stderr, "=== NF ===\n");
    h->NF_codec = cram_encoder_init(cram_stats_encoding(fd, c->NF_stats),
				    c->NF_stats, E_INT, NULL);

    //fprintf(stderr, "=== RL ===\n");
    h->RL_codec = cram_encoder_init(cram_stats_encoding(fd, c->RL_stats),
				    c->RL_stats, E_INT, NULL);

    //fprintf(stderr, "=== FN ===\n");
    h->FN_codec = cram_encoder_init(cram_stats_encoding(fd, c->FN_stats),
				    c->FN_stats, E_INT, NULL);

    //fprintf(stderr, "=== FC ===\n");
    h->FC_codec = cram_encoder_init(cram_stats_encoding(fd, c->FC_stats),
				    c->FC_stats, E_BYTE, NULL);

    //fprintf(stderr, "=== FP ===\n");
    h->FP_codec = cram_encoder_init(cram_stats_encoding(fd, c->FP_stats),
				    c->FP_stats, E_INT, NULL);

    //fprintf(stderr, "=== DL ===\n");
    h->DL_codec = cram_encoder_init(cram_stats_encoding(fd, c->DL_stats),
				    c->DL_stats, E_INT, NULL);

#ifdef BA_external
    h->BA_codec = cram_encoder_init(E_EXTERNAL, NULL, E_BYTE, (void *)5);
#else
    //fprintf(stderr, "=== BA ===\n");
    h->BA_codec = cram_encoder_init(cram_stats_encoding(fd, c->BA_stats),
				    c->BA_stats, E_BYTE, NULL);
#endif

    //fprintf(stderr, "=== BS ===\n");
    h->BS_codec = cram_encoder_init(cram_stats_encoding(fd, c->BS_stats),
				    c->BS_stats, E_BYTE, NULL);

    //fprintf(stderr, "=== TC ===\n");
    h->TC_codec = cram_encoder_init(cram_stats_encoding(fd, c->TC_stats),
				    c->TC_stats, E_BYTE, NULL);

    //fprintf(stderr, "=== TN ===\n");
    if (0) {
	//h->TN_codec = cram_encoder_init(cram_stats_encoding(fd, c->TN_stats),
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
	h->TN_codec = cram_encoder_init(E_EXTERNAL, NULL, E_INT, (void *)6);
#else
	h->TN_codec = cram_encoder_init(cram_stats_encoding(fd, c->TN_stats),
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
	if (fd->verbose)
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

	c_hdr = cram_encode_compression_header(fd, c, h);
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

    cram_free_compression_header(h);
    c->comp_hdr_block = c_hdr;

    return 0;
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
static int cram_add_feature(cram_container *c, cram_slice *s,
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

    s->features[s->nfeatures++] = *f;

    return 0;
}

static int cram_add_substitution(cram_fd *fd, cram_container *c,
				 cram_slice *s, cram_record *r,
				 int pos, char base, char qual, char ref) {
    cram_feature f;

    // seq=ACGTN vs ref=ACGT or seq=ACGT vs ref=ACGTN
    if (fd->L2[(uc)base]<4 || (fd->L2[(uc)base]<5 && fd->L2[(uc)ref]<4)) {
	f.X.pos = pos+1;
	f.X.code = 'X';
	f.X.base = fd->cram_sub_matrix[ref&0x1f][base&0x1f];
	cram_stats_add(c->BS_stats, f.X.base);
    } else {
	f.B.pos = pos+1;
	f.B.code = 'B';
	f.B.base = base;
	f.B.qual = qual;
	cram_stats_add(c->BA_stats, f.B.base);
	cram_stats_add(c->QS_stats, f.B.qual);
	BLOCK_APPEND_CHAR(s->qual_blk, qual);
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
    f.S.seq_idx = BLOCK_SIZE(s->base_blk);
    BLOCK_APPEND(s->base_blk, base, len);
    BLOCK_APPEND_CHAR(s->base_blk, '\0');
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
	f.S.seq_idx = BLOCK_SIZE(s->base_blk);
	BLOCK_APPEND(s->base_blk, base, len);
	BLOCK_APPEND_CHAR(s->base_blk, '\0');
    }
    return cram_add_feature(c, s, r, &f);
}

/*
 * Encodes auxiliary data.
 * Returns the read-group parsed out of the BAM aux fields on success
 *         NULL on failure or no rg present (FIXME)
 */
static char *cram_encode_aux(cram_fd *fd, bam_seq_t *b, cram_container *c,
			     cram_slice *s, cram_record *cr) {
    char *aux, *tmp, *rg = NULL, *tmp_tn;
    int aux_size = b->blk_size - ((char *)bam_aux(b) - (char *)&b->ref);
	
    /* Worst case is 1 nul char on every ??:Z: string, so +33% */
    BLOCK_GROW(s->aux_blk, aux_size*1.34+1);
    tmp = (char *)BLOCK_END(s->aux_blk);

#ifdef TN_external
    BLOCK_GROW(s->tn_blk, aux_size);
    tmp_tn = (char *)BLOCK_END(s->tn_blk);
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
	tmp_tn += itf8_put(tmp_tn, i32);
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

    cr->aux = BLOCK_SIZE(s->aux_blk);
    cr->aux_size = (uc *)tmp - (BLOCK_DATA(s->aux_blk) + cr->aux);
    BLOCK_SIZE(s->aux_blk) = (uc *)tmp - BLOCK_DATA(s->aux_blk);
    assert(s->aux_blk->byte <= s->aux_blk->alloc);

#ifdef TN_external
    cr->tn = BLOCK_SIZE(s->tn_blk);
    BLOCK_SIZE(s->tn_blk) = (uc *)tmp_tn - BLOCK_DATA(s->tn_blk);
    assert(s->tn_blk->byte <= s->tn_blk->alloc);
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
	if (fd->verbose)
	    fprintf(stderr, "Flush container %d/%d..%d\n",
		    c->ref_seq_id, c->ref_seq_start,
		    c->ref_seq_start + c->ref_seq_span -1);

	/* Encode slices */
	if (-1 == cram_flush_container(fd, c))
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
    cram_stats_add(c->BF_stats, fd->cram_flag_swap[cr->flags & 0x7ff]);

    cr->cram_flags  = CRAM_FLAG_PRESERVE_QUAL_SCORES; // FIXME
    //cram_stats_add(c->CF_stats, cr->cram_flags);

    cr->len         = bam_seq_len(b);  cram_stats_add(c->RL_stats, cr->len);
    cr->apos        = b->pos+1;
    cram_stats_add(c->AP_stats, cr->apos - s->last_apos);
    s->last_apos = cr->apos;

    cr->name        = BLOCK_SIZE(s->name_blk);
    cr->name_len    = bam_name_len(b); cram_stats_add(c->RN_stats, cr->name_len);
    BLOCK_APPEND(s->name_blk, bam_name(b), bam_name_len(b));


    /*
     * This seqs_ds is largely pointless and it could reuse the same memory
     * over and over.
     * s->base_ds is what we need for encoding.
     */
    cr->seq         = BLOCK_SIZE(s->seqs_blk);
    cr->qual        = BLOCK_SIZE(s->qual_blk);
    BLOCK_GROW(s->seqs_blk, cr->len);
    BLOCK_GROW(s->qual_blk, cr->len);
    seq = cp = (char *)BLOCK_END(s->seqs_blk);

    for (i = 0; i < cr->len; i++) {
	// FIXME: do 2 char at a time for efficiency
	cp[i] = bam_nt16_rev_table[bam_seqi(bam_seq(b), i)];
    }
    BLOCK_SIZE(s->seqs_blk) += cr->len;

    qual = cp = bam_qual(b);

    /* Copy and parse */
    if (!(cr->flags & BAM_FUNMAP)) {
	int32_t *cig_to, *cig_from;
	int apos = cr->apos-1, spos = 0;

	cr->cigar       = s->ncigar;
	cr->ncigar      = bam_cigar_len(b);
	while (cr->cigar + cr->ncigar >= s->cigar_alloc) {
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
			cram_add_substitution(fd, c, s, cr, spos,
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
		    cram_add_substitution(fd, c, s, cr, spos,
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
    BLOCK_GROW(s->qual_blk, cr->len);
    qual = cp = (char *)BLOCK_END(s->qual_blk);
    for (i = 0; i < cr->len; i++) {
	cp[i] = bam_qual(b)[i];
    }
    BLOCK_SIZE(s->qual_blk) += cr->len;

    /* Now we know apos and aend both, update mate-pair information */
    {
	int new;
	HashData hd;
	HashItem *hi;

	hd.i = c->curr_rec-1;
	//fprintf(stderr, "Checking %"PRId64"/%.*s\t", hd.i,
	//	cr->name_len, DSTRING_STR(s->name_ds)+cr->name);
	if (cr->flags & BAM_FPAIRED) {
	    hi = HashTableAdd(s->pair,
			      (char *)BLOCK_DATA(s->name_blk)+cr->name,
			      cr->name_len, hd, &new);
	} else {
	    new = 1;
	}

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

    if (fd->last_base < cr->aend)
	fd->last_base = cr->aend;

    return 0;
}



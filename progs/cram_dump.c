/*
 * A debugging program to dump out information on the layout of a CRAM file.
 * It's an abomination frankly, but isn't intended for production use.
 */

#include "io_lib_config.h"

#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>

#include <io_lib/cram.h>

void HashTableDumpMap(HashTable *h, FILE *fp, char *prefix, char *data) {
    int i, j, k;
    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi;
	cram_map *m;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    m = hi->data.p;
	    fprintf(fp, "%s%.*s => %16s {",
		    prefix ? prefix : "",
		    hi->key_len, hi->key,
		    cram_encoding2str(m->encoding));
	    for (k = m->offset, j = 0; j < m->size; j++, k++) {
		printf(j ? ", %d" : "%d", (unsigned char)data[k]);
	    }
	    printf("}\n");
	}
    }
}

void DumpMap2(cram_map **ma, FILE *fp, char *prefix, char *data) {
    int i, j, k;
    for (i = 0; i < CRAM_MAP_HASH; i++) {
	cram_map *m;
	for (m = ma[i]; m; m = m->next) {
	    fprintf(fp, "%s%c%c%c => %16s {",
		    prefix ? prefix : "",
		    (m->key>>16) & 0xff ? (m->key>>16) & 0xff : ' ',
		    (m->key>> 8) & 0xff,
		    (m->key>> 0) & 0xff,
		    cram_encoding2str(m->encoding));
	    for (k = m->offset, j = 0; j < m->size; j++, k++) {
		printf(j ? ", %d" : "%d", (unsigned char)data[k]);
	    }
	    printf("}\n");
	}
    }
}

static cram_map *map_find(cram_map **map, unsigned char *key, int id) {
    cram_map *m;

    m = map[CRAM_MAP(key[0],key[1])];
    while (m && m->key != id)
	m= m->next;

    assert(m);

    return m;
}

void dump_core_block(cram_block *b, int verbose) {
    int i;

    printf("Data = {");
    for (i = 0; (verbose || i < 100) && i < b->uncomp_size; i++) {
	printf(i ? ", %02x" : "%02x", (unsigned char)b->data[i]);
    }
    if (i < b->uncomp_size)
	printf(", ...}\n");
    else
	printf("}\n");
}

void dump_seq_block(cram_block *b, int verbose) {
    int i;
    for (i = 0; (verbose || i < 100) && i < b->uncomp_size; i++) {
	if (isprint(b->data[i]))
	    putchar(b->data[i]);
	else
	    printf("\\%03o", b->data[i]);
    }
}

void dump_quality_block(cram_block *b, int verbose) {
    int i;
    for (i = 0; (verbose || i < 100) && i < b->uncomp_size; i++) {
	putchar(b->data[i] + '!');
    }
    putchar('\n');
}

void dump_name_block(cram_block *b, int verbose) {
    int i;
    for (i = 0; (verbose || i < 100) && i < b->uncomp_size; i++) {
	if (isprint(b->data[i]))
	    putchar(b->data[i]);
	else
	    printf("\\%03o", b->data[i]);
    }
}

void dump_mate_info_block(cram_block *b, int verbose) {
    return dump_core_block(b, verbose);
}

void dump_tag_block(cram_block *b, int verbose) {
    return dump_core_block(b, verbose);
}

int main(int argc, char **argv) {
    cram_fd *fd;
    cram_container *c;
    off_t pos, pos2, hpos;
    int verbose = 0;

    static int bsize[100], bmax = 0;

    if (argc >= 2 && strcmp(argv[1], "-v") == 0) {
	argc--;
	argv++;
	verbose = 1;
    }

    if (argc < 2) {
	fprintf(stderr, "Usage: cram_dump [-v] filename.cram\n");
	return 1;
    }

    if (NULL == (fd = cram_open(argv[1], "rb"))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", argv[1]);
	return 1;
    }

    printf("File definition structure\n");
    printf("    Major vers: %d\n", fd->file_def->major_version);
    printf("    Minor vers: %d\n", fd->file_def->minor_version);
    printf("    file_id:    %.20s\n", fd->file_def->file_id);

    printf("\nBAM header:\n%.*s\n",
	   sam_header_length(fd->SAM_hdr),
	   sam_header_str(fd->SAM_hdr));

    pos = ftello(fd->fp);
    while ((c = cram_read_container(fd))) {
	int i, j;

	if (fd->err) {
	    perror("Cram container read");
	    return 1;
	}

	printf("\nContainer pos %"PRId64" size %d\n", (int64_t)pos, c->length);
	printf("    Ref id:            %d\n", c->ref_seq_id);
	printf("    Ref pos:           %d + %d\n", c->ref_seq_start, c->ref_seq_span);
	printf("    Rec counter:       %d\n", c->record_counter);
       	printf("    No. recs:          %d\n", c->num_records);
	printf("    No. bases          %"PRId64"\n", c->num_bases);
	printf("    No. blocks:        %d\n", c->num_blocks);
	printf("    No. landmarks:     %d\n", c->num_landmarks);

	printf("    {");
	for (i = 0; i < c->num_landmarks; i++) {
	    printf(i ? ", %d" : "%d", c->landmark[i]);
	}
	printf("}\n");

	hpos = ftello(fd->fp);

	printf("\n    Container_header block pos %"PRId64"\n", (int64_t)hpos);
	if (!(c->comp_hdr_block = cram_read_block(fd)))
	    return 1;
	assert(c->comp_hdr_block->content_type == COMPRESSION_HEADER);

	c->comp_hdr = cram_decode_compression_header(fd, c->comp_hdr_block);
	if (!c->comp_hdr)
	    return 1;

	printf("      Preservation map:\n");
	HashTableDump(c->comp_hdr->preservation_map, stdout, "\t");

	printf("\n      Record encoding map:\n");
	DumpMap2(c->comp_hdr->rec_encoding_map, stdout, "\t", 
		 (char *)c->comp_hdr_block->data);
	//HashTableDumpMap(c->comp_hdr->rec_encoding_map, stdout, "\t", c->comp_hdr_block->data);

	printf("\n      Tag encoding map:\n");
	DumpMap2(c->comp_hdr->tag_encoding_map, stdout, "\t",
		 (char *)c->comp_hdr_block->data);
	


	for (j = 0; j < c->num_landmarks; j++) {
	    cram_slice *s;
	    int id;
	    
	    pos2 = ftello(fd->fp);
	    assert(pos2 - pos - c->offset == c->landmark[j]);

	    s = cram_read_slice(fd);
	    printf("\n    Slice %d/%d, container offset %d, file offset %d\n", j+1, c->num_landmarks, (int)(pos2 - pos - c->offset), (int)pos2);
	    printf("\tSlice content type %s\n",
		   cram_content_type2str(s->hdr->content_type));

	    if (s->hdr->content_type == MAPPED_SLICE) {
		int i;
		printf("\tSlice ref seq    %d\n", s->hdr->ref_seq_id);
		printf("\tSlice ref start  %d\n", s->hdr->ref_seq_start);
		printf("\tSlice ref span   %d\n", s->hdr->ref_seq_span);
		printf("\tSlice MD5        ");
		for (i = 0; i < 16; i++)
		    printf("%02x", s->hdr->md5[i]);
		putchar('\n');
	    }
	    printf("\tRec counter      %d\n", s->hdr->record_counter);
	    printf("\tNo. records      %d\n", s->hdr->num_records);
	    printf("\tNo. blocks       %d\n", s->hdr->num_blocks);
	    printf("\tBlk IDS:         {");
	    for (id = 0; id < s->hdr->num_content_ids; id++) {
		printf(id ? ", %d" : "%d", s->hdr->block_content_ids[id]);
	    }
	    printf("}\n");
	    if (s->hdr->content_type == MAPPED_SLICE) {
		printf("\tRef base id:     %d\n", s->hdr->ref_base_id);
	    }
	
	    for (id = 0; id < s->hdr->num_blocks; id++)
		bsize[id] += s->block[id]->comp_size;
	    if (bmax < s->hdr->num_blocks)
		bmax = s->hdr->num_blocks;

	    for (id = 0; id < s->hdr->num_blocks; id++)
		cram_uncompress_block(s->block[id]);

	    /* Test decoding of 1st seq */
	    if (verbose) {
		cram_block *b = s->block[0];
		int32_t i32, bf, fn, prev_pos, rl;
		unsigned char cf;
		int out_sz, r, f;
		int rec;

		assert(b->content_type == CORE);

		for (rec = 0; rec < s->hdr->num_records; rec++) {
		    unsigned char ntags;

		    printf("Rec %d/%d at %d,%d\n", rec+1, s->hdr->num_records,
			   b->byte, b->bit);

		    out_sz = 1; /* decode 1 item */
		    r = c->comp_hdr->BF_codec->decode(s,c->comp_hdr->BF_codec, b, (char *)&bf, &out_sz);
		    printf("BF = %d => SAM 0x%x (ret %d, out_sz %d)\n", bf, fd->bam_flag_swap[bf], r, out_sz);
		    bf = fd->bam_flag_swap[bf];

		    r = c->comp_hdr->CF_codec->decode(s,c->comp_hdr->CF_codec, b, (char *)&cf, &out_sz);
		    printf("CF = %d (ret %d, out_sz %d)\n", cf, r, out_sz);

		    if (fd->version != CRAM_1_VERS) {
			int32_t ri;
			r |= c->comp_hdr->RI_codec->decode(s, c->comp_hdr->RI_codec, b, (char *)&ri, &out_sz);
			printf("RI = %d (ret %d, out_sz %d)\n", ri, r, out_sz);
		    }

		    r = c->comp_hdr->RL_codec->decode(s,c->comp_hdr->RL_codec, b, (char *)&rl, &out_sz);
		    printf("RL = %d (ret %d, out_sz %d)\n", rl, r, out_sz);

		    r = c->comp_hdr->AP_codec->decode(s,c->comp_hdr->AP_codec, b, (char *)&i32, &out_sz);
		    printf("AP = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

		    r = c->comp_hdr->RG_codec->decode(s,c->comp_hdr->RG_codec, b, (char *)&i32, &out_sz);
		    printf("RG = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

		    if (c->comp_hdr->read_names_included) {
			int32_t out_sz2 = 1;
			cram_block *dat = cram_new_block(EXTERNAL, 0);

			r = c->comp_hdr->RN_codec->decode(s,c->comp_hdr->RN_codec, b, (char *)dat, &out_sz2);
			printf("RN = %.*s (ret %d, out_sz %d)\n", out_sz2, BLOCK_DATA(dat), r, out_sz2);
			cram_free_block(dat);
		    }

		    if (cf & CRAM_FLAG_DETACHED) {
			char mf;
			puts("Detached");
			/* MF, RN if !captureReadNames, NS, NP, IS */
			r = c->comp_hdr->MF_codec->decode(s,c->comp_hdr->MF_codec, b, &mf, &out_sz);
			printf("MF = %d (ret %d, out_sz %d)\n", mf, r, out_sz);

			if (!c->comp_hdr->read_names_included) {
			    cram_block *dat = cram_new_block(EXTERNAL, 0);
			    int32_t out_sz2 = 1;

			    r = c->comp_hdr->RN_codec->decode(s,c->comp_hdr->RN_codec, b, (char *)dat, &out_sz2);
			    printf("RN = %.*s (ret %d, out_sz %d)\n", out_sz2, BLOCK_DATA(dat), r, out_sz2);
			    cram_free_block(dat);
			}
		    
			r = c->comp_hdr->NS_codec->decode(s,c->comp_hdr->NS_codec, b, (char *)&i32, &out_sz);
			printf("NS = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

			r = c->comp_hdr->NP_codec->decode(s,c->comp_hdr->NP_codec, b, (char *)&i32, &out_sz);
			printf("NP = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

			r = c->comp_hdr->TS_codec->decode(s,c->comp_hdr->TS_codec, b, (char *)&i32, &out_sz);
			printf("TS = %d (ret %d, out_sz %d)\n", i32, r, out_sz);
		    } else if (cf & CRAM_FLAG_MATE_DOWNSTREAM) {
			puts("Not detached, and mate is downstream");
			r = c->comp_hdr->NF_codec->decode(s,c->comp_hdr->NF_codec, b, (char *)&i32, &out_sz);
			printf("NF = %d+%d = %d (ret %d, out_sz %d)\n", i32, rec+1, i32+rec+1, r, out_sz);
		    }

		    if (fd->version == CRAM_1_VERS) {
			r = c->comp_hdr->TC_codec->decode(s,c->comp_hdr->TC_codec, b, (char *)&ntags, &out_sz);
			printf("TC = %d (ret %d, out_sz %d)\n", ntags, r, out_sz);

			for (f = 0; f < ntags; f++) {
			    int32_t id;
			    char key[3];
			    cram_map *m;
			    cram_block *tag = cram_new_block(EXTERNAL, 0);

			    r = c->comp_hdr->TN_codec->decode(s, c->comp_hdr->TN_codec,
							      b, (char *)&id, &out_sz);
			    key[0] = (id>>16)&0xff;
			    key[1] = (id>>8)&0xff;
			    key[2] = id&0xff;
			    printf("%3d: TN= %.3s\n", f, key);

			    printf("id=%d\n", id);
			    if ((m = map_find(c->comp_hdr->tag_encoding_map,
					      (unsigned char *)key, id))) {
				int i, out_sz;

				BLOCK_SIZE(tag) = 0;
				r = m->codec->decode(s, m->codec, b, (char *)tag, &out_sz);
				printf("%3d: Val", f);
				for(i = 0; i < out_sz; i++) {
				    printf(" %02x", BLOCK_DATA(tag)[i]);
				}
				printf("\n");
			    } else {
				fprintf(stderr, "*** ERROR: unrecognised aux key ***\n");
			    }

			    cram_free_block(tag);
			    // skip decoding of tag data itself and hope it's
			    // in an external block.
			}
		    } else {
			int32_t tl, ntags;
			char *tn;
			r = c->comp_hdr->TL_codec->decode(s,c->comp_hdr->TL_codec, b, (char *)&tl, &out_sz);
			printf("TL = %d (ret %d, out_sz %d)\n", tl, r, out_sz);

			tn = (char *)c->comp_hdr->TL[tl];
			ntags = strlen(tn)/3;

			for (f = 0; f < ntags; f++) {
			    int32_t id;
			    char key[3];
			    cram_map *m;
			    cram_block *tag = cram_new_block(EXTERNAL, 0);

			    key[0] = *tn++;
			    key[1] = *tn++;
			    key[2] = *tn++;
			    id = (key[0]<<16) | (key[1]<<8) | key[2];
			    printf("%3d: TN= %.3s\n", f, key);

			    printf("id=%d\n", id);
			    if ((m = map_find(c->comp_hdr->tag_encoding_map,
					      (unsigned char *)key, id))) {
				int i, out_sz;

				BLOCK_SIZE(tag) = 0;
				r = m->codec->decode(s, m->codec, b, (char *)tag, &out_sz);
				printf("%3d: Val", f);
				for(i = 0; i < out_sz; i++) {
				    printf(" %02x", BLOCK_DATA(tag)[i]);
				}
				printf("\n");
			    } else {
				fprintf(stderr, "*** ERROR: unrecognised aux key ***\n");
			    }

			    cram_free_block(tag);
			    // skip decoding of tag data itself and hope it's
			    // in an external block.
			}
		    }

		    if (!(bf & BAM_FUNMAP)) {
			r = c->comp_hdr->FN_codec->decode(s,c->comp_hdr->FN_codec, b, (char *)&fn, &out_sz);
			printf("FN = %d (ret %d, out_sz %d)\n", fn, r, out_sz);

			prev_pos = 0;
			for (f = 0; f < fn; f++) {
			    char op;
			    int32_t pos;

			    r = c->comp_hdr->FC_codec->decode(s,c->comp_hdr->FC_codec, b, &op, &out_sz);
			    printf("  %d: FC = %c (ret %d, out_sz %d)\n", f, op, r, out_sz);

			    r = c->comp_hdr->FP_codec->decode(s,c->comp_hdr->FP_codec, b, (char *)&pos, &out_sz);
			    printf("  %d: FP = %d+%d = %d (ret %d, out_sz %d)\n", f, pos, prev_pos, pos + prev_pos, r, out_sz);

			    pos += prev_pos;
			    prev_pos = pos;

			    switch(op) {
			    case 'S': { // soft clip: IN
				char dat[100];
				int32_t out_sz2 = 1;

				dat[0]='?';dat[1]=0;
				if (c->comp_hdr->IN_codec) {
				    r = c->comp_hdr->IN_codec->decode(s,c->comp_hdr->IN_codec, b, dat, &out_sz2);
				    printf("  %d: IN(S) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, dat, r, out_sz2);
				}
				break;
			    }

			    case 'X': { // Substitution; BS
				char bs;
				r = c->comp_hdr->BS_codec->decode(s,c->comp_hdr->BS_codec, b, &bs, &out_sz);
				printf("  %d: BS = %d (ret %d)\n", f, bs, r);
				break;
			    }

			    case 'D': { // Deletion; DL
				r = c->comp_hdr->DL_codec->decode(s,c->comp_hdr->DL_codec, b, (char *)&i32, &out_sz);
				printf("  %d: DL = %d (ret %d)\n", f, i32, r);
				break;
			    }

			    case 'I': { // Insertion (several bases); IN
				char dat[100];
				int32_t out_sz2 = 1;

				dat[0]='?';dat[1]=0;
				r = c->comp_hdr->IN_codec->decode(s,c->comp_hdr->IN_codec, b, dat, &out_sz2);
				printf("  %d: IN(I) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, dat, r, out_sz2);
				break;
			    }

			    case 'i': { // Insertion (single base); BA
				char cc;
				r = c->comp_hdr->BA_codec->decode(s,c->comp_hdr->BA_codec, b, &cc, &out_sz);
				printf("  %d: BA = %c (ret %d)\n", f, cc, r);
				break;
			    }

			    case 'B': { // Read base; BA, QS
				char cc, qc;
				r  = c->comp_hdr->BA_codec->decode(s,c->comp_hdr->BA_codec, b, &cc, &out_sz);
				r |= c->comp_hdr->QS_codec->decode(s,c->comp_hdr->QS_codec, b, &qc, &out_sz);
				printf("  %d: BA/QS(B) = %c/%d (ret %d)\n", f, cc, qc, r);
				break;
			    }

			    case 'Q': { // Quality score; QS
				char qc;
				r = c->comp_hdr->QS_codec->decode(s,c->comp_hdr->QS_codec, b, &qc, &out_sz);
				printf("  %d: QS = %d (ret %d)\n", f, qc, r);
			    }

			    default:
				abort();
			    }
			}

			r = c->comp_hdr->MQ_codec->decode(s,c->comp_hdr->MQ_codec, b, (char *)&i32, &out_sz);
			printf("MQ = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

			if (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES) {
			    char dat[1024];
			    int32_t out_sz2 = rl, i;

			    dat[0]='?';dat[1]=0;
			    r = c->comp_hdr->Qs_codec->decode(s,c->comp_hdr->Qs_codec, b, dat, &out_sz2);
			    for (i = 0; i < rl; i++)
				dat[i] += '!';
			    printf("Qs = %.*s (ret %d, out_sz %d)\n", out_sz2, dat, r, out_sz2);
			}
		    } else {
			puts("Unmapped");
			char dat[1024];
			int len = rl;

			do {
			    int32_t out_sz2 = len > 1024 ? 1024 : len;
			    r = c->comp_hdr->BA_codec->decode(s, c->comp_hdr->BA_codec, b, dat, &out_sz2);
			    printf("SQ = %.*s (out_sz %d)\n", out_sz2, dat, out_sz2);
			    len -= 1024;
			} while (len > 0);

			if (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES) {
			    int len = rl, i;

			    do {
				int32_t out_sz2 = len > 1024 ? 1024 : len;
				r = c->comp_hdr->Qs_codec->decode(s, c->comp_hdr->Qs_codec, b, dat, &out_sz2);
				for (i = 0; i < len; i++)
				    dat[i] += '!';
				printf("Qs = %.*s (out_sz %d)\n", out_sz2, dat, out_sz2);
				len -= 1024;
			    } while (len > 0);
			}
		    }
		}
	    }

	    for (id = 0; id < s->hdr->num_blocks; id++) {
		cram_block *b = s->block[id];
		printf("\n\tBlock %d/%d\n", id+1, s->hdr->num_blocks);
		printf("\t    Size:         %d comp / %d uncomp\n",
		       b->comp_size, b->uncomp_size);
		printf("\t    Method:       %s\n",
		       cram_block_method2str(b->orig_method));
		printf("\t    Content type: %s\n",
		       cram_content_type2str(b->content_type));
		printf("\t    Content id:   %d\n", b->content_id);

		if (b->method != RAW)
		    cram_uncompress_block(b);

		if (b->content_type == CORE) {
		    dump_core_block(b, verbose);
		} else {
		    switch (b->content_id) {
		    case 0:
			dump_seq_block(b, verbose);
			break;
			
		    case 1:
			dump_quality_block(b, verbose);
			break;
			
		    case 2:
			dump_name_block(b, verbose);
			break;
			
		    case 3:
			dump_mate_info_block(b, verbose);
			break;
			
		    case 4:
			dump_tag_block(b, verbose);
			break;

		    default:
			dump_core_block(b, verbose); // hex
			break;
		    }
		}
	    }

	    cram_free_slice(s);
	}

	pos = ftello(fd->fp);
	assert(pos == hpos + c->length);

	cram_free_container(c);
    }

    cram_close(fd);

    {
	int id;
	puts("");
	for (id = 0; id < bmax; id++)
	    printf("Block %d, total size %d\n", id, bsize[id]);
    }

    return 0;
}

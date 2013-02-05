/*
 * A debugging program to dump out information on the layout of a CRAM file.
 */

#include <stdio.h>
#include <assert.h>

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

void dump_core_block(cram_block *b) {
    int i;

    printf("Data = {");
    for (i = 0; i < 16 && i < b->uncomp_size; i++) {
	printf(i ? ", %02x" : "%02x", (unsigned char)b->data[i]);
    }
    if (i < b->uncomp_size)
	printf(", ...}\n");
    else
	printf("}\n");
}

void dump_seq_block(cram_block *b) {
    printf("%.*s\n", b->uncomp_size, b->data);
}

void dump_quality_block(cram_block *b) {
    int i;
    for (i = 0; i < b->uncomp_size; i++) {
	putchar(b->data[i] + '!');
    }
    putchar('\n');
}

void dump_name_block(cram_block *b) {
    printf("%.*s\n", b->uncomp_size, b->data);
}

void dump_mate_info_block(cram_block *b) {
    return dump_core_block(b);
}

void dump_tag_block(cram_block *b) {
    return dump_core_block(b);
}

int main(int argc, char **argv) {
    cram_fd *fd;
    cram_container *c;
    size_t pos, pos2;

    static int bsize[100], bmax = 0;

    if (argc != 2) {
	fprintf(stderr, "Usage: cram_dump filename.cram\n");
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
	   fd->SAM_hdr->length, fd->SAM_hdr->header);

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
	printf("    No. recs:          %d\n", c->num_records);
	printf("    No. blocks:        %d\n", c->num_blocks);
	printf("    No. landmarks:     %d\n", c->num_landmarks);

	printf("    {");
	for (i = 0; i < c->num_landmarks; i++) {
	    printf(i ? ", %d" : "%d", c->landmark[i]);
	}
	printf("}\n");

	printf("\n    Preservation map:\n");
	HashTableDump(c->comp_hdr->preservation_map, stdout, "\t");

	printf("\n    Record encoding map:\n");
	HashTableDumpMap(c->comp_hdr->rec_encoding_map, stdout, "\t", c->comp_hdr_block->data);

	printf("\n    Tag encoding map:\n");
	HashTableDumpMap(c->comp_hdr->tag_encoding_map, stdout, "\t",  c->comp_hdr_block->data);
	


	for (j = 0; j < c->num_landmarks; j++) {
	    cram_slice *s;
	    int id;
	    
	    pos2 = ftello(fd->fp);
	    assert(pos2 - pos - c->offset == c->landmark[j]);

	    s = cram_read_slice(fd);
	    printf("\n    Slice %d/%d, container offset %d\n", j+1, c->num_landmarks, (int)(pos2 - pos - c->offset));
	    printf("\tSlice content type %s\n",
		   cram_content_type2str(s->hdr->content_type));

	    if (s->hdr->content_type == MAPPED_SLICE) {
		printf("\tSlice ref seq    %d\n", s->hdr->ref_seq_id);
		printf("\tSlice ref start  %d\n", s->hdr->ref_seq_start);
		printf("\tSlice ref span   %d\n", s->hdr->ref_seq_span);
	    }
	    printf("\tNo. records      %d\n", s->hdr->num_records);
	    printf("\tNo. blocks       %d\n", s->hdr->num_blocks);
	    printf("\tBlk IDS:         {");
	    for (id = 0; id < s->hdr->num_blocks; id++) {
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
	    {
		cram_block *b = s->block[0];
		int32_t i32, bf, cf, fn, prev_pos, rl;
		int out_sz, r, f;
		block_t *blk = block_create((unsigned char *)b->data, b->uncomp_size);
		int rec;

		blk->bit = 7; // MSB first
		assert(b->content_type == CORE);

		for (rec = 0; rec < s->hdr->num_records; rec++) {
		    printf("Rec %d/%d\n", rec+1, s->hdr->num_records);

		    out_sz = 1; /* decode 1 item */
		    r = c->comp_hdr->BF_codec->decode(s,c->comp_hdr->BF_codec, blk, (char *)&bf, &out_sz);
		    printf("BF = %d (ret %d, out_sz %d)\n", bf, r, out_sz);

		    r = c->comp_hdr->CF_codec->decode(s,c->comp_hdr->CF_codec, blk, (char *)&cf, &out_sz);
		    printf("CF = %d (ret %d, out_sz %d)\n", cf, r, out_sz);

		    r = c->comp_hdr->RL_codec->decode(s,c->comp_hdr->RL_codec, blk, (char *)&rl, &out_sz);
		    printf("RL = %d (ret %d, out_sz %d)\n", rl, r, out_sz);

		    r = c->comp_hdr->AP_codec->decode(s,c->comp_hdr->AP_codec, blk, (char *)&i32, &out_sz);
		    printf("AP = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

		    r = c->comp_hdr->RG_codec->decode(s,c->comp_hdr->RG_codec, blk, (char *)&i32, &out_sz);
		    printf("RG = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

		    if (c->comp_hdr->read_names_included) {
			char dat[100];
			int32_t out_sz2 = 1;

			dat[0]='?';dat[1]=0;
			r = c->comp_hdr->RN_codec->decode(s,c->comp_hdr->RN_codec, blk, dat, &out_sz2);
			printf("RN = %.*s (ret %d, out_sz %d)\n", out_sz2, dat, r, out_sz2);
		    }

		    if (cf & CRAM_FLAG_DETACHED) {
			puts("Detached");
			/* MF, RN if !captureReadNames, NS, NP, IS */
			r = c->comp_hdr->MF_codec->decode(s,c->comp_hdr->MF_codec, blk, (char *)&i32, &out_sz);
			printf("MF = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

			if (!c->comp_hdr->read_names_included) {
			    char dat[100];
			    int32_t out_sz2 = 1;

			    dat[0]='?';dat[1]=0;
			    r = c->comp_hdr->RN_codec->decode(s,c->comp_hdr->RN_codec, blk, dat, &out_sz2);
			    printf("RN = %.*s (ret %d, out_sz %d)\n", out_sz2, dat, r, out_sz2);
			}
		    
			r = c->comp_hdr->NS_codec->decode(s,c->comp_hdr->NS_codec, blk, (char *)&i32, &out_sz);
			printf("NS = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

			r = c->comp_hdr->NP_codec->decode(s,c->comp_hdr->NP_codec, blk, (char *)&i32, &out_sz);
			printf("NP = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

			r = c->comp_hdr->TS_codec->decode(s,c->comp_hdr->TS_codec, blk, (char *)&i32, &out_sz);
			printf("TS = %d (ret %d, out_sz %d)\n", i32, r, out_sz);
		    } else if (cf & CRAM_FLAG_MATE_DOWNSTREAM) {
			puts("Not detached, and mate is downstream");
			r = c->comp_hdr->NF_codec->decode(s,c->comp_hdr->NF_codec, blk, (char *)&i32, &out_sz);
			printf("NF = %d (ret %d, out_sz %d)\n", i32, r, out_sz);
		    }

		    r = c->comp_hdr->TC_codec->decode(s,c->comp_hdr->TC_codec, blk, (char *)&i32, &out_sz);
		    printf("TC = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

		    if (!(bf & CRAM_FUNMAP)) {
			r = c->comp_hdr->FN_codec->decode(s,c->comp_hdr->FN_codec, blk, (char *)&fn, &out_sz);
			printf("FN = %d (ret %d, out_sz %d)\n", fn, r, out_sz);

			prev_pos = 0;
			for (f = 0; f < fn; f++) {
			    int32_t op, pos;

			    r = c->comp_hdr->FC_codec->decode(s,c->comp_hdr->FC_codec, blk, (char *)&op, &out_sz);
			    printf("  %d: FC = %c (ret %d, out_sz %d)\n", f, op, r, out_sz);

			    r = c->comp_hdr->FP_codec->decode(s,c->comp_hdr->FP_codec, blk, (char *)&pos, &out_sz);
			    printf("  %d: FP = %d+%d = %d (ret %d, out_sz %d)\n", f, pos, prev_pos, pos + prev_pos, r, out_sz);

			    pos += prev_pos;
			    prev_pos = pos;

			    switch(op) {
			    case 'S': { // soft clip: IN
				char dat[100];
				int32_t out_sz2 = 1;

				dat[0]='?';dat[1]=0;
				if (c->comp_hdr->IN_codec) {
				    r = c->comp_hdr->IN_codec->decode(s,c->comp_hdr->IN_codec, blk, dat, &out_sz2);
				    printf("  %d: IN(S) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, dat, r, out_sz2);
				}
				break;
			    }

			    case 'X': { // Substitution; BS
				r = c->comp_hdr->BS_codec->decode(s,c->comp_hdr->BS_codec, blk, (char *)&i32, &out_sz);
				printf("  %d: BS = %d (ret %d)\n", f, i32, r);
				break;
			    }

			    case 'D': { // Deletion; DL
				r = c->comp_hdr->DL_codec->decode(s,c->comp_hdr->DL_codec, blk, (char *)&i32, &out_sz);
				printf("  %d: DL = %d (ret %d)\n", f, i32, r);
				break;
			    }

			    case 'I': { // Insertion (several bases); IN
				char dat[100];
				int32_t out_sz2 = 1;

				dat[0]='?';dat[1]=0;
				r = c->comp_hdr->IN_codec->decode(s,c->comp_hdr->IN_codec, blk, dat, &out_sz2);
				printf("  %d: IN(I) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, dat, r, out_sz2);
				break;
			    }

			    case 'i': { // Insertion (single base); BA
				r = c->comp_hdr->BA_codec->decode(s,c->comp_hdr->BA_codec, blk, (char *)&i32, &out_sz);
				printf("  %d: BA = %c (ret %d)\n", f, i32, r);
				break;
			    }

			    case 'B': { // Read base; BA, QS
				int32_t qc;
				r  = c->comp_hdr->BA_codec->decode(s,c->comp_hdr->BA_codec, blk, (char *)&i32, &out_sz);
				r |= c->comp_hdr->QS_codec->decode(s,c->comp_hdr->QS_codec, blk, (char *)&qc, &out_sz);
				printf("  %d: BA/QS(B) = %c/%d (ret %d)\n", f, i32, qc, r);
				break;
			    }

			    case 'Q': { // Quality score; QS
				int32_t qc;
				r = c->comp_hdr->QS_codec->decode(s,c->comp_hdr->QS_codec, blk, (char *)&qc, &out_sz);
				printf("  %d: QS = %d (ret %d)\n", f, qc, r);
			    }

			    default:
				abort();
			    }
			}

			r = c->comp_hdr->MQ_codec->decode(s,c->comp_hdr->MQ_codec, blk, (char *)&i32, &out_sz);
			printf("MQ = %d (ret %d, out_sz %d)\n", i32, r, out_sz);

			if (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES) {
			    char dat[100];
			    int32_t out_sz2 = rl, i;

			    dat[0]='?';dat[1]=0;
			    r = c->comp_hdr->Qs_codec->decode(s,c->comp_hdr->Qs_codec, blk, dat, &out_sz2);
			    for (i = 0; i < rl; i++)
				dat[i] += '!';
			    printf("Qs = %.*s (ret %d, out_sz %d)\n", out_sz2, dat, r, out_sz2);
			}
		    } else {
			puts("Unmapped");
		    }
		}

		block_destroy(blk, 1);
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
		    dump_core_block(b);
		} else {
		    switch (b->content_id) {
		    case 0:
			dump_seq_block(b);
			break;
			
		    case 1:
			dump_quality_block(b);
			break;
			
		    case 2:
			dump_name_block(b);
			break;
			
		    case 3:
			dump_mate_info_block(b);
			break;
			
		    case 4:
			dump_tag_block(b);
			break;
		    }
		}
	    }

	    cram_free_slice(s);
	}

	cram_free_container(c);
	pos = ftello(fd->fp);
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

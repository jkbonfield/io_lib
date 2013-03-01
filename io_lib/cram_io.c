/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 *
 * CRAM I/O primitives.
 *
 * - ITF8 encoding and decoding.
 * - Block based I/O
 * - Zlib inflating and deflating (memory)
 * - CRAM basic data structure reading and writing
 * - File opening / closing
 * - Reference sequence handling
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

/* ----------------------------------------------------------------------
 * ITF8 encoding and decoding.
 *
 * Also see the itf8_get and itf8_put macros in cram_io.h
 */

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
 * Encodes and writes a single integer in ITF-8 format.
 * Returns 0 on success
 *        -1 on failure
 */
int itf8_encode(cram_fd *fd, int32_t val) {
    char buf[5];
    int len = itf8_put(buf, val);
    return fwrite(buf, 1, len, fd->fp) == len ? 0 : -1;
}

#ifndef ITF8_MACROS
/*
 * As above, but decoding from memory
 */
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
#endif

/*
 * Pushes a value in ITF8 format onto the end of a block.
 * This shouldn't be used for high-volume data as it is not the fastest
 * method.
 *
 * Returns the number of bytes written
 */
int itf8_put_blk(cram_block *blk, int val) {
    char buf[5];
    int sz;

    sz = itf8_put(buf, val);
    BLOCK_APPEND(blk, buf, sz);
    return sz;
}

/* ----------------------------------------------------------------------
 * zlib compression code - from Gap5's tg_iface_g.c
 * They're static here as they're only used within the cram_compress_block
 * and cram_uncompress_block functions, which are the external interface.
 */

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

/* ----------------------------------------------------------------------
 * CRAM blocks - the dynamically growable data block. We have code to
 * create, update, (un)compress and read/write.
 *
 * These are derived from the deflate_interlaced.c blocks, but with the
 * CRAM extension of content types and IDs.
 */

/*
 * Allocates a new cram_block structure with a specified content_type and
 * id.
 *
 * Returns block pointer on success
 *         NULL on failure
 */
cram_block *cram_new_block(enum cram_content_type content_type,
			   int content_id) {
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

/*
 * Frees a CRAM block, deallocating internal data too.
 */
void cram_free_block(cram_block *b) {
    if (!b)
	return;
    if (b->data)
	free(b->data);
    free(b);
}

/*
 * Uncompresses a CRAM block, if compressed.
 */
void cram_uncompress_block(cram_block *b) {
    char *uncomp;
    size_t uncomp_size = 0;

    switch (b->method) {
    case RAW:
	b->uncomp_size = b->comp_size;
	return;

    case GZIP:
	uncomp = zlib_mem_inflate((char *)b->data, b->comp_size, &uncomp_size);
	assert((int)uncomp_size == b->uncomp_size);
	free(b->data);
	b->data = (unsigned char *)uncomp;
	b->method = RAW;
	break;

    case BZIP2:
	fprintf(stderr, "Bzip2 compression not yet implemented\n");
	abort();
	break;
    }
}

/*
 * Compresses a block using one of two different zlib strategies. If we only
 * want one choice set strat2 to be -1.
 *
 * The logic here is that sometimes Z_RLE does a better job than Z_FILTERED
 * or Z_DEFAULT_STRATEGY on quality data. If so, we'd rather use it as it is
 * significantly faster.
 */
void cram_compress_block(cram_fd *fd, cram_block *b, cram_metrics *metrics,
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
	if (fd->verbose > 1)
	    fprintf(stderr, "metrics trial %d, next_trial %d, m1 %d, m2 %d\n",
		    metrics->trial, metrics->next_trial,
		    metrics->m1, metrics->m2);

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

	c1 = zlib_mem_deflate((char *)b->data, b->uncomp_size,
			      &s1, level, strat);
	c2 = zlib_mem_deflate((char *)b->data, b->uncomp_size,
			      &s2, level2, strat2);
	if (s1 < s2) {
	    if (fd->verbose > 1)
		fprintf(stderr, "M1 wins\n");
	    comp = c1; comp_size = s1;
	    free(c2);
	    metrics->m1++;
	} else {
	    if (fd->verbose > 1)
		fprintf(stderr, "M2 wins\n");
	    comp = c2; comp_size = s2;
	    free(c1);
	    metrics->m2++;
	}
    } else if (strat2 >= 0) {
	comp = zlib_mem_deflate((char *)b->data, b->uncomp_size, &comp_size,
				metrics->m1 > metrics->m2 ? level : level2,
				metrics->m1 > metrics->m2 ? strat : strat2);
    } else {
	comp = zlib_mem_deflate((char *)b->data, b->uncomp_size, &comp_size,
				level, strat);
    }

    free(b->data);
    b->data = (unsigned char *)comp;
    b->method = GZIP;
    b->comp_size = comp_size;

    if (fd->verbose)
	fprintf(stderr, "Compressed block ID %d from %d to %d\n",
		b->content_id, b->uncomp_size, b->comp_size);
}

cram_metrics *cram_new_metrics(void) {
    cram_metrics *m = malloc(sizeof(*m));
    m->m1 = m->m2 = 0;
    m->trial = 2;
    m->next_trial = 100;
    return m;
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

/* ----------------------------------------------------------------------
 * Reference sequence handling
 */

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

static void free_refs(refs *r) {
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

void cram_load_reference(cram_fd *fd, char *fn) {
    fd->refs = load_reference(fn);
    refs2id(fd->refs, fd->SAM_hdr);
}

/* ----------------------------------------------------------------------
 * Containers
 */

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

    c->comp_hdr_block = NULL;

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
 * Reads a container header.
 *
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

    c->slices = NULL;
    c->curr_slice = 0;
    c->max_slice = c->num_landmarks;
    c->curr_rec = 0;
    c->max_rec = 0;

    return c;
}

/*
 * Writes a container structure.
 *
 * Returns 0 on success
 *        -1 on failure
 */
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

/*
 * Flushes a completely or partially full container to disk, writing
 * container structure, header and blocks. This also calls the encoder
 * functions.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_flush_container(cram_fd *fd, cram_container *c) {
    int i, j;

    /* Encode the container blocks and generate compression header */
    if (0 != cram_encode_container(fd, c))
	return -1;

    /* Write the container struct itself */
    if (0 != cram_write_container(fd, c))
	return -1;

    /* And the compression header */
    cram_write_block(fd, c->comp_hdr_block);

    /* Followed by the slice blocks */
    for (i = 0; i < c->curr_slice; i++) {
	cram_slice *s = c->slices[i];

	cram_write_block(fd, s->hdr_block);

	for (j = 0; j < s->hdr->num_blocks; j++) {
	    if (-1 == cram_write_block(fd, s->block[j]))
		return -1;
	}
    }

    return fflush(fd->fp) == 0 ? 0 : -1;
}


/* ----------------------------------------------------------------------
 * Compression headers; the first part of the container
 */

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


/* ----------------------------------------------------------------------
 * Slices and slice headers
 */

void cram_free_slice_header(cram_block_slice_hdr *hdr) {
    if (!hdr)
	return;

    if (hdr->block_content_ids)
	free(hdr->block_content_ids);

    free(hdr);

    return;
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

    s->seqs_blk = cram_new_block(EXTERNAL, 0);
    s->qual_blk = cram_new_block(EXTERNAL, 1);
    s->name_blk = cram_new_block(EXTERNAL, 2);
    s->aux_blk  = cram_new_block(EXTERNAL, 4);
    s->base_blk = cram_new_block(EXTERNAL, 0);
#ifdef TN_external
    s->tn_blk   = cram_new_block(EXTERNAL, 6);
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


/* ----------------------------------------------------------------------
 * CRAM file definition (header)
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

    fd->first_container += 26;

    return def;
}

/*
 * Writes a cram_file_def structure to cram_fd.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_file_def(cram_fd *fd, cram_file_def *def) {
    return (fwrite(&def->magic[0], 1, 26, fd->fp) == 26) ? 0 : -1;
}

void cram_free_file_def(cram_file_def *def) {
    if (def) free(def);
}

/* ----------------------------------------------------------------------
 * SAM header I/O
 */

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
cram_SAM_hdr *cram_new_SAM_hdr(char *str, size_t len) {
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

    fd->first_container += 4 + hdr->header_len;

    /* Parse */
    bam_parse_header(hdr);
    
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

/* ----------------------------------------------------------------------
 * The top-level cram opening, closing and option handling
 */

/*
 * Initialises the lookup tables. These could be global statics, but they're
 * clumsy to setup in a multi-threaded environment unless we generate
 * verbatim code and include that.
 */
static void cram_init_tables(cram_fd *fd) {
    int i;

    memset(fd->L1, 4, 256);
    fd->L1['A'] = 0; fd->L1['a'] = 0;
    fd->L1['C'] = 1; fd->L1['c'] = 1;
    fd->L1['G'] = 2; fd->L1['g'] = 2;
    fd->L1['T'] = 3; fd->L1['t'] = 3;

    memset(fd->L2, 5, 256);
    fd->L2['A'] = 0; fd->L2['a'] = 0;
    fd->L2['C'] = 1; fd->L2['c'] = 1;
    fd->L2['G'] = 2; fd->L2['g'] = 2;
    fd->L2['T'] = 3; fd->L2['t'] = 3;
    fd->L2['N'] = 4; fd->L2['n'] = 4;

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

	fd->bam_flag_swap[i]  = f;
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

	fd->cram_flag_swap[i] = g;
    }

    memset(fd->cram_sub_matrix, 4, 32*32);
    for (i = 0; i < 32; i++) {
	fd->cram_sub_matrix[i]['A'&0x1f]=0;
	fd->cram_sub_matrix[i]['C'&0x1f]=1;
	fd->cram_sub_matrix[i]['G'&0x1f]=2;
	fd->cram_sub_matrix[i]['T'&0x1f]=3;
	fd->cram_sub_matrix[i]['N'&0x1f]=4;
    }
    for (i = 0; i < 20; i+=4) {
	int j;
	for (j = 0; j < 20; j++) {
	    fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
	    fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
	    fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
	    fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
	}
	fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+0]&0x1f]=0;
	fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+1]&0x1f]=1;
	fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+2]&0x1f]=2;
	fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+3]&0x1f]=3;
    }
}

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
    fd->first_container = 0;

    cram_init_tables(fd);

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
    fd->verbose = 0;
    fd->seqs_per_slice = SEQS_PER_SLICE;
    fd->slices_per_container = SLICE_PER_CNT;

    fd->index = NULL;

    for (i = 0; i < 7; i++)
	fd->m[i] = cram_new_metrics();

    fd->range.refid = -2; // no ref.

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
    int i;

    if (!fd)
	return -1;

    if (fd->mode == 'w' && fd->ctr) {
	if(fd->ctr->slice)
	    fd->ctr->curr_slice++;
	if (-1 == cram_flush_container(fd, fd->ctr))
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

    if (fd->index)
	cram_index_free(fd);

    free(fd);
    return 0;
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

    case CRAM_OPT_SEQS_PER_SLICE:
	fd->seqs_per_slice = val->i;
	break;

    case CRAM_OPT_SLICES_PER_CONTAINER:
	fd->slices_per_container = val->i;
	break;

    case CRAM_OPT_RANGE:
	fd->range = *(cram_range *)val->s;
	cram_seek_to_refpos(fd, &fd->range);
	break;

    default:
	fprintf(stderr, "Unknown CRAM option code %d\n", opt);
	return -1;
    }

    return 0;
}

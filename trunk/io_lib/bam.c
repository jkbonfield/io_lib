/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2010-3
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <fcntl.h>
#include <zlib.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stddef.h>

#include "io_lib/bam.h"
#include "io_lib/os.h"

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

/* Macros to store integers of various sizes in little endian byte order.
 * The value is put in the location pointed to by ucp, which should be 
 * an unsigned char pointer.  ucp is incremented by the size of the
 * stored value. */

#define STORE_UINT16(ucp, val) \
    *(ucp)++ = ((uint16_t) val)      & 0xff; \
    *(ucp)++ = ((uint16_t) val >> 8) & 0xff;

#define STORE_UINT32(ucp, val) \
    *(ucp)++ = ((uint32_t) (val))       & 0xff;	\
    *(ucp)++ = ((uint32_t) (val) >>  8) & 0xff;	\
    *(ucp)++ = ((uint32_t) (val) >> 16) & 0xff;	\
    *(ucp)++ = ((uint32_t) (val) >> 24) & 0xff;

#define STORE_UINT64(ucp, val) \
    *(ucp)++ = ((uint64_t) (val))       & 0xff;	\
    *(ucp)++ = ((uint64_t) (val) >>  8) & 0xff;	\
    *(ucp)++ = ((uint64_t) (val) >> 16) & 0xff;	\
    *(ucp)++ = ((uint64_t) (val) >> 24) & 0xff; \
    *(ucp)++ = ((uint64_t) (val) >> 32) & 0xff; \
    *(ucp)++ = ((uint64_t) (val) >> 40) & 0xff; \
    *(ucp)++ = ((uint64_t) (val) >> 48) & 0xff; \
    *(ucp)++ = ((uint64_t) (val) >> 56) & 0xff;

static int bam_more_input(bam_file_t *b);
static int bam_more_output(bam_file_t *b);
static int reg2bin(int start, int end);
static int bgzf_write(int fd, int level, const void *buf, size_t count);

/*
 * Reads len bytes from fp into data.
 *
 * Returns the number of bytes read.
 *         0 for eof.
 *        -1 for failure.
 */
static int bam_read(bam_file_t *b, void *data, size_t len) {
    int nb = 0, n;
    unsigned char *cdata = data;

    while (len) {
	/* Consume any available uncompressed output */
	if (b->out_sz) {
	    size_t l = MIN(b->out_sz, len);
	    memcpy(cdata, b->out_p, l);
	    b->out_p += l;
	    b->out_sz -= l;
	    cdata += l;
	    len -= l;
	    nb += l;

	    if (!len)
		return nb;
	}

	if (!b->gzip) {
	    /* Already uncompressed, so easy to deal with */
	    if (!b->in_sz)
		if (-1 == bam_more_input(b))
		    return nb ? nb : 0;
		    
	    b->out_p  = b->in_p;
	    b->out_sz = b->in_sz;
	    b->in_sz  = 0;
	    continue;
	}

	n = bam_more_output(b);
	if (n == -1)
	    return -1;
	if (n == 0)
	    return nb;
    }

    return nb;
}

/*
 * Reads a line of text of unknown length.
 * 'str' is both input and output. If *str == NULL then memory is allocated
 * for the line. If *str != NULL then it is expected to point to an existing
 * block of memory that we can write into and realloc as required.
 *
 * Similarly *len is both input and output. It is expected to hold the
 * current allocated size of *str. It is modified if we realloc it.
 *
 * Lines have the \n removed and will be null terminated.
 *
 * Returns actual line length used (note note the same as *len) on success
 *        -1 on failure
 */
int bam_get_line(bam_file_t *b, unsigned char **str, size_t *len) {
    unsigned char *buf = *str;
    int used_l = 0;
    size_t alloc_l = *len;
    int next_condition, r = 0;

    while (b->out_sz || (r=bam_more_output(b)) > 0) {
	int tmp;
	unsigned char *from = b->out_p;
	unsigned char *to   = &buf[used_l];

	/*
	 * Next condition is the number of loop iterations before something
	 * has to be done - either getting more uncompressed output or
	 * resizing the buffer. We don't care which, but it allows us to
	 * have just one check per loop instead of two. Once out of the loop
	 * we can then afford to determine which case is and deal with it.
	 */
	tmp = next_condition = MIN(b->out_sz, alloc_l-used_l);

	while (next_condition-- > 0) { /* these 3 lines are 50% of SAM cpu */
	    if (*from != '\n') {
		*to++ = *from++;
	    } else {
		b->out_p = from;
		used_l = to-buf;
		b->out_p++;
		buf[used_l] = 0;
		*str = buf;
		*len = alloc_l;
		b->out_sz -= (tmp - next_condition);
		return used_l;
	    }
	}

	used_l = to-buf;
	b->out_p = from;
	b->out_sz -= tmp;

	if (used_l >= alloc_l) {
	    alloc_l = alloc_l ? alloc_l * 2 : 1024;
	    if (NULL == (buf = realloc(buf, alloc_l)))
		return -1;
	}
    }

    if (r == -1)
	return -1;
    return b->out_sz ? -1 : 0;
}

static int load_bam_header(bam_file_t *b) {
    char magic[4], *header;
    int i;
    int32_t header_len, nref;

    if (4 != bam_read(b, magic, 4))
	return -1;
    if (memcmp(magic, "BAM\x01",4) != 0)
	return -1;
    if (4 != bam_read(b, &header_len, 4))
	return -1;
    header_len = le_int4(header_len);
    if (!(header = malloc(header_len+1)))
	return -1;
    *header = 0;
    if (header_len != bam_read(b, header, header_len))
	return -1;

    if (!(b->header = sam_hdr_parse(header, header_len)))
	return -1;
    free(header);

    /* Load the reference data and check it matches the parsed header */
    if (4 != bam_read(b, &nref, 4))
	return -1;
    nref = le_int4(nref);
    if (b->header->nref != nref && b->header->nref) {
	fprintf(stderr, "Error: @RG lines are at odds with "
		"binary encoded reference data\n");
	return -1;
    }

    for (i = 0; i < nref; i++) {
	uint32_t nlen, len;
	char name_a[1024], *name;;

	if (4 != bam_read(b, &nlen, 4))
	    return -1;
	name = nlen < 1024 ? name_a : malloc(nlen);
	if (!name)
	    return -1;
	nlen = le_int4(nlen);
	if (nlen != bam_read(b, name, nlen))
	    return -1;

	if (4 != bam_read(b, &len, 4))
	    return -1;
	len = le_int4(len);

	if (i < b->header->nref && b->header->ref[i].name) {
	    if (strcmp(b->header->ref[i].name, name)) {
		fprintf(stderr, "Error: @RG lines are at odds with "
			"binary encoded reference data\n");
		return -1;
	    }

	    if (b->header->ref[i].len != len) {
		fprintf(stderr, "Error: @RG lines are at odds with "
			"binary encoded reference data\n");
		return -1;
	    }
	} else {
	    char len_c[100];
	    sprintf(len_c, "%d", len);
	    if (sam_hdr_add(b->header, "SQ", "SN", name, "LN", len_c, NULL))
		return -1;
	}

	if (name != name_a)
	    free(name);
    }

    b->line = 0; // FIXME

    return 0;
}

static int load_sam_header(bam_file_t *b) {
    unsigned char *str = NULL;
    size_t alloc = 0, len;
    dstring_t *header = dstring_create(NULL);;
    int r = 0;

    while ((b->out_sz > 0 || (r=bam_more_output(b)) > 0) && *b->out_p == '@') {
	b->line++;
	if ((len = bam_get_line(b, &str, &alloc)) == -1)
	    return -1;

	if (-1 == dstring_nappend(header, (char *)str, len))
	    return -1;
	if (-1 == dstring_append_char(header, '\n'))
	    return -1;
    }
    if (r == -1)
	return -1;
    b->line = 0; // FIXME

    if (!(b->header = sam_hdr_parse((char *)dstring_str(header),
				    dstring_length(header))))
	return -1;

    dstring_destroy(header);
    free(str);

    return 0;
}

/* --------------------------------------------------------------------------
 * 
 */

#ifndef O_BINARY
#    define O_BINARY 0
#endif

/*! Opens a SAM or BAM file.
 *
 * The mode parameter indicates the file
 * type (if not auto-detecting) and whether it is for reading or
 * writing. Use "rb" or "wb" for reading or writing BAM and "r" or
 * "w" or reading or writing SAM. When writing BAM, the mode may end
 * with a digit from 0 to 9 to indicate the compression to use with 0
 * indicating uncompressed data.
 *
 * @param fn The filename to open or create.
 * @param mode The input/output mode, similar to fopen().
 *
 * @return
 * Returns a bam_file_t pointer on success;
 *         NULL on failure.
 */
bam_file_t *bam_open(const char *fn, const char *mode) {
    bam_file_t *b = calloc(1, sizeof *b);
    
    if (!b)
	return NULL;

    b->in_p     = b->in;
    b->in_sz    = 0;
    b->out_p    = b->out;
    b->out_sz   = 0;
    b->next_len = -1;
    b->bs       = NULL;
    b->bs_size  = 0;
    b->z_finish = 1;
    b->bgzf     = 0;
    b->no_aux   = 0;
    b->line     = 0;
    b->binary   = 0;
    b->level    = Z_DEFAULT_COMPRESSION;
    b->sam_str  = NULL;

    /* Creation */
    if (*mode == 'w') {
	b->mode = O_WRONLY | O_TRUNC | O_CREAT;
	if (mode[1] == 'b') {
	    b->mode |= O_BINARY;
	    b->binary = 1;
	}
	if (mode[2] >= '0' && mode[2] <= '9')
	    b->level = mode[2] - '0';

	if (strcmp(fn, "-") == 0) {
	    b->fd = 1; /* Stdout */
	} else {
	    if (-1 == (b->fd = open(fn, b->mode, 0666)))
		goto error;
	}

	return b;
    }

    if (*mode != 'r')
	return NULL;

    if (strcmp(mode, "rb") == 0) {
	b->mode = O_RDONLY | O_BINARY;
    } else {
	b->mode = O_RDONLY;
    }
    if (strcmp(fn, "-") == 0) {
	b->fd = 0; /* Stdin */
    } else {
	if (-1 == (b->fd = open(fn, b->mode, 0)))
	    goto error;
    }

    /* Load first block so we can check */
    bam_more_input(b);
    if (b->in_sz < 2)
	return NULL;
    if (b->in_p[0] == 31 && b->in_p[1] == 139)
	b->gzip = 1;
    else
	b->gzip = 0;

    if (b->gzip) {
	/* Set up zlib */
	b->s.zalloc    = NULL;
	b->s.zfree     = NULL;
	b->s.opaque    = NULL;
	inflateInit2(&b->s, -15);
    }

    if (-1 == bam_more_output(b))
	return NULL;
    /* Auto-correct open file type if we detect a BAM */
    if (b->out_sz >= 3 && strncmp("BAM", (char *)b->out_p, 3) != 0) {
	b->mode &= ~O_BINARY;
	mode = "r";
    } else if (b->out_sz >= 3 && strncmp("BAM", (char *)b->out_p, 3) == 0) {
	b->mode |= O_BINARY;
	mode = "rb";
    }

    /* Load header */
    if (strcmp(mode, "rb") == 0) {
	if (-1 == load_bam_header(b))
	    goto error;
	b->bam = 1;
    } else {
	if (-1 == load_sam_header(b))
	    goto error;
	b->bam = 0;
    }

    return b;

 error:
    if (b) {
	if (b->header)
	    free(b->header);
	free(b);
    }

    return NULL;
}

int bam_close(bam_file_t *b) {
    int r;

    if (!b)
	return 0;

    if (b->mode & O_WRONLY) {
	if (b->binary) {
	    if (bgzf_write(b->fd, b->level, b->out, b->out_p - b->out)) {
		fprintf(stderr, "Write failed in bam_close()\n");
	    }

	    /* Output a blank BGZF block too to mark EOF */
	    if (28 != write(b->fd, "\037\213\010\4\0\0\0\0\0\377\6\0\102\103"
			    "\2\0\033\0\3\0\0\0\0\0\0\0\0\0", 28)) {
		fprintf(stderr, "Write failed in bam_close()\n");
	    }
	} else {
	    if (b->out_p - b->out != write(b->fd, b->out, b->out_p - b->out)) {
		fprintf(stderr, "Write failed in bam_close()\n");
	    }
	}
    }

    if (b->bs)
	free(b->bs);
    if (b->header)
	sam_hdr_free(b->header);

    if (b->gzip)
	inflateEnd(&b->s);

    if (b->sam_str)
	free(b->sam_str);

    r = close(b->fd);

    free(b);

    return r;
}

/*
 * Loads more data into the input (compressed) buffer.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
static int bam_more_input(bam_file_t *b) {
    size_t l;

    if (b->in != b->in_p) {
	memmove(b->in, b->in_p, b->in_sz);
	b->in_p = b->in;
    }

    l = read(b->fd, &b->in[b->in_sz], Z_BUFF_SIZE - b->in_sz);
    if (l <= 0)
	return -1;
    
    b->in_sz += l;
    return 0;
}

/*
 * Converts compressed input to the uncompressed output buffer
 *
 * Returns number of additional output bytes on success
 *         0 on eof
 *        -1 on failure.
 */
static int bam_more_output(bam_file_t *b) {
    int err = Z_OK;
    unsigned char *bgzf;
    int xlen, bsize;

    assert(b->out_sz == 0);

    if (!b->gzip) {
	/* Already uncompressed, so easy to deal with */
	if (!b->in_sz)
	    if (-1 == bam_more_input(b))
		return 0;
		    
	b->out_p  = b->in_p;
	b->out_sz = b->in_sz;
	b->in_sz  = 0;
	return b->out_sz;
    }

    /* Uncompress another BGZF block */
    /* BGZF header */
    if (b->in_sz < 18) {
	if (-1 == bam_more_input(b))
	    return 0;
	if (b->in_sz < 18)
	    return -1;
	if (b->in_sz == 0)
	    return 0; /* eof */
    }
	
    if (b->z_finish) {
	/*
	 * BGZF header is gzip + extra fields.
	 */
	bgzf = b->in_p;
	b->in_p += 10; b->in_sz -= 10;

	if (bgzf[0] != 31 || bgzf[1] != 139)
	    return -1; /* magic number failure */
	if ((bgzf[3] & 4) == 4) {
	    /* has extra fields, eg BGZF */
	    xlen = bgzf[10] + bgzf[11]*256;
	    b->in_p += 2; b->in_sz -= 2;
	} else {
	    xlen = 0;
	}
    } else {
	/* Continuing with an existing data stream */
	xlen = 0;
    }

    /* BGZF */
    if (xlen == 6) {
	b->bgzf = 1;
	b->in_p += 6; b->in_sz -= 6;

	if (bgzf[12] != 'B' || bgzf[13] != 'C' ||
	    bgzf[14] !=  2  || bgzf[15] !=  0)
	    return -1;
	bsize = bgzf[16] + bgzf[17]*256;
	bsize -= 6+19;

	/* Inflate */
	if (b->in_sz < bsize + 8) {
	    do {
		if (bam_more_input(b) == -1)
		    return -1; /* Truncated */
	    } while (b->in_sz < bsize + 8);
	}
	b->s.avail_in  = bsize;
	b->s.next_in   = b->in_p;
	b->s.avail_out = Z_BUFF_SIZE;
	b->s.next_out  = b->out;
	b->s.total_out = 0;
	    
	inflateReset(&b->s);
	err = inflate(&b->s, Z_FINISH);
	if (err != Z_STREAM_END) {
	    fprintf(stderr, "Inflate returned error code %d\n", err);
	    return -1;
	}
	b->z_finish = 1;

	b->in_p   += bsize + 8; /* crc & isize */
	b->in_sz  -= bsize + 8;
	b->out_sz  = b->s.total_out;
	b->out_p   = b->out;
    } else {
	/* Some other gzip variant, but possibly still having xlen */
	while (xlen) {
	    int d = MIN(b->in_sz, xlen);
	    xlen     -= d;
	    b->in_p  += d;
	    b->in_sz -= d;
	    if (b->in_sz == 0)
		bam_more_input(b);
	    if (b->in_sz == 0)
		return -1; /* truncated file */
	}
	    
	b->s.avail_in  = b->in_sz;
	b->s.next_in   = b->in_p;
	b->s.avail_out = Z_BUFF_SIZE;
	b->s.next_out  = b->out;
	b->s.total_out = 0;
	if (b->z_finish)
	    inflateReset(&b->s);
	    
	do {
	    err = inflate(&b->s, Z_BLOCK);
	    //printf("err=%d\n", err);

	    if (err == Z_OK || err == Z_STREAM_END) {
		b->in_p  += b->in_sz - b->s.avail_in;
		b->in_sz  = b->s.avail_in;
		b->out_sz = b->s.total_out;
		b->out_p  = b->out;
	    }

	    if (err == Z_STREAM_END) {
		b->z_finish = 1;
		/* Consume (ignore) CRC & ISIZE */
		if (b->in_sz < 8)
		    bam_more_input(b);

		if (b->in_sz < 8)
		    return -1; /* truncated file */

		b->in_sz -= 8;
		b->in_p  += 8;
	    } else {
		b->z_finish = 0;
	    }
	} while (err != Z_OK && err != Z_STREAM_END);
    }

    /*
     * Zero length blocks may not actually be EOF, just bizarre. We return
     * 0 elsewhere for the EOF case, so if we got here and b->out_sz is 0
     * then go around again.
     */
    return b->out_sz ? b->out_sz : bam_more_output(b);
}

/*
 * Decodes the next line of SAM into a bam_seq_t struct.
 *
 * Returns 1 on success
 *         0 on eof
 *        -1 on error
 */
static int sam_next_seq(bam_file_t *b, bam_seq_t **bsp) {
    int used_l, n, sign;
    unsigned char *cpf, *cpt, *cp;
    int cigar_len;
    bam_seq_t *bs;
    HashItem *hi;
    int start, end;
    SAM_hdr *sh = b->header;
    static const int lookup[256] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 00 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 10 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 20 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 30 */
	0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0, /* 40 */
	0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0, /* 50 */
	0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0, /* 60 */
	0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0, /* 70 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 80 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 90 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* a0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* b0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* c0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* d0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* e0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};/* f0 */

    /* Fetch a single line */
    if ((used_l = bam_get_line(b, &b->sam_str, &b->alloc_l)) <= 0) {
	return used_l;
    }

    used_l *= 4; // FIXME

    /* Over sized memory, for worst case? FIXME: cigar can break this! */
    if (!*bsp || used_l + sizeof(*bs) > (*bsp)->alloc) {
	if (!(*bsp = realloc(*bsp, used_l + sizeof(*bs))))
	    return -1;
	(*bsp)->alloc = used_l + sizeof(*bs);
	(*bsp)->blk_size = 0; /* compute later */
    }

    bs = *bsp;
    bs->flag_packed = 0;
    bs->bin_packed = 0;
    
    /* Decode line */
    cpf = b->sam_str;
    cpt = (unsigned char *)&bs->data;
    
    /* Name */
    cp = cpf;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	*cpt++ = *cpf++;
    *cpt++ = 0;
    if (!*cpf++) return -1;
    if (cpf-cp > 255) {
	fprintf(stderr, "SAM name length >= 256 characters are invalid\n");
	return -1;
    }
    bam_set_name_len(bs, cpf-cp);

    /* flag */
    n = 0;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bam_set_flag(bs, n);

    /* ref */
    cp = cpf;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	cpf++;
    if (*cp == '*') {
	/* Unmapped */
	bs->ref = -1;
    } else {
	hi = HashTableSearch(b->header->ref_hash, (char *)cp, cpf-cp);
	if (!hi) {
	    SAM_hdr *sh = b->header;
	    HashData hd;

	    fprintf(stderr, "Reference seq %.*s unknown\n", (int)(cpf-cp), cp);

	    /* Fabricate it instead */
	    sh->ref = realloc(sh->ref, (sh->nref+1)*sizeof(*sh->ref));
	    if (!sh->ref)
		return -1;
	    sh->ref[sh->nref].len  = 0; /* Unknown value */
	    sh->ref[sh->nref].name = malloc(cpf-cp+1);
	    if (!sh->ref[sh->nref].name)
		return -1;
	    memcpy(sh->ref[sh->nref].name, cp, cpf-cp);
	    sh->ref[sh->nref].name[cpf-cp] = 0;

	    hd.i = sh->nref;
	    hi = HashTableAdd(sh->ref_hash, sh->ref[sh->nref].name, 0,
			      hd, NULL);
	    sh->nref++;
	}
	bs->ref = hi->data.i;
    }
    cpf++;

    /* Pos */
    n = 0;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bs->pos = n-1;
    start = end = bs->pos;

    /* map qual */
    n = 0;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bam_set_map_qual(bs, n);

    /* cigar */
    n = 0;
    cigar_len = 0;
    cpt = (unsigned char *)bam_cigar(bs);
    if (*cpf == '*') {
	cpf++;
    } else {
	//while (*cpf && *cpf != '\t') {
	while (*cpf > '\t') {
	    if (isdigit(*cpf)) {
		n = n*10 + *cpf++ - '0';
	    } else {
		unsigned char op;
		union {
		    unsigned char c[4];
		    uint32_t i;
		} c4i;

		switch (*cpf++) {
		case 'M': op=BAM_CMATCH;         end+=n; break;
		case 'I': op=BAM_CINS;                   break;
		case 'D': op=BAM_CDEL;           end+=n; break;
		case 'N': op=BAM_CREF_SKIP;      end+=n; break;
		case 'S': op=BAM_CSOFT_CLIP;             break;
		case 'H': op=BAM_CHARD_CLIP;             break;
		case 'P': op=BAM_CPAD;                   break;
		case '=': op=BAM_CBASE_MATCH;    end+=n; break;
		case 'X': op=BAM_CBASE_MISMATCH; end+=n; break;
		default:
		    fprintf(stderr, "Unknown cigar opcode '%c'\n", cpf[-1]);
		    return -1;
		}

		c4i.i = (n << 4) | op;
		*cpt++ = c4i.c[0];
		*cpt++ = c4i.c[1];
		*cpt++ = c4i.c[2];
		*cpt++ = c4i.c[3];

		n = 0;
		cigar_len++;
	    }
	}
    }
    bam_set_cigar_len(bs, cigar_len);
    //printf("pos %d, %d..%d => bin %d\n", bs->pos, start, end, reg2bin(start, end));
    bam_set_bin(bs, reg2bin(start,end));
    if (!*cpf++) return -1;
    
    /* mate ref name */
    cp = cpf;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	cpf++;
    if (*cp == '*' && cp[1] == '\t') {
	bs->mate_ref = -1;
    } else if (*cp == '=' && cp[1] == '\t') {
	bs->mate_ref = bs->ref;
    } else {
	hi = HashTableSearch(sh->ref_hash, (char *)cp, cpf-cp);
	if (!hi) {
	    HashData hd;

	    fprintf(stderr, "Mate ref seq \"%.*s\" unknown\n", (int)(cpf-cp), cp);

	    /* Fabricate it instead */
	    sh->ref = realloc(sh->ref, (sh->nref+1)*sizeof(*sh->ref));
	    if (!sh->ref)
		return -1;
	    sh->ref[sh->nref].len  = 0; /* Unknown value */
	    sh->ref[sh->nref].name = malloc(cpf-cp+1);
	    if (!sh->ref[sh->nref].name)
		return -1;
	    memcpy(sh->ref[sh->nref].name, cp, cpf-cp);
	    sh->ref[sh->nref].name[cpf-cp] = 0;

	    hd.i = sh->nref;
	    hi = HashTableAdd(sh->ref_hash, sh->ref[sh->nref].name, 0,
			      hd, NULL);
	    sh->nref++;
	}
	bs->mate_ref = hi->data.i;
    }
    if (!*cpf++) return -1;

    /* mate pos */
    n = 0;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bs->mate_pos = n-1;

    /* insert size */
    n = 0;
    if (*cpf == '-') {
	sign = -1;
	cpf++;
    } else {
	sign = 1;
    }
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bs->ins_size = n*sign;

    /* seq */
    cp = cpf;
    n = 0;
    //while (*cpf && *cpf != '\t') {
    if (cpf[0] == '*' && cpf[1] == '\t') {
	cpf++;
	bs->len = 0;
    } else {
	while (*cpf > '\t') {
	    if (n == 0) {
		*cpt = lookup[*cpf]<<4;
		n = 1;
	    } else {
		n = 0;
		*cpt++ |= lookup[*cpf];
	    }
	    cpf++;
	}
	if (n == 1)
	    cpt++;
	bs->len = cpf-cp;
    }
    if (!*cpf++) return -1;

    /* qual */
    cp = cpf;
    //while (*cpf && *cpf != '\t')
    if (cpf[0] == '*' && (cpf[1] == '\0' || cpf[1] == '\t')) {
	/* no qual */
	memset(cpt, '\xFF', bs->len);
	cpt += bs->len;
	cpf++;
    } else {
	while (*cpf > '\t')
	    *cpt++ = *cpf++ - '!';
    }

    assert((char *)cpt == (char *)(bam_aux(bs)));

    if (!*cpf++ || b->no_aux) goto skip_aux;

    /* aux */
    while (*cpf) {
	unsigned char *key = cpf, *value;
	if (!(key[0] && key[1] && key[2] == ':' && key[3] && key[4] == ':'))
	    return -1;
	cpf += 5;

	value = cpf;
	while (*cpf && *cpf != '\t')
	    cpf++;

	*cpt++ = key[0];
	*cpt++ = key[1];

	switch(key[3]) {
	case 'A':
	    *cpt++ = 'A';
	    *cpt++ = *value;
	    break;

	case 'i':
	    n = atoi((char *)value);
	    if (n >= 0) {
		if (n < 256) {
		    *cpt++ = 'C';
		    *cpt++ = n;
		} else if (n < 65536) {
		    *cpt++ = 'S';
		    STORE_UINT16(cpt, n);
		} else {
		    *cpt++ = 'I';
		    STORE_UINT32(cpt, n);
		}
	    } else {
		if (n >= -128 && n < 128) {
		    *cpt++ = 'c';
		    *cpt++ = n;
		} else if (n >= -32768 && n < 32768) {
		    *cpt++ = 's';
		    STORE_UINT16(cpt, n);
		} else {
		    *cpt++ = 'i';
		    STORE_UINT32(cpt, n);
		}
	    }
	    break;

	case 'f': {
	    union {
		float f;
		int i;
	    } u;
	    u.f = atof((char *)value);
	    *cpt++ = 'f';
	    STORE_UINT32(cpt, u.i);
	    break;
	}

	case 'Z':
	    *cpt++ = 'Z';
	    while (value != cpf)
		*cpt++=*value++;
	    *cpt++ = 0;
	    break;

	case 'H':
	    *cpt++ = 'H';
	    while (value != cpf)
		*cpt++=*value++;
	    *cpt++ = 0;
	    break;

	case 'B': {
	    char subtype = *value++;
	    unsigned char *sz;
	    int count = 0;

	    *cpt++ = 'B';
	    *cpt++ = subtype;
	    sz = cpt; cpt += 4; /* Fill out later */

	    while (*value == ',') {
		value++;
		switch (subtype) {
		case 'c': case 'C':
		    *cpt++ = strtol((char *)value, (char **)&value, 10);
		    break;

		case 's': case 'S':
		    n = strtol((char *)value, (char **)&value, 10);
		    STORE_UINT16(cpt, n);
		    break;
		    
		case 'i': case 'I':
		    n = strtoll((char *)value, (char **)&value, 10);
		    STORE_UINT32(cpt, n);
		    break;
		    
		case 'f': {
		    union {
			float f;
			int i;
		    } u;

		    u.f = strtod((char *)value, (char **)&value);
		    STORE_UINT32(cpt, u.i);
		    break;
		}
		}
		count++;
	    }
	    if (value != cpf) {
		fprintf(stderr, "Malformed %c%c:B:... auxiliary field\n",
			key[0], key[1]);
		value = cpf;
	    }
	    STORE_UINT32(sz, count);
	    break;
	}

	default:
	    fprintf(stderr, "Unknown aux format code '%c'\n", key[3]);
	    break;
	}

	if (*cpf == '\t')
	    cpf++;
    }

 skip_aux:
    *cpt++ = 0;

    bs->blk_size = (unsigned char *)cpt - (unsigned char *)&bs->ref - 1;
    if (bs->blk_size >= bs->alloc)
	abort();

    return 1;
}

/*
 * Fills out the next bam_seq_t struct.
 * bs must be non-null, but *bs may be NULL or an existing bam_seq_t pointer.
 * This function will alloc and/or grow the memory accordingly, allowing for
 * efficient reuse.
 *
 * Returns 1 on success
 *         0 on eof
 *        -1 on error
 */
#ifdef ALLOW_UAC
int bam_get_seq(bam_file_t *b, bam_seq_t **bsp) {
    int32_t blk_size, blk_ret;
    bam_seq_t *bs;
    uint32_t i32;

    b->line++;

    if (!b->bam)
	return sam_next_seq(b, bsp);

    if (b->next_len > 0) {
	blk_size = b->next_len;
    } else {
	if (4 != bam_read(b, &blk_size, 4))
	    return 0;
	blk_size = le_int4(blk_size);
    }

    if (!*bsp || blk_size+20 > (*bsp)->alloc) {
	/* 20 extra is for bs->alloc to bs->cigar_len plus next_len */
	if (!(*bsp = realloc(*bsp, blk_size+20)))
	    return -1;
	(*bsp)->alloc = blk_size+20;
	(*bsp)->blk_size = blk_size;
    }
    bs = *bsp;
    
    if ((blk_ret = bam_read(b, &bs->ref, blk_size+4)) == 0)
	return 0;

    if (blk_size+4 != blk_ret) {
	if (blk_size != blk_ret) {
	    return -1;
	} else {
	    b->next_len = 0;
	    ((char *)(&bs->ref))[blk_size] = 0;
	}
    } else {
	memcpy(&b->next_len, &((char *)(&bs->ref))[blk_size], 4);
	((char *)(&bs->ref))[blk_size] = 0;
    }
    b->next_len = le_int4(b->next_len);

    bs->blk_size  = blk_size;
    bs->ref       = le_int4(bs->ref);
    bs->pos       = le_int4(bs->pos);

    // order of bit-fields in struct is platform specific, so manually decode
    i32           = le_int4(bs->bin_packed);
    bs->bin      = i32 >> 16;
    bs->map_qual = (i32 >> 8) & 0xff;
    bs->name_len = i32 & 0xff;

    i32           = le_int4(bs->flag_packed);
    bs->flag      = i32 >> 16;
    bs->cigar_len = i32 & 0xffff;

    bs->len       = le_int4(bs->len);
    bs->mate_ref  = le_int4(bs->mate_ref);
    bs->mate_pos  = le_int4(bs->mate_pos);
    bs->ins_size  = le_int4(bs->ins_size);

    if (10 == be_int4(10)) {
	int i, cigar_len = bam_cigar_len(bs);
	uint32_t *cigar = bam_cigar(bs);
	for (i = 0; i < cigar_len; i++) {
	    cigar[i] = le_int4(cigar[i]);
	}
    }

    return 1;
}

#else

int bam_get_seq(bam_file_t *b, bam_seq_t **bsp) {
    int32_t blk_size, blk_ret;
    bam_seq_t *bs;
    uint32_t i32;

    b->line++;

    if (!b->bam)
	return sam_next_seq(b, bsp);

    if (b->next_len > 0) {
	blk_size = b->next_len;
    } else {
	if (4 != bam_read(b, &blk_size, 4))
	    return 0;
	blk_size = le_int4(blk_size);
    }

    if (!*bsp || blk_size+24 > (*bsp)->alloc) {
	if (!(*bsp = realloc(*bsp, blk_size+24)))
	    return -1;
	(*bsp)->alloc = blk_size+24;
	(*bsp)->blk_size = blk_size;
    }
    bs = *bsp;

    /* The fixed-sized fields */
    if ((blk_ret = bam_read(b, &bs->ref, 32)) == 0)
	return 0;

    if (blk_ret != 32)
	return -1;

    bs->blk_size  = blk_size;
    bs->ref       = le_int4(bs->ref);
    bs->pos       = le_int4(bs->pos);

    // order of bit-fields in struct is platform specific, so manually decode
    i32           = le_int4(bs->bin_packed);
    bs->bin      = i32 >> 16;
    bs->map_qual = (i32 >> 8) & 0xff;
    bs->name_len = i32 & 0xff;

    i32           = le_int4(bs->flag_packed);
    bs->flag      = i32 >> 16;
    bs->cigar_len = i32 & 0xffff;

    bs->len       = le_int4(bs->len);
    bs->mate_ref  = le_int4(bs->mate_ref);
    bs->mate_pos  = le_int4(bs->mate_pos);
    bs->ins_size  = le_int4(bs->ins_size);

    /* Name */
    if (bam_read(b, &bs->data, bam_name_len(bs)) != bam_name_len(bs))
	return -1;

    /* Pad name out to end on a word-aligned boundary */
    blk_ret = blk_size - 32 - bam_name_len(bs);
    //bam_set_name_len(bs, round4(bam_name_len(bs)));

    bs->blk_size += round4(bam_name_len(bs)) - bam_name_len(bs);

    /* The remainder, word aligned */
    blk_size = blk_ret;
    if ((blk_ret = bam_read(b, (char *)bam_cigar(bs), blk_size+4)) == 0)
	return 0;
    if (blk_size+4 != blk_ret) {
	if (blk_size != blk_ret) {
	    return -1;
	} else {
	    b->next_len = 0;
	    ((char *)bam_cigar(bs))[blk_size] = 0;
	}
    } else {
	memcpy(&b->next_len, &((char *)bam_cigar(bs))[blk_size], 4);
	((char *)bam_cigar(bs))[blk_size] = 0;
    }
    b->next_len = le_int4(b->next_len);

    if (10 == be_int4(10)) {
	int i, cigar_len = bam_cigar_len(bs);
	uint32_t *cigar = bam_cigar(bs);
	for (i = 0; i < cigar_len; i++) {
	    cigar[i] = le_int4(cigar[i]);
	}
    }

    return 1;
}
#endif

/* Old name */
int bam_next_seq(bam_file_t *b, bam_seq_t **bsp) {
    return bam_get_seq(b, bsp);
}

static int8_t aux_type_size[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 1, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 8, 0, 4, 0, 0, 4, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

/*
 * Looks for aux field 'key' and returns the type + value.
 * The type is the first char and the value is the 2nd character onwards.
 *
 * Returns NULL if not found.
 */
char *bam_aux_find(bam_seq_t *b, const char *key) {
    char *cp = bam_aux(b);

    while (*cp) {
	int sz;

	//printf("%c%c:%c:?\n", cp[0], cp[1], cp[2]);

	if (cp[0] == key[0] && cp[1] == key[1])
	    return cp+2;
	
	if ((sz = aux_type_size[(uint8_t) cp[2]])) {
	    /* Fixed length fields */
	    cp += sz + 3;
	} else {
	    switch (cp[2]) {
	    case 'Z':
	    case 'H': {  /* Variable length, null terminated */
		cp += 3;
		while (*cp++)
		    ;
		break;
	    }
	    case 'B': {  /* Array types */
		uint32_t count;
		if ((sz = aux_type_size[(uint8_t) cp[3]])) {
		    count = ((uint32_t)    cp[4]
			     + ((uint32_t) cp[5] << 8)
			     + ((uint32_t) cp[6] << 16)
			     + ((uint32_t) cp[7] << 24));
		    cp += 8 + count * sz;
		} else {
		    return NULL;
		}
	    }
	    default:
		return NULL;
	    }
	}
    }

    return NULL;
}

int32_t bam_aux_i(const uint8_t *dat) {
    switch(dat[0]) {
    case 'i':
	return (int32_t)(dat[1] + (dat[2]<<8) + (dat[3]<<16) + (dat[4]<<24));
    case 'I':
	return (uint32_t)(dat[1] + (dat[2]<<8) + (dat[3]<<16) + (dat[4]<<24));
	break;
    case 's':
	return (int16_t)(dat[1] + (dat[2]<<8));
    case 'S':
	return (uint16_t)(dat[1] + (dat[2]<<8));
    case 'c':
	return (int8_t)dat[1];
    case 'C':
	return (uint8_t)dat[1];
    }

    abort();
}

float bam_aux_f(const uint8_t *dat) {
    assert(dat[0] == 'f');
    return (float)((int32_t)((dat[1]<<0)+
			     (dat[2]<<8)+
			     (dat[3]<<16)+
			     (dat[4]<<24)));
}

double bam_aux_d(const uint8_t *dat) { 
    assert(dat[0] == 'd');
    return (double)((int64_t)((((uint64_t)dat[1])<<0)+
			      (((uint64_t)dat[2])<<8)+
			      (((uint64_t)dat[3])<<16)+
			      (((uint64_t)dat[4])<<24)+
			      (((uint64_t)dat[5])<<32)+
			      (((uint64_t)dat[6])<<40)+
			      (((uint64_t)dat[7])<<48)+
			      (((uint64_t)dat[8])<<54)));
}

char bam_aux_A(const uint8_t *dat) {
    assert(dat[0] == 'A');
    return dat[1];
}

char *bam_aux_Z(const uint8_t *dat) {
    assert(dat[0] == 'Z' || dat[0] == 'H');
    return (char *)(dat+1);
}

/*
 * An iterator on bam aux fields. NB: This code is not reentrant or multi-
 * thread capable. The values returned are valid until the next call to
 * this function.
 * key:  points to an array of 2 characters (eg "RG", "NM")
 * type: points to an address of 1 character (eg 'Z', 'i')
 * val:  points to an address of a bam_aux_t union.
 *
 * Pass in *iter_handle as NULL to initialise the search and then
 * pass in the modified value on each subsequent call to continue the search.
 *
 * Returns 0 if the next value is valid, setting key, type and val.
 *        -1 when no more found.
 */
int bam_aux_iter(bam_seq_t *b, char **iter_handle,
		 char *key, char *type, bam_aux_t *val) {
    char *s;

    if (!iter_handle || !*iter_handle) {
	s = (char *)bam_aux(b);
    } else {
	s = *iter_handle;
    }

    /* We null terminate our aux list for ease */
    if (s[0] == 0)
	return -1;

    key[0] = s[0];
    key[1] = s[1];
    
    switch (s[2]) {
    case 'A':
	if (type) *type = 'A';
	if (val) val->i = *(s+3);
	s+=4;
	break;

    case 'C':
	if (type) *type = 'i';
	if (val) val->i = *(uint8_t *)(s+3);
	s+=4;
	break;

    case 'c':
	if (type) *type = 'i';
	if (val) val->i = *(int8_t *)(s+3);
	s+=4;
	break;

    case 'S':
	if (type) *type = 'i';
	if (val)
	    val->i = (uint16_t)((((unsigned char *)s)[3]<< 0) +
				(((unsigned char *)s)[4]<< 8));
	s+=5;
	break;

    case 's':
	if (type) *type = 'i';
	if (val)
	    val->i = (int16_t)((((unsigned char *)s)[3]<< 0) +
			       (((unsigned char *)s)[4]<< 8));
	s+=5;
	break;

    case 'I':
	if (type) *type = 'i';
	if (val)
	    val->i = (uint32_t)((((unsigned char *)s)[3]<< 0) +
				(((unsigned char *)s)[4]<< 8) +
				(((unsigned char *)s)[5]<<16) +
				(((unsigned char *)s)[6]<<24));
	s+=7;
	break;

    case 'i':
	if (type) *type = 'i';
	if (val)
	    val->i = (int32_t)((((unsigned char *)s)[3]<< 0) +
			       (((unsigned char *)s)[4]<< 8) +
			       (((unsigned char *)s)[5]<<16) +
			       (((unsigned char *)s)[6]<<24));
	s+=7;
	break;

    case 'f':
	if (type) *type = 'f';
	if (val) /* Assume same endianness as integer */
	    val->i = (int32_t)((((unsigned char *)s)[3]<< 0) +
			       (((unsigned char *)s)[4]<< 8) +
			       (((unsigned char *)s)[5]<<16) +
			       (((unsigned char *)s)[6]<<24));
	s+=7;
	break;

    case 'd':
	if (type) *type = 'd';
	if (val) /* Assume same endianness as integer */
	    val->i64 = (uint64_t)(((uint64_t)(((unsigned char *)s)[ 3])<< 0) +
				  ((uint64_t)(((unsigned char *)s)[ 4])<< 8) +
				  ((uint64_t)(((unsigned char *)s)[ 5])<<16) +
				  ((uint64_t)(((unsigned char *)s)[ 6])<<24) +
				  ((uint64_t)(((unsigned char *)s)[ 7])<<32) +
				  ((uint64_t)(((unsigned char *)s)[ 8])<<40) +
				  ((uint64_t)(((unsigned char *)s)[ 9])<<48) +
				  ((uint64_t)(((unsigned char *)s)[10])<<54));
	s+=11;
	break;

    case 'Z': case 'H':
	if (type) *type = s[2];
	s+=3;
	if (val) val->s = s;
	while (*s++);
	break;

    case 'B': {
	uint32_t count;
	if (type) *type = 'B';
	count = (unsigned int)((((unsigned char *)s)[4]<< 0) +
			       (((unsigned char *)s)[5]<< 8) +
			       (((unsigned char *)s)[6]<<16) +
			       (((unsigned char *)s)[7]<<24));

	if (val) {
	    val->B.n = count;
	    val->B.t = s[3];
	    val->B.s = (unsigned char *)s+8;
	}
	s+=8;

	switch(val->B.t) {
	case 'c': case 'C': s +=   count; break;
	case 's': case 'S': s += 2*count; break;
	case 'i': case 'I': s += 4*count; break;
	case 'f':           s += 4*count; break;
	default:
	    fprintf(stderr, "Unknown sub-type '%c' for aux type 'B'\n",
		    val->B.t);
	    return -1;
	}
	break;
    }

    default:
	fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	return -1;
    }

    if (iter_handle)
	*iter_handle = s;

    return 0;
}

static int reg2bin(int start, int end) {
    if (end>start) end--;
    if ((start>>14) == (end>>14)) return ((1<<15)-1)/7 + (start>>14);
    if ((start>>17) == (end>>17)) return ((1<<12)-1)/7 + (start>>17);
    if ((start>>20) == (end>>20)) return ((1<<9 )-1)/7 + (start>>20);
    if ((start>>23) == (end>>23)) return ((1<<6 )-1)/7 + (start>>23);
    if ((start>>26) == (end>>26)) return ((1<<3 )-1)/7 + (start>>26);
    return 0;
}

/*
 * Constructs a bam_seq_t from components.
 * Ignores auxiliary tags for now.
 *
 * Returns -1 on error
 *          number of bytes written to bam_seq_t on success (ie tag offset)
 */
int bam_construct_seq(bam_seq_t **b, size_t extra_len,
		      const char *qname, size_t qname_len,
		      int flag,
		      int rname,      // Ref ID
		      int pos, // first aligned base (1-based)
		      int end, // last aligned base (to calculate bin)
		      int mapq,
		      uint32_t ncigar, const uint32_t *cigar,
		      int mrnm,       // Mate Ref ID
		      int mpos,
		      int isize,
		      int len,
		      const char *seq,
		      const char *qual) {
    size_t required;
    char *cp;
    int i;
    uint32_t *ip;

    /*
     * cp = "=ACMGRSVTWYHKDBN";
     * memset(L, 15, 256);
     * for (i = 0; i < 16; i++) {
     *     L[cp[i]] = L[tolower(cp[i])] = i;
     * }
     */
    static const char L[256] = {
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15, 0,15,15,
	15, 1,14, 2,13,15,15, 4,11,15,15,12,15, 3,15,15,
	15,15, 5, 6, 8,15, 7, 9,15,10,15,15,15,15,15,15,
	15, 1,14, 2,13,15,15, 4,11,15,15,12,15, 3,15,15,
	15,15, 5, 6, 8,15, 7, 9,15,10,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
	15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15
    };

    /* Sanity checks */
    if (NULL == b) return -1;
    if (len < 0) return -1;  /* not sure why the spec has it as an int */
    if (qname_len > 0 && NULL == qname) return -1;
    if (ncigar > 0 && NULL == cigar) return -1;
    if (len > 0 && NULL == seq) return -1;

    /* Reallocate if needed */
    required = (sizeof(**b)              /* the struct itself */
#ifdef ALLOW_UAC
		+ qname_len + 1         /* query name (unaligned) */
#else
		+ round4(qname_len + 1) /* query name (aligned) */
#endif
		+ 4 * ncigar            /* CIGAR string */
		+ (len + 1) / 2         /* Sequence, 2 bases per byte */
		+ len                   /* Quality */
		+ extra_len + 1);       /* Extra for optional tags */

    if (NULL == *b || (*b)->alloc < required) {
	bam_seq_t *new_bam = realloc(*b, required);
	if (NULL == new_bam) return -1;
	*b = new_bam;
	(*b)->alloc = required;
    }

    (*b)->ref = rname;
    (*b)->pos = pos-1;
    bam_set_map_qual(*b, mapq);
    bam_set_name_len(*b, qname_len+1);
    bam_set_flag(*b, flag);
    bam_set_cigar_len(*b, ncigar);
    (*b)->len = len;
    (*b)->mate_ref = mrnm;
    (*b)->mate_pos = mpos-1;
    (*b)->ins_size = isize;

    cp = bam_name(*b);
    memcpy(cp, qname, qname_len);
    cp[qname_len] = 0;

    /* Cigar */
    cp = (char *)bam_cigar(*b);
    ip = (uint32_t *)cp;
    for (i = 0; i < ncigar; i++) {
	ip[i] = cigar[i];
    }
    cp += ncigar*4;

    /* Bin */
    if (0 == end) { /* Calculate end from pos and cigar */
	end = pos;
	for (i = 0; i < ncigar; i++) {
	    int op = cigar[i] & BAM_CIGAR_MASK;
	    if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
		end += cigar[i] >> BAM_CIGAR_SHIFT;
	}
    }

    bam_set_bin(*b, reg2bin(pos-1,end-1));

    /* Seq */
    for (i = 0; i < len-1; i += 2) {
	*cp++ = (L[(uc)seq[i]]<<4) + L[(uc)seq[i+1]];
    }
    if (i < len)
	*cp++ = L[(uc)seq[i]]<<4;

    /* Qual */
    if (qual) {
	memcpy(cp, qual, len);
	cp += len;
    } else {
	for (i = 0; i < len; i++) {
	    *cp++ = '\xff';
	}
    }

    *cp = 0; /* terminate aux list, for ease of parsing later */

    /* cp now points to the auxiliary tags if required */
    (*b)->blk_size = (int)(cp-(char *)&(*b)->ref);
    return (int)(cp-(char *)(*b));
}

int bam_aux_add(bam_seq_t **b, const char tag[2], char type,
		uint32_t array_len, const void *data) {
    int tlen;
    size_t len;
    size_t used;
    uint8_t *cp;
#ifndef SP_LITTLE_ENDIAN
    uint32_t i;
#endif

    if (NULL == b || NULL == *b) return -1;

    /* Find size  of data type and how much space is needed */
    if (0 == (tlen = aux_type_size[(uint8_t) type])) {
	if (type == 'H' || type == 'Z') { /* Variable length types */
	    if (array_len != 0) return -1; /* No arrays for these allowed */
	    tlen = strlen((const char *) data) + 1;
	} else {
	    /* unknown type */
	    return -1;
	}
    }

    len = array_len > 0 ? 8 + tlen * array_len : 3 + tlen;
    
    /* Find the end of the existing tags and ensure there is enough space */
    cp = (uint8_t *)&(*b)->ref + (*b)->blk_size;
    used = cp - (uint8_t *)(*b);

    if ((*b)->alloc < used + len + 1) { /* + 1 for NUL terminator */
	size_t required = used + len + 1;
	bam_seq_t *new_bam = realloc((*b), required);
	if (NULL == new_bam) return -1;
	*b = new_bam;
	(*b)->alloc = required;
	cp = (uint8_t *)new_bam + used;
    }

    /* Append the data */
    *cp++ = tag[0];
    *cp++ = tag[1];
    if (array_len > 0) {  /* Array type */
	*cp++ = 'B';
	*cp++ = type;
	STORE_UINT32(cp, array_len);
    } else {
	*cp++ = type;
    }

    if (array_len == 0) array_len = 1;
#ifdef SP_LITTLE_ENDIAN
    memcpy(cp, data, array_len * tlen);
    cp += array_len * tlen;
#else
    switch (type) {
    case 'A': case 'c': case 'C':
	memcpy(cp, data, array_len);
	cp += array_len;
	break;
    case 's':
    case 'S': {
	uint16_t *sdata = (uint16_t *) data;
	for (i = 0; i < array_len; i++) {
	    STORE_UINT16(cp, sdata[i]);
	}
	break;
    }
    case 'i':
    case 'I':
    case 'f': {
	uint32_t *idata = (uint32_t *) data;
	for (i = 0; i < array_len; i++) {
	    STORE_UINT32(cp, idata[i]);
	}
	break;
    }
    case 'd': {
	uint64_t *ddata = (uint64_t *) data;
	for (i = 0; i < array_len; i++) {
	    STORE_UINT64(cp, ddata[i]);
	}
	break;
    }
    case 'H': case 'Z':
	memcpy(cp, data, tlen);
	cp += tlen;
	break;
    }
#endif
    
    /* Put a NUL at the end for bam_aux_iter */
    *cp = 0;

    /* Update block_size */
    (*b)->blk_size = (uint32_t)(cp - (uint8_t *)&(*b)->ref);
    
    return 0;
}

/* Calculate space needed to store tags */

ssize_t bam_aux_size_vec(uint32_t count, bam_aux_tag_t *tags) {
    uint32_t i;
    ssize_t len = 0;
    int sz;

    if (NULL == tags) return -1;
    
    for (i = 0; i < count; i++) {
	switch (tags[i].type) {
	case 'C': case 'S': case 'I':
	    if (tags[i].value.ui < 256) {
		sz = 1;
	    } else if (tags[i].value.ui < 65536) {
		sz = 2;
	    } else {
		sz = 4;
	    }
	    break;
	case 'c': case 's': case 'i':
	    if (tags[i].value.i >= -128 && tags[i].value.i < 128) {
		sz = 1;
	    } else if (tags[i].value.i >= -32768 && tags[i].value.i < 32768) {
		sz = 2;
	    } else {
		sz = 4;
	    }
	    break;
	case 'A':
	    sz = 1;
	    break;
	case 'f':
	    sz = 4;
	    break;
	case 'd':
	    sz = 8;
	    break;
	case 'H': case 'Z':
	    if (tags[i].array_len != 0) return -1;
	    sz = strlen(tags[i].value.z) + 1;
	    break;
	default:
	    return -1; /* bad data type */
	}
	len += tags[i].array_len == 0 ? sz + 3 : sz * tags[i].array_len + 8;
    }

    return len + 1;  /* + 1 for NUL byte at end */
}

int bam_aux_add_vec(bam_seq_t **b, uint32_t count, bam_aux_tag_t *tags) {
    ssize_t required = bam_aux_size_vec(count, tags);
    uint32_t i;
    size_t used;

    if (required < 0) return -1;
    if (NULL == b || NULL == *b) return -1;

    /* Find the end of the existing tags and ensure there is enough space.
       Do this once for the entire vector so we don't keep reallocing */
    used = (uint8_t *)&(*b)->ref + (*b)->blk_size - (uint8_t *)(*b);

    if ((*b)->alloc < used + required) {
	bam_seq_t *new_bam = realloc((*b), used + required);
	if (NULL == new_bam) return -1;
	*b = new_bam;
	(*b)->alloc = used + required;
    }

    /* Add the tags, storing integers in the most appropriate size */
    for (i = 0; i < count; i++) {
	if (tags[i].array_len > 0) {
	    /* Deal with array tags */
	    if (bam_aux_add(b, tags[i].tag, tags[i].type,
			    tags[i].array_len, tags[i].value.array)) return -1;
	    continue;
	}

	/* Non-array tags, storing integers as the most appropriate size */
	switch (tags[i].type) {
	case 'C': case 'S': case 'I':
	    if (tags[i].value.ui < 256) {
		uint8_t byte = tags[i].value.ui;
		if (bam_aux_add(b, tags[i].tag, 'C', 0, &byte)) return -1;
	    } else if (tags[i].value.ui < 65536) {
		uint16_t word = tags[i].value.ui;
		if (bam_aux_add(b, tags[i].tag, 'S', 0, &word)) return -1;
	    } else {
		if (bam_aux_add(b, tags[i].tag, 'I', 0, &tags[i].value.ui)) {
		    return -1;
		}
	    }
	    break;
	case 'c': case 's': case 'i':
	    if (tags[i].value.i >= -128 && tags[i].value.i < 128) {
		int8_t byte = tags[i].value.i;
		if (bam_aux_add(b, tags[i].tag, 'c', 0, &byte)) return -1;
	    } else if (tags[i].value.i >= -32768 && tags[i].value.i < 32768) {
		int16_t word = tags[i].value.i;
		if (bam_aux_add(b, tags[i].tag, 's', 0, &word)) return -1;
	    } else {
		if (bam_aux_add(b, tags[i].tag, 'i', 0, &tags[i].value.i))
		    return -1;
	    }
	    break;
	case 'A':
	    if (bam_aux_add(b, tags[i].tag, 'A', 0, &tags[i].value.a))
		return -1;
	    break;
	case 'f': case 'd':
	    if (bam_aux_add(b, tags[i].tag, tags[i].type, 0, &tags[i].value.f))
		return -1;
	    break;
	case 'H': case 'Z':
	    if (bam_aux_add(b, tags[i].tag, tags[i].type, 0, tags[i].value.z))
		return -1;
	    break;
	default:
	    return -1; /* unknown type */
	}	    
    }
    return 0;
}

/* Add SAM-formatted aux tags to a bam_seq_t struct.
   This is basically a copy of the code in sam_next_seq.  Unfortunately
   trying to get them to use a common version slows sam_next_seq down
   rather a lot, even when inlined.  Hence this extra copy. */

int bam_aux_add_from_sam(bam_seq_t **bsp, char *sam) {
    unsigned char *cpf = (unsigned char *) sam;
    unsigned char *cpt = (unsigned char *)&(*bsp)->ref + (*bsp)->blk_size;
    unsigned char *end = (unsigned char *)(*bsp) + (*bsp)->alloc;
    int n;

    while (*cpf) {
	unsigned char *key = cpf, *value;
	size_t max_len;

	if (!(key[0] && key[1] && key[2] == ':' && key[3] && key[4] == ':'))
	    return -1;
	cpf += 5;

	value = cpf;
	while (*cpf && *cpf != '\t')
	    cpf++;

	if (aux_type_size[key[3]]) {
            max_len = aux_type_size[key[3]] + 3;
        } else if (key[3] != 'B') {
            max_len = cpf - value + 4;
        } else {
	    /* Worst case */
            max_len = (cpf - value) * 4 + 8;
        }

	/* ensure we have enough room */
        if (end - cpt < max_len) {
            size_t used = cpt - (unsigned char *)(*bsp);
            bam_seq_t *new_bam = realloc(*bsp, used + max_len);
            if (NULL == new_bam) return -1;
            *bsp = new_bam;
            (*bsp)->alloc += used + max_len;
            cpt = (unsigned char *)(*bsp) + used;
            end = (unsigned char *)(*bsp) + (*bsp)->alloc;
        }

	*cpt++ = key[0];
	*cpt++ = key[1];

	switch(key[3]) {
	case 'A':
	    *cpt++ = 'A';
	    *cpt++ = *value;
	    break;

	case 'i':
	    n = atoi((char *)value);
	    if (n >= 0) {
		if (n < 256) {
		    *cpt++ = 'C';
		    *cpt++ = n;
		} else if (n < 65536) {
		    *cpt++ = 'S';
		    STORE_UINT16(cpt, n);
		} else {
		    *cpt++ = 'I';
		    STORE_UINT32(cpt, n);
		}
	    } else {
		if (n >= -128 && n < 128) {
		    *cpt++ = 'c';
		    *cpt++ = n;
		} else if (n >= -32768 && n < 32768) {
		    *cpt++ = 's';
		    STORE_UINT16(cpt, n);
		} else {
		    *cpt++ = 'i';
		    STORE_UINT32(cpt, n);
		}
	    }
	    break;

	case 'f': {
	    union {
		float f;
		int i;
	    } u;
	    u.f = atof((char *)value);
	    *cpt++ = 'f';
	    STORE_UINT32(cpt, u.i);
	    break;
	}

	case 'Z':
	    *cpt++ = 'Z';
	    while (value != cpf)
		*cpt++=*value++;
	    *cpt++ = 0;
	    break;

	case 'H':
	    *cpt++ = 'H';
	    while (value != cpf)
		*cpt++=*value++;
	    *cpt++ = 0;
	    break;

	case 'B': {
	    char subtype = *value++;
	    unsigned char *sz;
	    int count = 0;

	    *cpt++ = 'B';
	    *cpt++ = subtype;
	    sz = cpt; cpt += 4; /* Fill out later */

	    while (*value == ',') {
		value++;
		switch (subtype) {
		case 'c': case 'C':
		    *cpt++ = strtol((char *)value, (char **)&value, 10);
		    break;

		case 's': case 'S':
		    n = strtol((char *)value, (char **)&value, 10);
		    STORE_UINT16(cpt, n);
		    break;
		    
		case 'i': case 'I':
		    n = strtoll((char *)value, (char **)&value, 10);
		    STORE_UINT32(cpt, n);
		    break;
		    
		case 'f': {
		    union {
			float f;
			int i;
		    } u;

		    u.f = strtod((char *)value, (char **)&value);
		    STORE_UINT32(cpt, u.i);
		    break;
		}
		}
		count++;
	    }
	    if (value != cpf) {
		fprintf(stderr, "Malformed %c%c:B:... auxiliary field\n",
			key[0], key[1]);
		value = cpf;
	    }
	    STORE_UINT32(sz, count);
	    break;
	}

	default:
	    fprintf(stderr, "Unknown aux format code '%c'\n", key[3]);
	    break;
	}

	if (*cpf == '\t')
	    cpf++;
    }

    if (cpt == end) { /* Hopefully very unlikely */
        size_t used = cpt - (unsigned char *)(*bsp);
        bam_seq_t *new_bam = realloc(*bsp, used + 1);
        if (NULL == new_bam) return -1;
        *bsp = new_bam;
        (*bsp)->alloc += used + 1;
        cpt = (unsigned char *)(*bsp) + used;
    }

    *cpt = 0;
    (*bsp)->blk_size = cpt - (unsigned char *)&(*bsp)->ref;
    return 0;
}

/*! Add preformated raw aux data to the bam_seq.
 *
 * Consider using bam_aux_add instead if you have information in a more
 * integer or string form.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
int bam_aux_add_data(bam_seq_t **b, const char tag[2], char type,
		     size_t len, const uint8_t *data) {
    uint8_t *cp;
    size_t used;

    if (NULL == b || NULL == data) return -1;

    /* Find the end of the existing tags and ensure there is enough space */
    cp = (uint8_t *)&(*b)->ref + (*b)->blk_size;
    used = cp - (uint8_t *)(*b);

    if ((*b)->alloc < used + len + 4) {
	size_t required = used + len + 4;
	bam_seq_t *new_bam = realloc((*b), required);
	if (NULL == new_bam) return -1;
	*b = new_bam;
	(*b)->alloc = required;
	cp = (uint8_t *) new_bam + used;
    }

    *cp++ = tag[0];
    *cp++ = tag[1];
    *cp++ = type;
    memcpy(cp, data, len);
    *cp = 0;

    (*b)->blk_size = (uint32_t)(cp - (uint8_t *)&(*b)->ref);

    return 0;
}

static unsigned char *append_int(unsigned char *cp, int32_t i) {
    int32_t j;

    if (i < 0) {
	*cp++ = '-';
	if (i == INT_MIN) {
	    *cp++ = '2'; *cp++ = '1'; *cp++ = '4'; *cp++ = '7';
	    *cp++ = '4'; *cp++ = '8'; *cp++ = '3'; *cp++ = '6';
	    *cp++ = '4'; *cp++ = '8';
	    return cp;
	}

	i = -i;
    } else if (i == 0) {
	*cp++ = '0';
	return cp;
    }

    //if (i < 10)         goto b0;
    if (i < 100)        goto b1;
    //if (i < 1000)       goto b2;
    if (i < 10000)      goto b3;
    //if (i < 100000)     goto b4;
    if (i < 1000000)    goto b5;
    //if (i < 10000000)   goto b6;
    if (i < 100000000)  goto b7;

     if ((j = i / 1000000000)) {*cp++ = j + '0'; i -= j*1000000000; goto x8;}
     if ((j = i / 100000000))  {*cp++ = j + '0'; i -= j*100000000;  goto x7;}
 b7: if ((j = i / 10000000))   {*cp++ = j + '0'; i -= j*10000000;   goto x6;}
     if ((j = i / 1000000))    {*cp++ = j + '0', i -= j*1000000;    goto x5;}
 b5: if ((j = i / 100000))     {*cp++ = j + '0', i -= j*100000;     goto x4;}
     if ((j = i / 10000))      {*cp++ = j + '0', i -= j*10000;      goto x3;}
 b3: if ((j = i / 1000))       {*cp++ = j + '0', i -= j*1000;       goto x2;}
     if ((j = i / 100))        {*cp++ = j + '0', i -= j*100;        goto x1;}
 b1: if ((j = i / 10))         {*cp++ = j + '0', i -= j*10;         goto x0;}
     if (i)                     *cp++ = i + '0';
    return cp;

 x8: *cp++ = i / 100000000 + '0', i %= 100000000;
 x7: *cp++ = i / 10000000  + '0', i %= 10000000;
 x6: *cp++ = i / 1000000   + '0', i %= 1000000;
 x5: *cp++ = i / 100000    + '0', i %= 100000;
 x4: *cp++ = i / 10000     + '0', i %= 10000;
 x3: *cp++ = i / 1000      + '0', i %= 1000;
 x2: *cp++ = i / 100       + '0', i %= 100;
 x1: *cp++ = i / 10        + '0', i %= 10;
 x0: *cp++ = i             + '0';

    return cp;
}

/*
 * Unsigned version of above.
 * Only differs when the int has the top bit set (~2.15 billion and above).
 */
static unsigned char *append_uint(unsigned char *cp, uint32_t i) {
    uint32_t j;

    if (i == 0) {
	*cp++ = '0';
	return cp;
    }

    //if (i < 10)         goto b0;
    if (i < 100)        goto b1;
    //if (i < 1000)       goto b2;
    if (i < 10000)      goto b3;
    //if (i < 100000)     goto b4;
    if (i < 1000000)    goto b5;
    //if (i < 10000000)   goto b6;
    if (i < 100000000)  goto b7;

     if ((j = i / 1000000000)) {*cp++ = j + '0'; i -= j*1000000000; goto x8;}
     if ((j = i / 100000000))  {*cp++ = j + '0'; i -= j*100000000;  goto x7;}
 b7: if ((j = i / 10000000))   {*cp++ = j + '0'; i -= j*10000000;   goto x6;}
     if ((j = i / 1000000))    {*cp++ = j + '0', i -= j*1000000;    goto x5;}
 b5: if ((j = i / 100000))     {*cp++ = j + '0', i -= j*100000;     goto x4;}
     if ((j = i / 10000))      {*cp++ = j + '0', i -= j*10000;      goto x3;}
 b3: if ((j = i / 1000))       {*cp++ = j + '0', i -= j*1000;       goto x2;}
     if ((j = i / 100))        {*cp++ = j + '0', i -= j*100;        goto x1;}
 b1: if ((j = i / 10))         {*cp++ = j + '0', i -= j*10;         goto x0;}
     if (i)                     *cp++ = i + '0';
    return cp;

 x8: *cp++ = i / 100000000 + '0', i %= 100000000;
 x7: *cp++ = i / 10000000  + '0', i %= 10000000;
 x6: *cp++ = i / 1000000   + '0', i %= 1000000;
 x5: *cp++ = i / 100000    + '0', i %= 100000;
 x4: *cp++ = i / 10000     + '0', i %= 10000;
 x3: *cp++ = i / 1000      + '0', i %= 1000;
 x2: *cp++ = i / 100       + '0', i %= 100;
 x1: *cp++ = i / 10        + '0', i %= 10;
 x0: *cp++ = i             + '0';

    return cp;
}

/*
 * This is set up so that count should never be more than BGZF_BUFF_SIZE. 
 * This has been chosen to deliberately be small enough such that the
 * bgzf header/footer + worst-case expansion (deflateBound() func) of 'buf'
 * are <= 65536, thus ensuring BGZF BSIZE is always 16-bit.
 *
 * Returns 0 on success;
 *        -1 on error
 */
static int bgzf_write(int fd, int level, const void *buf, size_t count) {
    unsigned char blk[Z_BUFF_SIZE];
    z_stream s;
    int cdata_pos;
    int cdata_size;
    int cdata_alloc;
    int err;
    uint32_t crc;

    /* Initialise zlib stream */
    cdata_pos = 18;
    cdata_alloc = Z_BUFF_SIZE;
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)buf;
    s.avail_in = count;
    s.total_in = 0;
    s.next_out  = blk + cdata_pos;
    s.avail_out = cdata_alloc;
    s.total_out = 0;
    s.data_type = Z_BINARY;

    /* Compress it */
    //err = deflateInit2(&s, level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY);
    err = deflateInit2(&s, level, Z_DEFLATED, -15, 8, Z_FILTERED);

    if (err != Z_OK) {
	fprintf(stderr, "zlib deflateInit2 error: %s\n", s.msg);
	return -1;
    }

    /* Encode to 'cdata' array */
    for (;s.avail_in;) {
	s.next_out = blk + cdata_pos;
	s.avail_out = cdata_alloc - cdata_pos;
	if (cdata_alloc - cdata_pos <= 0) {
	    fprintf(stderr, "Deflate produced larger output than expected. Abort\n"); 
	    return -1;
	}
	err = deflate(&s, Z_NO_FLUSH); // or Z_FINISH?
	cdata_pos = cdata_alloc - s.avail_out;
	if (err != Z_OK) {
	    fprintf(stderr, "zlib deflate error: %s\n", s.msg);
	    break;
	}
    }
    if (deflate(&s, Z_FINISH) != Z_STREAM_END) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    cdata_size = s.total_out;

    if (deflateEnd(&s) != Z_OK) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }

    assert(cdata_size <= 65536);

    /* Fill out gzip header */
    blk[ 0] =  31; // ID1
    blk[ 1] = 139; // ID2
    blk[ 2] =   8; // CM (deflate)
    blk[ 3] =   4; // FLAGS (FEXTRA)
    blk[ 4] = blk[5] = blk[6] = blk[7] = 0; // MTIME
    blk[ 8] =   0; // XFL
    blk[ 9] = 255; // OS (unknown)
    blk[10] =   6; // XLEN
    blk[11] =   0; // XLEN

    /* Extra BGZF fields */
    blk[12] =  66; // SI1
    blk[13] =  67; // SI2
    blk[14] =   2; // SLEN
    blk[15] =   0; // SLEN
    blk[16] = ((cdata_size + 25) >> 0) & 0xff;
    blk[17] = ((cdata_size + 25) >> 8) & 0xff;

    crc = crc32(0L, NULL, 0L);
    crc = crc32(crc, (unsigned char *)buf, count);
    blk[18+cdata_size+0] = (crc >> 0) & 0xff;
    blk[18+cdata_size+1] = (crc >> 8) & 0xff;
    blk[18+cdata_size+2] = (crc >>16) & 0xff;
    blk[18+cdata_size+3] = (crc >>24) & 0xff;
    blk[18+cdata_size+4] = (count >> 0) & 0xff;
    blk[18+cdata_size+5] = (count >> 8) & 0xff;
    blk[18+cdata_size+6] = (count >>16) & 0xff;
    blk[18+cdata_size+7] = (count >>24) & 0xff;

    //printf("count=%d/%x, cdata_size=%d/%x\n", count, count, cdata_size, cdata_size);
    if (18+cdata_size+8 != write(fd, blk, 18+cdata_size+8))
	return -1;

    return 0;
}

/*
 * Writes a single bam sequence object.
 * Returns 0 on success
 *        -1 on failure
 */
int bam_put_seq(bam_file_t *fp, bam_seq_t *b) {
    char *auxh, aux_key[2], type;
    bam_aux_t val;

    /*
     * Thread safe version of:
     *
     *   static int init_done = 0;
     *   static uint16_t code2base[256];
     *   
     *   if (!init_done) {
     *       int i;
     *       char c[2];
     *       uint16_t s;
     *       for (i = 0; i < 256; i++) {
     *           c[0] = "=ACMGRSVTWYHKDBN"[i >> 4];
     *           c[1] = "=ACMGRSVTWYHKDBN"[i & 15];
     *           s = *(uint16_t *)c;
     *           code2base[i] = s;
     *	     
     *           //printf("%5d,%c", code2base[i],i%8==7?'\n':' ');
     *       }
     *       init_done = 1;
     *   }
     */
#ifdef ALLOW_UAC
    static const uint16_t code2base[256] = {
	15677, 16701, 17213, 19773, 18237, 21053, 21309, 22077,
	21565, 22333, 22845, 18493, 19261, 17469, 16957, 20029,
	15681, 16705, 17217, 19777, 18241, 21057, 21313, 22081,
	21569, 22337, 22849, 18497, 19265, 17473, 16961, 20033,
	15683, 16707, 17219, 19779, 18243, 21059, 21315, 22083,
	21571, 22339, 22851, 18499, 19267, 17475, 16963, 20035,
	15693, 16717, 17229, 19789, 18253, 21069, 21325, 22093,
	21581, 22349, 22861, 18509, 19277, 17485, 16973, 20045,
	15687, 16711, 17223, 19783, 18247, 21063, 21319, 22087,
	21575, 22343, 22855, 18503, 19271, 17479, 16967, 20039,
	15698, 16722, 17234, 19794, 18258, 21074, 21330, 22098,
	21586, 22354, 22866, 18514, 19282, 17490, 16978, 20050,
	15699, 16723, 17235, 19795, 18259, 21075, 21331, 22099,
	21587, 22355, 22867, 18515, 19283, 17491, 16979, 20051,
	15702, 16726, 17238, 19798, 18262, 21078, 21334, 22102,
	21590, 22358, 22870, 18518, 19286, 17494, 16982, 20054,
	15700, 16724, 17236, 19796, 18260, 21076, 21332, 22100,
	21588, 22356, 22868, 18516, 19284, 17492, 16980, 20052,
	15703, 16727, 17239, 19799, 18263, 21079, 21335, 22103,
	21591, 22359, 22871, 18519, 19287, 17495, 16983, 20055,
	15705, 16729, 17241, 19801, 18265, 21081, 21337, 22105,
	21593, 22361, 22873, 18521, 19289, 17497, 16985, 20057,
	15688, 16712, 17224, 19784, 18248, 21064, 21320, 22088,
	21576, 22344, 22856, 18504, 19272, 17480, 16968, 20040,
	15691, 16715, 17227, 19787, 18251, 21067, 21323, 22091,
	21579, 22347, 22859, 18507, 19275, 17483, 16971, 20043,
	15684, 16708, 17220, 19780, 18244, 21060, 21316, 22084,
	21572, 22340, 22852, 18500, 19268, 17476, 16964, 20036,
	15682, 16706, 17218, 19778, 18242, 21058, 21314, 22082,
	21570, 22338, 22850, 18498, 19266, 17474, 16962, 20034,
	15694, 16718, 17230, 19790, 18254, 21070, 21326, 22094,
	21582, 22350, 22862, 18510, 19278, 17486, 16974, 20046
    };
#endif

    if (!fp->binary) {
	/* SAM */
	unsigned char *end = fp->out + BGZF_BUFF_SIZE, *dat;
	int sz, i, n;

#define BF_FLUSH() \
	do {			  \
	    if (fp->out_p - fp->out != \
		write(fp->fd, fp->out, fp->out_p - fp->out)) \
		return -1;				     \
	    fp->out_p=fp->out;				\
	} while(0)

	/* QNAME */
	if (end - fp->out_p < (sz = bam_name_len(b))) BF_FLUSH();
	memcpy(fp->out_p, bam_name(b), sz-1); fp->out_p += sz-1;
	*fp->out_p++ = '\t';

	/* FLAG */
	if (end-fp->out_p < 5) BF_FLUSH();
	fp->out_p = append_int(fp->out_p, bam_flag(b));
	*fp->out_p++ = '\t';

	/* RNAME */
	if (b->ref < -1 || b->ref >= fp->header->nref)
	    return -1;

	if (b->ref != -1) {
	    size_t l = strlen(fp->header->ref[b->ref].name);
	    if (end-fp->out_p < l+1) BF_FLUSH();
	    memcpy(fp->out_p, fp->header->ref[b->ref].name, l);
	    fp->out_p += l;
	} else {
	    if (end-fp->out_p < 2) BF_FLUSH();
	    *fp->out_p++ = '*';
	}
	*fp->out_p++ = '\t';

	/* POS */
	if (b->pos < -1) return -1;
	if (end-fp->out_p < 12) BF_FLUSH();
	fp->out_p = append_int(fp->out_p, b->pos+1); *fp->out_p++ = '\t';

	/* MAPQ */
	if (end-fp->out_p < 5) BF_FLUSH();
	fp->out_p = append_int(fp->out_p, bam_map_qual(b)); *fp->out_p++ = '\t';

	/* CIGAR */
	n = bam_cigar_len(b);dat = (uc *)bam_cigar(b);
	if (n < 0 ||
	    dat - (uc *)b + n*4 > b->blk_size + offsetof(bam_seq_t, ref))
	    return -1;
	for (i = 0; i < n; i++, dat+=4) {
	    uint32_t c = *(uint32_t *)dat;
	    if (end-fp->out_p < 13) BF_FLUSH();
	    fp->out_p = append_int(fp->out_p, c>>4);
	    *fp->out_p++="MIDNSHP=X"[c&15];
	}
	if (n==0) {
	    if (end-fp->out_p < 2) BF_FLUSH();
	    *fp->out_p++='*';
	}
	*fp->out_p++='\t';

	/* NRNM */
	if (b->mate_ref < -1 || b->mate_ref >= fp->header->nref)
	    return -1;

	if (b->mate_ref != -1) {
	    if (b->mate_ref == b->ref) {
		if (end-fp->out_p < 2) BF_FLUSH();
		*fp->out_p++ = '=';
	    } else {
		size_t l = strlen(fp->header->ref[b->mate_ref].name);
		if (end-fp->out_p < l+1) BF_FLUSH();
		memcpy(fp->out_p, fp->header->ref[b->mate_ref].name, l);
		fp->out_p += l;
	    }
	} else {
	    if (end-fp->out_p < 2) BF_FLUSH();
	    *fp->out_p++ = '*';
	}
	*fp->out_p++ = '\t';

	/* MPOS */
	if (end-fp->out_p < 12) BF_FLUSH();
	fp->out_p = append_int(fp->out_p, b->mate_pos+1); *fp->out_p++ = '\t';

	/* ISIZE */
	if (end-fp->out_p < 12) BF_FLUSH();
	fp->out_p = append_int(fp->out_p, b->ins_size); *fp->out_p++ = '\t';

	/* SEQ */
	n = (b->len+1)/2;
	dat = (uc *)bam_seq(b);

	/* BAM encoding */
//	while (n) {
//	    int l = end-fp->out_p < n ? end-fp->out_p : n;
//	    memcpy(fp->out_p, dat, l); fp->out_p += l;
//	    n -= l; dat += l;
//	    if (end == fp->out_p) BF_FLUSH();
//	}
	if (b->len != 0) {
	    if (end - fp->out_p < b->len + 3) BF_FLUSH();
	    if (end - fp->out_p < b->len + 3) {
		/* Extra long seqs need more regular checks */
		for (i = 0; i < b->len-1; i+=2) {
		    if (end - fp->out_p < 3) BF_FLUSH();
		    *fp->out_p++ = "=ACMGRSVTWYHKDBN"[*dat >> 4];
		    *fp->out_p++ = "=ACMGRSVTWYHKDBN"[*dat++ & 15];
		}
		if (i < b->len) {
		    if (end - fp->out_p < 3) BF_FLUSH();
		    *fp->out_p++ = "=ACMGRSVTWYHKDBN"[*dat >> 4];
		}
	    } else {
		unsigned char *cp = fp->out_p;
		int n = b->len & ~1;
		for (i = 0; i < n; i+=2) {
#ifdef ALLOW_UAC
		    *(int16_t *)cp = le_int2(code2base[*dat++]);
		    cp += 2;
#else
		    cp[0] = "=ACMGRSVTWYHKDBN"[*dat >> 4];
		    cp[1] = "=ACMGRSVTWYHKDBN"[*dat++ & 15];
		    cp += 2;
#endif
		}
		if (i < b->len) {
		    *cp++ = "=ACMGRSVTWYHKDBN"[*dat >> 4];
		}
		fp->out_p = cp;
	    }
	} else {
	    if (end - fp->out_p < 2) BF_FLUSH();
	    *fp->out_p++ = '*';
	}
	*fp->out_p++ = '\t';

	/* QUAL */
	n = b->len;
	if (b->len < 0) return -1;
	dat = (uc *)bam_qual(b);
	if (dat - (uc *)b + b->len > b->blk_size + offsetof(bam_seq_t, ref))
	    return -1;
	/* BAM encoding */
//	while (n) {
//	    int l = end-fp->out_p < n ? end-fp->out_p : n;
//	    memcpy(fp->out_p, dat, l); fp->out_p += l;
//	    n -= l; dat += l;
//	    if (end == fp->out_p) BF_FLUSH();
//	}
	if (b->len != 0) {
	    if (*dat == 0xff) {
		if (end - fp->out_p < 2) BF_FLUSH();
		*fp->out_p++ = '*';
		dat += b->len;
	    } else {
		if (end - fp->out_p < b->len + 3) BF_FLUSH();
		if (end - fp->out_p < b->len + 3) {
		    /* Long seqs */
		    for (i = 0; i < b->len; i++) {
			if (end - fp->out_p < 3) BF_FLUSH();
			*fp->out_p++ = *dat++ + '!';
		    }
		} else {
		    unsigned char *cp = fp->out_p;
		    i = 0;
#ifdef ALLOW_UAC
		    int n = b->len & ~3;
		    for (; i < n; i+=4) {
			//*cp++ = *dat++ + '!';
			*(uint32_t *)cp = *(uint32_t *)dat + 0x21212121;
			cp  += 4;
			dat += 4;
		    }
#endif
		    for (; i < b->len; i++) {
			*cp++ = *dat++ + '!';
		    }
		    fp->out_p = cp;
		}
	    }
	} else {
	    if (end - fp->out_p < 2) BF_FLUSH();
	    *fp->out_p++ = '*';
	}

	/* Auxiliary tags */
	auxh = NULL;
	while (0 == bam_aux_iter(b, &auxh, aux_key, &type, &val)) {
	    if (end - fp->out_p < 20) BF_FLUSH();
	    *fp->out_p++ = '\t';
	    *fp->out_p++ = aux_key[0];
	    *fp->out_p++ = aux_key[1];
	    *fp->out_p++ = ':';
	    *fp->out_p++ = type;
	    *fp->out_p++ = ':';
	    switch(type) {
	    case 'A':
		*fp->out_p++ = val.i;
		break;

	    case 'C':
		fp->out_p = append_int(fp->out_p, (uint8_t)val.i);
		break;

	    case 'c':
		fp->out_p = append_int(fp->out_p, (int8_t)val.i);
		break;

	    case 'S':
		fp->out_p = append_int(fp->out_p, (uint16_t)val.i);
		break;

	    case 's':
		fp->out_p = append_int(fp->out_p, (int16_t)val.i);
		break;

	    case 'I':
		fp->out_p = append_int(fp->out_p, (uint32_t)val.i);
		break;

	    case 'i':
		fp->out_p = append_int(fp->out_p, (int32_t)val.i);
		break;

	    case 'f':
		fp->out_p += sprintf((char *)fp->out_p, "%g", val.f);
		break;

	    case 'd':
		fp->out_p += sprintf((char *)fp->out_p, "%g", val.d);
		break;

	    case 'Z':
	    case 'H': {
		size_t l = strlen(val.s), l2;
		char *dat = val.s;
		do {
		    if (end - fp->out_p < l+2) BF_FLUSH();
		    l2 = MIN(l, end-fp->out_p);
		    memcpy(fp->out_p, dat, l2);
		    fp->out_p += l2;
		    l   -= l2;
		    dat += l2;
		} while (l);
		break;
	    }

	    case 'B': {
		uint32_t count = val.B.n, sz, j;
		unsigned char *s = val.B.s;
		*fp->out_p++ = val.B.t;

		/*
		 * Chew through count items 4000 at a time.
		 * This is because 4000*14 (biggest %g output plus comma?)
		 * is just shy of 64k, so we avoid buffer overflows.
		 */
		switch (val.B.t) {
		case 'C': case 'c': sz = 4; break;
		case 'S': case 's': sz = 6; break;
		default:            sz = 14; break;
		}

		for (j = 0; j < count; j += 4000) {
		    int i_start = j;
		    int i_end = j + 4000 < count ? j + 4000 : count;

		    if (end - fp->out_p < 5+(i_end-i_start)*sz) BF_FLUSH();

		    switch (val.B.t) {
			int i;
		    case 'C':
			for (i = i_start; i < i_end; i++, s++) {
			    *fp->out_p++ = ',';
			    fp->out_p = append_int(fp->out_p, (uint8_t)s[0]);
			}
			break;

		    case 'c':
			for (i = i_start; i < i_end; i++, s++) {
			    *fp->out_p++ = ',';
			    fp->out_p = append_int(fp->out_p, (int8_t)s[0]);
			}
			break;

		    case 'S':
			for (i = i_start; i < i_end; i++, s+=2) {
			    *fp->out_p++ = ',';
			    fp->out_p = append_int(fp->out_p,
						   (uint16_t)((s[0] << 0) +
							      (s[1] << 8)));
			}
			break;

		    case 's':
			for (i = i_start; i < i_end; i++, s+=2) {
			    *fp->out_p++ = ',';
			    fp->out_p = append_int(fp->out_p,
						   (int16_t)((s[0] << 0) +
							     (s[1] << 8)));
			}
			break;

		    case 'I':
			for (i = i_start; i < i_end; i++, s+=4) {
			    *fp->out_p++ = ',';
			    fp->out_p = append_uint(fp->out_p,
						    (uint32_t)((s[0] << 0) +
							       (s[1] << 8) +
							       (s[2] <<16) +
							       (s[3] <<24)));
			}
			break;

		    case 'i':
			for (i = i_start; i < i_end; i++, s+=4) {
			    *fp->out_p++ = ',';
			    fp->out_p = append_int(fp->out_p,
						   (int32_t)((s[0] << 0) +
							     (s[1] << 8) +
							     (s[2] <<16) +
							     (s[3] <<24)));
			}
			break;

		    case 'f': {
			union {
			    float f;
			    unsigned char c[4];
			} u;
			for (i = i_start; i < i_end; i++, s+=4) {
			    *fp->out_p++ = ',';
			    u.c[0] = s[0];
			    u.c[1] = s[1];
			    u.c[2] = s[2];
			    u.c[3] = s[3];
			    fp->out_p += sprintf((char *)fp->out_p, "%g", u.f);
			}
			break;
		    }

		    default:
			fprintf(stderr, "Unhandled sub-type of aux type B\n");
		    }
		}
		break;
	    }

	    default:
		fprintf(stderr, "Unhandled auxilary type '%c' in "
			"bam_put_seq()\n", type);
	    }
	}

	*fp->out_p++ = '\n';

    } else {
	/* BAM */
	unsigned char *end = fp->out + BGZF_BUFF_SIZE, *ptr;
	size_t to_write;
	uint32_t i32;
#ifndef ALLOW_UAC
	int name_len = bam_name_len(b);
#endif
#if !defined(ALLOW_UAC) || defined(SP_BIG_ENDIAN)
	uint32_t *cigar = bam_cigar(b);
#endif
#if defined(SP_BIG_ENDIAN)
	int i, n = bam_cigar_len(b);
#endif

#define CF_FLUSH() \
	do {			  \
	    if (bgzf_write(fp->fd, fp->level, fp->out, fp->out_p - fp->out)) \
		return -1;				     \
	    fp->out_p=fp->out;				\
	} while(0)

	/* If big endian, byte swap inline, write it out, and byte swap back */
	b->bin_packed  = (b->bin  << 16) | (b->map_qual << 8) | b->name_len;
	b->flag_packed = (b->flag << 16) | b->cigar_len;

#ifdef SP_BIG_ENDIAN
	b->ref         = le_int4(b->ref);
	b->pos         = le_int4(b->pos);
	b->bin_packed  = le_int4(b->bin_packed);
	b->flag_packed = le_int4(b->flag_packed);
	b->len         = le_int4(b->len);
	b->mate_ref    = le_int4(b->mate_ref);
	b->mate_pos    = le_int4(b->mate_pos);
	b->ins_size    = le_int4(b->ins_size);
	    
	for (i = 0; i < n; i++) {
	    cigar[i] = le_int4(cigar[i]);
	}
#endif

#ifdef ALLOW_UAC
	/* Room for fixed size bits + name */
	if (end - fp->out_p < 4) CF_FLUSH();
	to_write = b->blk_size;
	STORE_UINT32(fp->out_p, to_write);

        ptr = (unsigned char *)&b->ref;
#else
	/* Room for fixed size bits + name */
	if (end - fp->out_p < 36+257) CF_FLUSH();
	to_write = b->blk_size - (round4(name_len) - name_len);
	//to_write = b->blk_size;
	STORE_UINT32(fp->out_p, to_write);

        ptr = (unsigned char *)&b->ref;

	/* Do fixed size bits + name first */
	memcpy(fp->out_p, ptr, 32 + name_len);
	fp->out_p += 32 + name_len;
	to_write  -= 32 + name_len;
	ptr        = (unsigned char *)cigar;
#endif

        do {
            size_t blk_len = MIN(to_write, end - fp->out_p);
            memcpy(fp->out_p, ptr, blk_len);
            fp->out_p += blk_len;
            to_write  -= blk_len;
            ptr       += blk_len;

            if (to_write) {
                //printf("flushing %d+%d\n",
                //       (int)(ptr-(unsigned char *)&b->ref),
                //       (int)(fp->out_p-fp->out));
                CF_FLUSH();
            }
        } while(to_write > 0);

#ifdef SP_BIG_ENDIAN
	b->ref         = le_int4(b->ref);
	b->pos         = le_int4(b->pos);
	b->bin_packed  = le_int4(b->bin_packed);
	b->flag_packed = le_int4(b->flag_packed);
	b->len         = le_int4(b->len);
	b->mate_ref    = le_int4(b->mate_ref);
	b->mate_pos    = le_int4(b->mate_pos);
	b->ins_size    = le_int4(b->ins_size);
	    
	for (i = 0; i < n; i++) {
	    cigar[i] = le_int4(cigar[i]);
	}
#endif

	i32 = b->bin_packed;
	b->bin       = i32 >> 16;
	b->map_qual  = (i32 >> 8) & 0xff;
	b->name_len  = i32 & 0xff;

	i32          = b->flag_packed;
	b->flag      = i32 >> 16;
	b->cigar_len = i32 & 0xffff;
    }

    return 0;
}

/*
 * Writes a header block.
 * Returns 0 for success
 *        -1 for failure
 */
int bam_write_header(bam_file_t *out) {
    char *header, *hp, *htext;
    size_t hdr_size;
    int i, htext_len;

    if (sam_hdr_rebuild(out->header))
	return -1;

    htext = sam_hdr_str(out->header);
    htext_len = sam_hdr_length(out->header);

    hdr_size = 12 + htext_len+1;
    for (i = 0; i < out->header->nref; i++) {
	hdr_size += strlen(out->header->ref[i].name)+1 + 8;
    }
    if (NULL == (hp = header = malloc(hdr_size)))
	return -1;

    if (out->binary) {
	*hp++ = 'B'; *hp++ = 'A'; *hp++ = 'M'; *hp++ = 1;
	STORE_UINT32(hp, htext_len);
    }
    memcpy(hp, htext, htext_len);
    hp += htext_len;

    if (out->binary) {
	int i;

	STORE_UINT32(hp, out->header->nref);

	for (i = 0; i < out->header->nref; i++) {
	    size_t l = strlen(out->header->ref[i].name)+1;
	    STORE_UINT32(hp, l);

	    strcpy(hp, out->header->ref[i].name);
	    hp += l;

	    l = out->header->ref[i].len;
	    STORE_UINT32(hp, l);
	}
    }

    if (out->binary) {
	int len = hp-header;
	char *cp = header;

	while (len) {
	    int sz = BGZF_BUFF_SIZE < len ? BGZF_BUFF_SIZE : len;
	    if (bgzf_write(out->fd, out->level, cp, sz))
		return -1;
	    cp  += sz;
	    len -= sz;
	}
    } else {
	if (hp-header != write(out->fd, header, hp-header))
	    return -1;
    }

    free(header);

    return 0;
}



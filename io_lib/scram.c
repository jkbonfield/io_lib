/*
 * Copyright (c) 2013 Genome Research Ltd.
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
 * Multi-threading.

Trying to multi-thread BAM input is problematic as the BGZF boundaries
share no common ground with the sequence boundaries. Header and
sequences may be larger than a BGZF block, and the header may also be
in the same block as the sequence data.

Therefore the granularity of multi-threading needs to be the
compression and uncompression code.

Threadable sections:

1) zlib portion bgzf_more_output?
2) decode (may merge with 1)
3) bgzf_write


Similary for cram writing:
1) Building crecs
2) Encode + basic huff
3) Block compression (may merge with 2)


Ie input is hard to multi-thread
Output can be distributed to multiple threads and aggregated together
before sending on to FILE *fp.
 */


/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 *
 * A wrapper around SAM, BAM and CRAM I/O to give a unified interface.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <string.h>
#include <assert.h>

#include "io_lib/scram.h"

#define SCRAM_BUF_SIZE (1024*1024)

/*
 * Expands the input buffer.
 * Returns 0 on sucess
 *        -1 on failure
 */
static int scram_more_input(scram_fd *fd) {
    size_t avail = fd->alloc - fd->used;
    size_t l;

    l = fread(&fd->buf[fd->used], 1, avail, fd->fp);
    if (l <= 0)
	return -1;
    
    fd->used += l;
    return 0;
}

/*
 * Consumes a block of data from the input stream and returns a malloced
 * copy of it. The input buffer is then copied down (FIXME: inefficient).
 */
static unsigned char *scram_input_next_block(scram_fd *fd, size_t max_size,
					     size_t *out_size) {
    ssize_t l = MIN(max_size, fd->used);
    ssize_t i;
    unsigned char *r = NULL;

    if (max_size > fd->used) {
	scram_more_input(fd);
	if (fd->used == 0)
	    return NULL;
    }

    if (fd->b->binary) {
	uint32_t bsize;

	if (l < 19)
	    return NULL;
	bsize = fd->buf[16] + 256*fd->buf[17] + 1;
	fprintf(stderr, "block_size=%d\n", bsize);
	
	l = MIN(bsize, l);
    } else {
	for (i = l-1; i >= 0; i--) {
	    while (fd->buf[i] != '\n')
		i--;
	}
	assert(i >= 0);

	l = i;
    }

    if (!(r = malloc(l)))
	return NULL;
    memcpy(r, fd->buf, l);
    memcpy(fd->buf, &fd->buf[l], fd->used - l);
    fd->used -= l;

    if (out_size)
	*out_size = l;

    return r;
}

int scram_input_bam_block(scram_fd *fd) {
    size_t sz;
    unsigned char *r;

    if (!fd->is_bam)
	return -1;
    r = scram_input_next_block(fd, Z_BUFF_SIZE, &sz);
    if (!r)
	return -1;

//    if (fd->b->comp_p && fd->b->comp_p != fd->b->comp)
//	free(fd->b->comp_p);
    fd->b->comp_p = r;
    fd->b->comp_sz = sz;
    
    return 0;
}

/*
 * Opens filename.
 * If reading we initially try cram first and then bam/sam if that fails.
 * The exception is when reading from stdin, where bam/sam is first.
 *
 * If writing we look at the mode parameter:
 *     w  => SAM
 *     ws => SAM
 *     wb => BAM
 *     wc => CRAM
 *
 * Returns scram pointer on success
 *         NULL on failure
 */
scram_fd *scram_open(const char *filename, const char *mode) {
    char mode2[10];
    scram_fd *fd = calloc(1, sizeof(*fd));
    if (!fd)
	return NULL;

    fd->eof = 0;

    /* I/O buffer */
    fd->fp = NULL;
    fd->buf = NULL;
    fd->alloc = fd->used = 0;
    fd->pool = NULL;

    if (strcmp(filename, "-") == 0 && mode[0] == 'r'
	&& mode[1] != 'b' && mode[1] != 'c' && mode[1] != 's') { 
	int c;
	/*
	 * Crude auto-detection.
	 * First char @ = sam, 0x1f = bam (gzip), C = cram
	 * Headerless SAM will need explicit mode setting.
	 */
	c = fgetc(stdin);
	ungetc(c, stdin);

	if (c == '@')
	    sprintf(mode2, "rs%.7s", mode+1), mode = mode2;
	else if (c == 0x1f)
	    sprintf(mode2, "rb%.7s", mode+1), mode = mode2;
	else if (c == 'C')
	    sprintf(mode2, "rc%.7s", mode+1), mode = mode2;
    }

    if (*mode == 'r') {
	if (mode[1] != 'b' && mode[1] != 's') {
	    if (!(fd->c = sam_open(filename, mode))) {
		fprintf(stderr, "Error opening \"%s\"\n", filename);
		return NULL;
	    }
	    fd->is_bam = 0;
	    if (!(fd->hdr = sam_hdr_read(fd->c))) {
		fprintf(stderr, "Failed to read header\n");
		return NULL;
	    }
	    return fd;
	}

	if ((fd->b = bam_open(filename, mode))) {
	    fd->is_bam = 1;
	    return fd;
	}
	
	free(fd);
	return NULL;
    }

    /* For writing we cannot auto detect, so create the file type based
     * on the format in the mode string.
     */
    if (strncmp(mode, "wc", 2) == 0) {
	if (!(fd->c = sam_open(filename, mode))) {
	    fprintf(stderr, "Error opening \"%s\"\n", filename);
	    return NULL;
	}
	fd->is_bam = 0;
	return fd;
    }

    /* Otherwise assume bam/sam */
    if (!(fd->b = bam_open(filename, mode))) {
	free(fd);
	return NULL;
    }
    fd->is_bam = 1;
    return fd;
}

#if defined(CRAM_IO_CUSTOM_BUFFERING)
/*
 * Open CRAM file for reading via callbacks
 *
 * Returns scram pointer on success
 *         NULL on failure
 */
scram_fd *scram_open_cram_via_callbacks(
    char const * filename,
    cram_io_allocate_read_input_t   callback_allocate_function,
    cram_io_deallocate_read_input_t callback_deallocate_function,
    size_t const bufsize            
)
{
    scram_fd *fd = calloc(1, sizeof(*fd));
    if (!fd)
	return NULL;

    fd->eof = 0;

    /* I/O buffer */
    fd->fp = NULL;
    fd->buf = NULL;
    fd->alloc = fd->used = 0;
    fd->pool = NULL;

    if ((fd->c = cram_open_by_callbacks(filename,
					callback_allocate_function,
					callback_deallocate_function,
					bufsize))) 
    {
	cram_load_reference(fd->c, NULL);
	fd->is_bam = 0;
	return fd;
    }

    return NULL;
}
#endif

int scram_close(scram_fd *fd) {
    int r;

    if (fd->is_bam) {
	r = bam_close(fd->b);
    } else {
	r = sam_close(fd->c);
    }

    if (fd->pool)
	t_pool_destroy(fd->pool, 0);

    if (fd->hdr)
	sam_hdr_destroy(fd->hdr);

    if (fd->bc)
	bam_destroy1(fd->bc);

    free(fd);
    return r;
}

SAM_hdr *scram_get_header(scram_fd *fd) {
#ifdef __INTEL_COMPILER
    // avoids cmovne generation from icc 2015 (bug)
    return fd->is_bam && fd->b ? fd->b->header : fd->hdr;
#else
    return fd->is_bam ? fd->b->header : fd->hdr;
#endif
}

refs_t *scram_get_refs(scram_fd *fd) {
    return NULL;
    //return fd->is_bam ? NULL : fd->c->refs;
}

void scram_set_refs(scram_fd *fd, refs_t *refs) {
    return;
//    if (fd->is_bam)
//	return;
//    if (fd->c->refs)
//	refs_free(fd->c->refs);
//    fd->c->refs = refs;
//    if (refs)
//	refs->count++;
}

void scram_set_header(scram_fd *fd, SAM_hdr *sh) {
    if (fd->is_bam) {
	fd->b->header = sh;
    } else {
	fd->c->bam_header = sh;
    }
    sam_hdr_incr_ref(sh);

    fd->hdr = sh;
    sam_hdr_incr_ref(fd->hdr);
}

int scram_write_header(scram_fd *fd) {
    return fd->is_bam
	? bam_write_header(fd->b)
	: sam_hdr_write(fd->c, fd->c->bam_header);
}

int bam1_to_bam_seq(bam1_t *b, bam_seq_t **bsp_p) {
    bam_seq_t *bsp = *bsp_p;
    if (!bsp) {
	bsp = calloc(1, sizeof(*bsp));
	if (!bsp)
	    return -1;
	*bsp_p = bsp;
    }
    if (bsp->alloc < sizeof(*bsp) + b->l_data+1) {
	bsp->alloc = sizeof(*bsp) + b->l_data + 8;
	bam_seq_t *n = realloc(bsp, bsp->alloc);
	if (!n)
	    return -1;
	bsp = *bsp_p = n;
    }
    bsp->blk_size    = b->l_data + 32;
    bsp->pos         = b->core.pos;
    bsp->mate_pos    = b->core.mpos;
    bsp->ins_size    = b->core.isize;
    // As per raw BAM block below
    bsp->ref         = b->core.tid;
    bsp->pos_32      = b->core.pos; // bottom 32-bits
    bsp->name_len    = b->core.l_qname;
    bsp->map_qual    = b->core.qual;
    bsp->bin         = b->core.bin;
    bsp->cigar_len   = b->core.n_cigar;
    bsp->flag        = b->core.flag;
    bsp->len         = b->core.l_qseq;
    bsp->mate_ref    = b->core.mtid;
    bsp->mate_pos_32 = b->core.mpos;
    bsp->ins_size_32 = b->core.isize;

    memcpy(&bsp->data, b->data, b->l_data);
    (&bsp->data)[b->l_data] = 0; // io_lib's AUX end of tag marker
    return 0;
}

int bam_seq_to_bam1(bam_seq_t *bsp, bam1_t *b) {
    // NB: bam_set1 doesn't work as bam_seq(bsp) is 4-bit encoding while
    // bam_set1 uses ASCII.
    if (b->m_data < bsp->blk_size) {
	b->m_data = bsp->blk_size + 8;
	uint8_t *n = realloc(b->data, b->m_data);
	if (!n)
	    return -1;
	b->data = n;
    }
    b->l_data = bsp->blk_size - 32;
    b->core.pos = bsp->pos;
    b->core.mpos = bsp->mate_pos;
    b->core.isize = bsp->ins_size;
    b->core.tid = bsp->ref;
    b->core.l_qname = bsp->name_len;
    // SADLY this is costly in the main thread
    b->core.l_extranul = bsp->name_len - (strlen(bam_name(bsp))+1);
    b->core.qual = bsp->map_qual;
    b->core.bin = bsp->bin;
    b->core.n_cigar = bsp->cigar_len;
    b->core.flag = bsp->flag;
    b->core.l_qseq = bsp->len;
    b->core.mtid = bsp->mate_ref;
    
    memcpy(b->data, &bsp->data, b->l_data);
    return 0;
}

int scram_get_seq(scram_fd *fd, bam_seq_t **bsp) {
    if (fd->is_bam) {
	switch (bam_get_seq(fd->b, bsp)) {
	case 1:
	    return 0;

	case 0:
	    // FIXME: if we ever implement range queries for BAM this will
	    // need amendments to not claim a sub-range is invalid EOF.
	    fd->eof = fd->b->eof_block ? 1 : 2;
	    return -1;

	default:
	    fd->eof = -1; // err
	    return -1;
	}
    }

    // CRAM
    if (!fd->bc)
	fd->bc = bam_init1();
    int ret;
    if ((ret = sam_read1(fd->c, fd->hdr, fd->bc)) < 0) {
	fd->eof = ret == -1;
	return -1;
    }

    // convert bam1_t to bam_seq
    return bam1_to_bam_seq(fd->bc, bsp);
}

int scram_next_seq(scram_fd *fd, bam_seq_t **bsp) {
    return scram_get_seq(fd, bsp);
}

int scram_put_seq(scram_fd *fd, bam_seq_t *s) {
    if (fd->is_bam)
	return bam_put_seq(fd->b, s);

    if (!fd->bc)
	fd->bc = bam_init1();
    if (bam_seq_to_bam1(s, fd->bc) < 0)
	return -1;
    return sam_write1(fd->c, fd->hdr, fd->bc);
}

int scram_set_option(scram_fd *fd, enum cram_option opt, ...) {
    int r = 0;
    va_list args;

    va_start(args, opt);

    if (opt == CRAM_OPT_THREAD_POOL) {
	t_pool *p = va_arg(args, t_pool *);
	if (fd->is_bam) {
	    return bam_set_option(fd->b, BAM_OPT_THREAD_POOL, p);
	} else {
	    htsThreadPool tp = {
		.pool = p,
		.qsize = hts_tpool_size(p)*2
	    };
	    return hts_set_thread_pool(fd->c, &tp);
	}
    } else if (opt == CRAM_OPT_NTHREADS) {
	int nthreads = va_arg(args, int);
	if (nthreads > 1) {
	    if (fd->is_bam) {
		if (!(fd->pool = t_pool_init(nthreads*2, nthreads)))
		    return -1;

		return bam_set_option(fd->b, BAM_OPT_THREAD_POOL, fd->pool);
	    } else {
		return hts_set_threads(fd->c, nthreads);
	    }
	} else {
	    fd->pool = NULL;
	    return 0;
	}
//    // unsupported in HTSlib
//    } else if (opt == CRAM_OPT_BINNING) {
//	int bin = va_arg(args, int);
//
//	return fd->is_bam
//	    ? bam_set_option (fd->b,  BAM_OPT_BINNING, bin)
//	    : cram_set_option(fd->c, CRAM_OPT_BINNING, bin);
    } else if (opt == CRAM_OPT_IGNORE_CHKSUM) {
	int chk = va_arg(args, int);

	return fd->is_bam
	    ? bam_set_option(fd->b,  BAM_OPT_IGNORE_CHKSUM, chk)
	    : hts_set_opt(fd->c, CRAM_OPT_IGNORE_CHKSUM, chk);
    } else if (opt == CRAM_OPT_WITH_BGZIP_INDEX) {
        bgzi *idx = va_arg(args, bgzi *);
        if (fd->is_bam)
	    return bam_set_option(fd->b,  BAM_OPT_WITH_BGZIP_IDX, idx);
    } else if (opt == CRAM_OPT_OUTPUT_BGZIP_IDX) {
        char *idx_fn = va_arg(args, char *);
        if (fd->is_bam)
	    return bam_set_option(fd->b,  BAM_OPT_OUTPUT_BGZIP_IDX, idx_fn);
    } else if (opt == CRAM_OPT_EMBED_CONS) {
	return hts_set_opt(fd->c, CRAM_OPT_EMBED_REF, 2);
    }

    if (!fd->is_bam) {
	r = cram_set_voption(fd->c->fp.cram, opt, args);
    }

    va_end(args);

    return r;
}

/*! Returns the line number when processing a SAM file
 *
 * @return
 * Returns line number if input is SAM;
 *         0 for CRAM / BAM input.
 */
int scram_line(scram_fd *fd) {
    if (fd->is_bam)
	return fd->b->line;
    else
	return 0;
}

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

/*! Advises the memory allocator of CRAM usage patterns
 *
 * CRAM decoding will typically allocate & deallocate blocks for each
 * slice.  Under certain conditions this can cause a large number of
 * page faults where malloc gives a page back to the OS (free) and
 * then requests it again (the next malloc).  We could write our own
 * memory cache layer on top of malloc to keep track of previously
 * freed blocks, but it is complex in a multi-threaded environment and
 * arguably this is what malloc does anyway.
 *
 * Under GNU malloc we can simply tune it to avoid too many page faults.
 */
void scram_init(void) {
#if defined(HAVE_MALLOPT) && defined(M_MMAP_MAX)
    mallopt(M_MMAP_MAX, 0);
#endif
#if defined(HAVE_MALLOPT) && defined(M_TRIM_THRESHOLD)
    mallopt(M_TRIM_THRESHOLD, 100000000);
#endif
}


/* ---------------------------------------------------------------------------
 * A generalisation of CRAM_IO_CUSTOM_BUFFERING mode, used by libmaus.
 *
 * We use htslib's hfile plugin system to register a scramio type, and
 * then open the file via this new plugin backend.
 * This permits us to hijack all the read/write calls and call the
 * CRAM_IO_CUSTOM_BUFFERING function pointers.
 */

#include "hfile_internal.h"

typedef struct hFILE_scram {
    hFILE base;

    // CRAM_IO_CUSTOM_BUFFERING compatability
    // cram_io_{input,output}_t are tables of function pointers, used by
    // libmause2 to redirect the underlying I/O to read or write to buffers
    // instead of file streams.  See io_lib/cram.h for the definition.

    cram_io_input_t  *in_callbacks;
    cram_io_output_t *out_callbacks;

    //cram_io_deallocate_read_input_t  in_deallocate_func;
    //cram_io_deallocate_read_output_t out_deallocate_func;

    size_t bufsize;
} hFILE_scram;

static ssize_t scramio_read(hFILE *fpv, void *buffer, size_t nbytes) {
    hFILE_scram *fp = (hFILE_scram *)fpv;
    return fp->in_callbacks->fread_callback(buffer, 1, nbytes,
					    fp->in_callbacks->user_data);
    //return hread(&fp->base, buffer, nbytes);
}

static ssize_t scramio_write(hFILE *fpv, const void *buffer, size_t nbytes) {
    hFILE_scram *fp = (hFILE_scram *)fpv;
    return fp->out_callbacks->fwrite_callback((void *)buffer, 1, nbytes,
					      fp->out_callbacks->user_data);
    //return hwrite(&fp->base, buffer, nbytes);
}

static off_t scramio_seek(hFILE *fpv, off_t offset, int whence) {
    hFILE_scram *fp = (hFILE_scram *)fpv;
    // in or out?  How do I know?  add to struct maybe
    return fp->in_callbacks->fseek_callback(fp->in_callbacks->user_data,
					    offset, whence);
    //return hseek(&fp->base, offset, whence);
}

// FIXME: no ftell_callback?  It's in cram.h, but htslib lacks it.
// Let's hope we don't need this!
// I think it was only used in old io_lib to update buf_start so the tell
// coord matches seek coord, but I'm not 100% sure.

static int scramio_flush(hFILE *fpv) {
    return 0;
    //hFILE_scram *fp = (hFILE_scram *)fpv;
    //return hflush(&fp->base);
}

static int scramio_close(hFILE *fpv) {
    hFILE_scram *fp = (hFILE_scram *)fpv;
    free(fp);
    return 0;

    //int ret = hclose(&fp->base);
    //return ret;
}


static const struct hFILE_backend scramio_backend =
{
    scramio_read, scramio_write, scramio_seek, scramio_flush,
    scramio_close
};

static hFILE *hopen_scramio(const char *filename, const char *modestr) {
    // For now, just make this a NOP operation
    hFILE_scram *fp = calloc(1, sizeof(*fp));
    if (!fp)
	return NULL;

    // hfile_init fills out hFILE_scram->base only.  The rest is ours
    fp = (hFILE_scram *)hfile_init(sizeof(*fp), modestr, 1024*1024);
    if (!fp) {
	free(fp);
	return NULL;
    }
    fp->base.backend = &scramio_backend;

    return (hFILE *)fp;
}

/*
 * Open CRAM file for reading via callbacks
 *
 * Returns scram pointer on success
 *         NULL on failure
 */
scram_fd *scram_open_cram_via_callbacks(
    char const *filename,
    cram_io_allocate_read_input_t   callback_allocate_function,
    cram_io_deallocate_read_input_t callback_deallocate_function,
    size_t const bufsize            
)
{
    scram_fd *fd = calloc(1, sizeof(*fd));
    if (!fd)
	return NULL;

    // Register scheme handler.  It doesn't matter if it's done many times
    hfile_has_plugin("load-me");
    static const struct hFILE_scheme_handler handler =
        { hopen_scramio, hfile_always_local, "scramio", 10 };
    hfile_add_scheme_handler("scramio", &handler);

    // Open the filename with scheme scramio to get an hFILE_scram
    char fn[1024];
    snprintf(fn, 1024, "scramio:%s", filename);
    hFILE *hf = hopen(fn, "r");
    if (!hf)
	return NULL;

    // Extend the hFILE with the function callbacks given to us here.
    hFILE_scram *sio = (hFILE_scram *)hf;
    sio->in_callbacks = callback_allocate_function(filename, 0);
    sio->bufsize = bufsize;
    // TODO: close needs to call callback_deallocate_function, which means we
    // need to cache this somewhere to call it later.
    //fd->sio->in_deallocate_func = callback_deallocate_function;

    // Now expand the hFILE into an hts_file.
    fd->c = hts_hopen(hf, filename, "r");

    fd->is_bam = 0;
    if (!(fd->hdr = sam_hdr_read(fd->c))) {
	// FIXME: mem leak
	fprintf(stderr, "Failed to read header\n");
	return NULL;
    }

    return fd;
}

// Experiment to test scram_open_cram_via_callbacks()
cram_io_input_t *alloc_read_funcs(const char *fn, const int decompress) {
    cram_io_input_t *io = malloc(sizeof(*io));
    if (!io)
	return NULL;
    io->user_data = (void *)fopen(fn, "r");
    io->fread_callback  = (cram_io_C_FILE_fread_t)fread;
    io->fseek_callback  = (cram_io_C_FILE_fseek_t)fseeko;
    io->ftell_callback  = (cram_io_C_FILE_ftell_t)ftello;

    return io;
}

cram_io_input_t *free_read_funcs(cram_io_input_t *io) {
    fclose((FILE *)io->user_data);
    free(io);
    return NULL;
}

scram_fd *scram_open_(const char *filename, const char *mode) {
    return scram_open_cram_via_callbacks(filename,
					 alloc_read_funcs,
					 free_read_funcs,
					 1024*1024);
}

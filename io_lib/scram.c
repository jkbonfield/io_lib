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
    fd->line = 0;

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

    fd->is_bam = mode[1] == 'c' ? 0 : 1;

    if (*mode == 'r') {
	if (!(fd->sc= sam_open(filename, mode))) {
	    fprintf(stderr, "Error opening \"%s\"\n", filename);
	    return NULL;
	}
	sam_hdr_t *sh = sam_hdr_read(fd->sc);
	if (!(fd->hdr = sam_hdr_convert(sh))) {
	    fprintf(stderr, "Failed to read header\n");
	    return NULL;
	}
	sam_hdr_destroy(sh);
	fd->c = calloc(1, sizeof(*fd->c));
	if (!fd->c)
	    return NULL;
	fd->c->sc = fd->sc;
	fd->c->header = fd->hdr;

	// Count lines for SAM header
	if (fd->sc->format.format == sam) {
	    char *cp = fd->hdr->text->str;
	    char *cp_end = cp + fd->hdr->text->length;
	    while (cp < cp_end && (cp = strchr(cp, '\n'))) {
		fd->line++;
		cp++;
	    }
	}
	return fd;
    }

    /* For writing we cannot auto detect, so create the file type based
     * on the format in the mode string.
     */
    //if (strncmp(mode, "wc", 2) == 0) {
    if (*mode == 'w') {
	// Mode s is for socket in htslib, or SAM for iolib
	char mode2[100];
	snprintf(mode2, 100, "%c%s", mode[0],
		 mode[1] == 's' ? mode+2 : mode+1);
	if (!(fd->sc = sam_open(filename, mode2))) {
	    fprintf(stderr, "Error opening \"%s\"\n", filename);
	    return NULL;
	}
	fd->c = calloc(1, sizeof(*fd->c));
	if (!fd->c)
	    return NULL;
	fd->c->sc = fd->sc;
	return fd;
    }

    /* Otherwise unknown mode */
    return NULL;
}

int scram_close(scram_fd *fd) {
    int r = sam_close(fd->sc);

    if (fd->pool)
	t_pool_destroy(fd->pool, 0);

    if (fd->hdr)
	sam_hdr_free(fd->hdr);

    if (fd->bc)
	bam_destroy1(fd->bc);

    if (fd->c) {
	if (fd->c->index)
	    hts_idx_destroy(fd->c->index);
	free(fd->c);
    }

    free(fd);
    return r;
}


//// Use htslib's sam_hdr_parse and create a shadow struct matching the
//// io_lib name.  This is to permit some fields to be exposed.
//SAM_hdr *sam_hdr_parse_(const char *hdr, int len) {
//    SAM_hdr *h = sam_hdr_convert(sam_hdr_parse_htslib(hdr, len));
//    if (!h)
//	return NULL;
//    return h;
//}

SAM_hdr *scram_get_header(scram_fd *fd) {
    return fd->hdr;
}

refs_t *scram_get_refs(scram_fd *fd) {
    return fd->c->refs;
}

void scram_set_refs(scram_fd *fd, refs_t *refs) {
    fd->c->refs = refs;
}

void scram_set_header(scram_fd *fd, SAM_hdr *sh) {
    fd->sc->bam_header = sam_hdr_parse_htslib(sh->text->length,
					      sh->text->str);
    fd->c->header = sh;
    sam_hdr_incr_ref(sh);

    fd->hdr = sh;
}

int scram_write_header(scram_fd *fd) {
    // TODO: keep fd->hdr->hdr updated on the fly?
    if (fd->hdr->hdr)
	sam_hdr_destroy(fd->hdr->hdr);
    fd->hdr->hdr = sam_hdr_convert_to_htslib(fd->hdr);
    return sam_hdr_write(fd->sc, fd->hdr->hdr);
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

    // libmaus2 validates read names length matches strlen(name)+1.
    // So we must remove the padding bytes sadly.
    // NB: This meant libmaus2 also didn't work on systems without ALLOC_UAC.

    bsp->blk_size    = b->l_data + 32;
    bsp->pos         = b->core.pos;
    bsp->mate_pos    = b->core.mpos;
    bsp->ins_size    = b->core.isize;
    // As per raw BAM block below
    bsp->ref         = b->core.tid;
    bsp->pos_32      = b->core.pos; // bottom 32-bits
    bsp->name_len    = b->core.l_qname - b->core.l_extranul;
    bsp->map_qual    = b->core.qual;
    bsp->bin         = b->core.bin;
    bsp->cigar_len   = b->core.n_cigar;
    bsp->flag        = b->core.flag;
    bsp->len         = b->core.l_qseq;
    bsp->mate_ref    = b->core.mtid;
    bsp->mate_pos_32 = b->core.mpos;
    bsp->ins_size_32 = b->core.isize;

    //memcpy(&bsp->data, b->data, bsp->name_len);
    //memcpy(&bsp->data + bsp->name_len, b->data+b->core.l_qname,
    //	   b->l_data - b->core.l_qname);
    //bsp->blk_size -= b->core.l_qname - bsp->name_len;
    //(&bsp->data)[bsp->blk_size - 32] = 0; // io_lib's AUX end of tag marker

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
    b->core.l_qname = round4(bsp->name_len);
    b->core.l_extranul = (4-(bam_name_len(bsp)&3))&3;
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
    fd->line++;
    if (!fd->bc)
	fd->bc = bam_init1();
    int ret;
    if ((ret = sam_read1(fd->sc, fd->hdr->hdr, fd->bc)) < 0) {
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
    fd->line++;
    if (!fd->bc)
	fd->bc = bam_init1();
    if (bam_seq_to_bam1(s, fd->bc) < 0)
	return -1;

    if (fd->do_binning) {
	uint8_t *qual = bam_get_qual(fd->bc);
	uint32_t len = fd->bc->core.l_qseq;
	for (uint32_t i = 0; i < len; i++)
	    qual[i] = illumina_bin[qual[i]];
    }

    return sam_write1(fd->sc, fd->hdr->hdr, fd->bc);
}

int scram_set_option(scram_fd *fd, enum cram_option opt, ...) {
    int r = 0;
    va_list args;

    va_start(args, opt);

    if (opt == CRAM_OPT_THREAD_POOL) {
	t_pool *p = va_arg(args, t_pool *);
	htsThreadPool tp = {
	    .pool = p,
	    .qsize = hts_tpool_size(p)*2
	};
	return hts_set_thread_pool(fd->sc, &tp);
    } else if (opt == CRAM_OPT_NTHREADS) {
	int nthreads = va_arg(args, int);
	if (nthreads > 1) {
	    return hts_set_threads(fd->sc, nthreads);
	} else {
	    fd->pool = NULL;
	    return 0;
	}
    } else if (opt == CRAM_OPT_BINNING) {
	int bin = va_arg(args, int);
	fd->do_binning = bin;
    } else if (opt == CRAM_OPT_IGNORE_CHKSUM) {
	int chk = va_arg(args, int);
	hts_set_opt(fd->sc, CRAM_OPT_IGNORE_CHKSUM, chk);
    } else if (opt == CRAM_OPT_EMBED_REF) {
	return hts_set_opt(fd->sc, CRAM_OPT_EMBED_REF, 1);
    } else if (opt == CRAM_OPT_EMBED_CONS) {
	return hts_set_opt(fd->sc, CRAM_OPT_EMBED_REF, 2);
    } else if (opt == CRAM_OPT_PROFILE) {
	char *prof = va_arg(args, char *);
	int iprof = HTS_PROFILE_NORMAL;
	if (strcasecmp(prof, "fast") == 0) {
	    iprof = HTS_PROFILE_FAST;
	} else if (strcasecmp(prof, "normal") == 0) {
	    iprof = HTS_PROFILE_NORMAL;
	} else if (strcasecmp(prof, "small") == 0) {
	    iprof = HTS_PROFILE_SMALL;
	} else if (strcasecmp(prof, "archive") == 0) {
	    iprof = HTS_PROFILE_ARCHIVE;
	} else {
	    fprintf(stderr, "Unknown profile '%s', assuming 'normal'\n",
		    prof);
	}
	return hts_set_opt(fd->sc, HTS_OPT_PROFILE, iprof);
    }

    if (!fd->is_bam)
	r = cram_set_voption(fd->sc->fp.cram, opt, args);

    va_end(args);

    return r;
}

// Copied from htslib
int cram_set_option(cram_fd *fd, enum hts_fmt_option opt, ...) {
    int r;
    va_list args;

    va_start(args, opt);
    r = cram_set_voption(fd, opt, args);
    va_end(args);

    return r;
}

/*! Returns the line number when processing a SAM file
 *
 * @return
 * Returns line number in SAM (counting from 1),
 *         record number in BAM/CRAM (counting from 1)
 */
uint64_t scram_line(scram_fd *fd) {
    return fd->line;
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

int cram_index_load(cram_fd *fd, char const *fn) {
    fd->index = sam_index_load(fd->sc, fn);
    return fd->index ? 0 : -1;
}



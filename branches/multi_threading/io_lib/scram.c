/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 *
 * A wrapper around SAM, BAM and CRAM I/O to give a unified interface.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <string.h>

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
	    if ((fd->c = cram_open(filename, mode))) {
		cram_load_reference(fd->c, NULL);
		fd->is_bam = 0;
		return fd;
	    }
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
	if (!(fd->c = cram_open(filename, mode))) {
	    free(fd);
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

int scram_close(scram_fd *fd) {
    int r;

    if (fd->is_bam) {
	r = bam_close(fd->b);
	free(fd);
	return r;
    }

    r = cram_close(fd->c);
    free(fd);
    return r;
}

SAM_hdr *scram_get_header(scram_fd *fd) {
    return fd->is_bam ? fd->b->header : fd->c->header;
}

refs_t *scram_get_refs(scram_fd *fd) {
    return fd->is_bam ? NULL : fd->c->refs;
}

void scram_set_refs(scram_fd *fd, refs_t *refs) {
    if (fd->is_bam)
	return;
    if (fd->c->refs)
	refs_free(fd->c->refs);
    fd->c->refs = refs;
    if (refs)
	refs->count++;
}

void scram_set_header(scram_fd *fd, SAM_hdr *sh) {
    if (fd->is_bam) {
	fd->b->header = sh;
    } else {
	fd->c->header = sh;
    }
}

int scram_write_header(scram_fd *fd) {
    return fd->is_bam
	? bam_write_header(fd->b)
	: cram_write_SAM_hdr(fd->c, fd->c->header);
}

int scram_get_seq(scram_fd *fd, bam_seq_t **bsp) {
    if (fd->is_bam) {
	switch (bam_get_seq(fd->b, bsp)) {
	case 1:
	    return 0;

	case 0:
	    fd->eof = 1;

	default:
	    return -1;
	}
    }

    if (-1 == cram_get_bam_seq(fd->c, bsp)) {
	if (cram_eof(fd->c))
	    fd->eof = 1;
	return -1;
    }
    return 0;
}

int scram_next_seq(scram_fd *fd, bam_seq_t **bsp) {
    return scram_get_seq(fd, bsp);
}

int scram_put_seq(scram_fd *fd, bam_seq_t *s) {
    return fd->is_bam
	? bam_put_seq(fd->b, s)
	: cram_put_bam_seq(fd->c, s);
}

int scram_set_option(scram_fd *fd, enum cram_option opt, ...) {
    int r;
    va_list args;

    if (fd->is_bam)
	return 0;

    va_start(args, opt);
    r = cram_set_voption(fd->c, opt, args);
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

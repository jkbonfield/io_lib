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
 * If reading we initially try bam and then cram if that fails.
 *
 * If writing we look at the mode parameter:
 *     w  => SAM
 *     wb => BAM
 *     wc => CRAM
 *
 * Returns scram pointer on success
 *         NULL on failure
 */
scram_fd *scram_open(char *filename, char *mode) {
//    char mode2[10];
    scram_fd *fd = calloc(1, sizeof(*fd));
    if (!fd)
	return NULL;

    fd->eof = 0;

#if 0
    // Non functioning at present as bam IO is read() while cram is fread().
    if (strcmp(filename, "-") == 0 &&
	mode[0] == 'r' && mode[1] != 'b' && mode[1] != 'c') { 
	int c;
	/*
	 * Crude auto-detection.
	 * First char @ = sam, 0x1f = bam (gzip), C = cram
	 */
	c = fgetc(stdin);
	ungetc(c, stdin);

	if (c == 0x1f)
	    sprintf(mode2, "rb%.7s", mode+1), mode = mode2;
	else if (c == 'C')
	    sprintf(mode2, "rc%.7s", mode+1), mode = mode2;
    }
#endif

    if (*mode == 'r') {
	if (mode[1] == 'c') {
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
	
	if ((fd->c = cram_open(filename, mode))){ 
	    cram_load_reference(fd->c, NULL);
	    fd->is_bam = 0;
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

    /* Otherwise assume bam */
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
	bam_close(fd->b);
	free(fd);
	return 0;
    }

    r = cram_close(fd->c);
    free(fd);
    return r;
}

SAM_hdr *scram_get_header(scram_fd *fd) {
    return fd->is_bam ? fd->b->header : fd->c->SAM_hdr;
}

refs *scram_get_refs(scram_fd *fd) {
    return fd->is_bam ? NULL : fd->c->refs;
}

void scram_set_refs(scram_fd *fd, refs *refs) {
    if (fd->is_bam)
	return;
    fd->c->refs = refs;
}

void scram_set_header(scram_fd *fd, SAM_hdr *sh) {
    if (fd->is_bam) {
	fd->b->header = sh;
    } else {
	fd->c->SAM_hdr = sh;
    }
}

int scram_write_header(scram_fd *fd) {
    return fd->is_bam
	? bam_write_header(fd->b)
	: cram_write_SAM_hdr(fd->c, fd->c->SAM_hdr);
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


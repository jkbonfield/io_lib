#ifndef _CRAM_IO_H_
#define _CRAM_IO_H_

/*
 * Implements the low level CRAM I/O primitives.
 * This includes basic data types such as byte, int, ITF-8,
 * maps, bitwise I/O, etc.
 */

#include <stdint.h>
#include "io_lib/deflate_interlaced.h"

/* ----------------------------------------------------------------------------
 * Utility / debugging functions
 */
char *cram_block_method2str(enum cram_block_method m);
char *cram_content_type2str(enum cram_content_type t);

/* ----------------------------------------------------------------------------
 * Low level I/O functions for basic encoding and decoding of bits and bytes
 */

/*
 * Reads an integer in ITF-8 encoding from 'cp' and stores it in
 * *val.
 *
 * Returns the number of bytes read on success
 *        -1 on failure
 */
int itf8_decode(cram_fd *fd, int32_t *val);

/*
 * As above, but decoding from memory
 */
int itf8_get(char *cp, int32_t *val_p);

/*
 * Stores a value to memory in ITF-8 format.
 *
 * Returns the number of bytes required to store the number.
 * This is a maximum of 5 bytes.
 */
int itf8_put(char *cp, int32_t val);

/* ----------------------------------------------------------------------------
 * Mid level I/O functions for manipulating CRAM file structures:
 * Headers, containers, blocks, etc
 */

/*
 * Statistics gathering functions.
 */
cram_stats *cram_stats_create(void);
void cram_stats_add(cram_stats *st, int32_t val);
enum cram_encoding cram_stats_encoding(cram_stats *st);
void cram_stats_free(cram_stats *st);

/*
 * Reads a CRAM file definition structure.
 * Returns file_def ptr on success
 *         NULL on failure
 */
cram_file_def *cram_read_file_def(cram_fd *fd);
int cram_write_file_def(cram_fd *fd, cram_file_def *def);
void cram_free_file_def(cram_file_def *def);

/*
 * Reads the SAM header from the first CRAM data block.
 * Also performs minimal parsing to extract read-group
 * and sample information.
 *
 * Returns SAM hdr ptr on success
 *         NULL on failure
 */
cram_SAM_hdr *cram_read_SAM_hdr(cram_fd *fd);
cram_SAM_hdr *cram_create_SAM_hdr(char *str, size_t len);
int cram_write_SAM_hdr(cram_fd *fd, cram_SAM_hdr *hdr);
void cram_free_SAM_hdr(cram_SAM_hdr *hdr);


/*
 * Reads a container header.
 * Returns cram_container on success
 *         NULL on failure or no container left (fd->err == 0).
 */
cram_container *cram_read_container(cram_fd *fd);
void cram_free_container(cram_container *c);

/*
 * Reads a block from a cram file.
 * Returns cram_block pointer on success.
 *         NULL on failure
 */
cram_block *cram_read_block(cram_fd *fd);
void cram_free_block(cram_block *b);

/*
 * Decodes a CRAM block compression header.
 *         Returns header ptr on success
 *         NULL on failure
 */
#define CRAM_KEY(a,b) (((a)<<8)|((b)))
cram_block_compression_hdr *cram_decode_compression_header(cram_block *b);
void cram_free_compression_header(cram_block_compression_hdr *hdr);

/*
 * Decodes a CRAM mapped slice header block.
 * Returns slice header ptr on success
 *         NULL on failure
 */
cram_block_slice_hdr *cram_decode_slice_header(cram_block *b);
void cram_free_slice_header(cram_block_slice_hdr *hdr);


/*
 * Loads an entire slice.
 * Returns cram_slice ptr on success
 *         NULL on failure
 */
cram_slice *cram_read_slice(cram_fd *fd);
void cram_free_slice(cram_slice *s);

int cram_decode_slice(cram_fd *fd, cram_container *c, cram_slice *s,
		      bam_file_t *bfd);

/*
 * Loads a reference.
 * FIXME: use the sam_comp.cpp get_ref_base() equivalent to allow
 * random access within an indexed reference instead?
 *
 * Returns a ref_seq structure on success
 *         NULL on failure
 */
refs *load_reference(char *fn);
void free_refs(refs *r);
void refs2id(refs *r, bam_file_t *bfd);


/* ----------------------------------------------------------------------------
 * High level I/O functions for opening, getting or putting sequences and
 * closing CRAM files.
 */

/*
 * Opens a CRAM file for read (mode "rb") or write ("wb").
 * Returns file handle on success
 *         NULL on failure.
 */
cram_fd *cram_open(char *filename, char *mode);

/*
 * Closes a CRAM file.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_close(cram_fd *fd);

/* 
 * Sets options on the cram_fd. See CRAM_OPT_* definitions in cram_structs.h.
 * Use this immediately after opening.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_set_option(cram_fd *fd, enum cram_option opt, cram_opt *val);

void cram_uncompress_block(cram_block *b);

/*
 * Returns the next 'size' bytes from a block, or NULL if insufficient
 * data left. This is just a pointer into the block data and not an
 * allocated object, so do not free the result.
 */
char *cram_extract_block(cram_block *b, int size);

/*
 * Converts a cram in-memory record into a bam in-memory record. We
 * pass a pointer to a bam_seq_t pointer along with the a pointer to
 * the allocated size. These can initially be pointers to NULL and zero.
 *
 * This function will reallocate the bam buffer as required and update
 * bam_alloc accordingly, allowing it to be used within a loop
 * efficiently without needing to allocate new bam objects over and
 * over again.
 *
 * Returns the used size of the bam record on success
 *         -1 on failure.
 */
int cram_to_bam(bam_file_t *bfd, cram_fd *fd, cram_slice *s, cram_record *cr,
		int rec, bam_seq_t **bam, size_t *bam_alloc);

/*
 * Write iterator: put BAM format sequences into a CRAM file.
 * We buffer up a containers worth of data at a time.
 *
 * FIXME: break this into smaller pieces.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_put_bam_seq(cram_fd *fd, bam_seq_t *b);

#endif /* _CRAM_IO_H_ */

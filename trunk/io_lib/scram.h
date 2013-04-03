/*! \file
 * Generic SAM/BAM/CRAM interface.
 *
 * This file implements a higher level scram_*() API for programs that
 * wish to be file format agnostic.
 */

#ifndef _SCRAM_H_
#define _SCRAM_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include "io_lib/bam.h"
#include "io_lib/cram.h"

/*! The primary file handle for reading and writing. */
typedef struct {
    int is_bam;
    int eof;
    union {
	bam_file_t *b;
	cram_fd    *c;
    };
} scram_fd;

/*!@return
 * Returns 1 if we hit the end of file; 0 otherwise.
 */
#define scram_eof(fd) ((fd)->eof)


/*! Opens a file.
 *
 * If reading we look for the following mode parameters:
 * -    r  => Try SAM/BAM first, if fail try CRAM
 * -    rb => BAM
 * -    rc => CRAM
 *
 * If writing we look at the mode parameter:
 * -    w  => SAM
 * -    wb => BAM
 * -    wc => CRAM
 *
 * Additionally we can specify the compression level when writing
 * after the file type character, as 0 to 9. Eg "wb9" for maximum
 * compression of BAM or "wc0" for uncompressed CRAM.
 *
 * @return
 * Returns scram pointer on success
 *         NULL on failure
 */
scram_fd *scram_open(char *filename, char *mode);


/*! Closes a scram_fd handle
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int scram_close(scram_fd *fd);


/*! Returns the SAM_hdr struct.
 *
 * @return
 * The SAM_hdr struct on success; NULL on failure.
 */
SAM_hdr *scram_get_header(scram_fd *fd);


/*! Sets the SAM_hdr struct.
 *
 * Note that this sets the raw pointer and does not take an internal
 * copy of it. If you need to do this call sam_header_dup() first.
 */
void scram_set_header(scram_fd *fd, SAM_hdr *sh);


/*! Writes the SAM hdr.
 *
 * This calls the appropriate SAM, BAM or CRAM I/O function to write
 * out the SAM_hdr currently associated with this fd.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int scram_write_header(scram_fd *fd);


/*! Returns the reference sequence array.
 *
 * Note: this only works for CRAM files.
 *
 * @return
 * Returns the refs structure on success;
 *         NULL on failure.
 * 
 * After failure, check with scram_eof(fd) to see whether an genuine
 * error occurred or whether we hit the end of file.
 */
refs_t *scram_get_refs(scram_fd *fd);


/*! Sets the reference sequence array.
 *
 * Note: this only works for CRAM files.
 */
void scram_set_refs(scram_fd *fd, refs_t *refs);


/*! Fetches the next sequence and returns it in BAM format.
 *
 * This reads a new sequence line from fd and returns it in the BAM
 * in-memory format, regardless of whether the input file was SAM, BAM
 * or CRAM.
 *
 * @param bsp bsp is a pointer to a bam_seq_t*, as our usual bam_seq_t
 * structure pointer may be reallocated internally by this
 * function. It is permitted to pass in the address of a bam_seq_t*
 * that points to NULL. This behaviour differs to the Samtools API due
 * to the bam_seq_t structure being a single contiguous block of
 * memory instead of in two halves; the static and variable "data"
 * component.
 *
 * Note: For maximum speed of CRAM I/O you may wish to use the cram
 * specific layer and return cram_record objects instead.
 *
 * @return
 * Returns 0 on success and fills out bsp;
 *        -1 on failure
 */
int scram_get_seq(scram_fd *fd, bam_seq_t **bsp);

/*! Deprecated: please use scram_get_seq() instead */
int scram_next_seq(scram_fd *fd, bam_seq_t **bsp);


/*! Writes a BAM encoded bam_seq_t to fd.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int scram_put_seq(scram_fd *fd, bam_seq_t *s);


/*! Sets a CRAM option on fd.
 *
 * This is only supported for CRAM files currently.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int scram_set_option(scram_fd *fd, enum cram_option opt, ...);

#ifdef __cplusplus
}
#endif

#endif /* _SCRAM_H_ */

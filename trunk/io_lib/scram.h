#ifndef _SCRAM_H_
#define _SCRAM_H_

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include "io_lib/bam.h"
#include "io_lib/cram.h"

typedef struct {
    int is_bam;
    int eof;
    union {
	bam_file_t *b;
	cram_fd    *c;
    };
} scram_fd;

#define scram_eof(fd) ((fd)->eof)

scram_fd *scram_open(char *filename, char *mode);
int scram_close(scram_fd *fd);


SAM_hdr *scram_get_header(scram_fd *fd);
void scram_set_header(scram_fd *fd, SAM_hdr *sh);
int scram_write_header(scram_fd *fd);

refs *scram_get_refs(scram_fd *fd);
void scram_set_refs(scram_fd *fd, refs *refs);

int scram_next_seq(scram_fd *fd, bam_seq_t **bsp);
int scram_put_seq(scram_fd *fd, bam_seq_t *s);

int scram_set_option(scram_fd *fd, enum cram_option opt, cram_opt *val);

#endif /* _SCRAM_H_ */

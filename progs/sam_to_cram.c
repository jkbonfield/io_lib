#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

#include <io_lib/cram.h>


int main(int argc, char **argv) {
    cram_fd *out;
    cram_container *c;
    cram_SAM_hdr *hdr;
    size_t pos, pos2;
    bam_file_t *in;
    bam_seq_t *s = NULL;
    refs *refs;
    size_t bam_alloc = 0;

    /* opening */
    if (NULL == (in = bam_open(argv[1], "rb"))) {
	perror(argv[1]);
	return 1;
    }

    if (argc == 3) {
	refs = load_reference(argv[2]);
    } else {
	refs = NULL;
    }
    refs2id(refs, in);

    if (NULL == (out = cram_open("-", "wb"))) {
	fprintf(stderr, "Error opening CRAM file '%s'.\n", argv[1]);
	return 1;
    }
    out->refs = refs;


    /* SAM Header */
    if (NULL == (hdr = cram_create_SAM_hdr(in->header, in->header_len)))
	return 1;
    if (-1 == cram_write_SAM_hdr(out, hdr))
	return 1;
    out->SAM_hdr = hdr;
    //cram_free_SAM_hdr(hdr);
    
    /* Sequence iterators */
    while (bam_next_seq(in, &s) > 0) {
	//if (-1 == bam_put_seq(out, s))
	//return 1;
	if (-1 == cram_put_bam_seq(out, s)) {
	    fprintf(stderr, "Failed in cram_put_bam_seq()\n");
	    return 1;
	}
    }

    bam_close(in);
    cram_close(out);

    if (refs)
	free_refs(refs);

    if (s)
	free(s);

    return 0;
}

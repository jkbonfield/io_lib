/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 *
 * Support for CRAM index format: foo.cram.crai
 */

/*
 * The index is a gzipped tab-delimited text file with one line per slice.
 * The columns are:
 * 1: reference number (0 to N-1, as per BAM ref_id)
 * 2: reference position of 1st read in slice (1..?)
 * 3: number of reads in slice
 * 4: offset of container start (relative to end of SAM header, so 1st
 *    container is offset 0).
 * 5: slice number within container (ie which landmark).
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
#include "io_lib/zfio.h"

/*
 * Loads a CRAM .crai index into memory.
 * Returns 0 for success
 *        -1 for failure
 */
int cram_index_load(cram_fd *fd, char *fn) {
    char line[1024], fn2[PATH_MAX];
    zfp *fp;
    int nalloc = 0;
    int nind = 0;
    cram_index_entry *e = NULL;

    /* Check if already loaded */
    if (fd->index)
	return 0;

    fd->index = malloc(sizeof(*fd->index));
    if (!fd->index)
	return -1;

    sprintf(fn2, "%s.crai", fn);
    if (!(fp = zfopen(fn2, "r"))) {
	perror(fn2);
	return -1; 
    }

    while (zfgets(line, 1024, fp)) {
	if (nind >= nalloc) {
	    nalloc = nalloc ? nalloc*2 : 1024;
	    e = realloc(e, nalloc * sizeof(*e));
	    if (!e)
		return -1;
	}

	if (5 != sscanf(line, "%d\t%d\t%d\t%"PRId64"\t%d",
			&e[nind].refid,
			&e[nind].start,
			&e[nind].nseq,
			&e[nind].offset,
			&e[nind].slice)) {
	    fprintf(stderr, "Malformed index line\n");
	    return -1;
	}

	nind++;
    }
    zfclose(fp);

    fd->index->nslice = nind;
    fd->index->e = e;

    return 0;
}

void cram_index_free(cram_fd *fd) {
    if (!fd->index)
	return;

    if (fd->index->e)
	free(fd->index->e);
    free(fd->index);

    fd->index = NULL;
}

/*
 * Searches the index for the first slice overlapping a reference ID
 * and position, or one immediately preceeding it if none is found in
 * the index to overlap this position. (Our index may have missing
 * entries.)
 *
 * Returns the cram_index_entry pointer on sucess
 *         NULL on failure
 */
cram_index_entry *cram_index_query(cram_fd *fd, int refid, int pos) {
    int i, j, k;
    cram_index_entry *e;

    i = 0, j = fd->index->nslice-1;

    if (refid == -1) {
	i = fd->index->nslice-1;
	while (i > 0 && fd->index->e[i-1].refid == -1)
	    i--;
	goto skip;
    }

    for (k = j/2; k != i; k = (j-i)/2 + i) {
//	printf("%d/%d\t%d..%d..%d %d..%d..%d %d..%d..%d\n",
//	       refid, pos, i, j, k,
//	       fd->index->e[i].refid,
//	       fd->index->e[j].refid,
//	       fd->index->e[k].refid,
//	       fd->index->e[i].start,
//	       fd->index->e[j].start,
//	       fd->index->e[k].start);
	if (fd->index->e[k].refid > refid) {
	    j = k;
	    continue;
	}

	if (fd->index->e[k].refid < refid) {
	    i = k;
	    continue;
	}

	if (fd->index->e[k].start >= pos) {
	    j = k;
	    continue;
	}

	if (fd->index->e[k].start < pos) {
	    i = k;
	    continue;
	}
    }

    /* Special case for matching a start pos */
    if (i+1 < fd->index->nslice &&
	fd->index->e[i+1].start == pos &&
	fd->index->e[i+1].refid == refid)
	i++;

 skip:
    e = &fd->index->e[i];

    return e;
}

/*
 * Seek within a cram file.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek(cram_fd *fd, off_t offset, int whence) {
    char buf[65536];

    if (fseeko(fd->fp, offset, whence) == 0)
	return 0;

    if (!(whence == SEEK_CUR && offset >= 0))
	return -1;

    /* Couldn't fseek, but we're in SEEK_CUR mode so read instead */
    while (offset > 0) {
	int len = MIN(65536, offset);
	if (len != fread(buf, 1, len, fd->fp))
	    return -1;
	offset -= len;
    }

    return 0;
}


/*
 * Skips to a container overlapping the start coordinate listed in
 * cram_range.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek_to_refpos(cram_fd *fd, cram_range *r) {
    cram_index_entry *e;

    // Ideally use an index, so see if we have one.
    if ((e = cram_index_query(fd, r->refid, r->start))) {
	if (0 != cram_seek(fd, e->offset + fd->first_container, SEEK_SET))
	    return -1;
    } else {
	// Skip manually. NB: only works for first seek
	if (0 != cram_seek(fd, e->offset, SEEK_CUR))
	    return -1;
    }

    // FIXME: ignores slice number currently

    return 0;
}

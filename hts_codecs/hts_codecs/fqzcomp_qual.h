#ifndef FQZ_COMP_QUAL_H
#define FQZ_COMP_QUAL_H

/* Bit flags, deliberately mirroring BAM ones */
#define FQZ_FREVERSE 16
#define FQZ_FREAD2 128

/* Maximum length of an individual quality string.
 * If longer than this, simply break it up into
 * smaller portions.
 */
#ifndef MAX_SEQ
#  define MAX_SEQ 100000
#endif

/* A single record.
 * To compress we need to know the junction from one quality string to
 * the next (len, qual), whether it is first/second read and whether it is
 * reverse complemented (flags).
 */
typedef struct {
    int len;
    int qual; // FIXME: merge len and qual.  Artificial
    int flags;
} fqz_rec;

typedef struct {
    int num_records;
    fqz_rec *crecs;
} fqz_slice;

char *fqz_compress(int vers, fqz_slice *s, char *in, size_t uncomp_size,
                   size_t *comp_size, int level);
char *fqz_decompress(char *in, size_t comp_size, size_t *uncomp_size, int *lengths);

#endif

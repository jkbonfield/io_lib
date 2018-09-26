#ifndef FQZ_H
#define FQZ_H

char *fqz_compress(int vers, cram_slice *s, char *in, size_t uncomp_size,
		   size_t *comp_size, int level);
char *fqz_decompress(char *in, size_t comp_size, size_t *uncomp_size);

#endif /* FQZ_H */

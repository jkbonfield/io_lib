#ifndef ARITH_DYNAMIC_H
#define ARITH_DYNAMIC_H

unsigned char *arith_compress(unsigned char *in, unsigned int in_size,
			      unsigned int *out_size, int order);

unsigned char *arith_uncompress(unsigned char *in, unsigned int in_size,
				unsigned int *out_size);

unsigned char *arith_compress_to(unsigned char *in,  unsigned int in_size,
				 unsigned char *out, unsigned int *out_size,
				 int order);

unsigned char *arith_uncompress_to(unsigned char *in, unsigned int in_size,
				   unsigned char *out, unsigned int *out_sz);

unsigned int arith_compress_bound(unsigned int size, int order);

#endif /* ARITH_DYNAMIC_H */

/* Fuzz testing target. */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <ctype.h>

#include "htscodecs/fqzcomp_qual.h"
#include "htscodecs/fqzcomp_qual.c"
#undef NSYM
#define MODEL_256 // Prevent double definition
#include "htscodecs/arith_dynamic.c"

int LLVMFuzzerTestOneInput(uint8_t *in, size_t in_size) {
    size_t uncomp_size;
    char *uncomp = fqz_decompress((char *)in, in_size, &uncomp_size, NULL);
    if (uncomp)
	free(uncomp);
    
    return 0;
}

#ifdef NOFUZZ
#include <stdio.h>
#include <fcntl.h>
#include <errno.h>

#define BS 1024*1024
static unsigned char *load(char *fn, uint64_t *lenp) {
    unsigned char *data = NULL;
    uint64_t dsize = 0;
    uint64_t dcurr = 0;
    signed int len;

    int fd = open(fn, O_RDONLY);
    if (!fd) {
	perror(fn);
	return NULL;
    }

    do {
	if (dsize - dcurr < BS) {
	    dsize = dsize ? dsize * 2 : BS;
	    data = realloc(data, dsize);
	}

	len = read(fd, data + dcurr, BS);
	if (len > 0)
	    dcurr += len;
    } while (len > 0);

    if (len == -1) {
	perror("read");
    }
    close(fd);

    *lenp = dcurr;
    return data;
}

int main(int argc, char **argv) {
    uint64_t in_size;
    unsigned char *in = load(argv[1], &in_size);

    LLVMFuzzerTestOneInput(in, in_size);

    free(in);
    return 0;
}
#endif

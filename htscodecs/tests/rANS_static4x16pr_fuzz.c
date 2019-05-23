/*
Fuzz testing target.

Local instructions: compile, from a build subdir, with
/software/badger/opt/llvm/7.0.0/bin/clang -O3 -g ../../tests/rANS_static4x16pr_fuzz.c -I../.. ../../htscodecs/rANS_static4x16pr.c  -pthread -fsanitize=fuzzer,address /software/badger/opt/gcc/8.1.0/lib64/libstdc++.a

(This bizarrity is because our local clang install wasn't built with
C++ support.)

Run with:
    export ASAN_OPTIONS=allow_addr2line=true
    ./a.out corpus
or 
    ./a.out -detect_leaks=0 corpus

I generated corpus as a whole bunch of precompressed tiny inputs from
tests/dat/q4 for different compression modes.

For debugging purposes, we can compile a non-fuzzer non-ASAN build using
-DNOFUZZ which creates a binary we can debug on any libfuzzer generated
output using valgrind.  (The rans4x16 command line test won't quite work as
it's a slightly different input format with explicit sizes in the binary
stream.)
*/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>

#include "htscodecs/rANS_static4x16.h"
#include "htscodecs/rANS_static4x16pr.c"

int LLVMFuzzerTestOneInput(uint8_t *in, size_t in_size) {
    unsigned int uncomp_size;
    unsigned char *uncomp = rans_uncompress_4x16(in, in_size, &uncomp_size);
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

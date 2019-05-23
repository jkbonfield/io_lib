Htscodecs
=========

This repository implements the custom CRAM codecs used for "EXTERNAL"
block types.  These consist of two variants of the rANS codec (8-bit
and 16-bit renormalisation, with run-length encoding and bit-packing
also supported in the latter), a dynamic arithmetic coder, and custom
codecs for name/ID compression and quality score compression derived
from fqzcomp.


They come with small command line test tools to act as both
compression exploration programs and as part of the test harness.


Usage
-----

The library can be used as a git sub-module or as a completely
separate entity.  If you are attempting to make use of these codecs
within your own library, such as we do within Staden io_lib, it may be
useful to configure this with `--disable-shared --with-pic'.


API
---

Many functions just take an input buffer and size and return an output
buffer, setting *out_size with the decoded size.  NULL is returned for
error.  This buffer is malloced and is expected to be freed by the
caller.  These are the *`compress` and *`uncompress` functions.

A second variant sometimes exists where the output buffer is
optionally allocated by the caller (it may be NULL in which case it
has the same operation as above).  If specified, `*out_size` must also
be set to the allocated size of `out`.  These are the `compress_to`
and `uncompress_to` functions.

The compress size sometimes needs additional options.  For the rANS
and arithmetic coder this is the "order".  Values of 0 and 1 are
simple order-0 and order-1 entropy encoder, but this is a bit field
and the more advanced codecs have additional options to pass in order
(so it should really be renamed to flags).  See below.  Fqzcomp
requires more input data - also see below.  In all cases, sufficient
information is stored in the compressed byte stream such that the
decompression will work without needing these input paramaters.

Finally the various `compress_bound` functions give the size of buffer
needed to be allocated when compressing a block of data.


### Static rANS 4x8 (introduced in CRAM v3.0)

```
#include "htscodecs/rANS_static.h"

unsigned char *rans_compress(unsigned char *in, unsigned int in_size,
                             unsigned int *out_size, int order);
unsigned char *rans_uncompress(unsigned char *in, unsigned int in_size,
                               unsigned int *out_size);
```

No (un)compress_to functions exist for this older codec.


### Static rANS 4x16 with bit-pack/RLE (CRAM v3.1):

```
#include "htscodecs/rANS_static4x16.h"

unsigned int rans_compress_bound_4x16(unsigned int size, int order);
unsigned char *rans_compress_to_4x16(unsigned char *in,  unsigned int in_size,
                                     unsigned char *out, unsigned int *out_size,
                                     int order);
unsigned char *rans_compress_4x16(unsigned char *in, unsigned int in_size,
                                  unsigned int *out_size, int order);
unsigned char *rans_uncompress_to_4x16(unsigned char *in,  unsigned int in_size,
                                       unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_4x16(unsigned char *in, unsigned int in_size,
                                    unsigned int *out_size);
```

### Adaptive arithmetic coding (CRAM v3.1):

```
#include "htscodecs/arith_dynamic.h"

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
```

### Name tokeniser (CRAM v3.1):

```
#include "htscodecs/tokenise_name3.h"

uint8_t *encode_names(char *blk, int len, int level, int use_arith,
                      int *out_len, int *last_start_p);

uint8_t *decode_names(uint8_t *in, uint32_t sz, uint32_t *out_len);
```

This differs to the general purpose entropy encoders as it takes a
specific type of data.  The names should be newline or nul separated
for `encode_names`.  `decode_names` will alway return nul terminated
names, so you may need to swap these to newlines if you do round-trip
tests.

The compression level controls how hard it tries to find the optimum
compression method per internal token column.  By default it'll use
the rANS 4x16 codec, but with non-zero `use_arith` it'll use the
adaptive arithmetic coder instead.

If non-NULL, last_start_p can be used to point to a partial name if an
arbitrary block of names were supplied that don't end of a whole read
name. (Is this useful?  Probably not.)


### FQZComp Qual (CRAM v3.1):


```
#include "htscodecs/fqzcomp_qual.h"

#define FQZ_FREVERSE 16
#define FQZ_FREAD2 128

#define MAX_SEQ 100000

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
```

This is derived from the quality compression in fqzcomp.  The input
buffer is a concatenated block of quality strings, without any
separator.  In order to achieve maximum compression it needs to know
where these separators are, so they must be passed in via the
`fqz_rec` struct.

For now there are two length fields (just set them the same?), but
this will be changed.  The reason is one holds the actual length and
the other the quality string length as encoded in CRAM; this may be
subtly different depending on the CRAM features being encoded, but
it's this latter value that matters.

It can also be beneficial to supply per-record flags so fqzcomp can
determine whether orientation (complement strand) helps and whether
the READ1 vs READ2 quality distributions differ.  These are just
sub-fields from BAM FLAG.

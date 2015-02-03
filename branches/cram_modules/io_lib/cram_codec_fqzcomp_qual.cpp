#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "io_lib/cram_block_compression.h"
#include "clr.cdr"
#include "simple_model.h"

static const char *name(void) {
    return "fqzcomp-qual";
}

#define QMAX 64
#define QBITS 12
#define QSIZE (1<<QBITS)

unsigned char *compress_block(int level,
                              unsigned char *in,
                              size_t in_size,
                              size_t *out_size) {
    unsigned char *comp = (unsigned char *)malloc(in_size*1.1+1000);
    size_t i;
    unsigned int last = 0;
    RangeCoder rc;
    unsigned char q1 = 0, q2 = 0, q3 = 0;
    SIMPLE_MODEL<QMAX> *model_qual;
    model_qual = new SIMPLE_MODEL<QMAX>[QSIZE];
    int rle = 0;

    rc.output((char *)comp);
    rc.StartEncode();

    for (i = 0; i < in_size; i++) {
        unsigned char q = in[i] & (QMAX-1);
	assert(in[i] < QMAX && in[i] >= 0);
        model_qual[last].encodeSymbol(&rc, q);
        //last = ((MAX(q1, q2)<<6) + q) & (QSIZE-1);
        last = (((q1>>2)<<6) + q) & (QSIZE-1); // quicker learning

        //last += (q2==q3) << QBITS;
        //q3 = q2;

        //_mm_prefetch((const char *)&model_qual[last], _MM_HINT_T0);
        q2 = q1; q1 = q;
    }

    model_qual[last].encodeSymbol(&rc, QMAX-1); /* terminator */
    rc.FinishEncode();

    *out_size = rc.size_out();

    return comp;
}

unsigned char *uncompress_block(unsigned char *in,
                                size_t in_size,
                                size_t *out_size) {
    int block_size, data_size;
    unsigned char *uncomp = NULL;
    RangeCoder rc;
    size_t i, len = *out_size;
    unsigned char q1 = 0, q2 = 0;
    unsigned int last = 0;
    SIMPLE_MODEL<QMAX> *model_qual;
    model_qual = new SIMPLE_MODEL<QMAX>[QSIZE];

    uncomp = (unsigned char *)malloc(*out_size);

    rc.input((char *)in);
    rc.StartDecode();

    for (i = 0; i < len; i++) {
        unsigned char q = model_qual[last].decodeSymbol(&rc);
        uncomp[i] = q;

        //last = ((MAX(q1, q2)<<6) + q) & (QSIZE-1);
        last = (((q1>>2)<<6) + q) & (QSIZE-1); // quicker learning

        q2 = q1; q1 = q;
    }

    rc.FinishDecode();

    return uncomp;
}

static cram_compressor c = {
    FOUR_CC("FQZq"),
    1<<DS_QS, // quality only
    1.1,
    name,
    compress_block,
    uncompress_block,
};

extern "C" {
cram_compressor *cram_compressor_init(void) {
    return &c;
}
}

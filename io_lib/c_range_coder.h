/*
 * Note it is up to the calling code to ensure that no overruns on input and
 * output buffers occur.
 *
 * Call the input() and output() functions to set and query the current
 * buffer locations.
 *
 * Derived from Eugene Shelwien's work: http://ctxmodel.net/rem-37.html.
 */

#ifndef C_RANGER_CODER_H
#define C_RANGER_CODER_H

#define  DO(n)	   int _;for (_=0; _<n; _++)
#define  TOP	   (1<<24)

typedef unsigned char uc;

typedef struct {
    uint64_t low;
    uint32_t range, code;
    uc *in_buf;
    uc *out_buf;
} RangeCoder;

static inline void RC_SetInput(RangeCoder *rc, char *in) { rc->out_buf = rc->in_buf = (uc *)in; }
static inline void RC_SetOutput(RangeCoder *rc, char *out) { rc->in_buf = rc->out_buf = (uc *)out; }
static inline char *RC_GetInput(RangeCoder *rc) { return (char *)rc->in_buf; }
static inline char *RC_GetOutput(RangeCoder *rc) { return (char *)rc->out_buf; }
static inline size_t RC_OutSize(RangeCoder *rc) { return rc->out_buf - rc->in_buf; }
static inline size_t RC_InSize(RangeCoder *rc) { return rc->in_buf - rc->out_buf; }

static inline void RC_StartEncode(RangeCoder *rc)
{ 
    rc->low=0;  
    rc->range=(uint32_t)-1; 
}

static inline void RC_StartDecode(RangeCoder *rc)
{ 
    rc->code = rc->low = 0;  
    rc->range = (uint32_t)-1;
    DO(8) rc->code = (rc->code<<8) | *rc->in_buf++;
}

static inline void RC_FinishEncode(RangeCoder *rc) 
{ 
    DO(8) (*rc->out_buf++ = rc->low>>56), rc->low<<=8; 
}

static inline void RC_FinishDecode(RangeCoder *rc) {}

#if DIV_RCP
static uint64_t rcp[65536];
static inline void build_rcp_freq(void) {
    int i;
    for (i = 1; i < 65536; i++)
	rcp[i] = (((uint64_t)1<<32)) / i;
}
#else
static inline void build_rcp_freq(void) {}
#endif

static inline void RC_Encode (RangeCoder *rc, uint32_t cumFreq, uint32_t freq, uint32_t totFreq) 
{
    //fprintf(stderr, "                       RC %d+%d of %d\n", cumFreq, freq, totFreq);

    // division-less doesn't help much in this case as only a single division anyway
    //uint32_t d = (rcp[totFreq] * rc->range)>>32;
    //d += (rc->range - totFreq * d >= totFreq);
    //rc->low  += cumFreq * (rc->range = d);

    rc->low  += cumFreq * (rc->range/= totFreq);
    rc->range*= freq;

    if (cumFreq + freq > totFreq) {
	fprintf(stderr, "cumFreq %d + freq %d > tot %d\n", cumFreq, freq, totFreq);
	abort();
    }

    while( rc->range<TOP ) {
	// range = 0x00ffffff..
	// low/high may be matching
	//       eg 88332211/88342211 (range 00010000)
	// or differing
	//       eg 88ff2211/89002211 (range 00010000)
	//
	// If the latter, we need to reduce range down
	// such that high=88ffffff.
	// Eg. top-1      == 00ffffff
	//     low|top-1  == 88ffffff
	//     ...-low    == 0000ddee
	if ( (uc)((rc->low^(rc->low+rc->range))>>56) ) 
	    rc->range = (((uint32_t)(rc->low)|(TOP-1))-(uint32_t)(rc->low));
	*rc->out_buf++ = rc->low>>56, rc->range<<=8, rc->low<<=8;
    }
}

static inline uint32_t RC_GetFreq (RangeCoder *rc, uint32_t totFreq) {
#if DIV_RCP
    // Replacing two divisions by one is beneficial as they get stuck waiting.
    // 2.53s
    uint32_t d = (rcp[totFreq] * rc->range)>>32;
    d += (rc->range - totFreq * d >= totFreq);
    return rc->code/(rc->range=d);
#else
    // 2.67s
    return rc->code/(rc->range/=totFreq);
#endif
}

static inline void RC_Decode (RangeCoder *rc, uint32_t cumFreq, uint32_t freq, uint32_t totFreq) 
{
    //fprintf(stderr, "                       RC %d+%d of %d\n", cumFreq, freq, totFreq);

    uint32_t temp = cumFreq*rc->range;
    rc->low  += temp;
    rc->code -= temp;
    rc->range*= freq;
 
    while( rc->range<TOP ) {
	if ( (uc)((rc->low^(rc->low+rc->range))>>56) ) 
	    rc->range = (((uint32_t)(rc->low)|(TOP-1))-(uint32_t)(rc->low));
	rc->code = (rc->code<<8) | *rc->in_buf++, rc->range<<=8, rc->low<<=8;
    }
}

#endif /* C_RANGER_CODER_H */

#include <stdint.h>
#include "c_range_coder.h"

/*
 *--------------------------------------------------------------------------
 * A simple frequency model.
 *
 * Define NSYM to be an integer value before including this file.
 * It will then generate types and functions specific to that
 * maximum number of symbols.
 *
 * This keeps a list of symbols and their frequencies, approximately
 * sorted by symbol frequency. We allow for a single symbol to periodically
 * move up the list when emitted, effectively doing a single step of
 * bubble sort periodically. This means it's largely the same complexity
 * irrespective of alphabet size.
 * It's more efficient on strongly biased distributions than random data.
 *
 * There is no escape symbol, so the model is tailored to relatively
 * stationary samples (although we do have occasional normalisation to
 * avoid frequency counters getting too high).
 *--------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
// Bits we want included once only - constants, types, etc
#ifndef C_SIMPLE_MODEL_H
#define C_SIMPLE_MODEL_H

#define MAX_FREQ (1<<16)-32
#define PASTE3(a,b,c) a##b##c
#define SIMPLE_MODEL(a,b) PASTE3(SIMPLE_MODEL,a,b)
#define STEP 8
typedef struct {
    uint16_t Freq;
    uint16_t Symbol;
} SymFreqs;
#endif /* C_SIMPLE_MODEL_H */


//-----------------------------------------------------------------------------
// Bits we regenerate for each NSYM value.

typedef struct {
    uint32_t TotFreq;  // Total frequency
    uint32_t BubCnt;   // Periodic counter for bubble sort step

    // Array of Symbols approximately sorted by Freq. 
    SymFreqs sentinel, F[NSYM+1];
} SIMPLE_MODEL(NSYM,_);


static inline void SIMPLE_MODEL(NSYM,_init)(SIMPLE_MODEL(NSYM,_) *m, int max_sym) {
    int i;
    
    for (i=0; i<max_sym; i++) {
	m->F[i].Symbol = i;
	m->F[i].Freq   = 1;
    }
    for (; i<NSYM; i++) {
	m->F[i].Symbol = i;
	m->F[i].Freq   = 0;
    }

    m->TotFreq         = max_sym;
    m->sentinel.Symbol = 0;
    m->sentinel.Freq   = MAX_FREQ; // Always first; simplifies sorting.
    m->BubCnt          = 0;

    m->F[NSYM].Freq    = 0; // terminates normalize() loop. See below.
}


static inline void SIMPLE_MODEL(NSYM,_normalize)(SIMPLE_MODEL(NSYM,_) *m) {
    SymFreqs *s;

    /* Faster than F[i].Freq for 0 <= i < NSYM */
    m->TotFreq=0;
    for (s = m->F; s->Freq; s++) {
	s->Freq -= s->Freq>>1;
	m->TotFreq += s->Freq;
    }
}

static inline void SIMPLE_MODEL(NSYM,_encodeSymbol)(SIMPLE_MODEL(NSYM,_) *m,
                                                    RangeCoder *rc, uint16_t sym) {
    SymFreqs *s = m->F;
    uint32_t AccFreq  = 0;

    while (s->Symbol != sym)
	AccFreq += s++->Freq;

//    if (s->Freq == 0) {
//        fprintf(stderr, "sym=%d, s->Freq=%d, idx=%d\n", sym, s->Freq, s-m->F);
//        abort();
//    }

    RC_Encode(rc, AccFreq, s->Freq, m->TotFreq);
    s->Freq    += STEP;
    m->TotFreq += STEP;

    if (m->TotFreq > MAX_FREQ)
	SIMPLE_MODEL(NSYM,_normalize)(m);

    /* Keep approx sorted */
    if (((++m->BubCnt&15)==0) && s[0].Freq > s[-1].Freq) {
	SymFreqs t = s[0];
	s[0] = s[-1];
	s[-1] = t;
    }
}

static inline uint16_t SIMPLE_MODEL(NSYM,_decodeSymbol)(SIMPLE_MODEL(NSYM,_) *m, RangeCoder *rc) {
    SymFreqs* s = m->F;
    uint32_t freq = RC_GetFreq(rc, m->TotFreq);
    uint32_t AccFreq;

    for (AccFreq = 0; (AccFreq += s->Freq) <= freq; s++);
    AccFreq -= s->Freq;

    RC_Decode(rc, AccFreq, s->Freq, m->TotFreq);
    s->Freq    += STEP;
    m->TotFreq += STEP;

    if (m->TotFreq > MAX_FREQ)
	SIMPLE_MODEL(NSYM,_normalize)(m);

    /* Keep approx sorted */
    if (((++m->BubCnt&15)==0) && s[0].Freq > s[-1].Freq) {
	SymFreqs t = s[0];
	s[0] = s[-1];
	s[-1] = t;
	return t.Symbol;
    }

    return s->Symbol;
}

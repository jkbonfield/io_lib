#ifndef _CRAM_STATS_H_
#define _CRAM_STATS_H_

#include "io_lib/hash_table.h"

#define MAX_STAT_VAL 1024
//#define MAX_STAT_VAL 16
typedef struct {
    int freqs[MAX_STAT_VAL];
    HashTable *h;
    int nsamp; // total number of values added
    int nvals; // total number of unique values added
} cram_stats;

cram_stats *cram_stats_create(void);
void cram_stats_add(cram_stats *st, int32_t val);
void cram_stats_del(cram_stats *st, int32_t val);
void cram_stats_dump(cram_stats *st);

#endif

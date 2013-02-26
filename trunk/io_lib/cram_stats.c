#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <ctype.h>

#include "io_lib/cram.h"
#include "io_lib/os.h"

cram_stats *cram_stats_create(void) {
    cram_stats *st = calloc(1, sizeof(cram_stats));
    return st;
}

void cram_stats_add(cram_stats *st, int32_t val) {
    st->nsamp++;

    //assert(val >= 0);

    if (val < MAX_STAT_VAL && val >= 0) {
	st->freqs[val]++;
    } else {
	HashItem *hi;

	if (!st->h) {
	    st->h = HashTableCreate(2048, HASH_DYNAMIC_SIZE|HASH_NONVOLATILE_KEYS|HASH_INT_KEYS);
	}

	if ((hi = HashTableSearch(st->h, (char *)val, 4))) {
	    hi->data.i++;
	} else {
	    HashData hd;
	    hd.i = 1;
	    HashTableAdd(st->h, (char *)val, 4, hd, NULL);
	}
    }
}

void cram_stats_del(cram_stats *st, int32_t val) {
    st->nsamp--;

    //assert(val >= 0);

    if (val < MAX_STAT_VAL && val >= 0) {
	st->freqs[val]--;
    } else if (st->h) {
	HashItem *hi;

	if ((hi = HashTableSearch(st->h, (char *)val, 4))) {
	    if (--hi->data.i == 0)
		HashTableDel(st->h, hi, 0);
	} else {
	    fprintf(stderr, "Failed to remove val %d from cram_stats\n", val);
	    st->nsamp++;
	}
    } else {
	fprintf(stderr, "Failed to remove val %d from cram_stats\n", val);
	st->nsamp++;
    }
}

void cram_stats_dump(cram_stats *st) {
    int i;
    fprintf(stderr, "cram_stats:\n");
    for (i = 0; i < MAX_STAT_VAL; i++) {
	if (!st->freqs[i])
	    continue;
	fprintf(stderr, "\t%d\t%d\n", i, st->freqs[i]);
    }
    if (st->h) {
	HashIter *iter=  HashTableIterCreate();
	HashItem *hi;

	while ((hi = HashTableIterNext(st->h, iter))) {
	    fprintf(stderr, "\t%d\t%d\n", (int)hi->key, (int)hi->data.i);
	}
	HashTableIterDestroy(iter);
    }
}

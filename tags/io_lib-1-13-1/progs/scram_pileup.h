#ifndef _PILEUP_H_
#define _PILEUP_H_

#include "io_lib/scram.h"

typedef struct pileup {
    struct pileup *next;  // A link list, for active seqs
    void *cd;		  // General purpose per-seq client-data

    bam_seq_t *b;	  // Bam entry associated with struct
    unsigned char *b_qual;// cached bam_qual
    unsigned char *b_seq; // cached bam_seq
    uint32_t *b_cigar;    // cached bam_cigar
    int  b_strand;        // 0 => fwd, 1 => rev

    int  pos;             // Current unpadded position in seq
    int  nth;		  // nth base at unpadded position 'pos'
    int  seq_offset;      // Current base position in s->seq[] array.

    int  cigar_ind;       // Current location in s->alignment cigar str
    int  cigar_op;        // Current cigar operation
    int  cigar_len;       // Remaining length of this cigar op

    int  first_del;       // Used when first base is a deletion

    int  eof;		  // True if this sequence has finished
    int  qual;            // Current qual (for active seq only)
    char base;		  // Current base (for active seq only)
    char start;		  // True if this is a new sequence
    char ref_skip;        // True if the cause of eof or start is cigar N
    char padding;         // True if the base was added due to another seq
} pileup_t;

int pileup_loop(scram_fd *fp,
		int (*seq_init)(void *client_data,
				scram_fd *fp,
				pileup_t *p),
		int (*seq_add)(void *client_data,
			       scram_fd *fp,
			       pileup_t *p,
			       int depth,
			       int pos,
			       int nth,
			       int is_insert),
		void *client_data);

#endif

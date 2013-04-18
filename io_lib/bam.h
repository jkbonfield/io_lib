/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2010-3
 */

/*! \file
 * The primary SAM/BAM API.
 *
 * Consider using scram.h if you wish to also have support for CRAM.
 */ 

#ifndef _BAM_H_
#define _BAM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include <zlib.h>

#include "io_lib/os.h"
#include "io_lib/hash_table.h"
#include "io_lib/sam_header.h"

/* BAM header structs */
typedef struct tag_list {
    char *value; /* NULL => end of tags */
    int key;
    int length;  
} tag_list_t;

/* The main bam sequence struct */
typedef struct bam_seq_s {
    uint32_t alloc; /* total size of this struct + 'data' onwards */
    uint32_t blk_size;

    /* The raw bam block follows, in same order as on the disk */
    /* This is the analogue of a bam1_core_t in samtools */
    int32_t  ref;
    int32_t  pos;

    union {
	struct {
	    uint32_t name_len:8, map_qual:8, bin:16;
	};
	uint32_t bin_packed;
    };
    union {
	struct {
	    uint32_t cigar_len:16, flag:16;
	};
	uint32_t flag_packed;
    };

    int32_t  len;
    int32_t  mate_ref;
    int32_t  mate_pos;
    int32_t  ins_size;

    /* Followed by arbitrary bytes of packed data */
    unsigned char data; /* unknown size */
} bam_seq_t;


/* Auxillary field handling */
typedef union {
    char    *s;
    int      i;
    uint64_t i64;
    float    f;
    double   d;
    struct {
	int n, t;
	unsigned char *s;
    } B;
} bam_aux_t;


/*
 * Our bam stream consists of a zlib gzFile stream and a buffer for it to
 * output to. This allows us to call small bam_read requests while
 * translating these to fewer, larger, gzread calls. The overhead is
 * therefore minimal.
 */
#define Z_BUFF_SIZE 65536    /* Max size of a zlib block */
#define BGZF_BUFF_SIZE 65400 /* Max size of a BGZF block, 65477 actual */
typedef struct {
    int fd, mode, binary, level;
    z_stream s;

    unsigned char in[Z_BUFF_SIZE];
    unsigned char *in_p;
    size_t in_sz;

    unsigned char out[Z_BUFF_SIZE];
    unsigned char *out_p;
    size_t out_sz;

    /* BAM specifics */
    int32_t next_len;

    SAM_hdr *header;      /* Parsed SAM header */

    /* Cached bam_seq_t, to avoid excessive mallocs */
    bam_seq_t *bs;
    int bs_size;

    /* Boolean to indicate if we've finished the most recent z stream */
    int z_finish;

    /* Indicates whether gzipped, and if so with bgzf extra fields */
    int gzip;
    int bgzf; 

    /* Whether BAM or SAM format */
    int bam;

    /* If true, skip auxillary field parsing while reading SAM */
    int no_aux;

    /* line number (when in SAM mode) */
    int line;

    /* Static avoidance: used in sam_next_seq() */
    unsigned char *sam_str;
    size_t alloc_l;
} bam_file_t;

/* Decoding the above struct */
#define bam_blk_size(b)       ((b)->blk_size)
#define bam_set_blk_size(b,v) ((b)->blk_size = (v))

#define bam_ref(b)       ((b)->ref)
#define bam_pos(b)       ((b)->pos)
#define bam_mate_ref(b)  ((b)->mate_ref)
#define bam_mate_pos(b)  ((b)->mate_pos)
#define bam_ins_size(b)  ((b)->ins_size)
#define bam_cigar_len(b) ((b)->cigar_len)
#define bam_name_len(b)  ((b)->name_len)
#define bam_map_qual(b)  ((b)->map_qual)
#define bam_bin(b)       ((b)->bin)
#define bam_flag(b)      ((b)->flag)
#define bam_seq_len(b)   ((b)->len)
#define bam_strand(b)    ((bam_flag((b)) & BAM_FREVERSE) != 0)
#define bam_mstrand(b)   ((bam_flag((b)) & BAM_FMREVERSE) != 0)

#define bam_set_ref(b,v)        ((b)->ref = (v))
#define bam_set_pos(b,v)        ((b)->pos = (v))
#define bam_set_mate_ref(b,v)   ((b)->mate_ref = (v))
#define bam_set_mate_pos(b,v)   ((b)->mate_pos = (v))
#define bam_set_ins_size(b,v)   ((b)->ins_size = (v))
#define bam_set_cigar_len(b, v) ((b)->cigar_len = (v))
#define bam_set_name_len(b,v)   ((b)->name_len = (v))
#define bam_set_map_qual(b,v)   ((b)->map_qual = (v))
#define bam_set_bin(b,v)        ((b)->bin = (v));
#define bam_set_flag(b,v)       ((b)->flag = (v));
#define bam_set_seq_len(b,v)    ((b)->len = (v));

#define bam_name(b)      ((char *)(&(b)->data))
#ifdef ALLOW_UAC
#define bam_cigar(b)     ((uint32_t *)(bam_name((b)) + bam_name_len((b))))
#else
#define bam_cigar(b)     ((uint32_t *)(bam_name((b)) + round4(bam_name_len((b)))))
#endif
#define bam_seq(b)       (((char *)bam_cigar((b))) + 4*bam_cigar_len(b))
#define bam_qual(b)      (bam_seq(b) + (int)(((b)->len+1)/2))
#define bam_aux(b)       (bam_qual(b) + (b)->len)

/* Rounds up to the next multiple of 4 */
#define round4(v) (((v-1)&~3)+4)

/* BAM flags */
#define BAM_FPAIRED        1
#define BAM_FPROPER_PAIR   2
#define BAM_FUNMAP         4
#define BAM_FMUNMAP        8
#define BAM_FREVERSE      16
#define BAM_FMREVERSE     32
#define BAM_FREAD1        64
#define BAM_FREAD2       128
#define BAM_FSECONDARY   256
#define BAM_FQCFAIL      512
#define BAM_FDUP        1024

/* CIGAR operations, taken from samtools bam.h */
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

enum cigar_op {
    BAM_CMATCH=0,
    BAM_CINS=1,
    BAM_CDEL=2,
    BAM_CREF_SKIP=3,
    BAM_CSOFT_CLIP=4,
    BAM_CHARD_CLIP=5,
    BAM_CPAD=6,
    BAM_CBASE_MATCH=7,
    BAM_CBASE_MISMATCH=8
};

/* ----------------------------------------------------------------------
 * Function prototypes
 * We only support reading, so basically we have open, read, close along
 * with some utility functions for querying aux records.
 */

/*! Opens a SAM or BAM file.
 *
 * The mode parameter indicates the file
 * type (if not auto-detecting) and whether it is for reading or
 * writing. Use "rb" or "wb" for reading or writing BAM and "r" or
 * "w" or reading or writing SAM. When writing BAM, the mode may end
 * with a digit from 0 to 9 to indicate the compression to use with 0
 * indicating uncompressed data.
 *
 * @param fn The filename to open or create.
 * @param mode The input/output mode, similar to fopen().
 *
 * @return
 * Returns a bam_file_t pointer on success;
 *         NULL on failure.
 */
bam_file_t *bam_open(char *fn, char *mode);

/*! Closes a SAM or BAM file.
 * 
 * @param b The file to close.
 *
 * @return
 * Retrurns 0 on success;
 *         -1 on failure.
 */
int bam_close(bam_file_t *b);

/*! Deprecated: please use bam_get_seq() instead.
 */
int bam_next_seq(bam_file_t *b, bam_seq_t **bsp);

/*! Reads the next sequence.
 *
 * Fills out the next bam_seq_t struct.
 * This function will alloc and/or grow the memory accordingly, allowing for
 * efficient reuse.
 *
 * @param bsp Must be non-null, but *bsp may be NULL or an existing
 * bam_seq_t pointer.
 *
 * @return
 * Returns 1 on success;
 *         0 on eof;
 *        -1 on error.
 */
int bam_get_seq(bam_file_t *b, bam_seq_t **bsp);

/*!Looks for aux field 'key' and returns the value.
 * The type is the first char and the value is the 2nd character onwards.
 *
 * @return
 * Returns the value for key; NULL if not found.
 */
char *bam_aux_find(bam_seq_t *b, char *key);

//!Converts an encoded integer value return by bam_aux_find to an integer.
int32_t bam_aux2i(const uint8_t *dat);
float bam_aux2f(const uint8_t *s);
double bam_aux2d(const uint8_t *s);
char bam_aux2A(const uint8_t *s);
char *bam_aux2Z(const uint8_t *s);
int bam_aux_del(bam_seq_t *b, uint8_t *s);
void bam_aux_append(bam_seq_t *b, const char tag[2],
		    char type, int len, uint8_t *data);

/*! An iterator on bam_aux_t fields.
 *
 * NB: This code is not reentrant or multi-thread capable. The values
 * returned are valid until the next call to this function.
 *
 * @param key  points to an array of 2 characters (eg "RG", "NM")
 * @param type points to an address of 1 character (eg 'Z', 'i')
 * @paran val  points to an address of a bam_aux_t union.
 * @param iter_handle NULL to initialise the search, and then the
 * returned (modified) iter_handle on each subsequent call to continue
 * the iteration.
 *
 * @return
 * Returns 0 if the next value is valid, setting key, type and val;
 *        -1 when no more found.
 */
int bam_aux_iter(bam_seq_t *b, char **iter_handle,
		 char *key, char *type, bam_aux_t *val);

/* Taken from samtools/bam.h */
#define bam_seqi(s, i) ((s)[(i)/2] >> 4*(1-(i)%2) & 0xf)
#define bam_nt16_rev_table "=ACMGRSVTWYHKDBN"

/* Output code */

/*! Writes a single bam sequence object.
 *
 * @param fp The SAM/BAM file handle.
 * @param b  The bam_seq_t pointer
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int bam_put_seq(bam_file_t *fp, bam_seq_t *b);

/*! Constructs a bam_seq_t from separate components.
 *
 * Note: ignores auxiliary tags for now. These need to be appended
 * manually by the calling function.
 *
 * @param b         Points to the location of a bam_seq_t *.  If *b is NULL
 *                  or (*b)->alloc is too small, the bam_seq_t struct will
 *                  be reallocated.
 * @param extra_len Extra space to allocate for auxiliary tags
 * @param qname     Query name
 * @param qname_len Query name length
 * @param flag      BAM flags
 * @param rname     Reference ID
 * @param pos       Mapped position (N.B. 1-based)
 * @param end       Last aligned base
 * @param mapq      Mapping quality
 * @param ncigar    Number of CIGAR elements
 * @param cigar     CIGAR alignment information
 * @param mrnm      Mate reference ID
 * @param mpos      Mate position (N.B. 1-based)
 * @param isize     Insert size
 * @param len       Sequence length
 * @param seq       Sequence (ASCII format)
 * @param qual      Quality values (phred scale, 8-bit binary, no offset)
 *                  Passing in NULL to qual will cause all quality values
 *                  to be treated as absent (i.e. set to 0xff).
 *
 * @return
 * Returns number of bytes written to bam_seq_t on success (ie tag offset);
 *         -1 on error.
 */
int bam_construct_seq(bam_seq_t **b, size_t extra_len,
		      const char *qname, size_t qname_len,
		      int flag,
		      int rname,      // Ref ID
		      int pos,
		      int end, // aligned start/end coords
		      int mapq,
		      uint32_t ncigar, const uint32_t *cigar,
		      int mrnm,       // Mate Ref ID
		      int mpos,
		      int isize,
		      int len,
		      const char *seq,
		      const char *qual);

/*! Writes a SAM header block.
 *
 * @return
 * Returns 0 for success;
 *        -1 for failure
 */
int bam_write_header(bam_file_t *out);

#ifdef __cplusplus
}
#endif

#endif /* _BAM_H_ */

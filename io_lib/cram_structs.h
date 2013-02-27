#ifndef _CRAM_STRUCTS_H_
#define _CRAM_STRUCTS_H_

/*
 * Defines in-memory structs for the basic file-format objects in the
 * CRAM format.
 *
 * The basic file format is:
 *     File-def SAM-hdr Container Container ...
 *
 * Container:
 *     Service-block data-block data-block ...
 *
 * Multiple blocks in a container are grouped together as slices,
 * also sometimes referred to as landmarks in the spec.
 */


#include <stdint.h>

#include "io_lib/hash_table.h"       // From io_lib aka staden-read
#include "io_lib/bam.h"              // For BAM header parsing

#define MAX_NAME_LEN 1024

#define SEQS_PER_SLICE 10000
#define SLICE_PER_CNT  1

#define CRAM_SUBST_MATRIX "CGTNAGTNACTNACGNACGT"

#define TN_external
//#define NS_external
#define TS_external
//#define BA_external

#define MAX_STAT_VAL 1024
//#define MAX_STAT_VAL 16
typedef struct {
    int freqs[MAX_STAT_VAL];
    HashTable *h;
    int nsamp; // total number of values added
    int nvals; // total number of unique values added
} cram_stats;

/* NB: matches java impl, not the spec */
enum cram_encoding {
    E_NULL               = 0,
    E_EXTERNAL           = 1,
    E_GOLOMB             = 2,
    E_HUFFMAN            = 3,
    E_BYTE_ARRAY_LEN     = 4,
    E_BYTE_ARRAY_STOP    = 5,
    E_BETA               = 6,
    E_SUBEXP             = 7,
    E_GOLOMB_RICE        = 8,
    E_GAMMA              = 9
};

enum cram_external_type {
    E_INT                = 1,
    E_LONG               = 2,
    E_BYTE               = 3,
    E_BYTE_ARRAY         = 4,
    E_BYTE_ARRAY_BLOCK   = 5,
};

/* "File Definition Structure" */
typedef struct {
    char    magic[4];
    uint8_t major_version;
    uint8_t minor_version;
    char    file_id[20];      // Filename or SHA1 checksum
} cram_file_def;


struct cram_slice;

/* SAM header */
//typedef struct {
//    uint32_t   length;
//    char      *header;
//    HashTable *ref_hash; /* Keyed on SN, value = numeric id */
//    HashTable *rg_hash;  /* Keyed on ID, value = taglist ptr */
//    uint32_t   nref;
//    bam_ref_t *ref;
//} cram_SAM_hdr;

// Reuse a bam structure for our cram_SAM_hdr as it contains all the
// necessary functions for parsing and manipulating the header already.
// In time, merge these to be the same structure.
typedef bam_file_t cram_SAM_hdr;

enum cram_block_method {
    RAW   = 0,
    GZIP  = 1,
    BZIP2 = 2,
};

enum cram_content_type {
    FILE_HEADER        = 0,
    COMPRESSION_HEADER = 1,
    MAPPED_SLICE       = 2,
    UNMAPPED_SLICE     = 3,
    EXTERNAL           = 4,
    CORE               = 5,
};

/* Compression metrics */
typedef struct {
    int m1;
    int m2;
    int trial;
    int next_trial;
} cram_metrics;

/* Block */
typedef struct {
    enum cram_block_method  method, orig_method;
    enum cram_content_type  content_type;
    int32_t  content_id;
    int32_t  comp_size;
    int32_t  uncomp_size;
    int32_t  idx; /* offset into data */
    unsigned char    *data;

    // For bit I/O
    size_t alloc;
    size_t byte;
    int bit;
} cram_block;

struct cram_codec; /* defined in cram_codecs.h */
struct cram_map;

#define CRAM_MAP_HASH 32
#define CRAM_MAP(a,b) (((a)*3+(b))&(CRAM_MAP_HASH-1))

/* Compression header block */
typedef struct {
    int32_t ref_seq_id;
    int32_t ref_seq_start;
    int32_t ref_seq_span;
    int32_t num_records;
    int32_t num_landmarks;
    int32_t *landmark;

    /* Flags from preservation map */
    int mapped_qs_included;
    int unmapped_qs_included;
    int unmapped_placed;
    int qs_included;
    int read_names_included;
    // indexed by ref-base and subst. code
    char substitution_matrix[5][5];
    
    HashTable *preservation_map;
    struct cram_map *rec_encoding_map[CRAM_MAP_HASH];
    struct cram_map *tag_encoding_map[CRAM_MAP_HASH];

    struct cram_codec *BF_codec; // bam bit flags
    struct cram_codec *CF_codec; // compression flags
    struct cram_codec *RL_codec; // read length
    struct cram_codec *AP_codec; // alignment pos
    struct cram_codec *RG_codec; // read group
    struct cram_codec *MF_codec; // mate flags
    struct cram_codec *NS_codec; // next frag ref ID
    struct cram_codec *NP_codec; // next frag pos
    struct cram_codec *TS_codec; // template size
    struct cram_codec *NF_codec; // next frag distance
    struct cram_codec *TC_codec; // tag count
    struct cram_codec *TN_codec; // tag name/type
    struct cram_codec *FN_codec; // no. features
    struct cram_codec *FC_codec; // feature code
    struct cram_codec *FP_codec; // feature pos
    struct cram_codec *BS_codec; // base subst feature
    struct cram_codec *IN_codec; // insertion feature
    struct cram_codec *DL_codec; // deletion len feature
    struct cram_codec *BA_codec; // base feature
    struct cram_codec *MQ_codec; // mapping quality
    struct cram_codec *RN_codec; // read names
    struct cram_codec *QS_codec; // quality value (single)
    struct cram_codec *Qs_codec; // quality values (string)
    struct cram_codec *TM_codec; // ?
    struct cram_codec *TV_codec; // ?

    char *uncomp; // A single block of uncompressed data
    size_t uncomp_size, uncomp_alloc;
} cram_block_compression_hdr;

typedef struct cram_map {
    int key;    /* 0xe0 + 3 bytes */
    enum cram_encoding encoding;
    int offset; /* Offset into a single block of memory */
    int size;   /* Size */
    struct cram_codec *codec;
    struct cram_map *next; // for noddy internal hash
} cram_map;

/* Mapped or unmapped slice header block */
typedef struct {
    enum cram_content_type content_type;
    int32_t ref_seq_id;     /* if content_type == MAPPED_SLICE */
    int32_t ref_seq_start;  /* if content_type == MAPPED_SLICE */
    int32_t ref_seq_span;   /* if content_type == MAPPED_SLICE */
    int32_t num_records;
    int32_t num_blocks;
    int32_t num_content_ids;
    int32_t *block_content_ids;
    int32_t ref_base_id;    /* if content_type == MAPPED_SLICE */
} cram_block_slice_hdr;

/*
 * Container.
 *
 * Conceptually a container is split into slices, and slices into blocks.
 * However on disk it's just a list of blocks and we need to query the
 * block types to identify the start/end points of the slices.
 *
 * OR... are landmarks the start/end points of slices?
 */
typedef struct {
    int32_t  length;
    int32_t  ref_seq_id;
    int32_t  ref_seq_start;
    int32_t  ref_seq_span;
    int32_t  num_records;
    int32_t  num_blocks;
    int32_t  num_landmarks;
    int32_t *landmark;

    /* Size of container header above */
    size_t   offset;

    /* Compression header is always the first block? */
    cram_block_compression_hdr *comp_hdr;
    cram_block *comp_hdr_block;

    /* For construction purposes */
    int max_slice, curr_slice;   // maximum number of slices
    int max_rec, curr_rec;       // current and max recs per slice
    int curr_ref;                // current ref ID. -2 for no previous
    int curr_ctr_rec;            // current record number in total container
    int last_pos;                // last record position
    struct cram_slice **slices, *slice;

    /* Statistics for encoding */
    cram_stats *TS_stats;
    cram_stats *RG_stats;
    cram_stats *FP_stats;
    cram_stats *NS_stats;
    cram_stats *RN_stats;
    cram_stats *CF_stats;
    cram_stats *TN_stats;
    cram_stats *BA_stats;
    cram_stats *TV_stats;
    cram_stats *BS_stats;
    cram_stats *FC_stats;
    cram_stats *BF_stats;
    cram_stats *AP_stats;
    cram_stats *NF_stats;
    cram_stats *MF_stats;
    cram_stats *FN_stats;
    cram_stats *RL_stats;
    cram_stats *DL_stats;
    cram_stats *TC_stats;
    cram_stats *MQ_stats;
    cram_stats *TM_stats;
    cram_stats *IN_stats;
    cram_stats *QS_stats;
    cram_stats *NP_stats;
    //cram_stats *aux_B_stats;

    HashTable *tags_used; // hash of tag types in use, for tag encoding map
} cram_container;

/*
 * A single cram record
 */
typedef struct {
    struct cram_slice *s; // Filled out by cram_decode only

    int32_t ref_id;       // fixed for all recs in slice?
    int32_t flags;        // BF
    int32_t cram_flags;   // CF
    int32_t len;          // RL
    int32_t apos;         // AP
    int32_t rg;           // RG
    int32_t name;         // RN; idx to s->names_blk
    int32_t name_len;
    int32_t mate_line;    // index to another cram_record
    int32_t mate_ref_id;
    int32_t mate_pos;     // NP
    int32_t tlen;         // TS

    // Auxiliary data
    int32_t ntags;        // TC
    int32_t aux;          // idx to s->aux_blk
    int32_t aux_size;     // total size of packed ntags in aux_blk
#ifndef TN_external
    int32_t TN_idx;       // TN; idx to s->TN;
#else
    int32_t tn;           // idx to s->tn_blk
#endif

    int32_t seq;          // idx to s->seqs_blk
    int32_t qual;         // idx to s->qual_blk
    int32_t cigar;        // idx to s->cigar
    int32_t ncigar;
    int32_t aend;         // alignment end
    int32_t mqual;        // MQ

    int32_t feature;      // idx to s->feature
    int32_t nfeature;     // number of features
    int32_t mate_flags;   // MF
} cram_record;

// Accessor macros as an analogue of the bam ones
#define cram_qname(c)    (&(c)->s->name_blk->data[(c)->name])
#define cram_seq(c)      (&(c)->s->seqs_blk->data[(c)->seq])
#define cram_qual(c)     (&(c)->s->qual_blk->data[(c)->qual])
#define cram_aux(c)      (&(c)->s->aux_blk->data[(c)->aux])
#define cram_seqi(c,i)   (cram_seq((c))[(i)])
#define cram_name_len(c) ((c)->name_len)
#define cram_strand(c)   (((c)->flags & BAM_FREVERSE) != 0)
#define cram_mstrand(c)  (((c)->flags & BAM_FMREVERSE) != 0)
#define cram_cigar(c)    (&((cr)->s->cigar)[(c)->cigar])

/*
 * A feature is a base difference, used for the sequence reference encoding.
 * (We generate these internally when writing CRAM.)
 */
typedef struct {
    union {
	struct {
	    int pos;
	    int code;
	    int base;    // substitution code
	} X;
	struct {
	    int pos;
	    int code;
	    int base;    // actual base & qual
	    int qual;
	} B;
	struct {
	    int pos;
	    int code;
	    int len;
	    int seq_idx; // soft-clip multiple bases
	} S;
	struct {
	    int pos;
	    int code;
	    int len;
	    int seq_idx; // insertion multiple bases
	} I;
	struct {
	    int pos;
	    int code;
	    int base; // insertion single base
	} i;
	struct {
	    int pos;
	    int code;
	    int len;
	} D;
    };
} cram_feature;

/*
 * A slice is really just a set of blocks, but it
 * is the logical unit for decoding a number of
 * sequences.
 */
typedef struct cram_slice {
    cram_block_slice_hdr *hdr;
    cram_block *hdr_block;
    cram_block **block;
    cram_block **block_by_id;

    /* State used during encoding/decoding */
    int last_apos;

    /* Identifier used for auto-assigning read names */
    uint64_t id;

    /* Array of decoded cram records */
    cram_record *crecs;

    /* An dynamically growing buffers for data pointed
     * to by crecs[] array.
     */
    uint32_t  *cigar;
    uint32_t   cigar_alloc;
    uint32_t   ncigar;
    cram_block *name_blk;
    cram_block *seqs_blk;
    cram_block *qual_blk;
    cram_block *aux_blk;
    cram_block *base_blk; // substitutions, soft-clips

    cram_feature *features;
    int           nfeatures;
    int           afeatures; // allocated size of features

#ifndef TN_external
    // TN field (Tag Name)
    uint32_t      *TN;
    int           nTN, aTN;  // used and allocated size for TN[]
#else
    cram_block *tn_blk;
    int tn_id;
#endif

    HashTable *pair;         // for identifying read-pairs in this slice.

    char *ref;               // slice of current reference
    int ref_start;           // start position of current reference;

#ifdef BA_external
    int BA_len;
    int ba_id;
#endif
} cram_slice;

/*-----------------------------------------------------------------------------
 * Consider moving reference handling to cram_refs.[ch]
 */
// from fa.fai / samtools faidx files
typedef struct {
    char name[256];
    int64_t length;
    int64_t offset;
    int bases_per_line;
    int line_length;
} ref_entry;

// References structure.
typedef struct {
    HashTable *h_seq;
    HashTable *h_meta;
    ref_entry **ref_id;
    FILE *fp;
} refs;

/*-----------------------------------------------------------------------------
 */
/* CRAM File handle */
typedef struct {
    FILE          *fp;
    int            mode;     // 'r' or 'w'
    cram_file_def *file_def;
    cram_SAM_hdr  *SAM_hdr;
    
    char          *prefix;
    int            slice_num;
    int            err;

    // Most recent compression header decoded
    //cram_block_compression_hdr *comp_hdr;
    //cram_block_slice_hdr       *slice_hdr;

    // Current container being processed.
    cram_container *ctr;

    // positions for encoding or decoding
    int first_base, last_base;

    // cached reference portion
    refs *refs;     // ref meta-data structure
    char *ref;      // current portion held in memory
    int   ref_id;
    int   ref_start;
    int   ref_end;

    // compression level and metrics
    int level;
    cram_metrics *m[7];

    int decode_md; // Whether to export MD and NM tags
    int verbose;

    // lookup tables, stored here so we can be trivially multi-threaded
    unsigned int bam_flag_swap[0x200];  // cram -> bam flags
    unsigned int cram_flag_swap[0x800]; // bam -> cram flags
    unsigned char L1[256];              // ACGT{*} ->0123{4}
    unsigned char L2[256];              // ACGTN{*}->01234{5}
    char cram_sub_matrix[32][32];	// base substituion codes
} cram_fd;

enum cram_option {
    CRAM_OPT_DECODE_MD,
    CRAM_OPT_PREFIX,
    CRAM_OPT_VERBOSITY
};

typedef union {
    int   i;
    char *s;
} cram_opt;


/* BF bitfields */
/* FIXME: still not correct? */
#define CRAM_FPAIRED      256
#define CRAM_FPROPER_PAIR 128
#define CRAM_FUNMAP        64
#define CRAM_FREVERSE      32
#define CRAM_FREAD1        16
#define CRAM_FREAD2         8
#define CRAM_FSECONDARY     4
#define CRAM_FQCFAIL        2
#define CRAM_FDUP           1

#define CRAM_M_REVERSE  1
#define CRAM_M_UNMAP    2


/* CF bitfields */
#define CRAM_FLAG_PRESERVE_QUAL_SCORES (1<<0)
#define CRAM_FLAG_DETACHED             (1<<1)
#define CRAM_FLAG_MATE_DOWNSTREAM      (1<<2)

#endif /* _CRAM_STRUCTS_H_ */

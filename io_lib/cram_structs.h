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
#include "io_lib/dstring.h"

#define MAX_NAME_LEN 1024

#define SEQS_PER_SLICE 10000
#define SLICE_PER_CNT  1

#define CRAM_SUBST_MATRIX "CGTNAGTNACTNACGNACGT"

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
    E_INT,
    E_LONG,
    E_BYTE,
    E_BYTE_ARRAY,
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
typedef struct {
    uint32_t   length;
    char      *header;
    HashTable *ref_hash; /* Keyed on SN, value = numeric id */
    HashTable *rg_hash;  /* Keyed on ID, value = taglist ptr */
    uint32_t   nref;
    bam_ref_t *ref;
} cram_SAM_hdr;

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

/* Block */
typedef struct {
    enum cram_block_method  method, orig_method;
    enum cram_content_type  content_type;
    int32_t  content_id;
    int32_t  comp_size;
    int32_t  uncomp_size;
    int32_t  idx; /* offset into data */
    char    *data;
} cram_block;

struct cram_codec; /* defined in cram_encodings.h */

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
    char substitution_matrix[5][4];
    
    HashTable *preservation_map;
    HashTable *rec_encoding_map; /* value is cram_map */
    HashTable *tag_encoding_map; /* value is cram_map */

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

typedef struct {
    //int key;    /* 2 or 3 bytes */
    enum cram_encoding encoding;
    int offset; /* Offset into a single block of memory */
    int size;   /* Size */
    struct cram_codec *codec;
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

#define MAX_STAT_VAL 1024
//#define MAX_STAT_VAL 16
typedef struct {
    int freqs[MAX_STAT_VAL];
    HashTable *h;
    int nsamp; // total number of values added
    int nvals; // total number of unique values added
} cram_stats;

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
} cram_container;

/*
 * A single cram record
 */
typedef struct {
    int32_t ref_id;       // fixed for all recs in slice?
    int32_t flags;        // BF
    int32_t cram_flags;   // CF
    int32_t len;          // RL
    int32_t apos;         // AP
    int32_t rg;           // RG
    int32_t name;         // RN; idx to s->names_ds
    int32_t name_len;
    int32_t mate_line;    // index to another cram_record
    int32_t mate_flags;   // MF
    int32_t mate_pos;     // NP
    int32_t tlen;         // TS
    int32_t ntags;        // TC
    int32_t seq;          // idx to s->seqs_ds
    int32_t cigar;        // idx to s->cigar
    int32_t ncigar;
    int32_t aend;         // alignment end
    int32_t mqual;        // MQ
    int32_t mate_ref_id;
    int32_t qual;         // idx to s->qual_ds
    int32_t aux;          // idx to s->aux_ds
    int32_t aux_size;     // total size of packed ntags in aux_ds
    int32_t feature;      // idx to s->feature
    int32_t nfeature;     // number of features
} cram_record;

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
    dstring_t *name_ds;
    dstring_t *seqs_ds;
    dstring_t *qual_ds;
    dstring_t *aux_ds;

    dstring_t    *base_ds; // substitutions, soft-clips
    cram_feature *features;
    int           nfeatures;
    int           afeatures; // allocated size of features
} cram_slice;

typedef struct {
    HashTable *h;
    char **ref_id;
    char *file;
} refs;

/* CRAM File handle */
typedef struct {
    FILE          *fp;
    int            mode;     // 'r' or 'w'
    cram_file_def *file_def;
    cram_SAM_hdr  *SAM_hdr;
    
    char          *prefix;
    int            slice_num;
    int            err;

    refs	  *refs; // FIXME: move load ref into here and free()s too

    // Most recent compression header decoded
    //cram_block_compression_hdr *comp_hdr;
    //cram_block_slice_hdr       *slice_hdr;

    // Current container being processed.
    cram_container *ctr;

    // positions for encoding
    int first_base, last_base;
} cram_fd;


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

/* CF bitfields */
#define CRAM_FLAG_PRESERVE_QUAL_SCORES (1<<0)
#define CRAM_FLAG_DETACHED             (1<<1)
#define CRAM_FLAG_MATE_DOWNSTREAM      (1<<2)

#endif /* _CRAM_STRUCTS_H_ */

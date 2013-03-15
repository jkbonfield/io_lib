/*
 * TODO.
 *
 * - Sort order (parse to struct, enum type, updating funcs)
 * - Removal of lines.
 * - Updating of lines
 * - header_dup() ? Just convert to text, copy text, reparse.
 */

#ifndef _SAM_HEADER_H_
#define _SAM_HEADER_H_

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdarg.h>

#include "io_lib/dstring.h"
#include "io_lib/hash_table.h"
#include "io_lib/string_alloc.h"

/*
 * Proposed new SAM header parsing

1 @SQ ID:foo LN:100
2 @SQ ID:bar LN:200
3 @SQ ID:ram LN:300 UR:xyz
4 @RG ID:r ...
5 @RG ID:s ...

Hash table for 2-char keys without dup entries.
If dup lines, we form a circular linked list. Ie hash keys = {RG, SQ}.

HASH("SQ")--\
            |
    (3) <-> 1 <-> 2 <-> 3 <-> (1)

HASH("RG")--\
            |
    (5) <-> 4 <-> 5 <-> (4)

Items stored in the hash values also form their own linked lists:
Ie SQ->ID(foo)->LN(100)
   SQ->ID(bar)->LN(200)
   SQ->ID(ram)->LN(300)->UR(xyz)
   RG->ID(r)
 */

typedef struct SAM_hdr_tag_s {
    struct SAM_hdr_tag_s *next;
    char *str;
    int   len;
} SAM_hdr_tag;

typedef struct SAM_hdr_item_s {
    struct SAM_hdr_item_s *next; // cirular
    struct SAM_hdr_item_s *prev;
    SAM_hdr_tag *tag;            // first tag
    int order;                   // 0 upwards
} SAM_hdr_type;

typedef struct {
    char *name;
    uint32_t len;
    SAM_hdr_tag *tag;
} SAM_SQ;

typedef struct {
    char *name;
    SAM_hdr_tag *tag;
    int name_len;
    int id;           // numerical ID
} SAM_RG;

typedef struct {
    char *name;
    SAM_hdr_tag *tag;
    int name_len;
    int id;           // numerical ID
    int prev_id;      // -1 if none
} SAM_PG;


typedef struct {
    dstring_t *text;      // concatenated text, indexed by SAM_hdr_tag
    HashTable *h;         // 2-char IDs, values are SAM_hdr_type.
    string_alloc_t *str_pool;
    pool_alloc_t   *type_pool;
    pool_alloc_t   *tag_pool;

    // @SQ lines / references
    int nref;
    SAM_SQ *ref;
    HashTable *ref_hash;

    // @RG lines / read-groups
    int nrg;
    SAM_RG *rg;
    HashTable *rg_hash;

    // @PG lines / programs
    int npg, npg_end, npg_end_alloc;
    SAM_PG *pg;
    HashTable *pg_hash;
    int *pg_end;        // @PG chain termination IDs.

    char ID_buf[1024];  // temporary buffer
    int ID_cnt;
} SAM_hdr;


/*
 * Tokenises a SAM header into a hash table.
 * Also extracts a few bits on specific data types, such as @RG lines.
 *
 * Returns a SAM_hdr struct on success (free with sam_header_free())
 *         NULL on failure
 */
SAM_hdr *sam_header_parse(char *hdr, int len);

/*
 * Deallocates all storage used by a SAM_hdr struct.
 */
void sam_header_free(SAM_hdr *hdr);

int sam_header_length(SAM_hdr *hdr);
char *sam_header_str(SAM_hdr *hdr);

/*
 * Appends a formatted line to an existing SAM header.
 * Line is a full SAM header record, eg "@SQ\tSN:foo\tLN:100", with
 * optional new-line. If it contains more than 1 line then multiple lines
 * will be added in order.
 *
 * Len is the length of the text data, or 0 if unknown (in which case
 * it should be null terminated).
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sam_header_add_lines(SAM_hdr *sh, char *lines, int len);

/*
 * Adds a single line to a SAM header.
 * Specify type and one or more key,value pairs, ending with the NULL key.
 * Eg. sam_header_add(h, "SQ", "ID", "foo", "LN", "100", NULL).
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sam_header_add(SAM_hdr *sh, char *type, ...);
int sam_header_vadd(SAM_hdr *sh, char *type, va_list ap, ...);

/*
 * Returns the first header item matching 'type'. If ID is non-NULL it checks
 * for the tag ID: and compares against the specified ID.
 *
 * Returns NULL if no type/ID is found
 */
SAM_hdr_type *sam_header_find(SAM_hdr *hdr, char *type,
			      char *ID_key, char *ID_value);

/*
 * As per SAM_hdr_type, but returns a complete line of formatted text
 * for a specific head type/ID combination. If ID is NULL then it returns
 * the first line of the specified type.
 *
 * The returned string is malloced and should be freed by the calling
 * function with free().
 *
 * Returns NULL if no type/ID is found.
 */
char *sam_header_find_line(SAM_hdr *hdr, char *type,
			   char *ID_key, char *ID_value);

/*
 * Looks for a specific key in a single sam header line.
 * If prev is non-NULL it also fills this out with the previous tag, to
 * permit use in key removal. *prev is set to NULL when the tag is the first
 * key in the list. When a tag isn't found, prev (if non NULL) will be the last
 * tag in the existing list.
 *
 * Returns the tag pointer on success
 *         NULL on failure
 */
SAM_hdr_tag *sam_header_find_key(SAM_hdr *sh,
				 SAM_hdr_type *type,
				 char *key,
				 SAM_hdr_tag **prev);

/*
 * Adds or updates tag key,value pairs in a header line.
 * Eg for adding M5 tags to @SQ lines or updating sort order for the
 * @HD line (although use the sam_header_sort_order() function for
 * HD manipulation, which is a wrapper around this funuction).
 *
 * Specify multiple key,value pairs ending in NULL.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sam_header_update(SAM_hdr *hdr, SAM_hdr_type *type, ...);

/*
 * Reconstructs the dstring from the header hash table.
 * Returns 0 on success
 *        -1 on failure
 */
int sam_header_rebuild(SAM_hdr *hdr);

/*
 * Looks up a reference sequence by name and returns the numerical ID.
 * Returns -1 if unknown reference.
 */
int sam_header_name2ref(SAM_hdr *hdr, char *ref);

/*
 * Looks up a read-group by name and returns a pointer to the start of the
 * associated tag list.
 *
 * Returns NULL on failure
 */
SAM_RG *sam_header_find_rg(SAM_hdr *hdr, char *rg);

/*
 * Fixes any PP links in @PG headers.
 * If the entries are in order then this doesn't need doing, but incase
 * our header is out of order this goes through the sh->pg[] array
 * setting the prev_id field.
 *
 * Returns 0 on sucess
 *        -1 on failure (indicating broken PG/PP records)
 */
int sam_header_link_pg(SAM_hdr *hdr);


/*
 * Add an @PG line.
 *
 * If we wish complete control over this use sam_header_add() directly. This
 * function uses that, but attempts to do a lot of tedious house work for
 * you too.
 *
 * - It will generate a suitable ID if the supplied one clashes.
 * - It will generate multiple @PG records if we have multiple PG chains.
 *
 * Call it as per sam_header_add() with a series of key,value pairs ending
 * in NULL.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sam_header_add_PG(SAM_hdr *sh, char *name, ...);

/*
 * A function to help with construction of CL tags in @PG records.
 * Takes an argc, argv pair and returns a single space-separated string.
 * This string should be deallocated by the calling function.
 * 
 * Returns malloced char * on success
 *         NULL on failure
 */
char *stringify_argv(int argc, char *argv[]);

#endif /* _SAM_HEADER_H_ */

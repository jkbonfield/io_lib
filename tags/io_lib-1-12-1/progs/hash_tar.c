#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <fcntl.h>
#include <io_lib/tar_format.h>
#include <io_lib/hash_table.h>

typedef struct {
    char member[256];
    uint64_t pos;
    uint32_t size;
} tar_file;

void seek_forward(FILE *fp, int size) {
    if (fp != stdin) {
	fseek(fp, size, SEEK_CUR);
    } else {
	/* Seeking on a pipe isn't supported, even for fwd seeks */
	char buf[8192];
	while (size) {
	    size -= fread(buf, 1, size > 8192 ? 8192 : size, fp);
	}
    }
}

int main(int argc, char **argv) {
    int directories = 0;
    FILE *fp;
    tar_block blk;
    char member[256];
    size_t size, extra;
    int LongLink = 0;
    size_t offset = 0;
    int verbose = 0;
    HashFile *hf;
    tar_file *files = NULL;
    int nfiles = 1024;
    int fileno = 0;
    int i;
    char *header = NULL, *footer = NULL;
    int found_header, found_footer;
    int basename = 0;
    char *archive = NULL;
    int append_mode = 0;
    int prepend_mode = 0;

    files = (tar_file *)malloc(nfiles * sizeof(tar_file));

    hf = HashFileCreate(0, HASH_DYNAMIC_SIZE);

    /* process command line arguments of the form -arg */
    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (**argv != '-' || strcmp(*argv, "--") == 0)
	    break;

	if (strcmp(*argv, "-a") == 0 && argc > 1) {
	    archive = argv[1];
	    argv++;
	    argc--;
	}

	if (strcmp(*argv, "-A") == 0)
	    append_mode = 1;

	if (strcmp(*argv, "-O") == 0)
	    prepend_mode = 1;

	if (strcmp(*argv, "-d") == 0)
	    directories = 1;

	if (strcmp(*argv, "-v") == 0)
	    verbose = 1;

	if (strcmp(*argv, "-b") == 0)
	    basename = 1;

	if (strcmp(*argv, "-h") == 0 && argc > 1) {
	    /* Common header */
	    hf->headers = (HashFileSection *)
		realloc(hf->headers, (hf->nheaders+1) *
			sizeof(HashFileSection));
	    header = argv[1];
	    hf->nheaders++;
	    argv++;
	    argc--;
	}

	if (strcmp(*argv, "-f") == 0 && argc > 1) {
	    /* Common footer */
	    hf->footers = (HashFileSection *)
		realloc(hf->footers, (hf->nfooters+1) *
			sizeof(HashFileSection));
	    footer = argv[1];
	    hf->nfooters++;
	    argv++;
	    argc--;
	}
    }

    if (argc != 1 && !archive) {
	fprintf(stderr, "Usage: hash_tar [options] [tarfile] > tarfile.hash\n");
	fprintf(stderr, "    -a fname  Tar archive filename: use if reading from stdin\n");
	fprintf(stderr, "    -A        Force no archive name (eg will concat to archive itself)\n");
	fprintf(stderr, "    -O        Set arc. offset to size of hash (use when prepending)\n");
	fprintf(stderr, "    -v        Verbose mode\n");
	fprintf(stderr, "    -d        Index directory names (useless?)\n");
	fprintf(stderr, "    -h name   Set tar entry 'name' to be a file header\n");
	fprintf(stderr, "    -f name   Set tar entry 'name' to be a file footer\n");
	fprintf(stderr, "    -b        Use only the filename portion of a pathname\n");
	return 1;
    }

    /* open the tarfile */
    if (argc >= 1) {
	archive = argv[0];
	if (NULL == (fp = fopen(archive, "rb"))) {
	    perror(archive);
	    return 1;
	}
    } else {
	fp = stdin;
	if (!archive) {
	    fprintf(stderr, "If reading from stdin you must use the "
		    "\"-a archivename\" option\n");
	    return 1;
	}
    }

    /* Fill out the files[] array with the offsets, size and names */
    while(fread(&blk, sizeof(blk), 1, fp) == 1) {
	/*
	 * If a directory is too large to fit in the name (>100) but short
	 * enough to fit in the prefix the name field will be empty, this is
	 * not the cas for ordinary files where the name field is always
	 * non-empty
	 */
	if (!blk.header.name[0] && !blk.header.prefix[0])
	    break;

        /* get size of member, rounded to a multiple of TBLOCK */
	size = strtoul(blk.header.size, NULL, 8);
        extra = TBLOCK*((size+TBLOCK-1)/TBLOCK) - size;

        /* skip directories unless requested */
        if (directories || blk.header.typeflag != DIRTYPE) {

            /*
	     * extract member name (prefix + name), unless last member
	     * was ././@LongLink
	     */
            if (LongLink == 0) {
                (void) strncpy(member, blk.header.prefix, 155);
	        if (strlen(blk.header.prefix) > 0 && blk.header.name[0])
		    (void) strcat(member, "/");
    	        (void) strncat(member, blk.header.name, 100);
            }
            
            /* account for gtar ././@LongLink */
            if (strcmp(member, "././@LongLink") == 0) {
                /* still expect filenames to fit into 256 bytes */
                if (size > 256) {
                    fread(member, 1, size > 256 ? 256 : size, fp);
                    fprintf(stderr,"././@LongLink too long size=%ld\n",
			    (long)size);
                    fprintf(stderr,"%s...\n", member);
                    exit(1);
                }
                /*
		 * extract full name of next member then rewind to start
		 * of header
		 */
                fread(member, 1, size > 256 ? 256 : size, fp);
                fseek(fp, -size, SEEK_CUR);
                LongLink = 1;
            } else {
                /* output offset, member name */
                /* printf("%lu %.256s\n", (long)offset, member); */
                LongLink = 0;

		if (fileno >= nfiles) {
		    nfiles *= 2;
		    files = (tar_file *)realloc(files,nfiles*sizeof(tar_file));
		}
		if (basename) {
		    char *cp = strrchr(member, '/');
		    strcpy(files[fileno].member, cp ? cp+1 : member);
		} else {
		    strcpy(files[fileno].member, member);
		}
		files[fileno].pos = offset+sizeof(blk);
		files[fileno].size = size;
		if (verbose)
		    fprintf(stderr, "File %d: pos %010ld+%06d: %s\n",
			    fileno,
			    (long)files[fileno].pos,
			    files[fileno].size,
			    files[fileno].member);

		fileno++;
            }
        }

        /* increment offset */
        size += extra;
	seek_forward(fp, size);
        offset += sizeof(blk) + size;
    }
   
    /*
     * Find the header/footer if specified. For now we only support one of
     * each.
     */
    found_header = found_footer = 0;
    for (i = 0; i < fileno; i++) {
	if (header && strncmp(header, files[i].member, 256) == 0) {
	    hf->headers[0].pos  = files[i].pos;
	    hf->headers[0].size = files[i].size;
	    hf->headers[0].cached_data = NULL;
	    found_header++;
	}
	if (footer && strncmp(footer, files[i].member, 256) == 0) {
	    hf->footers[0].pos  = files[i].pos;
	    hf->footers[0].size = files[i].size;
	    hf->footers[0].cached_data = NULL;
	    found_footer++;
	}
    }
    if (header && !found_header) {
	fprintf(stderr, "Warning: could not find header '%s' in file\n",
		header);
	hf->nheaders = 0;
    }
    if (footer && !found_footer) {
	fprintf(stderr, "Warning: could not find footer '%s' in file\n",
		footer);
	hf->nfooters = 0;
    }

    /*
     * Construct the hash
     */
    for (i = 0; i < fileno; i++) {
	HashData hd;
	HashFileItem *hfi = (HashFileItem *)calloc(1, sizeof(*hfi));

	/* Just use the last head/foot defined as we only allow 1 at the mo. */
	hfi->header = hf->nheaders;
	hfi->footer = hf->nfooters;
	hfi->pos = files[i].pos;
	hfi->size = files[i].size;
	hd.p = hfi;
	HashTableAdd(hf->h, files[i].member, strlen(files[i].member),
		     hd, NULL);
    }

    fclose(fp);
   
    HashTableStats(hf->h, stderr);
    if (!append_mode)
	hf->archive = strdup(archive);
	
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    HashFileSave(hf, stdout, prepend_mode ? HASHFILE_PREPEND : 0);
    HashFileDestroy(hf);

    free(files);

    return 0;
}

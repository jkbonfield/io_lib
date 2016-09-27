#ifndef _BGZIP_H_
#define _BGZIP_H_

struct gzi;
typedef struct gzi gzi;

gzi *gzi_index_load(const char *fn);
void gzi_index_free(gzi *idx);
uint64_t gzi_load(FILE *fp, gzi *idx, uint64_t ustart, uint64_t uend, char *out);

struct bzi_FILE;
typedef struct bzi_FILE bzi_FILE;

bzi_FILE *bzi_open(const char *path, const char *mode);
void bzi_close(bzi_FILE *zp);
size_t bzi_read(void *ptr, size_t size, size_t nmemb, bzi_FILE *zp);
int bzi_seek(bzi_FILE *zp, off_t offset, int whence);

#endif /* _BGZIP_H_ */

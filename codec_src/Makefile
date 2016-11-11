CODECS=libcram_codec_names3.so libcram_codec_bsc.so libcram_codec_qnib1.so libcram_codec_fqzcomp_qual.so libcram_codec_zstd.so libcram_codec_rans_dedup.so

all: $(CODECS)

CC=gcc
CXX=g++
CFLAGS=-O3 -g -Wall -fPIC
CXXFLAGS=-O3 -g -Wall -fPIC
INCLUDES=-I..
# Temporary libtool hack!
LDFLAGS=-L../io_lib/.libs
LIBS=-lstaden-read

libcram_codec_names3.so: cram_codec_names3.c
	$(CC) $(CFLAGS) $(INCLUDES) $(DEFINES) -shared -o $@ $< $(LDFLAGS) $(LIBS)

libcram_codec_bsc.so: cram_codec_bsc.cpp
	$(CXX) -std=gnu++0x $(CFLAGS) $(INCLUDES) $(DEFINES) -shared -fopenmp -o $@ $< -I/nfs/users/nfs_j/jkb/ftp/compression/libbsc/ -L/nfs/users/nfs_j/jkb/ftp/compression/libbsc/ -lbsc $(LDFLAGS) $(LIBS)

# Waaaay too slow in current usage
libcram_codec_squash-brotli.so: cram_codec_squash.c
	$(CC) $(CFLAGS) $(INCLUDES) $(DEFINES) -I/nfs/sam_scratch/jkb/opt/squash/include/squash-0.7 -shared -o $@ $< $(LDFLAGS) -L/nfs/sam_scratch/jkb/opt/squash/lib/ -Wl,-rpath,/nfs/sam_scratch/jkb/opt/squash/lib/ -lsquash0.7 $(LIBS)

libcram_codec_qnib1.so: cram_codec_qnib1.c
	$(CC) $(CFLAGS) $(INCLUDES) $(DEFINES) -shared -o $@ $< $(LDFLAGS) $(LIBS)

libcram_codec_fqzcomp_qual.so: cram_codec_fqzcomp_qual.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) $(DEFINES) -shared -fopenmp -o $@ $< $(LDFLAGS) $(LIBS)

libcram_codec_zstd.so: cram_codec_zstd.c
	$(CC) $(CFLAGS) $(INCLUDES) -I/nfs/sam_scratch/jkb/opt/zstd/include $(DEFINES) -shared -o $@ $< $(LDFLAGS) -L/nfs/sam_scratch/jkb/opt/zstd/lib -lzstd -Wl,-rpath,/nfs/sam_scratch/jkb/opt/zstd/lib $(LIBS) 

libcram_codec_rans_dedup.so: cram_codec_rans_dedup.c
	$(CC) $(CFLAGS) $(INCLUDES) $(DEFINES) -shared -o $@ $< $(LDFLAGS) $(LIBS)

clean: 
	-rm *.o
	-rm $(CODECS)

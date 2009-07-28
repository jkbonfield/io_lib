#-----------------------------------------------------------------------------
# Makefile for Read library.
#

LIBS	= sp_read
PROGS	= $(LIBS)
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/../..
include $(SRCROOT)/global.mk
include ../../system.mk

CFLAGS += $(SHLIB_CFLAGS) $(CURL_INC) -I$(SRC)/..

include $(SRCROOT)/io_lib/options.mk

#-----------------------------------------------------------------------------
# Mandatory object files
OBJS=\
	Read.o \
	translate.o \
	scf_extras.o \
	find.o \
	mach-io.o \
	traceType.o \
	read_alloc.o \
	compress.o \
	open_trace_file.o \
	hash_table.o \
	jenkins_lookup3.o \
	mFILE.o\
	vlen.o\
	srf.o

#-----------------------------------------------------------------------------
# Optional objects, depending on above IOLIB_* definitions
ifdef IOLIB_SCF
OBJS	+=\
	read_scf.o \
	write_scf.o \
	misc_scf.o
endif

ifdef IOLIB_EXP
OBJS	+=\
	expFileIO.o
endif

ifdef IOLIB_PLN
OBJS	+=\
	seqIOPlain.o
endif

ifdef IOLIB_ABI
OBJS	+=\
	fpoint.o \
	seqIOABI.o
endif

ifdef IOLIB_ALF
OBJS	+=\
	seqIOALF.o
endif

ifdef IOLIB_CTF
OBJS	+=\
	ctfCompress.o \
	seqIOCTF.o
endif

ifdef IOLIB_ZTR
OBJS	+=\
	compression.o\
	ztr_translate.o\
	deflate_interlaced.o\
	ztr.o
endif

ifdef IOLIB_SFF
OBJS	+=\
	sff.o
endif

#-----------------------------------------------------------------------------
# Build rules

RLIBS += $(CURL_LIB)

$(LIBS): $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(RLIBS)

# Rule used when $(DEF_FILE) defined - currently only for windows
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)


install:
	cp $(PROGLIBS) $(INSTALLLIB)

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

include dependencies

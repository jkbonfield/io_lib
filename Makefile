#-----------------------------------------------------------------------------
# Makefile for Read library.
#

LIBS	= read
PROGS	= $(LIBS) progs
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk
INCLUDES_S += -Iinclude
CFLAGS += $(SHLIB_CFLAGS) $(shell curl-config --cflags)
#COPTDEBUG=-O -g3

%/$(O)/.dir:
	-@mkdir $(@D) 2>/dev/null
	touch $@

dirs:	io_lib/$(O)/.dir \
	progs/$(O)/.dir

%.o:	../%.c
	-@mkdir $(@D) 2>/dev/null
	$(CC) $(CFLAGS) -DHAVE_LIBCURL $(COBJFLAG)$@ -c $<
	$(MKDEFC) $(MKFLAGS) $@

include options.mk

#-----------------------------------------------------------------------------
# Mandatory object files
OBJS=\
	io_lib/$(O)/Read.o \
	io_lib/$(O)/translate.o \
	io_lib/$(O)/scf_extras.o \
	io_lib/$(O)/find.o \
	io_lib/$(O)/mach-io.o \
	io_lib/$(O)/traceType.o \
	io_lib/$(O)/read_alloc.o \
	io_lib/$(O)/compress.o \
	io_lib/$(O)/open_trace_file.o \
	io_lib/$(O)/hash_table.o \
	io_lib/$(O)/jenkins_lookup3.o \
	io_lib/$(O)/mFILE.o\
	io_lib/$(O)/vlen.o\
	io_lib/$(O)/srf.o

#-----------------------------------------------------------------------------
# Optional objects, depending on above IOLIB_* definitions
ifdef IOLIB_SCF
OBJS	+=\
	io_lib/$(O)/read_scf.o \
	io_lib/$(O)/write_scf.o \
	io_lib/$(O)/misc_scf.o
endif

ifdef IOLIB_EXP
OBJS	+=\
	io_lib/$(O)/expFileIO.o
endif

ifdef IOLIB_PLN
OBJS	+=\
	io_lib/$(O)/seqIOPlain.o
endif

ifdef IOLIB_ABI
OBJS	+=\
	io_lib/$(O)/fpoint.o \
	io_lib/$(O)/seqIOABI.o
endif

ifdef IOLIB_ALF
OBJS	+=\
	io_lib/$(O)/seqIOALF.o
endif

ifdef IOLIB_CTF
OBJS	+=\
	io_lib/$(O)/ctfCompress.o \
	io_lib/$(O)/seqIOCTF.o
endif

ifdef IOLIB_ZTR
OBJS	+=\
	io_lib/$(O)/compression.o\
	io_lib/$(O)/ztr_translate.o\
	io_lib/$(O)/deflate_interlaced.o\
	io_lib/$(O)/ztr.o
endif

ifdef IOLIB_SFF
OBJS	+=\
	io_lib/$(O)/sff.o
endif

#-----------------------------------------------------------------------------
# Build rules

#RLIBS += $(MISC_LIB) $(shell curl-config --libs)
RLIBS += $(MISC_LIB) -lcurl

$(LIBS): $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): dirs $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(RLIBS)

# Rule used when $(DEF_FILE) defined - currently only for windows
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)


# Platform specific clean rules
ifeq ($(MACHINE),windows)
CLEANCOMMAND = -rm -f */$(O)/*.o */$(O)/*.def */$(O)/*.ilk */$(O)/*.pdb
else
CLEANCOMMAND = -rm -f $(OBJS); cd progs && $(MAKE) clean
endif


# To rebuild the dependencies, type "gmake depend"
DEPEND_OBJ = $(OBJS)
include dependencies

.PHONY:	progs io_lib

io_lib:
	cd $@ && $(MAKE)

progs:
	cd $@ && $(MAKE)

#clean:
#	$(CLEANRULE)


spotless:	clean
	-rm -f $(PROGLIBS) $(L)/so_locations
	cd progs && $(MAKE) spotless

install:
	-mkdir -p $(INSTALLBIN)
	-mkdir -p $(INSTALLLIB)
	-mkdir -p $(INSTALLLIB)/$(O)
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)
	cd progs && $(MAKE) install

FORCE:
distsrc: FORCE
	echo DIRNAME=$(DIRNAME)
	-rmdir $(DIRNAME)/alpha-binaries
	-rmdir $(DIRNAME)/solaris-binaries
	-rmdir $(DIRNAME)/sgi-binaries
	-rmdir $(DIRNAME)/linux-binaries
	cpd() { \
	    mkdir $(DIRNAME)/$$1; \
	    mkdir $(DIRNAME)/$$1/alpha-binaries; \
	    mkdir $(DIRNAME)/$$1/solaris-binaries; \
	    mkdir $(DIRNAME)/$$1/sgi-binaries; \
	    mkdir $(DIRNAME)/$$1/linux-binaries; \
	    cp -R $$1/*.[ch] $(DIRNAME)/$$1; \
	}; \
	cpd progs; \
	cpd io_lib; \
	cpd docs
	-cp -R CHANGES COPYRIGHT Makefile README options.mk include man mk \
	    $(DIRNAME)
	-cp -R progs/Makefile $(DIRNAME)/progs

depend:
	-DEPEND_SRC=`echo $(DEPEND_OBJ:.o=.c) $(DEPEND_OBJ:.o=.cpp) \
	| sed 's/$(O)\///g'`; \
	touch ./dependencies.tmp; \
	makedepend -f ./dependencies.tmp -- $(CFLAGS) -- $$DEPEND_SRC 2>&-
	sort < ./dependencies.tmp | uniq | sed -e 's; /usr/[^ ]*;;g' \
	  -e "s;`echo $(SRCROOT)|sed 's/\\./\\\\./g'`;\$$(SRCROOT);g" \
	  -e '/:/s;^[^/]*;&/$$O;' | \
	  grep -v '^[^:]*:[     ]*$$' > ./dependencies
	-rm ./dependencies.tmp*

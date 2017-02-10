Io_lib:  Version 1.14.9
=======================

Io_lib is a library of file reading and writing code to provide a general
purpose SAM/BAM/CRAM, trace file (and Experiment File) reading
interface.  Programmatically {S,B,CR}AM can be manipulated using the
scram_*() API functions while DNA Chromatogram ("trace") files  can be
read using the read_reading() function.

It has been compiled and tested on a variety of unix systems, MacOS X
and MS Windows.

The directories below here contain the io_lib code. These support the
following file formats:

	SAM/BAM sequence files
	CRAM sequence files
	SCF trace files
	ABI trace files
	ALF trace files
	ZTR trace files
	SFF trace archives
	SRF trace archives
	Experiment files
	Plain text files

These link together to form a single "libstaden-read" library supporting
all the file formats via a single read_reading (or fread_reading or
mfread_reading) function call and analogous write_reading functions
too. See the file include/Read.h for the generic 'Read' structure.

See the CHANGES for a summary of older updates or git logs for the
full details.

Version 1.14.9 (9th February 2017)
--------------

Updates:

* BAM: Added CRC checking.  Bizarrely this was absent here and in most
  other BAM implementations too.  Pure BAM decode of an uncompressed
  BAM is around 9% slower and compressed BAM to compressed BAM is
  almost identical.  The most significant hit is reading uncompressed
  BAM (and doing nothing else) which is 120% slower as CRC dominates.
  Options are available to disable the CRC checking incase this is an
  issue (scramble -!).

* CRAM: Now supports bgziped fasta references.

* CRAM/SAM: Headers are now kept in the same basic type order while
  transcoding. (Eg all @PG before all @SQ, or vice versa, depending on
  input ordering.)

* CRAM: Compression level 1 is now faster but larger. (The old -1 and
  -2 were too similar.)

* CRAM: Improved compression efficiency in some files, when switching
  from sorted to unsorted data.

* CRAM: Various speedups relating to memory handling,
  multi-threaded performance and the rANS codec.

* CRAM: Block CRC checks are now only done when the block is used,
  speeding up multi-threading and tools that do not decode all blocks
  (eg flagstat).

* Scramble -g and -G options to generate and reuse bgzip indices when
  reading and writing BAM files.

* Scramble -q option to omit updating the @PG header records.

* Experimental cram_filter tool has been added, to rapidly produce
  cram subsets.

* Migrated code base to git.  Use github for primary repository.

Bug fixes:

* BAM: Fixed the bin value calculation for placed but unmapped reads.

* CRAM: Fixed file descriptor leak in refs_load_fai().

* CRAM: Fixed a crash in MD5 calculation for sequences beyond the
  reference end.

* CRAM: Bug fixes when encoding malformed @SQ records.

* CRAM: Fixed a rare renormalisation bug in rANS codec.

* Fixed tests so make -j worked.

* Removed ancient, broken and unused popen() code.


Building
========


Zlib
----

This code makes heavy use of the Deflate algorithm, assuming a Zlib
interface.  The native Zlib bundled with most systems is now rather
old and better optimised versions exist for certain platforms
(e.g. using the SSE instructions on Intel and AMD CPUs).

Therefore the --with-zlib=/path/to/zlib configure option may be used
to point to a different Zlib.  I have tested it with the vanilla zlib,
Intel's zlib and CloudFlare's Zlib.  Of the three it appears the
CloudFlare one has the quickest implementation, but mileage may vary
depending on OS and CPU.  

CloudFlare: https://github.com/cloudflare/zlib
Intel:      https://github.com/jtkukunas/zlib
Zlib-ng:    https://github.com/Dead2/zlib-ng

The Zlib-ng one needs configuring with --zlib-compat and when you
build Io_lib you will need to define -DWITH_GZFILEOP too.  It also
doesn't work well when used in conjunction with LD_PRELOAD. Therefore
I wouldn't recommend it for now.

If you are using the CloudFlare implementation, you may also want to
disable the CRC implementation in this code if your CloudFlare zlib
was built with PCLMUL support as their implementation is faster.
Otherwise the CRC here is quicker than Zlib's own version.
Building io_lib with the internal CRC code disabled is done
with ./configure --disable-own-crc (or CFLAGS=-UIOLIB_CRC).


Linux
-----

We use the GNU autoconf build mechanism.

To build:

1. ./configure

"./configure --help" will give a list of the options for GNU autoconf. For
modifying the compiler options or flags you may wish to redefine the CC or
CFLAGS variable.

Eg (in sh or bash):
   CC=cc CFLAGS=-g ./configure

2. make (or gmake)

This will build the sources.

CFLAGS may also be changed a build time using (eg):
    make 'CFLAGS=-g ...'

3. make install

The default installation location is /usr/local/bin and /usr/local/lib. These
can be changed with the --prefix option to "configure".


Windows
-------

Under Microsoft Windows we recommend the use of MSYS and MINGW as a
build environment.

These contain enough tools to build using the configure script as per
Linux. Visit http://sourceforge.net/projects/mingw/files/ and
download/install Automated MinGW Installer (eg MinGW-5.1.4.exe), MSYS
Base System (eg MSYS-1.0.11.exe) and MSYS Supplementary Tools (eg
msysDTK-1.0.1.exe).

If you wish to use Microsoft Visual Studio you may need to add the
MSVC_includes subdirectory to your C include search path.  This
adds several missing header files (eg unistd.h and sys/time.h) needed
to build this software.  We do not have a MSVC project file available.

In this case you will also need to copy io_lib/os.h.in to io_lib/os.h
and either remove the @SET_ENDIAN@ and adjacent @ lines (as these are
normally filled out for you by autoconf) or add -DNO_AUTOCONF to your
compiler options.

The code should also build cleanly under a cross-compiler.  We use

    ./configure \
            --host=x86_64-w64-mingw32 \
            --prefix=$DIST \
            --with-io_lib=$DIST \
            --with-tcl=$DIST/lib \
            --with-tk=$DIST/lib \
            --with-tklib=$DIST/lib/tklib0.5 \
            --with-zlib=$DIST \
            LDFLAGS=-L$DIST/lib

with $DIST being pre-populated with already built and installed 3rd
party dependencies, some from MSYS mentioned above.



MacOS X
-------

The configure script should work by default, but if you are attempting
to build FAT binaries to work on both i386 and ppc targets you'll need
to disable dependency tracking. Ie:

    CFLAGS="-arch i386 -arch ppc" LDFLAGS="-arch i386 -arch ppc" \
      ../configure --disable-dependency-tracking

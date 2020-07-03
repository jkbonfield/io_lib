Io_lib:  Version 1.14.13
========================

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


Version 1.14.13 (3rd July 2020)
---------------

This release has a mixture of on-going CRAM 4 work (not compatible
with previous CRAM 4) and some more general quality of life
improvements for all CRAM versions including speed-ups and better
multi-threading.

Note both CRAM 3.1 and 4.0 are still to be considered an unofficial
CRAM extensions.

Updates:

* Scramble can now filter-in or filter-out aux tags during
  transcoding.  This is done using -d and -D options.  For example:

      scramble -D OQ,BI,BD in.bam out.cram

  removes the GATK added OQ, BI and BD aux tags.
  Requested by @jhaezebrouck in issue #24.

* The Scramble -X <profile> options are now implemented using a
  CRAM_OPT_PROFILE option.  This simplifies the scramble code and
  makes it easier to call from a library.  This also fixes a number of
  bugs in the order of argument parsing.

* Improved CRAM writing speeds.

  The bam_copy function now only copies the number of used bytes
  rather than the number of allocated bytes, which can sometimes be
  substantially smaller.  As this was done in the main thread it may
  have a significant benefit when multi-threading.

* Added libdeflate support into CRAM too (in addition to the existing
  support in BAM).  This isn't a huge change to CRAM speeds except at
  high levels (-8 and -9) which are now slower, but also better
  compression ratio.  A modest 2-3% speed gain is visible are low and
  mid levels, and at -1/-2 to -4 the compression ratio is also
  improved.

* CRAM 3.1 compression level -1 is now 25% faster, but 4% larger.
  This is achieved by difference choice of compression codecs, most
  notably disabling the name tokeniser for level 1.  Use level 2 for
  something comparable to the old behaviour.

* Added an io_lib/version.h to make it easier to detect the version
  being compiled against using IOLIB_VERSION macros.
  Requested by German Tischler in issue #25.

* Refactored the cram encoding interface used by biobambam.
  Implemented by German Tischler in PR#27.

* CRAM 4 now uses E_CONST instead of a uni-value version of
  E_HUFFMAN.  Also added offset field to VARINT_SIGNED and
  VARINT_UNSIGNED which helps for data series that have values from -1
  to MAXINT.

* CRAM 4 container structure has changed so that all values are
  variable sized integers instead of fixed size.

* Further improvements with CRAM 4's use of signed values.
  - Ref_seq_id is container and slice headers are now signed.
  - RI (ref ID) data series and NS (mate ref ID) are also now signed
    as -1 is a valid value.
  - Embedded ref id is now 0 for unusued instead of -1.

* Reversed the use of CRAM 4 delta encoding for the B array.  It only
  helps at the moment for ONT signal data, so it needs more work to
  make it auto-detect when delta makes sense. (Enabling it globally
  for CRAM4 B aux tags was accidental.)

* Htscodecs submodule has gained support for big-endian platforms
  Other big-endian improvements to parts of CRAM4 too.

Bug fixes:

* Fixed CRAM MD tag generatin when using the "b" feature code
  (NB: unused by known CRAM encoders).
  Also see https://github.com/samtools/htslib/pull/1086 for more details.

* Fixed CRAM quality string when using "q" feature code (unused by
  encoders?) and in lossy-quality mode (maybe utilised in old
  Cramtools).
  Also see https://github.com/samtools/htslib/pull/1094 for more details.

* Fixed some minor memory leaks.

* "Scramble -X archive -1" enabled lzma, which should only have
  arrived at level 7 and above. (It compared integer 7 vs ASCII '1'.)

* Removed minor compilation warning in printf debugging.

* Fixed a 7 year old bug in scram_pileup which couldn't cope with
  soft-clips being followed by hard-clips.


Technology Demo: CRAM 3.1 and 4.0
=================================

The current official GA4GH CRAM version is 3.0.

For purposes of *EVALUATION ONLY* this release of io_lib includes CRAM
version 3.1, with new compression codecs (but is otherwise identical
file layout to 3.0), and 4.0 with a few additional format
modifications, such as 64-bit sizes, deduplication of read names,
orientation changes of quality strings and a revised variable sized
integer encoding.

They can be turned on using e.g. scramble -V3.1 or scramble -V4.0.
It is likely CRAM v4.0 will be official significantly later, but we
plan on v3.1 being a recognised GA4GH standard this year.

By default enabling either of these will also enable the new codecs.
Which codecs are used also depends on the profile specified (eg via
"-X small").  Some of the new codecs are considerably slower,
especially at decompression, but by default CRAM 3.1 aims to be
comparable speed to 3.0.  Note the profiles also change the
granularity of random access (1k, 10k, 25k, 100k for fast, normal,
small and archive respectively).

Here are some example file sizes and timings with different codecs and
levels on 10 million 150bp NovaSeq reads, single threaded.  Decode
timing is checked using "scram_flagstat -b".  Tests were performed
on an Intel i5-4570 processor at 3.2GHz.

|Scramble opts.      |Size(MB) |Enc(s)|Dec(s)|Codecs used                |
|--------------------|--------:|-----:|-----:|---------------------------|
|-O bam              |    531.9|  92.3|   7.5|bgzf(zlib)                 |
|-O bam -1           |    611.4|  26.4|   5.4|bgzf(libdeflate)           |
|-O bam (default)    |    539.5|  45.0|   4.9|bgzf(libdeflate)           |
|-O bam -9           |    499.5| 920.2|   4.9|bgzf(libdeflate)           |
||||||
|-V2.0 -X fast       |    317.7|  38.8|  11.8|(default, level 1)         |
|-V2.0 (default)     |    267.6|  47.0|  10.5|(default)                  |
|-V2.0 -X small      |    218.0| 124.6|  33.1|bzip2                      |
||||||
|-V3.0 -X fast       |    264.9|  31.3|  10.8|(default, level 1)         |
|-V3.0 (default)     |    223.7|  34.7|  10.3|(default)                  |
|-V3.0 -X small      |    212.3|  88.3|  18.2|bzip2                      |
|-V3.0 -X archive    |    209.4|  98.7|  18.2|bzip2                      |
||||||
|-V3.1 -X fast       |    262.4|  29.1|   9.3|rANS++                     |
|-V3.1 (default)     |    186.4|  33.7|   8.3|rANS++,tok3                |
|-V3.1 -X small      |    176.8|  74.0|  35.2|rANS++,tok3,fqz            |
|-V3.1 -X archive    |    171.9| 127.9|  34.9|rANS++,tok3,fqz,bzip2,arith|
||||||
|-V4.0 -X fast       |    251.2|  28.9|   9.6|rANS++                     |
|-V4.0 (default)     |    182.1|  32.9|   8.2|rANS++,tok3                |
|-V4.0 -X small      |    170.9|  70.9|  35.0|rANS++,tok3,fqz            |
|-V4.0 -X archive    |    166.9| 116.4|  34.2|rANS++,tok3,fqz,bzip2,arith|

We also tested on a small human aligned HiSeq run (ERR317482)
representing older Illumina data with pre-binning era quality values.
This dataset shows less impressive gains with 4.0 over 3.0 in the
default profile, but major gains in small profile once fqzcomp quality
encoding is enabled.

Note for this file, the file sizes are larger meaning less disk
caching is possible (the test machine wasn't a memory stressed
desktop).  Threading was also enabled, albeit with just 4 threads,
which further exacerbates I/O bottlenecks.  The previous test
demonstrated BAM being faster to read than CRAM, but with large files
in a more I/O stressed situation this test demonstrates the default
profile of CRAM is faster to read than BAM, due to the smaller I/O
footprint.

NB: the table below was produced with 1.14.12.

|Scramble opts.         |Size(MB) |Enc(s)|Dec(s)|Codecs used                     |
|--------------------   |--------:|-----:|-----:|--------------------------------|
|-t4 -O bam (default)   |    6526 | 115.4|  44.7|bgzf(libdeflate)                |
||||||
|-t4 -V2.0 -X fast      |    3674 |  87.4|  31.4|(default, level 1)              |
|-t4 -V2.0 (default)    |    3435 |  91.4|  30.7|(default)                       |
|-t4 -V2.0 -X small     |    3373 | 145.5|  47.8|bzip2                           |
|-t4 -V2.0 -X archive   |    3377 | 166.3|  49.7|bzip2                           |
|-t4 -V2.0 -X archive -9|    3125 |1900.6|  76.9|bzip2                           |
||||||
|-t4 -V3.0 -X fast      |    3620 |  88.3|  29.3|(default, level 1)              |
|-t4 -V3.0 (default)    |    3287 |  90.5|  29.5|(default)                       |
|-t4 -V3.0 -X small     |    3238 | 128.5|  40.3|bzip2                           |
|-t4 -V3.0 -X archive   |    3220 | 164.9|  50.0|bzip2                           |
|-t4 -V3.0 -X archive -9|    3115 |1866.6|  75.2|bzip2, lzma                     |
||||||
|-t4 -V3.1 -X fast      |    3611 |  87.9|  29.2|rANS++                          |
|-t4 -V3.1 (default)    |    3161 |  88.8|  29.7|rANS++,tok3                     |
|-t4 -V3.1 -X small     |    2249 | 192.2| 146.1|rANS++,tok3,fqz                 |
|-t4 -V3.1 -X archive   |    2157 | 235.2| 127.5|rANS++,tok3,fqz,bzip2,arith     |
|-t4 -V3.1 -X archive   |    2145 | 480.3| 128.9|rANS++,tok3,fqz,bzip2,arith,lzma|
||||||
|-t4 -V4.0 -X fast      |    3551 |  87.8|  29.5|rANS++                          |
|-t4 -V4.0 (default)    |    3148 |  88.9|  30.0|rANS++,tok3                     |
|-t4 -V4.0 -X small     |    2236 | 189.7| 142.6|rANS++,tok3,fqz                 |
|-t4 -V4.0 -X archive   |    2139 | 226.7| 127.5|rANS++,tok3,fqz,bzip2,arith     |
|-t4 -V4.0 -X archive -9|    2132 | 453.5| 128.2|rANS++,tok3,fqz,bzip2,arith,lzma|


Building
========

Prerequisites
-------------

You will need a C compiler, a Unix "make" program plus zlib, bzip2 and
lzma libraries and associated development packages (including C header
files).  The appropriate operating system package names and comands
differ per system.  On Debian Linux derived systems use the command
below (or build and install your own copies from source):

  sudo apt-get install make zlib1g-dev libbz2-dev liblzma-dev

On RedHat derived systems the package names differ:

  sudo yum install make zlib-devel bzip2-devel xz-devel


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


Libdeflate
----------

The BAM reading and writing also has optional support for the
libdeflate library (https://github.com/ebiggers/libdeflate).  This can
be used instead of an optimised zlib (see above), and generally is
slightly faster.  Build using:

    ./configure --with-libdeflate=/path


Git clone
---------

We recommend building from a release tarball, which has the configure
script already created for you.  However if you wish to build from the
latest code then use "git clone -r" to clone recursively to get the
htscodecs submodule (or follow up a normal "git clone" with
"git submodule update --init --recursive").  You will then need to
create the configure script using autoreconf in both io_lib and
htscodecs directories.  This is easiest achieved using the supply
bootstrap script.

  ./bootstrap

The autotools programs may not be on your system.  If it fails, then
install autoconf, automake and libtool packages; see above for example
OS-specific installation commands.


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
Linux. The latest msys can be downloaded here:

   http://repo.msys2.org/distrib/msys2-x86_64-latest.exe

Once installed and setup ("pacman -Syu"; close window & relaunch msys;
"pacman -Syu" again), install mingw64 compilers via "pacman -S
--needed man base-devel git mingw-w64-x86_64-toolchain".

This should then be sufficient to configure and compile.  However note
that you may need to use "./configure --disable-shared" for the test
harness to work due to deficiences in the libtool wrapper script.

If you wish to use Microsoft Visual Studio you may need to add the
MSVC_includes subdirectory to your C include search path.  This
adds several missing header files (eg unistd.h and sys/time.h) needed
to build this software.  We do not have a MSVC project file available
and have not tested the build under this environment for a number of
years.

In this case you will also need to copy io_lib/os.h.in to io_lib/os.h
and either remove the @SET_ENDIAN@ and adjacent @ lines (as these are
normally filled out for you by autoconf) or add -DNO_AUTOCONF to your
compiler options.

The code should also build cleanly under a cross-compiler.  This has
not been tested recently, but a past successful invocation was:

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


Libbsc
------

This is experimental, just to see what we can get with a high quality
compression engine in CRAM.  It's hard to build right now, especially
given it's a C++ library and our code is C.  The hacky solution now
is (linux) e.g.:

  ../configure \
    CPPFLAGS=-I$HOME/ftp/compression/libbsc \
    LDFLAGS="-L$HOME/ftp/compression/libbsc -fopenmp" \
    LIBS=-lstdc++

Enable it using scramble -J, but note this requires experimental CRAM
versions 3.1 or 4.0.

** Neither of these should be used for production data. **


MacOS X
-------

The configure script should work by default, but if you are attempting
to build FAT binaries to work on both i386 and ppc targets you'll need
to disable dependency tracking. Ie:

    CFLAGS="-arch i386 -arch ppc" LDFLAGS="-arch i386 -arch ppc" \
      ../configure --disable-dependency-tracking

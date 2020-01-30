Io_lib:  Version 1.14.12
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


Version 1.14.12 (30th January 2020)
---------------

This is primarily a change to CRAM, focusing mainly on the unofficial
CRAM 3.1 and 4.0 file formats.  Note these newer experimental formats
are INCOMPATIBLE with the 1.14.11 output!

Some changes also affect CRAM 3.0 (current) though.  Main updates are:

* Added compression profiles to scramble: fast, normal (default),
  small and archive.  Specify using scramble -X profile-name.  These
  change compression codecs permitted as well as the granularity of
  random access ("fast" profile is 1/10th the size per block than
  normal).

* NM and MD tags are now checked during encode to validate
  auto-generation during decode.  If they differ they are stored
  verbatim.

* CRAM behaves better when many small chromosomes occur in the middle
  of larger ones (as it can switch out of multi-ref mode again).

* Numerous improvements to CRAM 4.0 compression ratios.

* Some speed improvements to CRAM 3.1 and 4.0 decoding.

* Fixes to github issues/bugs #12, #14-15, #17-22.

See CHANGES for more details.


Version 1.14.11 (16th October 2018)
---------------

Updates:

* CRAM: http(s) queries now honour redirects.
  The User-Agent header is also set, which is necessary in some
  proxies.

Bug fixes:

* CRAM: fix to major range query bug introduced in 1.14.10.

* CRAM: more bug fixing on range queries when multi-threading (EOF
  detection).

* The test harness now works correctly in bourne shell, without
  using bashisms.


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
|-V2.0 -X fast       |    302.6|  33.5|  12.7|(default, level 1)         |
|-V2.0 (default)     |    257.0|  39.7|  11.5|(default)                  |
|-V2.0 -X small      |    216.3| 123.8|  32.0|bzip2                      |
||||||
|-V3.0 -X fast       |    274.0|  30.8|  11.0|(default, level 1)         |
|-V3.0 (default)     |    223.7|  36.7|  10.4|(default)                  |
|-V3.0 -X small      |    212.2|  90.3|  18.2|bzip2                      |
|-V3.0 -X archive    |    209.3| 103.5|  18.2|bzip2, lzma                |
||||||
|-V3.1 -X fast       |    275.1|  28.6|  11.3|rANS++                     |
|-V3.1 (default)     |    186.2|  36.4|   8.5|rANS++,tok3                |
|-V3.1 -X small      |    176.8|  77.9|  34.9|rANS++,tok3,fqz            |
|-V3.1 -X archive    |    172.0| 134.7|  34.0|rANS++,tok3,fqz,bzip2,arith|
||||||
|-V4.0 -X fast       |    258.4|  29.9|  11.2|rANS++                     |
|-V4.0 (default)     |    181.9|  34.3|   8.3|rANS++,tok3                |
|-V4.0 -X small      |    170.8|  74.7|  34.4|rANS++,tok3,fqz            |
|-V4.0 -X archive    |    166.8| 122.0|  33.7|rANS++,tok3,fqz,bzip2,arith|

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
|-t4 -V3.0 -X archive   |    3220 | 164.9|  50.0|bzip2, lzma                     |
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

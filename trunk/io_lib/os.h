/*
 * File: os.h
 *
 * Author: 
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: operating system specific type definitions
 *
 */

#ifndef _OS_H_
#define _OS_H_

/*
 * The autoconf setup doesn't use options.mk as the Makefile is rebuilt.
 * I haven't worked out how to do the analogous via autoconf, so for now
 * I take the easy approach and hard-code all formats to be enabled
 */
#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#define IOLIB_ABI
#define IOLIB_SCF
#define IOLIB_ALF
#define IOLIB_PLN
#define IOLIB_CTF
#define IOLIB_EXP
#define IOLIB_SFF
#endif

#include <limits.h>
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Emulation of config.h stuff when not using autoconf
 */
#ifndef HAVE_CONFIG_H
#    ifndef HAVE_ZLIB
#        define HAVE_ZLIB
#    endif
#    ifndef HAVE_SYS_WAIT_H
#        define HAVE_SYS_WAIT_H
#    endif
#endif

/*
 *-----------------------------------------------------------------------------
 * <inttypes.h> definitions, incase they're not present
 *-----------------------------------------------------------------------------
 */
#ifndef PRId64
#define __PRI64__ "l"
#define PRId64 __PRI64__ "d"
#define PRId32 "d"
#define PRId16 "d"
#define PRId8  "d"
#define PRIu64 __PRI64__ "u"
#define PRIu32 "u"
#define PRIu16 "u"
#define PRIu8  "u"
#endif


/*
 *-----------------------------------------------------------------------------
 * Operating system specifics
 *-----------------------------------------------------------------------------
 */
/*
 * SunOS 4.x
 * Even though we use the ANSI gcc, we make use the the standard SunOS 4.x
 * libraries and include files, which are non-ansi
 */
#if defined(__sun__) && !defined(__svr4__)
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#define FOPEN_MAX 64
#define FILENAME_MAX 1024

/* Missing functions, defined in strings.c */
#define NOMEMMOVE
#define NOSTRERROR
#define BUGGY_SSCANF
#endif

/* 
 * SunOS 5.x - gcc or Sun's cc 
 */ 
#if (defined(__sun__) || defined(__sun)) && (defined(__svr4__) || defined(__SVR4))
#  define IMAGEDISPLAY 
#  define NOSTRDUP 
#endif 

/*
 * Microsoft Visual C++
 * Windows
 */
#if defined(_MSC_VER)
#define popen _popen
#define pclose _pclose
typedef int mode_t;
#define ftruncate(fd,len) _chsize(fd,len)
#define sysconf(x) 512
#define NOPIPE
#define NOLOCKF
#define NOSTRCASECMP
#define NO_STRPTIME
#undef HAVE_SYS_WAIT_H
#endif


/*
 * Microsoft Windows running MinGW
 */
#if defined(__MINGW32__)
/* #define mkdir(filename,mode) mkdir((filename)) */
#define NOPIPE
#define NOLOCKF
#define NO_STRPTIME
#define sysconf(x) 512
#define ftruncate(fd,len) _chsize(fd,len)
#undef HAVE_SYS_WAIT_H
#endif

/* Generic WIN32 API issues */
#ifdef _WIN32
#  ifndef HAVE_FSEEKO
#    if __MSVCRT_VERSION__ >= 0x800
       /* if you have MSVCR80 installed then you can use these definitions: */
#      define off_t __int64
#      define fseeko _fseeki64
#      define ftello _ftelli64
#    else
       /* otherwise we're stuck with 32-bit file support */
#      define off_t long
#      define fseeko fseek
#      define ftello ftell
#    endif
#  endif /* !HAVE_FSEEKO */
#endif /* _WIN32 */

/*
 * Linux on AMD64 also needs to use va_copy()
 */
#if defined(__linux__) && defined(__amd64__)
#  ifndef NEED_VA_COPY
#    define NEED_VA_COPY
#  endif
#endif

/*
 * DEC Alpha's running Digital UNIX
 */
#if defined(__alpha)
/* Nothing as yet */
#endif

/*
 * Silicon Graphics - Irix
 */
#if defined(__sgi)
#define NOSTRDUP
#define NO_STRPTIME
#endif

/*
 * Macs (<= OS 9) - yuk!
 */
#if defined(MAC)
#define NOSTRDUP
#endif

#if defined(__APPLE__) && defined(__ppc__)
#define NO_STRPTIME
#define NOLOCKF
#endif

#if defined(__APPLE__) && defined(__i386__)
/* nothing untoward as yet */
#endif

/*
 *-----------------------------------------------------------------------------
 * Machine endianness.
 * These are overridden at the end when using autoconf, but are here for
 * when being built as part of the Staden Package.
 *-----------------------------------------------------------------------------
 */

/*
 * First guess via autoconf. This is here as endianness can only be computed
 * by autoconf for the machine we're building on. Hence cross-compilations
 * or MacOSX "fat binaries" cannot work via autoconf in this manner.
 *
 * These provide a default for when the machine specific ifdefs below do you
 * catch your system type.
 */
#ifdef HAVE_CONFIG_H
#    if defined(WORDS_BIGENDIAN)
#        ifdef SP_LITTLE_ENDIAN
#            undef SP_LITTLE_ENDIAN
#        endif
#        define SP_BIG_ENDIAN
#    else
#        ifdef SP_BIG_ENDIAN
#            undef SP_BIG_ENDIAN
#        endif
#        define SP_LITTLE_ENDIAN
#    endif
#endif

/*
 * x86 equivalents
 */
#if defined(__i386__) || defined(__i386) || defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64) || defined(__i686__) || defined(__i686)
#  if defined(SP_BIG_ENDIAN)
#    undef SP_BIG_ENDIAN
#  endif
#  define SP_LITTLE_ENDIAN
#endif

/*
 * DEC Alpha
 */
#if defined(__alpha__) || defined(__alpha)
#  if defined(SP_LITTLE_ENDIAN)
#    undef SP_LITTLE_ENDIAN
#  endif
#  define SP_BIG_ENDIAN
#endif

/*
 * SUN Sparc
 */
#if defined(__sparc__) || defined(__sparc)
#  if defined(SP_LITTLE_ENDIAN)
#    undef SP_LITTLE_ENDIAN
#  endif
#  define SP_BIG_ENDIAN
#endif

/*
 * PowerPC
 */
#if defined(__ppc__) || defined(__ppc)
#  if defined(SP_LITTLE_ENDIAN)
#    undef SP_LITTLE_ENDIAN
#  endif
#  define SP_BIG_ENDIAN
#endif

/* Some catch-alls */
#if defined(__LITTLE_ENDIAN__) || defined(__LITTLEENDIAN__)
#    define SP_LITTLE_ENDIAN
#endif

#if defined(__BIG_ENDIAN__) || defined(__BIGENDIAN__)
#    define SP_BIG_ENDIAN
#endif

#if defined(SP_BIG_ENDIAN) && defined(SP_LITTLE_ENDIAN)
#    error Both BIG and LITTLE endian defined. Fix os.h and/or Makefile
#endif

#if !defined(SP_BIG_ENDIAN) && !defined(SP_LITTLE_ENDIAN)
#    error Neither BIG nor LITTLE endian defined. Fix os.h and/or Makefile
#endif


/*
 *-----------------------------------------------------------------------------
 * Typedefs for data sizes. Note there's umpteen versions of typedefs here
 * due to old code being supported. The ones that should be used everywhere
 * are {u,}int[124].
 *
 * C9X defines its own versions of these in inttypes.h so new code should
 * ideally be using the (for example) int32_t style of types.
 *-----------------------------------------------------------------------------
 */

/*
 * One byte integers
 */ 
typedef unsigned char	int1;
typedef unsigned char	uint1;

/*
 * Two byte integers
 */
typedef signed short	int2;
typedef unsigned short	uint2;

/*
 * Four byte integers
 */
typedef signed int	int4;
typedef unsigned int	uint4;

typedef int64_t		int8;
typedef uint64_t	uint8;

/*
 * Backwards compatibility
 */
typedef signed char	int_1;
typedef unsigned char	uint_1;
typedef signed short	int_2;
typedef unsigned short	uint_2;
typedef signed int	int_4;
typedef unsigned int	uint_4;


/*
 *-----------------------------------------------------------------------------
 * The FORTRAN interface.
 *-----------------------------------------------------------------------------
 */

typedef int4 f_int;
typedef int4 f_implicit;
typedef void f_proc_ret;	/* procedure return value */

/* James Bonfield compatability mode */
typedef int4 int_f;		/* f_int */
typedef int4 int_fl;		/* f_implicit */

#define f_proc_return() return /* (f_proc_ret) 0 */

/*
 * Use when calling/defining a Fortran function from C.
 */
#ifdef VMS
#    define FORT(symbol) (symbol)
#else
#    define FORT(symbol) (_symbol)
#endif


/*
 *-----------------------------------------------------------------------------
 * Some handy definitions.
 *-----------------------------------------------------------------------------
 */

#define MAXINT8 (9223372036854775807LL)
#define MAXINT4 (INT_MAX)
#define MAXINT2 (SHRT_MAX)


/*
 *=============================================================================
 * Anything below here should not be changed.
 *=============================================================================
 */

#define False 0
#define True 1


/*
 * Our new swap runs at the same speed on Ultrix, but substantially faster
 * (300% for swap_int4, ~50% for swap_int2) on an Alpha (due to the lack of
 * decent 'char' support).
 *
 * They also have the ability to swap in situ (src == dst). Newer code now
 * relies on this so don't change back!
 */
#define iswap_int8(x) \
    (((x & 0x00000000000000ffLL) << 56) + \
     ((x & 0x000000000000ff00LL) << 40) + \
     ((x & 0x0000000000ff0000LL) << 24) + \
     ((x & 0x00000000ff000000LL) <<  8) + \
     ((x & 0x000000ff00000000LL) >>  8) + \
     ((x & 0x0000ff0000000000LL) >> 24) + \
     ((x & 0x00ff000000000000LL) >> 40) + \
     ((x & 0xff00000000000000LL) >> 56))

#define iswap_int4(x) \
    (((x & 0x000000ff) << 24) + \
     ((x & 0x0000ff00) <<  8) + \
     ((x & 0x00ff0000) >>  8) + \
     ((x & 0xff000000) >> 24))

#define iswap_int2(x) \
    (((x & 0x00ff) << 8) + \
     ((x & 0xff00) >> 8))

#define swap_int8(src, dst) ((dst) = iswap_int8(src))
#define swap_int4(src, dst) ((dst) = iswap_int4(src))
#define swap_int2(src, dst) ((dst) = iswap_int2(src))


/*
 * Linux systems may use byteswap.h to get assembly versions of byte-swap
 * on intel systems. This can be as trivial as the bswap opcode, which works
 * out at over 2-times faster than iswap_int4 above.
 */
#if 0
#if defined(__linux__)
#    include <byteswap.h>
#    undef iswap_int8
#    undef iswap_int4
#    undef iswap_int2
#    define iswap_int8 bswap_64
#    define iswap_int4 bswap_32
#    define iswap_int2 bswap_16
#endif
#endif


/*
 * Macros to specify that data read in is of a particular endianness.
 * The macros here swap to the appropriate order for the particular machine
 * running the macro and return the new answer. These may also be used when
 * writing to a file to specify that we wish to write in (eg) big endian
 * format.
 *
 * This leads to efficient code as most of the time these macros are
 * trivial.
 */
#ifdef SP_BIG_ENDIAN
#define be_int8(x) (x)
#define be_int4(x) (x)
#define be_int2(x) (x)
#define be_int1(x) (x)

#define le_int8(x) iswap_int8((x))
#define le_int4(x) iswap_int4((x))
#define le_int2(x) iswap_int2((x))
#define le_int1(x) (x)
#endif

#ifdef SP_LITTLE_ENDIAN
#define be_int8(x) iswap_int8((x))
#define be_int4(x) iswap_int4((x))
#define be_int2(x) iswap_int2((x))
#define be_int1(x) (x)

#define le_int8(x) (x)
#define le_int4(x) (x)
#define le_int2(x) (x)
#define le_int1(x) (x)
#endif

#ifdef __cplusplus
}
#endif

#endif /*_OS_H_*/

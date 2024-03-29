dnl Process this file with autoconf to produce a configure script.
AC_INIT(io_lib, 1.15.0)
IOLIB_VERSION=$PACKAGE_VERSION
IOLIB_VERSION_MAJOR=`expr "$PACKAGE_VERSION" : '\([[0-9]]*\)'`
IOLIB_VERSION_MINOR=`expr "$PACKAGE_VERSION" : '[[0-9]]*\.\([[0-9]]*\)'`
IOLIB_VERSION_PATCH=`expr "$PACKAGE_VERSION" : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`

AC_CONFIG_HEADERS([io_lib_config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([serial-tests])
AM_MAINTAINER_MODE

# Htscodec subdirectory.
#
# For simplicity, we use --disable-shared and --with-pic here so we can
# compile these directly into libstaden-read itself.

dnl AC_CONFIG_SUBDIRS([htscodecs])
AX_SUBDIRS_CONFIGURE([htscodecs],[[--disable-shared],[--with-pic]])


# For my own benefit, here's a summary of how versioning should work.
#
# Library versioning, not to be confused with the release number:
# Update VERS_CURRENT when an interface changes in an incompatible manner.
#     eg. change of struct, removing a function, changing a prototype.
#     Also increment this when we add new functions in a backwards compatible
#     manner, but update VERS_AGE too (see below).
#
# Update VERS_REVISION whenever VERS_CURRENT hasn't changed, but the source
#     code has. (Otherwise reset revision to 0).
#
# Update VERS_AGE if new interfaces have been added, causing VERS_CURRENT to
#        be incremented too. The idea is that if the previous values are
#        VERS_CURRENT=2, VERS_AGE=0 then after adding a new function we
#        get VERS_CURRENT=3, VERS_AGE=1 (ie ABI compat with vers 2 or 3).
#
# Also see:
# http://www.gnu.org/software/libtool/manual/libtool.html#Updating-version-info

# Example: We have io_lib 1.11.6 implementing ABI vers 0, rev 6, age 0.
# Our library has an soname of libstaden-read.so.0 and files as follows:
#       (SONAME      libstaden-read.so.0)
# 	libstaden-read.so -> libstaden-read.so.0.0.6
#       libstaden-read.so.0 -> libstaden-read.so.0.0.6
#       libstaden-read.so.0.0.6
#
# We then release io_lib 1.12.0 with modified structs, so it's incompatible.
# We bump to ABI vers 1, rev 0, age 0 - giving:
#       (SONAME      libstaden-read.so.1)
# 	libstaden-read.so -> libstaden-read.so.1.0.0
#       libstaden-read.so.1 -> libstaden-read.so.1.0.0
#       libstaden-read.so.1.0.0
#
# If in 1.12.1 to 1.12.5 we then perform several successive minor tweaks (bug
# fix, optimisation, etc) that do not change the ABI at all, we may
# end up with ABI vers 1, rev 5, age 0:
#       (SONAME      libstaden-read.so.1)
# 	libstaden-read.so -> libstaden-read.so.1.0.5
#       libstaden-read.so.1 -> libstaden-read.so.1.0.5
#       libstaden-read.so.1.0.5
#
# And assume in 1.13.0 we add an extra public function. Code linked
# against 1.12.x still works (aka library version 1) meaning our
# soname doesn't change, but we do update ABI version: vers 2, rev 0, age 1.
#       (SONAME      libstaden-read.so.1)
# 	libstaden-read.so -> libstaden-read.so.1.1.0
#       libstaden-read.so.1 -> libstaden-read.so.1.1.0
#       libstaden-read.so.1.1.0

VERS_CURRENT=15
VERS_REVISION=3
VERS_AGE=1
AC_SUBST(VERS_CURRENT)
AC_SUBST(VERS_REVISION)
AC_SUBST(VERS_AGE)

dnl Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LIBTOOL

# # by Marcelo Magallon <mmagallo@efis.ucr.ac.cr>
# # Turn around -rpath problem with libtool 1.0c
# # This define should be improbable enough to not conflict with anything.
# echo host=$host
# case ${host} in
#   *-linux-gnu | *-apple-darwin*)
#     AC_MSG_RESULT([Fixing libtool for -rpath problems.])
#     sed < libtool > libtool-2 \
#     's/^hardcode_libdir_flag_spec.*$'/'hardcode_libdir_flag_spec=""/'
#     mv libtool-2 libtool
#     chmod 755 libtool
#   ;;
# esac

dnl Turn on C compiler warnings.
VL_PROG_CC_WARNINGS

dnl Check for libraries
LIBCURL_CHECK_CONFIG(,,[have_curl=yes], [have_curl=no])
ZLIB_CHECK_CONFIG(,,[have_zlib=yes], [have_zlib=no])

dnl Socket support for Solaris (et al)
AC_SEARCH_LIBS([gethostbyname], [nsl])
AC_SEARCH_LIBS([socket], [socket], [],
    [AC_CHECK_LIB([socket], [socket],
                  [LIBS="-lsocket -lnsl $LIBS"], [], [-lnsl])])

AC_ARG_ENABLE([bz2],
    [AS_HELP_STRING([--disable-bz2],
                    [omit support for BZ2-compressed CRAM files])],
    [], [enable_bz2=yes])

if test "$enable_bz2" != no
then
    bz2_devel=ok
    AC_CHECK_LIB([bz2], [BZ2_bzBuffToBuffCompress], [
        LIBS="-lbz2 $LIBS"
        AC_DEFINE([HAVE_LIBBZ2],1,[Define to 1 if you have the libbz2 library.])], [bz2_devel=missing])
    if test "$bz2_devel" != ok
    then
        AC_MSG_ERROR([libbzip2 development files not found.

This is used within CRAM.  If you wish to compile io_lib without bzip2
support, rerun configure with the --disable-bz2 option.
])
    fi
fi

AC_ARG_ENABLE([lzma],
    [AS_HELP_STRING([--disable-lzma],
                    [omit support for LZMA-compressed CRAM files])],
    [], [enable_lzma=yes])

if test "$enable_lzma" != no
then
    lzma_devel=ok
    AC_CHECK_LIB([lzma], [lzma_easy_buffer_encode], [
        LIBS="-llzma $LIBS"
        AC_DEFINE([HAVE_LIBLZMA],1,[Define to 1 if you have the liblzma library.])], [lzma_devel=missing])
    if test "$lzma_devel" != ok
    then
        AC_MSG_ERROR([liblzma development files not found.

This is used within CRAM.  If you wish to compile io_lib without lzma
support, rerun configure with the --disable-lzma option.
])
    fi
fi

AC_CHECK_LIB([bsc], [bsc_compress], [
	LIBS="-lbsc $LIBS"
	AC_DEFINE([HAVE_LIBBSC],1,[Define to 1 if you have the libbsc library.])])

AX_LIBDEFLATE
AX_ZSTD

dnl Check host endian-ness
AC_C_BIGENDIAN([SET_ENDIAN="#define SP_BIG_ENDIAN"],
	       [SET_ENDIAN="#define SP_LITTLE_ENDIAN"],
	       [SET_ENDIAN=""],
	       [SET_ENDIAN=""])

AC_SUBST(SET_ENDIAN)

dnl Do we have the va_copy function?
#
# Checks whether we the ability to use va_copy().
# AC_DEFINEs HAVE_VA_COPY if you do.
#
AC_DEFUN([AX_FUNC_VA_COPY],
[AC_LANG_PUSH(C)
AC_LINK_IFELSE([AC_LANG_SOURCE(
[[#include <stdarg.h>
void va_test(va_list ap) {
    va_list ap_local;
    va_copy(ap_local, ap);
    va_end(ap_local);
}
int main(void) { return 0; }
]])],
AC_DEFINE([HAVE_VA_COPY],1,[Define to 1 if you have the va_copy() function.]),
)
AC_LANG_POP(C)])

AX_FUNC_VA_COPY

dnl Checks for header files.
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS(fcntl.h limits.h unistd.h malloc.h)
AC_CHECK_HEADER(zlib.h)
AC_CHECK_HEADER(stdio_ext.h,[SET_STDIO_EXT="#define HAVE_STDIO_EXT_H"],[SET_STDIO_EXT=])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_CHECK_SIZEOF(short,2)
AC_CHECK_SIZEOF(int,4)
AC_CHECK_SIZEOF(long,4)
dnl AC_C_CONST
AC_C_INLINE
AC_TYPE_MODE_T
dnl AC_TYPE_OFF_T
dnl AC_TYPE_SIZE_T
dnl AC_STRUCT_TM

AC_ARG_ENABLE(custom-buffering,
        AS_HELP_STRING([--enable-custom-buffering],[enable custom buffering code for CRAM [default yes]]),
        [custombuffering=${enableval}],[custombuffering=yes])

SET_CRAM_IO_CUSTOM_BUFFERING=
if test ${custombuffering} = "yes" ; then
    SET_CRAM_IO_CUSTOM_BUFFERING="#define CRAM_IO_CUSTOM_BUFFERING"
fi

AC_ARG_ENABLE(own-crc,
	AS_HELP_STRING([--disable-own-crc],[Disable use of own CRC32 function (instead of Zlib's)]))

if test "$enable_own_crc" != "no"
then
    own_crc=yes
    AC_DEFINE([IOLIB_CRC],1,[Define to 1 to use own CRC function instead of from Zlib.])
else
    own_crc=no
fi

AC_CHECK_FUNCS(fileno)
AC_CHECK_FUNCS(fstat)
AC_CHECK_FUNCS(mallopt)

dnl Large file support
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

dnl Mmap is used for CRAM reference files
AC_FUNC_MMAP

dnl Checks for library functions.
AC_SEARCH_LIBS([pthread_join], [pthread])
AC_SEARCH_LIBS(cos, m)
dnl AC_FUNC_MEMCMP
dnl AC_FUNC_STRFTIME
dnl AC_FUNC_VPRINTF
AC_CHECK_FUNCS(strdup)
dnl AC_CHECK_FUNCS(mktime strspn strstr strtol)
AC_CHECK_FUNCS(fsync)

AC_SUBST([SET_STDIO_EXT])
AC_SUBST([SET_CRAM_IO_CUSTOM_BUFFERING])
AC_SUBST([IOLIB_VERSION])
AC_SUBST([IOLIB_VERSION_MAJOR])
AC_SUBST([IOLIB_VERSION_MINOR])
AC_SUBST([IOLIB_VERSION_PATCH])
AC_OUTPUT(Makefile io_lib/Makefile progs/Makefile tests/Makefile io_lib-config io_lib/os.h io_lib/version.h)


dnl Report summary
echo
echo "--------------------------------------------------"
echo CONFIGURE SUMMARY
echo "--------------------------------------------------"
echo Bzip2 enabled: $enable_bz2
echo lzma enabled: $enable_lzma
echo custom-buffering enabled: $custombuffering
echo Using own crc: $own_crc
echo Libcurl enabled: $have_curl
echo Libdeflate enabled: $ac_libdeflate_with
echo have libbsc: $ac_cv_lib_bsc_bsc_compress
echo Have zstd: $zstd_ok
echo
echo CC: $CC
echo CFLAGS: $CFLAGS
echo LDFLAGS: $LDFLAGS
echo LIBS: $LIBS

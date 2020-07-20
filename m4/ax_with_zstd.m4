# SYNOPSIS
#
#   AX_WITH_ZSTD([ACTION-IF-TRUE], [ACTION-IF-FALSE])
#
# DESCRIPTION
#
#   This macro checks whether ZSTD library is installed and adds a
# --with-zstd=DIR option to override the search path.
# See https://github.com/facebook/zstd for the library itself.
#
#   The following output variables are amended by this macro:
#
#     CPPFLAGS     Preprocessor flags for compiling
#     LDFLAGS      Linker flags for linking against the library
#     LIBS         Library list
#
#   It also sets ZSTD_LDFLAGS variable, to aid creation of
#   pkg-config files.
#
#   The HAVE_ZSTD cpp variable will be defined in a working
#   zstd was found.
#
# LICENSE
#
#   Copyright (C) 2020 Genome Research Ltd
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.
AC_DEFUN([AX_ZSTD],
[
  AC_ARG_WITH(zstd,
	      AC_HELP_STRING([--with-zstd=DIR],[look for zstd in DIR]),
	      [ac_zstd_with=$withval],[ac_zstd_with="no"])

  # Check if it's a working library
  zstd_ok=no
  _cppflags=$CPPFLAGS
  _ldflags=$LDFLAGS
  if test "x$ac_zstd_with" != "xno"
  then
    if test "$ac_zstd_with" != "yes"
    then
      if test -f "${ac_zstd_with}/include/zstd.h"
      then
        CPPFLAGS="$CPPFLAGS -I${ac_zstd_with}/include"
      else
        CPPFLAGS="$CPPFLAGS -I${ac_zstd_with}"
      fi
      if test -f "${ac_zstd_with}/lib/libzstd.a" -o -f "${ac_zstd_with}/lib/libzstd.so"
      then
        ZSTD_LDFLAGS="-L${ac_zstd_with}/lib"
      else
        ZSTD_LDFLAGS="-L${ac_zstd_with}"
      fi
      LDFLAGS="$LDFLAGS ${ZSTD_LDFLAGS}"
    fi
    AC_SEARCH_LIBS([ZSTD_compress], [zstd],
	[AC_CHECK_HEADER(zstd.h, [zstd_ok=yes LIBS="$LIBS -lzstd"], zstd_ok=no)])
    if test "$zstd_ok" != "yes"
    then
        AC_MSG_WARN("--with-zstd specified, but non functioning")
    fi

    # perform substitutions
    if test "$zstd_ok" = "yes"
    then
        AC_DEFINE(HAVE_ZSTD, 1,
           [Define to 1 if you have a functional libz.])
	ZSTD_LDFLAGS="$ZSTD_LDFLAGS $ac_cv_search_ZSTD_compress"
    else
      AC_MSG_WARN("No functioning zstd found")
      CPPFLAGS=$_cppflags
      LDFLAGS=$_ldflags
    fi
  fi

  AH_TEMPLATE([HAVE_ZSTD], [Define if zstd is installed])
  AM_CONDITIONAL(HAVE_ZSTD, test "$zstd_ok" = "yes")

  # Execute the conditional expressions
  if test "$zstd_ok" = "yes"
  then
     # This is the IF-YES path
     ifelse([$1],,:,[$1])
  else
     # This is the IF-NO path
     ifelse([$2],,:,[$2])
  fi

  # Tidy up
  unset _cppflags
  unset _ldflags
])

AC_PREREQ([2.61])
AC_INIT([sgp4], 1.0, [info@danrw.com])
AM_INIT_AUTOMAKE
AM_SILENT_RULES([yes])

AC_CONFIG_SRCDIR([libsgp4/Tle.cpp])

CFLAGS=" $CFLAGS"
CXXFLAGS=" $CXXFLAGS"

AC_PROG_RANLIB
AC_PROG_CXX

AC_ARG_ENABLE(debug,
              AC_HELP_STRING([--enable-debug],
                             [build debug library (default=no)])
              ,
              ,
              enable_debug=no)

if test x$enable_debug = xyes; then
  AM_CXXFLAGS="-g -O0 -Wall -Werror"
else
  AM_CXXFLAGS="-O2 -Wall -Werror -fomit-frame-pointer"
fi

AC_SUBST(AM_CXXFLAGS)

AC_CONFIG_FILES([Makefile
                 libsgp4/Makefile
                 passpredict/Makefile
                 runtest/Makefile
                 sattrack/Makefile])

AC_OUTPUT

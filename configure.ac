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
  AM_CXXFLAGS="-g       -O0                      -Werror -Wextra -W -Wall -Wconversion"
else
  AM_CXXFLAGS="-DNDEBUG -O2 -fomit-frame-pointer -Wextra -Werror -W -Wall -Wconversion"
fi

AC_SEARCH_LIBS([clock_gettime],
               [rt],
               [AC_DEFINE(HAVE_CLOCK_GETTIME,
                          [1],
                          [Define if clock_gettime is available.])])

AC_SUBST(AM_CXXFLAGS)

AC_CONFIG_FILES([Makefile
                 libsgp4/Makefile
                 passpredict/Makefile
                 runtest/Makefile
                 sattrack/Makefile])

AC_OUTPUT
echo "---"
echo "Configuration summary for $PACKAGE_NAME version $VERSION"
echo ""
echo " * Installation prefix: $prefix"
echo " * System type: $build_vendor-$build_os"
echo " * Host CPU: $build_cpu"
echo "        C Compiler: ${CC}"
echo "      C++ Compiler: ${CXX}"
echo "    Compiler flags: ${AM_CXXFLAGS}"
echo "      Linker flags: ${AM_LDFLAGS}"
echo "         Libraries: ${LIBS}"
echo "     Debug enabled: $enable_debug"
echo ""
echo "---"

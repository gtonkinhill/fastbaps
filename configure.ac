AC_PREREQ(2.61)

AC_INIT(your_package_name_here, m4_esyscmd_s([awk -e '/^Version:/ {print $2}' DESCRIPTION]))
AC_COPYRIGHT(Copyright (C) 2017 your name?)


## Determine Install Location of R
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])
fi

## Setup RBin
RBIN="${R_HOME}/bin/R"
CXX=`"${RBIN}" CMD config CXX`
CPPFLAGS=`"${RBIN}" CMD config CPPFLAGS`
CXXFLAGS=`"${RBIN}" CMD config CXXFLAGS`

## Package Requires C++
AC_LANG(C++)
AC_REQUIRE_CPP

## Compiler flag check
AC_PROG_CXX

## Borrowed from BHC/imager/icd/randomForest
# Check for OpenMP
AC_OPENMP

# since some systems have broken OMP libraries
# we also check that the actual package will work
ac_pkg_openmp=no
if test -n "${OPENMP_CFLAGS}"; then
  AC_MSG_CHECKING([OpenMP detected, checking if viable for package use])
  AC_LANG_CONFTEST([AC_LANG_PROGRAM([[#include <omp.h>]], [[ return omp_get_num_threads (); ]])])
  "$RBIN" CMD SHLIB conftest.c 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp=yes
  AC_MSG_RESULT([${ac_pkg_openmp}])
fi

# if ${ac_pkg_openmp} = "yes" then we have OMP, otherwise it will be "no"
if test "${ac_pkg_openmp}" = no; then
  AC_MSG_WARN([No OpenMP support. If using GCC, upgrade to >= 4.2. If using clang, upgrade to >= 3.8.0])
  AC_MSG_ERROR([Please use a different compiler.])
fi

# Fin
AC_OUTPUT

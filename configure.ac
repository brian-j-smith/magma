AC_INIT([magma], 1.3.0)
AC_CONFIG_AUX_DIR([tools])

#################### GCC Compiler ####################

AC_PROG_CC(gcc)


#################### CUDA Toolkit ####################

AC_ARG_WITH([cuda-home],
	[AC_HELP_STRING([--with-cuda-home=DIR],
		[full path to CUDA installation @<:@/usr/local/cuda@:>@])],
	[CUDA_HOME="$withval"], [CUDA_HOME=/usr/local/cuda])
AC_MSG_RESULT([setting CUDA home directory to ${CUDA_HOME}])

NVCC=${CUDA_HOME}/bin/nvcc
AC_CHECK_FILE([${NVCC}],,
	[AC_MSG_ERROR([nvcc not found!])])

NCFLAGS="-arch=sm_13 -O2 -DADD_"
CUDA_INC=${CUDA_HOME}/include
AC_MSG_CHECKING([for CUDA headers])
AC_CHECK_FILE([${CUDA_INC}/cublas.h],
	[NCFLAGS="${NCFLAGS} -I${CUDA_INC}"],
	[AC_MSG_ERROR([CUDA headers not found!])]
)

AC_MSG_CHECKING([whether this is a 64-bit version of CUDA])
AC_CHECK_FILE([${CUDA_HOME}/lib64/libcublas.so], [check_bit=64-bit])
if test "${check_bit}" == "64-bit"; then
	LDFLAGS="-L${CUDA_HOME}/lib64"
else
	LDFLAGS="-L${CUDA_HOME}/lib"
fi

LIBS="-lcuda -lcublas"


#################### LAPACK Library ####################

AC_ARG_WITH([lapack],
	[AC_HELP_STRING([--with-lapack=LIBS],
		[LAPACK shared library @<:@-llapack@:>@])],
   [LIBS="${LIBS} $withval"], [LIBS="${LIBS} -llapack"]
)


#################### MAGMA Libraries ####################

NCFLAGS="${NCFLAGS} -I./include"

AC_ARG_WITH([magma-lib],
	[AC_HELP_STRING([--with-magma-lib=DIR],
		[full path to MAGMA shared libraries @<:@@:>@])],
	[LDFLAGS="${LDFLAGS} -L$withval"]
)

LIBS="${LIBS} -lmagmablas -lmagma"


#################### R Software ####################

: ${R_HOME=`R RHOME`}
AC_ARG_WITH([r-home],
	[AC_HELP_STRING([--with-r-home=DIR],
		[full path to R installation @<:@@:>@])],
	[R_HOME=$withval]
)
AC_MSG_RESULT([setting R home directory to ${R_HOME}])
AC_CHECK_PROG([check_R], [R], [yes],, [${R_HOME}/bin])
if test -z "${check_R}"; then
	AC_MSG_ERROR([R not found!])
fi 

AC_MSG_RESULT([setting R header files to ${R_INCLUDE_DIR}])
R_CPICFLAGS=`"${R_HOME}/bin/R" CMD config CPICFLAGS`
NCFLAGS="${NCFLAGS} -Xcompiler \"-O2 -DADD_ ${R_CPICFLAGS} -I${R_INCLUDE_DIR}\""


AC_SUBST([NVCC])
AC_SUBST([NCFLAGS])
AC_SUBST([LDFLAGS])
AC_SUBST([LIBS])

AC_CONFIG_FILES([src/Makefile])
AC_OUTPUT

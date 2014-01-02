/********************************************************************************
 * File:   magFactors.c
 * Author: Brian J Smith <brian-j-smith@uiowa.edu>
 *
 * Created on June 18, 2010, 11:36 AM
 *
 * This file is part of the magma R package.
 *
 * magma is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * magma is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with magma.  If not, see <http://www.gnu.org/licenses/>.
 ********************************************************************************/

#include "magUtils.h"
#include "magFactors.h"

#include <cublas.h>
#include <cuda.h>
#include <magma.h>
#include <magma_lapack.h>


SEXP magChol(SEXP a)
{
   SEXP gpu = GET_SLOT(a, install("gpu")),
        b = PROTECT(NEW_OBJECT(MAKE_CLASS("magma")));
   int *DIMA = INTEGER(GET_DIM(a)), N = DIMA[0], N2 = N * N, LDB = N, info;
   double *B;

   if(DIMA[1] != N) error("non-square matrix");

   b = SET_SLOT(b, install(".Data"), AS_NUMERIC(a));
   SET_SLOT(b, install("gpu"), duplicate(gpu));
   B = REAL(b);
   
   if(LOGICAL_VALUE(gpu)) {
      double *dB;

      magma_malloc((void**)&dB, N2*sizeof(double));

      magma_dsetmatrix(N, N, B, LDB, dB, LDB);
      magma_dpotrf_gpu('U', N, dB, LDB, &info);
      magma_dgetmatrix(N, N, dB, LDB, B, LDB);

      magma_free(dB);
   } else {
      double *hB;

      magma_malloc_pinned((void**)&hB, N2*sizeof(double));
      lapackf77_dlacpy(MagmaUpperStr, &N, &N, B, &LDB, hB, &LDB);
      magma_dpotrf('U', N, hB, N, &info);
      lapackf77_dlacpy(MagmaUpperStr, &N, &N, hB, &LDB, B, &LDB);

      magma_free_pinned(hB);
   }

   if(info < 0) error("illegal argument %d in 'magChol", -1 * info);
   else if(info > 0) error("leading minor of order %d is not positive definite", info);

   int i, j;
   for(j = 0; j < N; j++) {
      for(i = j + 1; i < N; i++) {
         B[i + j * N] = 0.0;
      }
   }

   UNPROTECT(1);

   return b;
}


SEXP magLU(SEXP a)
{
   SEXP gpu = GET_SLOT(a, install("gpu")),
        b = PROTECT(NEW_OBJECT(MAKE_CLASS("magmaLU")));
   int *DIMA = INTEGER(GET_DIM(a)), M = DIMA[0], N = DIMA[1], LDA = M,
       MIN_MN = M < N ? M : N, *ipiv, info;
   double *A = REAL(PROTECT(AS_NUMERIC(a)));

   b = SET_SLOT(b, install(".Data"), AS_NUMERIC(a));
   SET_SLOT(b, install("pivot"), NEW_INTEGER(MIN_MN));
   ipiv = INTEGER(GET_SLOT(b, install("pivot")));
   SET_SLOT(b, install("gpu"), duplicate(gpu));

   if(LOGICAL_VALUE(gpu)) {
      double *dA;

      magma_malloc((void**)&dA, (M*N)*sizeof(double));

      magma_dsetmatrix(M, N, A, LDA, dA, LDA);
      magma_dgetrf_gpu(M, N, dA, LDA, ipiv, &info);
      magma_dgetmatrix(M, N, dA, LDA, REAL(b), LDA);

      magma_free(dA);
   } else {
      double *hA;

      magma_malloc_pinned((void**)&hA, (M*N)*sizeof(double));

      lapackf77_dlacpy(MagmaUpperLowerStr, &M, &N, A, &LDA, hA, &LDA);
      magma_dgetrf(M, N, hA, LDA, ipiv, &info);
      lapackf77_dlacpy(MagmaUpperLowerStr, &M, &N, hA, &LDA, REAL(b), &LDA);

      magma_free_pinned(hA);
   }

   if(info < 0) error("illegal argument %d in 'magLU'", -1 * info);
   else if(info > 0) error("factor U is singular");

   UNPROTECT(2);

   return b;
}


SEXP magQR(SEXP a)
{
   SEXP gpu = GET_SLOT(a, install("gpu")),
        b = PROTECT(NEW_OBJECT(MAKE_CLASS("magmaQR")));
   int *DIMA = INTEGER(GET_DIM(a)), M = DIMA[0], N = DIMA[1],
       MIN_MN = (M < N ? M : N), NB = magma_get_dgeqrf_nb(M), *pivot, info;
   double *A, *tau;

   A = REAL(SET_VECTOR_ELT(b, 0, AS_NUMERIC(duplicate(a))));
   SET_VECTOR_ELT(b, 1, ScalarInteger(MIN_MN));
   tau = REAL(SET_VECTOR_ELT(b, 2, NEW_NUMERIC(MIN_MN)));
   pivot = INTEGER(SET_VECTOR_ELT(b, 3, NEW_INTEGER(N)));

   int i;
   for(i = 1; i <= N; i++) *pivot++ = i;

   if(LOGICAL_VALUE(gpu)) {
      int LENT = (2*MIN_MN + (N+31)/32*32)*NB;
      double *dA, *dT, *work;

      SET_SLOT(b, install("work"), NEW_NUMERIC(LENT));
      work = REAL(GET_SLOT(b, install("work")));

      magma_malloc((void**)&dA, (M*N)*sizeof(double));
      magma_malloc((void**)&dT, LENT*sizeof(double));

      magma_dsetmatrix(M, N, A, M, dA, M);
      magma_dgeqrf_gpu(M, N, dA, M, tau, dT, &info);
      magma_dgetmatrix(M, N, dA, M, A, M);
      magma_dgetvector(LENT, dT, 1, work, 1);

      magma_free(dA);
      magma_free(dT);
   } else {
      int LWORK = N * NB;
      double *hA, *hwork;

      magma_malloc_pinned((void**)&hA, (M*N)*sizeof(double));
      magma_malloc_pinned((void**)&hwork, LWORK*sizeof(double));

      lapackf77_dlacpy(MagmaUpperLowerStr, &M, &N, A, &M, hA, &M);
      magma_dgeqrf(M, N, hA, M, tau, hwork, LWORK, &info);
      lapackf77_dlacpy(MagmaUpperLowerStr, &M, &N, hA, &M, A, &M);

      magma_free_pinned(hA);
      magma_free_pinned(hwork);
   }

   if(info < 0) error("illegal argument %d in 'magQR'", -1 * info);

   UNPROTECT(1);

   return b;
}


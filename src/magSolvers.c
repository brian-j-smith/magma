/********************************************************************************
 * File:   magSolvers.c
 * Author: Brian J Smith <brian-j-smith@uiowa.edu>
 *
 * Created on June 20, 2010, 4:19 AM
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
#include "magSolvers.h"

#include <cublas.h>
#include <magmablas.h>
#include <magma.h>
#include <magma_lapack.h>

SEXP magCholSolve(SEXP a, SEXP b)
{
   SEXP gpu = magGetGPU(a, b),
        c = PROTECT(NEW_OBJECT(MAKE_CLASS("magma")));
   int *DIMA = INTEGER(GET_DIM(a)), *DIMB = INTEGER(GET_DIM(b)),
       N = DIMA[0], NRHS = DIMB[1], info;
   double *A = REAL(PROTECT(AS_NUMERIC(a)));

   if(DIMA[1] != N) error("non-square matrix");
   if(DIMB[0] != N) error("non-conformable matrices");
   
   c = SET_SLOT(c, install(".Data"), AS_NUMERIC(b));
   SET_SLOT(c, install("gpu"), duplicate(gpu));

   if(LOGICAL_VALUE(gpu)) {
      double *dA, *dB;

      magma_malloc((void**)&dA, (N*N)*sizeof(double));
      magma_malloc((void**)&dB, (N*NRHS)*sizeof(double));

      magma_dsetmatrix(N, N, A, N, dA, N);
      magma_dsetmatrix(N, NRHS, REAL(c), N, dB, N);
      magma_dpotrs_gpu('U', N, NRHS, dA, N, dB, N, &info);
      magma_dgetmatrix(N, NRHS, dB, N, REAL(c), N);

      magma_free(dA);
      magma_free(dB);
   } else {
      lapackf77_dpotrs("U", &N, &NRHS, A, &N, REAL(c), &N, &info);
   }

   if(info < 0) error("Illegal argument %d in 'magCholSolve'", -1 * info);

   UNPROTECT(2);

   return c;
}


SEXP magLUSolve(SEXP a, SEXP b)
{
   SEXP gpu = magGetGPU(a, b),
        c = PROTECT(NEW_OBJECT(MAKE_CLASS("magma")));
   int *DIMA = INTEGER(GET_DIM(a)), *DIMB = INTEGER(GET_DIM(b)),
       N = DIMA[0], NRHS = DIMB[1], LDA = N, LDB = N,
       *ipiv = INTEGER(GET_SLOT(a, install("pivot"))), info;
   double *A = REAL(PROTECT(AS_NUMERIC(a)));

   if(DIMA[1] != N) error("non-square matrix");
   if(DIMB[0] != N) error("non-conformable matrices");
   
   c = SET_SLOT(c, install(".Data"), AS_NUMERIC(b));
   SET_SLOT(c, install("gpu"), duplicate(gpu));

   if(LOGICAL_VALUE(gpu)) {
      double *dA, *dB;

      magma_malloc((void**)&dA, (N*N)*sizeof(double));
      magma_malloc((void**)&dB, (N*NRHS)*sizeof(double));

      magma_dsetmatrix(N, N, A, N, dA, LDA);
      magma_dsetmatrix(N, NRHS, REAL(c), N, dB, LDB);
      magma_dgetrs_gpu('N', N, NRHS, dA, LDA, ipiv, dB, LDB, &info);
      magma_dgetmatrix(N, NRHS, dB, LDB, REAL(c), N);

      magma_free(dA);
      magma_free(dB);
   } else {
      dgetrs_("N", &N, &NRHS, A, &N, ipiv, REAL(c), &N, &info);
   }

   if(info < 0) error("illegal argument %d in 'magQRSolve'", -1 * info);

   UNPROTECT(2);

   return c;
}


SEXP magQRSolve(SEXP a, SEXP b)
{
   SEXP qr = VECTOR_ELT(a, 0),
        gpu = magGetGPU(qr, b),
        c = PROTECT(NEW_OBJECT(MAKE_CLASS("magma")));
   int *DIMA = INTEGER(GET_DIM(qr)), *DIMB = INTEGER(GET_DIM(b)),
       M = DIMA[0], N = DIMA[1], NRHS = DIMB[1],
       NB = magma_get_dgeqrf_nb(M), LWORK = (M-N+NB)*(NRHS + 2*NB),
       info;
   double *A = REAL(qr), *B = REAL(PROTECT(AS_NUMERIC(duplicate(b)))),
          *tau = REAL(VECTOR_ELT(a, 2)), *hwork;

   if(M < N) error("indeterminate linear system");
   if(DIMB[0] != M) error("non-conformable matrices");

   c = SET_SLOT(c, install(".Data"), allocMatrix(REALSXP, N, NRHS));
   SET_SLOT(c, install("gpu"), duplicate(gpu));

   magma_malloc_pinned((void**)&hwork, LWORK*sizeof(double));

   if(LOGICAL_VALUE(gpu) && 0) {
      SEXP workS = GET_SLOT(a, install("work"));
      int LENT = LENGTH(workS);
      double *dA, *dB, *dT, *work = REAL(workS);

      magma_malloc((void**)&dA, (M*N)*sizeof(double));
      magma_malloc((void**)&dB, (M*NRHS)*sizeof(double));
      magma_malloc((void**)&dT, LENT*sizeof(double));

      magma_dsetmatrix(M, N, A, M, dA, M);
      magma_dsetmatrix(M, NRHS, B, M, dB, M);
      magma_dsetvector(LENT, work, 1, dT, 1);

      magma_dgeqrs_gpu(M, N, NRHS, dA, M, tau, dB, dT, M, hwork, LWORK, &info);

      magma_dgetmatrix(N, NRHS, dA, M, REAL(c), N);

      magma_free(dA);
      magma_free(dB);
      magma_free(dT);
   } else {
      double ALPHA = 1.0;

      dormqr_("L", "T", &M, &NRHS, &N, A, &M, tau, B, &M,
              hwork, &LWORK, &info);
      dtrsm_("L", "U", "N", "N", &M, &NRHS, &ALPHA, A, &M, B, &M);

      magCopyMatrix(N, NRHS, REAL(c), N, B, M);
   }

   if(info < 0) error("illegal argument %d in 'magQRSolve'", -1 * info);

   magma_free_pinned(hwork);
   UNPROTECT(2);

   return c;
}


SEXP magSolve(SEXP a, SEXP b)
{
   SEXP gpu = magGetGPU(a, b),
        c = PROTECT(NEW_OBJECT(MAKE_CLASS("magma")));
   int *DIMA = INTEGER(GET_DIM(a)), *DIMB = INTEGER(GET_DIM(b)),
       N = DIMA[0], NRHS = DIMB[1], ipiv[N], info;

   if(DIMA[1] != N) error("non-square matrix");
   if(DIMB[0] != N) error("non-conformable matrices");

   c = SET_SLOT(c, install(".Data"), AS_NUMERIC(b));
   SET_SLOT(c, install("gpu"), duplicate(gpu));

   if(LOGICAL_VALUE(gpu)) {
      double *A = REAL(PROTECT(AS_NUMERIC(a))), *dA, *dB;

      magma_malloc((void**)&dA, (N*N)*sizeof(double));
      magma_malloc((void**)&dB, (N*NRHS)*sizeof(double));

      magma_dsetmatrix(N, N, A, N, dA, N);
      magma_dsetmatrix(N, NRHS, REAL(c), N, dB, N);

      magma_dgetrf_gpu(N, N, dA, N, ipiv, &info);
      if(info < 0) error("illegal argument %d in 'magSolve'", -1 * info);
      else if(info > 0) error("non-singular matrix");

      magma_dgetrs_gpu('N', N, NRHS, dA, N, ipiv, dB, N, &info);
      if(info < 0) error("illegal argument %d in 'magSolve'", -1 * info);

      magma_dgetmatrix(N, NRHS, dB, N, REAL(c), N);

      magma_free(dA);
      magma_free(dB);

      UNPROTECT(1);
   } else {
      double *A = REAL(PROTECT(AS_NUMERIC(duplicate(a))));

      dgesv_(&N, &NRHS, A, &N, ipiv, REAL(c), &N, &info);
      if(info < 0) error("illegal argument %d in 'magSolve'", -1 * info);
      else if(info > 0) error("non-singular matrix");
   
      UNPROTECT(1);
   }

   UNPROTECT(1);

   return c;
}


SEXP magTriSolve(SEXP a, SEXP b, SEXP k, SEXP uprtri, SEXP transa)
{
   SEXP gpu = magGetGPU(a, b),
        c = PROTECT(NEW_OBJECT(MAKE_CLASS("magma")));
   int *DIMA = INTEGER(GET_DIM(a)), *DIMB = INTEGER(GET_DIM(b)),
       M = DIMB[0], N = DIMB[1], K = INTEGER_VALUE(k);
   char UPLO = (LOGICAL_VALUE(uprtri) ? 'U' : 'L'),
        TRANSA = (LOGICAL_VALUE(transa) ? 'T' : 'N');
   double *A = REAL(PROTECT(AS_NUMERIC(a))), *B = REAL(PROTECT(AS_NUMERIC(b))),
          *dA, *dB;

   if((K <= 0) || (K > M)) error("invalid number of equations");

   c = SET_SLOT(c, install(".Data"), allocMatrix(REALSXP, K, N));
   SET_SLOT(c, install("gpu"), duplicate(gpu));

   magma_malloc((void**)&dA, (M*M)*sizeof(double));
   magma_malloc((void**)&dB, (M*N)*sizeof(double));

   magma_dsetmatrix(M, M, A, M, dA, M);
   magma_dsetmatrix(M, N, B, M, dB, M);

   if(LOGICAL_VALUE(gpu))
      magma_dtrsm('L', UPLO, TRANSA, 'N', K, N, 1.0, dA, M, dB, M);
   else
      cublasDtrsm('L', UPLO, TRANSA, 'N', K, N, 1.0, dA, M, dB, M);

   magma_dgetmatrix(K, N, dB, M, REAL(c), K);

   magma_free(dA);
   magma_free(dB);
   UNPROTECT(3);
   
   return c;
}


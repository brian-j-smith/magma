/********************************************************************************
 * File:   magAriths.c
 * Author: Brian J Smith <brian-j-smith@uiowa.edu>
 *
 * Created on June 18, 2010, 7:37 PM
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

#include "magAriths.h"
#include "magUtils.h"

#include <cublas.h>
#include <magmablas.h>


SEXP magMultmm(SEXP a, SEXP transa, SEXP b, SEXP transb)
{
   SEXP gpu = magGetGPU(a, b),
        c = PROTECT(NEW_OBJECT(MAKE_CLASS("magma")));
   int TA = LOGICAL_VALUE(transa), TB = LOGICAL_VALUE(transb),
       *DIMA = INTEGER(GET_DIM(a)), *DIMB = INTEGER(GET_DIM(b)),
       M = DIMA[TA], N = DIMB[!TB], K = DIMA[!TA],
       LDA = DIMA[0], LDB = DIMB[0], LDC = M;
   char TRANSA = (TA ? 'T' : 'N'), TRANSB = (TB ? 'T' : 'N');
   double *A = REAL(PROTECT(AS_NUMERIC(a))), *B = REAL(PROTECT(AS_NUMERIC(b))),
          *dA, *dB, *dC;
 
   if(DIMB[TB] != K) error("non-conformable matrices");

   c = SET_SLOT(c, install(".Data"), allocMatrix(REALSXP, M, N));
   SET_SLOT(c, install("gpu"), duplicate(gpu));
   
   magma_malloc((void**)&dA, (M*K)*sizeof(double));
   magma_malloc((void**)&dB, (K*N)*sizeof(double));
   magma_malloc((void**)&dC, (M*N)*sizeof(double));

   magma_dsetmatrix(DIMA[0], DIMA[1], A, LDA, dA, LDA);
   magma_dsetmatrix(DIMB[0], DIMB[1], B, LDB, dB, LDB);

   if(LOGICAL_VALUE(gpu))
      magmablas_dgemm(TRANSA, TRANSB, M, N, K, 1.0, dA, LDA, dB, LDB, 0.0, dC, LDC);
   else
      cublasDgemm(TRANSA, TRANSB, M, N, K, 1.0, dA, LDA, dB, LDB, 0.0, dC, LDC);

   magma_dgetmatrix(M, N, dC, LDC, REAL(c), LDC);

   magma_free(dA);
   magma_free(dB);
   magma_free(dC);

   UNPROTECT(3);

   return c;
}


SEXP magMultmv(SEXP a, SEXP transa, SEXP x, SEXP right)
{
   SEXP gpu = magGetGPU(a, x),
        y = PROTECT(NEW_OBJECT(MAKE_CLASS("magma")));
   int RHS = LOGICAL_VALUE(right), TA = (LOGICAL_VALUE(transa) ^ !RHS),
       *DIMA = INTEGER(GET_DIM(a)),       
       M = DIMA[0], N = DIMA[1], LENX = LENGTH(x), LENY = DIMA[TA], LDA=M;
   char TRANSA = (TA ? 'T' : 'N');
   double *A = REAL(PROTECT(AS_NUMERIC(a))), *X = REAL(PROTECT(AS_NUMERIC(x))),
          *dA, *dX, *dY;

   if(DIMA[!TA] != LENX) error("non-conformable matrices");

   y = SET_SLOT(y, install(".Data"),
                allocMatrix(REALSXP, (RHS ? LENY : 1), (RHS ? 1 : LENY)));
   SET_SLOT(y, install("gpu"), duplicate(gpu));

   magma_malloc((void**)&dA, (M*N)*sizeof(double));
   magma_malloc((void**)&dX, LENX*sizeof(double));
   magma_malloc((void**)&dY, LENY*sizeof(double));

   magma_dsetmatrix(M, N, A, LDA, dA, LDA);
   magma_dsetvector(LENX, X, 1, dX, 1);

   if(LOGICAL_VALUE(gpu)) {
      magmablas_dgemv(TRANSA, M, N, 1.0, dA, LDA, dX, 1, 0.0, dY, 1);
   } else {
      cublasDgemv(TRANSA, M, N, 1.0, dA, LDA, dX, 1, 0.0, dY, 1);
   }

   magma_dgetvector(LENY, dY, 1, REAL(y), 1);

   magma_free(dA);
   magma_free(dX);
   magma_free(dY);

   UNPROTECT(3);

   return y;
}


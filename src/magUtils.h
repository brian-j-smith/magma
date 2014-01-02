/********************************************************************************
 * File:   magUtils.h
 * Author: Brian J Smith <brian-j-smith@uiowa.edu>
 *
 * Created on June 21, 2010, 9:19 AM
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

#ifndef _MAGUTILS_H
#define	_MAGUTILS_H

#include <Rdefines.h>

#ifdef	__cplusplus
extern "C" {
#endif

   void magCopyMatrix(int m, int n, double *a, int lda, double *b, int ldb);
   SEXP magGetGPU(SEXP a, SEXP b);

   void magLoad();
   void magUnload();

   int checkCudaError(const char * msg);
   int checkCublasError(const char * msg);


#ifdef	__cplusplus
}
#endif

#endif	/* _MAGUTILS_H */

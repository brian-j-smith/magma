/********************************************************************************
 * File:   magSolvers.h
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

#ifndef _MAGSOLVERS_H
#define	_MAGSOLVERS_H

#include <Rdefines.h>

#ifdef	__cplusplus
extern "C" {
#endif


   SEXP magCholSolve(SEXP a, SEXP b);
   SEXP magLUSolve(SEXP a, SEXP b);
   SEXP magQRSolve(SEXP a, SEXP b);
   SEXP magSolve(SEXP a, SEXP b);
   SEXP magTriSolve(SEXP a, SEXP b, SEXP k, SEXP uprtri, SEXP transa);


#ifdef	__cplusplus
}
#endif

#endif	/* _MAGSOLVERS_H */

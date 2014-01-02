/********************************************************************************
 * File:   magAriths.h
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

#ifndef _MAGARITHS_H
#define	_MAGARITHS_H

#include <Rdefines.h>

#ifdef	__cplusplus
extern "C" {
#endif


   SEXP magMultmm(SEXP a, SEXP transa, SEXP b, SEXP transb);
   SEXP magMultmv(SEXP a, SEXP transa, SEXP x, SEXP right);


#ifdef	__cplusplus
}
#endif

#endif	/* _MAGARITHS_H */


################################################################################
## File:   magSolvers.r
## Author: Brian J Smith <brian-j-smith@uiowa.edu>
##
## This file is part of the magma R package.
##
## magma is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## any later version.
##
## magma is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with magma.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

################################################################################
# Linear Solvers
################################################################################


#################### Inverse from Cholesky Decomposition ####################

setMethod("chol2inv", signature(x = "magma"),
   function(x, ...) .Call("magCholSolve", x, diag(nrow(x)))
)


#################### Solve using LU decomposition ####################

setMethod("solve", signature(a = "magmaLU", b = "matrix"),
   function(a, b) .Call("magLUSolve", a, b)
)

setMethod("solve", signature(a = "magmaLU", b = "missing"),
   function(a, b) callGeneric(a, diag(nrow(a)))
)

setMethod("solve", signature(a = "magmaLU", b = "numeric"),
   function(a, b) as.numeric(callGeneric(a, as.matrix(b)))
)


#################### Solve using QR decomposition ####################

setMethod("solve", signature(a = "magmaQR", b = "matrix"),
   function(a, b) .Call("magQRSolve", a, b)
)

setMethod("solve", signature(a = "magmaQR", b = "missing"),
   function(a, b) {
      n <- nrow(a$qr)
      if(n != ncol(a$qr))
         stop("only square matrices can be inverted in 'solve'", call.=FALSE)
      callGeneric(a, diag(n))
   }
)

setMethod("solve", signature(a = "magmaQR", b = "numeric"),
   function(a, b) as.numeric(callGeneric(a, as.matrix(b)))
)


#################### Solve a System of Equations ####################

setMethod("solve", signature(a = "magma", b = "magma"),
   function(a, b, ...) .Call("magSolve", a, b)
)

setMethod("solve", signature(a = "magma", b = "matrix"),
   function(a, b, ...) .Call("magSolve", a, b)
)

setMethod("solve", signature(a = "magma", b = "missing"),
   function(a, b, ...) .Call("magSolve", a, diag(nrow(a)))
)

setMethod("solve", signature(a = "magma", b = "numeric"),
   function(a, b, ...) as.numeric(.Call("magSolve", a, as.matrix(b)))
)

setMethod("solve", signature(a = "matrix", b = "magma"),
   function(a, b, ...) .Call("magSolve", a, b)
)


#################### Solve an Upper or Lower Triangular System ####################

setMethod("backsolve", signature(r = "magma", x = "magma"),
   function(r, x, k, upper.tri, transpose) {
      .Call("magTriSolve", r, x, k, upper.tri, transpose)
   }
)

setMethod("backsolve", signature(r = "magma", x = "matrix"),
   function(r, x, k, upper.tri, transpose)
      .Call("magTriSolve", r, x, k, upper.tri, transpose)
)

setMethod("backsolve", signature(r = "magma", x = "numeric"),
   function(r, x, k, upper.tri, transpose)
      as.vector(.Call("magTriSolve", r, as.matrix(x), k, upper.tri, transpose))
)

setMethod("backsolve", signature(r = "matrix", x = "magma"),
   function(r, x, k, upper.tri, transpose)
      .Call("magTriSolve", r, x, k, upper.tri, transpose)
)

forwardsolve <- function(l, x, k = ncol(l), upper.tri=FALSE, transpose=FALSE)
{
   backsolve(l, x, k, upper.tri, transpose)
}


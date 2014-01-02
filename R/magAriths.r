################################################################################
## File:   margAriths.r
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
# Arithmetic Operatorations
################################################################################


#################### Generic Arithmetic Operators ####################

setMethod("Arith", signature(e1 = "magma", e2 = "magma"),
   function(e1, e2)
      magma(callGeneric(as(e1, "matrix"), as(e2, "matrix")), gpu=gpu(e1))
)

setMethod("Arith", signature(e1 = "magma", e2 = "matrix"),
   function(e1, e2) magma(callGeneric(as(e1, "matrix"), e2), gpu=gpu(e1))
)

setMethod("Arith", signature(e1 = "magma", e2 = "numeric"),
   function(e1, e2) magma(callGeneric(as(e1, "matrix"), e2), gpu=gpu(e1))
)

setMethod("Arith", signature(e1 = "matrix", e2 = "magma"),
   function(e1, e2) magma(callGeneric(e1, as(e2, "matrix")), gpu=gpu(e2))
)

setMethod("Arith", signature(e1 = "numeric", e2 = "magma"),
   function(e1, e2) magma(callGeneric(e1, as(e2, "matrix")), gpu=gpu(e2))
)


#################### Multiplication Operators ####################

setMethod("%*%", signature(x = "magma", y = "magma"),
   function(x, y) .Call("magMultmm", x, FALSE, y, FALSE)
)

setMethod("%*%", signature(x = "magma", y = "matrix"),
   function(x, y) .Call("magMultmm", x, FALSE, y, FALSE)
)

setMethod("%*%", signature(x = "magma", y = "numeric"),
   function(x, y) {
      if(ncol(x) == 1) .Call("magMultmm", x, FALSE, as.matrix(y), TRUE)
      else .Call("magMultmv", x, FALSE, y, TRUE)
   }
)

setMethod("%*%", signature(x = "matrix", y = "magma"),
   function(x, y) .Call("magMultmm", x, FALSE, y, FALSE)
)

setMethod("%*%", signature(x = "numeric", y = "magma"),
   function(x, y) {
      if(nrow(y) == 1) .Call("magMultmm", as.matrix(x), FALSE, y, FALSE)
      else .Call("magMultmv", y, FALSE, x, FALSE)
   }
)


#################### Crossproduct Functions: t(A) %*% B ####################

setMethod("crossprod", signature(x = "magma", y = "magma"),
   function(x, y) .Call("magMultmm", x, TRUE, y, FALSE)
)

setMethod("crossprod", signature(x = "magma", y = "matrix"),
   function(x, y) .Call("magMultmm", x, TRUE, y, FALSE)
)

setMethod("crossprod", signature(x = "magma", y = "missing"),
   function(x, y) .Call("magMultmm", x, TRUE, x, FALSE)
)

setMethod("crossprod", signature(x = "magma", y = "numeric"),
   function(x, y) .Call("magMultmv", x, TRUE, y, TRUE)
)

setMethod("crossprod", signature(x = "matrix", y = "magma"),
   function(x, y) .Call("magMultmm", x, TRUE, y, FALSE)
)

setMethod("crossprod", signature(x = "numeric", y = "magma"),
   function(x, y) .Call("magMultmv", y, FALSE, x, FALSE)
)


#################### tCrossproduct Functions: A %*% t(B) ####################

setMethod("tcrossprod", signature(x = "magma", y = "magma"),
   function(x, y) .Call("magMultmm", x, FALSE, y, TRUE)
)

setMethod("tcrossprod", signature(x = "magma", y = "matrix"),
   function(x, y) .Call("magMultmm", x, FALSE, y, TRUE)
)

setMethod("tcrossprod", signature(x = "magma", y = "missing"),
   function(x, y) .Call("magMultmm", x, FALSE, x, TRUE)
)

setMethod("tcrossprod", signature(x = "magma", y = "numeric"),
   function(x, y) .Call("magMultmm", x, FALSE, as.matrix(y), TRUE)
)

setMethod("tcrossprod", signature(x = "matrix", y = "magma"),
   function(x, y) .Call("magMultmm", x, FALSE, y, TRUE)
)

setMethod("tcrossprod", signature(x = "numeric", y = "magma"),
   function(x, y) {
      if(ncol(y) == 1) .Call("magMultmm", as.matrix(x), FALSE, y, TRUE)
      else .Call("magMultmv", y, TRUE, x, FALSE)
   }
)


\name{magma-class}
\Rdversion{1.1}
\docType{class}
\alias{magma-class}
\alias{Arith,magma,magma-method}
\alias{Arith,magma,matrix-method}
\alias{Arith,magma,numeric-method}
\alias{Arith,matrix,magma-method}
\alias{Arith,numeric,magma-method}
\alias{[,magma,ANY,ANY-method}
\alias{[,magma,ANY,missing-method}
\alias{[,magma,missing,ANY-method}
\alias{[,magma,missing,missing-method}
\alias{\%*\%,magma,magma-method}
\alias{\%*\%,magma,matrix-method}
\alias{\%*\%,magma,numeric-method}
\alias{\%*\%,matrix,magma-method}
\alias{\%*\%,numeric,magma-method}
\alias{chol2inv,magma-method}
\alias{chol,magma-method}
\alias{crossprod,magma,magma-method}
\alias{crossprod,magma,matrix-method}
\alias{crossprod,magma,missing-method}
\alias{crossprod,magma,numeric-method}
\alias{crossprod,matrix,magma-method}
\alias{crossprod,numeric,magma-method}
\alias{qr,magma-method}
\alias{show,magma-method}
\alias{solve,magma,magma-method}
\alias{solve,magma,matrix-method}
\alias{solve,magma,missing-method}
\alias{solve,magma,numeric-method}
\alias{solve,matrix,magma-method}
\alias{tcrossprod,magma,magma-method}
\alias{tcrossprod,magma,matrix-method}
\alias{tcrossprod,magma,missing-method}
\alias{tcrossprod,magma,numeric-method}
\alias{tcrossprod,matrix,magma-method}
\alias{tcrossprod,numeric,magma-method}

\title{Matrix Class "magma"}

\description{
The "\code{magma}" class is provided by the \pkg{magma} package and provides method functions for parallel matrix and linear algebra computations on hybrid multicore CPU and GPU architectures.  It inherits from the R \pkg{base} \linkS4class{matrix} class.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("magma", data, nrow, ncol, byrow, dimnames, ...)} or to the constructor function \code{\link[magma]{magma}}.
  \describe{
    \item{\code{data}}{a numeric data vector or matrix.}
    \item{\code{nrow}}{when \code{data} is not a matrix, the desired number of rows}
    \item{\code{ncol}}{when \code{data} is not a matrix, the desired number of columns}
    \item{\code{byrow}}{logical.  If \code{FALSE} (the default) the matrix is filled by
    columns, otherwise the matrix is filled by rows.}
    \item{\code{dimnames}}{a \code{\link{dimnames}} attribute for the matrix: a
    \code{list} of two character components.  They are set if not
    \code{\link{NULL}} (as per default).}
    \item{\code{gpu}}{logical.  If \code{TRUE}, the MAGMA library GPU interface is used
    for matrix operations, when available (default).  Otherwise, the CPU interface
    is used.}
  }
}

\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"}.}
    \item{\code{gpu}:}{Object of class \code{"logical"}.}
  }
}

\section{Extends}{
Class \code{"\linkS4class{matrix}"}, from data part.
Class \code{"\linkS4class{array}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{structure}"}, by class "matrix", distance 3.
Class \code{"\linkS4class{vector}"}, by class "matrix", distance 4, with explicit coerce.
}

\section{Methods}{
  \describe{
    \item{Arith}{\code{signature(e1 = "magma", e2 = "magma")}: ... }
    \item{Arith}{\code{signature(e1 = "magma", e2 = "matrix")}: ... }
    \item{Arith}{\code{signature(e1 = "magma", e2 = "numeric")}: ... }
    \item{Arith}{\code{signature(e1 = "matrix", e2 = "magma")}: ... }
    \item{Arith}{\code{signature(e1 = "numeric", e2 = "magma")}: ... }
    \item{backsolve}{\code{signature(r = "magma", x = "magma")}: ... }
    \item{backsolve}{\code{signature(r = "magma", x = "matrix")}: ... }
    \item{backsolve}{\code{signature(r = "magma", x = "numeric")}: ... }
    \item{backsolve}{\code{signature(r = "matrix", x = "magma")}: ... }
    \item{[}{\code{signature(x = "magma", i = "ANY", j = "ANY")}: ... }
    \item{[}{\code{signature(x = "magma", i = "ANY", j = "missing")}: ... }
    \item{[}{\code{signature(x = "magma", i = "missing", j = "ANY")}: ... }
    \item{[}{\code{signature(x = "magma", i = "missing", j = "missing")}: ... }
    \item{\%*\%}{\code{signature(x = "magma", y = "magma")}: ... }
    \item{\%*\%}{\code{signature(x = "magma", y = "matrix")}: ... }
    \item{\%*\%}{\code{signature(x = "magma", y = "numeric")}: ... }
    \item{\%*\%}{\code{signature(x = "matrix", y = "magma")}: ... }
    \item{\%*\%}{\code{signature(x = "numeric", y = "magma")}: ... }
    \item{chol2inv}{\code{signature(x = "magma")}: ... }
    \item{chol}{\code{signature(x = "magma")}: ... }
    \item{crossprod}{\code{signature(x = "magma", y = "magma")}: ... }
    \item{crossprod}{\code{signature(x = "magma", y = "matrix")}: ... }
    \item{crossprod}{\code{signature(x = "magma", y = "missing")}: ... }
    \item{crossprod}{\code{signature(x = "magma", y = "numeric")}: ... }
    \item{crossprod}{\code{signature(x = "matrix", y = "magma")}: ... }
    \item{crossprod}{\code{signature(x = "numeric", y = "magma")}: ... }
    \item{gpu<-}{\code{signature(object = "magma")}: ... }
    \item{gpu}{\code{signature(object = "magma")}: ... }
    \item{lu}{\code{signature(x = "magma")}: ... }
    \item{qr}{\code{signature(x = "magma")}: ... }
    \item{show}{\code{signature(object = "magma")}: ... }
    \item{solve}{\code{signature(a = "magma", b = "magma")}: ... }
    \item{solve}{\code{signature(a = "magma", b = "matrix")}: ... }
    \item{solve}{\code{signature(a = "magma", b = "missing")}: ... }
    \item{solve}{\code{signature(a = "magma", b = "numeric")}: ... }
    \item{solve}{\code{signature(a = "matrix", b = "magma")}: ... }
    \item{tcrossprod}{\code{signature(x = "magma", y = "magma")}: ... }
    \item{tcrossprod}{\code{signature(x = "magma", y = "matrix")}: ... }
    \item{tcrossprod}{\code{signature(x = "magma", y = "missing")}: ... }
    \item{tcrossprod}{\code{signature(x = "magma", y = "numeric")}: ... }
    \item{tcrossprod}{\code{signature(x = "matrix", y = "magma")}: ... }
    \item{tcrossprod}{\code{signature(x = "numeric", y = "magma")}: ... }
  }
}

\details{
It is possible to mix R \pkg{base} scalars, matrices, and vectors, with \code{"magma"} class objects.  In general, magma matrices can be used in operations and functions that support \code{"matrix"} class objects.

The MAGMA library of C functions offers two LAPACK-style interfaces, referred to
as the \bold{GPU interface} and the \bold{CPU interface}.  The GPU interface
takes input and produces results in the GPU's memory, whereas the CPU interface
produces results in the CPU's memory.  The GPU and CPU interfaces, although
similar, are not derivatives of each other, but instead have different
communication patters.  Magma objects have a \code{gpu} slot that specifies
the interface to be used for matrix and linear algebra operations.
}

\author{
Brian J. Smith <brian-j-smith@uiowa.edu>
}

\references{
Stanimire Tomov, Rajib Nath, Hatem Ltaief, and Jack Dongarra (2010)
\emph{Dense Linear Algebra Solvers for Multicore with {GPU} Accelerators},
Proceedings of IPDPS 2010: 24th IEEE International Parallel and Distributed Processing Symposium,
Atlanta, GA, April 2010
(\url{http://www.netlib.org/netlib/utk/people/JackDongarra/PAPERS/lawn225.pdf}).
}

\seealso{
\code{\link[magma]{magma}},
\code{\link[magma]{gpu}},
\code{\link{crossprod}},
\code{\link{tcrossprod}},
\code{\link[magma]{backsolve}},
\code{\link[magma]{forwardsolve}},
\code{\link{chol}},
\code{\link{chol2inv}},
\code{\link[magma]{lu}},
\code{\link{qr}},
\code{\link{solve}}
}

\examples{
## magma matrix and base matrix and vector
mA <- magma(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), 3, 3)
B <- matrix(1, 3, 4)
x <- c(1, 2, 3)

## element extraction
mA[1,2]
mA[1,]
mA[1,,drop=FALSE]
mA[,2]
mA[,2,drop=FALSE]

## matrix algebra
x + 10 * mA \%*\% B

## crossproduct
crossprod(B, mA)
crossprod(x, mA)

## general linear solver
solve(mA)
solve(mA, x)

## triangular solver
mB <- mA
mB[lower.tri(mB)] <- 0
backsolve(mB, diag(nrow(mB)))
backsolve(mB, x)

## Cholesky decomposition and inverse
U <- chol(mA)
print(U)
chol2inv(U)

## LU decomposition and solver
LU <- lu(mA)
solve(LU)
solve(LU, x)

## QR decomposition and solver
QR <- qr(mA)
solve(QR)
solve(QR, x)
}

\keyword{classes}
\keyword{algebra}

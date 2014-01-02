\name{backsolve}
\alias{backsolve}
\alias{backsolve,magma,magma-method}
\alias{backsolve,magma,matrix-method}
\alias{backsolve,magma,numeric-method}
\alias{backsolve,matrix,magma-method}
\alias{forwardsolve}

\title{Solve an Upper or Lower Triangular System}

\description{
Solves a system of linear equations where the coefficient matrix is upper (or \sQuote{right}, \sQuote{R}) or lower (\sQuote{left}, \sQuote{L}) triangular using parallelized computational algorithms.

\code{x <- backsolve(R, b)} solves \eqn{R x = b}{R * x = b}, and \cr
\code{x <- forwardsolve(L, b)} solve \eqn{L x = b}{L * x = b}, respectively.
}

\usage{
backsolve(r, x, k = ncol(r), upper.tri=TRUE, transpose=FALSE)
forwardsolve(l, x, k = ncol(l), upper.tri=FALSE, transpose=FALSE)
}

\arguments{
  \item{r, l}{an upper (or lower) triangular matrix-like object giving the
  coefficients for the system to be solved. Values below (above) the diagonal
  are ignored.}
  \item{x}{a matrix-like object whose columns give the right-hand sides for the
  equations.}
  \item{k}{The number of columns of \code{r} and rows of \code{x} to use.}
  \item{upper.tri}{logical; if \code{TRUE} (default), the upper triangular part
  of \code{r} is used. Otherwise, the lower one.}
  \item{transpose}{logical; if \code{TRUE}, solve
  \eqn{r^\top y = x}{r' * y = x} for \code{y}, i.e. \code{t(r) \%*\% y == x}.}
}

\section{Methods}{
  \describe{
    \item{\code{signature(r = "magma", x = "magma")}}{}
    \item{\code{signature(r = "magma", x = "matrix")}}{}
    \item{\code{signature(r = "magma", x = "numeric")}}{}
    \item{\code{signature(r = "matrix", x = "magma")}}{}
  }
}

\value{
The solution of the triangular system.  The result will be a vector if \code{x} is a vector and a \linkS4class{magma} class matrix if \code{x} is a magma or matrix object.

Note that \code{forwardsolve(L, b)} is just a wrapper for \code{backsolve(L, b, upper.tri=FALSE)}.
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
\code{\linkS4class{magma}}
}

\examples{
## upper triangular magma matrix and base vector
mA <- magma(c(1, 2, 3, 0, 1, 1, 0, 0, 2), nrow=3, ncol=3, byrow=TRUE)
x <- c(8, 4, 2)

## triangular linear solver
y <- backsolve(mA, x)
mA \%*\% y
}

\keyword{algebra}
\keyword{array}


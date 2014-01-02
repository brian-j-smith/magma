\name{lu}
\alias{lu}
\alias{lu,magma-method}


\title{The LU Decomposition}

\description{
Parallelized computation of an LU factorization of an \eqn{m \times n}{m * n} general matrix \eqn{A} using partial pivoting with row interchanges.
}

\usage{
\S4method{lu}{magma}(x, \dots)
}

\arguments{
  \item{x}{a \linkS4class{magma} matrix object.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
The decomposition has the form
\deqn{A = P L U}{A = P * L * U}
where \eqn{P} is a permutation matrix, \eqn{L} is lower triangular with unit diagonal elements (lower trapezoidal if \eqn{m > n}), and \eqn{U} is upper triangular (upper trapezoidal if \eqn{m < n}).
}

\value{
A \linkS4class{magmaLU} class object storing the factors \eqn{L} and \eqn{U} from the decomposition.  The unit diagonal elements of \eqn{L} are not stored.
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
\code{\linkS4class{magma}},
\code{\linkS4class{magmaLU}}
}

\examples{
mA <- magma(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), 3, 3)
lu(mA)
}

\keyword{array}
\keyword{algebra}


\name{magmaLU-class}
\Rdversion{1.1}
\docType{class}
\alias{magmaLU-class}
\alias{show,magmaLU-method}
\alias{solve,magmaLU,matrix-method}
\alias{solve,magmaLU,missing-method}
\alias{solve,magmaLU,numeric-method}

\title{Class "magmaLU"}
\description{
Storage for an LU decomposition as computed for a magma matrix using the \code{\link[magma]{lu}} method function.
}


\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("magmaLU", data, nrow, ncol, byrow, dimnames, ...)}.  More commonly the objects are created explicitly from calls of the form \code{\link[magma]{lu}(x)} where \code{x} is an object that inherits from the \code{"magma"} class or as a side-effect of other functions applied to \code{"magma"} objects.
}

\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"}}
    \item{\code{pivot}:}{Object of class \code{"numeric"}}
    \item{\code{gpu}:}{Object of class \code{"logical"}}
  }
}

\section{Extends}{
Class \code{"\linkS4class{magma}"}, directly.
Class \code{"\linkS4class{matrix}"}, by class "magma", distance 2.
Class \code{"\linkS4class{array}"}, by class "magma", distance 3.
Class \code{"\linkS4class{structure}"}, by class "magma", distance 4.
Class \code{"\linkS4class{vector}"}, by class "magma", distance 5, with explicit coerce.
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "magmaLU")}: ... }
    \item{solve}{\code{signature(a = "magmaLU", b = "matrix")}: ... }
    \item{solve}{\code{signature(a = "magmaLU", b = "missing")}: ... }
    \item{solve}{\code{signature(a = "magmaLU", b = "numeric")}: ... }
  }
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
\code{\link[magma]{lu}}
}

\examples{
mA <- magma(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), 3, 3)
y <- c(1, 2, 3)

## magmaLU object
LU <- lu(mA)

## solution to A \%*\% x = y
x <- solve(LU, y)

## check solution
val <- mA \%*\% x
all.equal(as.numeric(val), as.numeric(y))
}

\keyword{classes}
\keyword{algebra}


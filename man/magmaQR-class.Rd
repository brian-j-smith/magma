\name{magmaQR-class}
\Rdversion{1.1}
\docType{class}
\alias{magmaQR-class}
\alias{show,magmaQR-method}
\alias{solve,magmaQR,matrix-method}
\alias{solve,magmaQR,missing-method}
\alias{solve,magmaQR,numeric-method}

\title{Class "magmaQR"}
\description{
Storage for a QR decomposition as computed for a magma matrix using the \code{\link{qr}} function.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("magmaQR", ...)}.  More commonly the objects are created explicitly from calls of the form \code{\link{qr}(x)} where \code{x} is an object that inherits from the \code{"magma"} class or as a side-effect of other functions applied to \code{"magma"} objects.
}

\section{Slots}{
  \describe{
    \item{\code{work}:}{Object of class \code{"numeric"}.}
    \item{\code{.S3Class}:}{Object of class \code{"qr"}.}
  }
}

\section{Extends}{
Class \code{"\link{qr}"}, directly.
Class \code{"\linkS4class{oldClass}"}, by class "qr", distance 2.
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "magmaQR")}: ... }
    \item{solve}{\code{signature(a = "magmaQR", b = "matrix")}: ... }
    \item{solve}{\code{signature(a = "magmaQR", b = "missing")}: ... }
    \item{solve}{\code{signature(a = "magmaQR", b = "numeric")}: ... }
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
\code{\linkS4class{magma}},
\code{\link{qr}}
}

\examples{
mA <- magma(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), 3, 3)
y <- c(1, 2, 3)

## magmaQR object
QR <- qr(mA)

## solution to A \%*\% x = y
x <- solve(QR, y)

## check solution
val <- mA \%*\% x
all.equal(as.numeric(val), as.numeric(y))
}

\keyword{classes}
\keyword{algebra}


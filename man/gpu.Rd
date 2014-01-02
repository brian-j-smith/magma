\name{gpu}
\alias{gpu}
\alias{gpu<-}
\alias{gpu,magma-method}
\alias{gpu<-,magma-method}

\title{GPU Optimization for magma Matrices}

\description{
Retrieve or set the MAGMA interface (GPU or CPU) for an \linkS4class{magma} matrix.
}

\usage{
\S4method{gpu}{magma}(object)
\S4method{gpu}{magma}(object) <- value
}

\arguments{
  \item{object}{a \linkS4class{magma} class object.}
  \item{value}{logical.  If \code{TRUE}, the MAGMA library GPU interface is used
  for matrix operations, when available.  Otherwise, the CPU interface is used.}
}

\details{
The MAGMA library of C functions offers two LAPACK-style interfaces, referred to
as the \bold{GPU interface} and the \bold{CPU interface}.  The GPU interface
takes input and produces results in the GPU's memory, whereas the CPU interface
produces results in the CPU's memory.  The GPU and CPU interfaces, although
similar, are not derivatives of each other, but instead have different
communication patters.  The \code{gpu} function returns or sets the logical
flag indicating whether the GPU interface is to be used for caluculations
involving the specified magma matrix.
}

\author{
Brian J. Smith <brian-j-smith@uiowa.edu>
}

\seealso{
\code{\linkS4class{magma}}
}

\examples{
x <- magma(1:6, nrow=2, ncol=3)
gpu(x)

gpu(x) <- !gpu(x)
gpu(x)
}

\keyword{classes}


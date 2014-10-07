\name{magma}
\alias{magma}

\title{Matrix Class "magma" Constructor}

\description{
Constructs a \linkS4class{magma} class object that inherits from \linkS4class{matrix}.
}

\usage{
magma(data = NA, nrow = 1, ncol = 1, byrow = FALSE,
      dimnames = NULL, gpu = TRUE)
}

\arguments{
  \item{data}{a numeric data vector or matrix.}
  \item{nrow}{when \code{data} is not a matrix, the desired number of rows}
  \item{ncol}{when \code{data} is not a matrix, the desired number of columns}
  \item{byrow}{logical.  If \code{FALSE} (the default) the matrix is filled by
  columns, otherwise the matrix is filled by rows.}
  \item{dimnames}{a \code{\link{dimnames}} attribute for the matrix: a
  \code{list} of two character components.  They are set if not
  \code{\link{NULL}} (as per default).}
  \item{gpu}{logical.  If \code{TRUE}, the MAGMA library GPU interface is used
  for matrix operations, when available (default).  Otherwise, the CPU interface
  is used.}
}

\value{
Returns a \linkS4class{magma} class matrix.  The arguments \code{nrow}, \code{ncol}, \code{byrow}, and \code{dimnames} are made use of only if \code{data} is not a matrix or magma object.
}

\details{
If either of \code{nrow} or \code{ncol} is not given, an attempt is made to infer its value from the length of \code{data} and the other parameter.

The MAGMA library of C functions offers two LAPACK-style interfaces, referred to
as the \bold{GPU interface} and the \bold{CPU interface}.  The GPU interface
takes input and produces results in the GPU's memory, whereas the CPU interface
produces results in the CPU's memory.  The GPU and CPU interfaces, although
similar, are not derivatives of each other, but instead have different
communication patters.  The \code{gpu} argument allows specification of the
interface to be used for caluculations involving the resulting magma matrix.
}

\author{
Brian J. Smith <brian-j-smith@uiowa.edu>
}

\seealso{
\code{\link[magma]{gpu}},
\code{\linkS4class{magma}},
\code{\linkS4class{matrix}}
}

\examples{
magma(1:6, nrow = 2, ncol = 3, byrow = TRUE,
      dimnames = list(c("row1", "row2"), c("col1", "col2", "col3")),
      gpu = TRUE)
}

\keyword{array}
\keyword{algebra}

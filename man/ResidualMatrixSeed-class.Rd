% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seed.R
\docType{class}
\name{ResidualMatrixSeed-class}
\alias{ResidualMatrixSeed-class}
\alias{ResidualMatrixSeed}
\alias{dim,ResidualMatrixSeed-method}
\alias{dimnames,ResidualMatrixSeed-method}
\alias{extract_array,ResidualMatrixSeed-method}
\alias{DelayedArray,ResidualMatrixSeed-method}
\alias{show,ResidualMatrixSeed-method}
\title{The ResidualMatrixSeed class}
\description{
This is a seed class that powers the \pkg{DelayedArray} machinery underlying the \linkS4class{ResidualMatrix}.
}
\section{Construction}{

\code{ResidualMatrixSeed(x, design=NULL, keep=NULL)} returns a ResidualMatrixSeed object, given:
\itemize{
\item \code{x}, a matrix-like object.
This can alternatively be a ResidualMatrixSeed, in which case \code{design} is ignored.
\item \code{design}, a numeric matrix containing the experimental design,
to be used for linear model fitting on each \emph{column} of \code{x}.
This defaults to an intercept-only matrix.
\item \code{keep}, an integer vector specifying the columns of \code{design} to \emph{not} regress out.
By default, all columns of \code{design} are regressed out.
\item \code{restrict}, an integer or logical vector specifying the rows of \code{x} to use for model fitting.
If \code{NULL}, all rows of \code{x} are used.
}
}

\section{Methods}{

ResidualMatrixSeed objects are implemented as \linkS4class{DelayedMatrix} backends.
They support standard operations like \code{dim}, \code{dimnames} and \code{extract_array}.

Passing a ResidualMatrixSeed object to the \code{\link{DelayedArray}} or \code{\link{ResidualMatrix}} constructors
will create a \linkS4class{ResidualMatrix} (which is what most users should be working with, anyway).
}

\examples{
design <- model.matrix(~gl(5, 50))

library(Matrix)
y0 <- rsparsematrix(nrow(design), 200, 0.1)
s <- ResidualMatrixSeed(y0, design)
s

ResidualMatrix(s)

DelayedArray(s)
}
\author{
Aaron Lun
}

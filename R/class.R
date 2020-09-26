#' The ResidualMatrix class
#'
#' The ResidualMatrix class supports delayed calculation of the residuals from a linear model fit.
#' This serves as a light-weight representation of what would otherwise be a large dense matrix in memory.
#' It also enables efficient matrix multiplication based on features of the the original matrix (e.g., sparsity).
#' 
#' @section Construction:
#' \code{ResidualMatrix(x, design=NULL, keep=NULL)} returns a ResidualMatrix object, given:
#' \itemize{
#' \item \code{x}, a matrix-like object.
#' This can alternatively be a ResidualMatrixSeed, in which case \code{design} and \code{keep} are ignored.
#' \item \code{design}, a numeric matrix containing the experimental design,
#' to be used for linear model fitting on each \emph{column} of \code{x}.
#' This defaults to an intercept-only matrix.
#' \item \code{keep}, an integer vector specifying the columns of \code{design} to \emph{not} regress out.
#' By default, all columns of \code{design} are regressed out.
#' \item \code{restrict}, an integer or logical vector specifying the rows of \code{x} to use for model fitting.
#' If \code{NULL}, all rows of \code{x} are used.
#' }
#' 
#' When \code{keep=NULL}, the ResidualMatrix contains values equivalent to \code{lm.fit(x=design, y=x)$residuals}.
#' 
#' @section Methods:
#' In the following code chunks, \code{x} is a ResidualMatrix object:
#' \itemize{ 
#' \item \code{x[i, j, .., drop=FALSE]} will return a ResidualMatrix object for the specified row and column subsets,
#' or a numeric vector if either \code{i} or \code{j} are of length 1.
#' \item \code{t(x)} will return a ResidualMatrix object with transposed contents.
#' \item \code{dimnames(x) <- value} will return a ResidualMatrix object where the rows and columns are renamed by \code{value},
#' a list of two character vectors (or \code{NULL}).
#' }
#'
#' \code{\link{colSums}(x)}, \code{\link{colMeans}(x)}, \code{\link{rowSums}(x)} and \code{\link{rowMeans}(x)}
#' will return the relevant statistics for a ResidualMatrix \code{x}.
#' 
#' \code{\%*\%}, \code{\link{crossprod}} and \code{\link{tcrossprod}} can also be applied
#' where one or both of the arguments are ResidualMatrix objects.
#'
#' ResidualMatrix objects are derived from \linkS4class{DelayedMatrix} objects and support all of valid operations on the latter.
#' All operations not listed here will use the underlying \pkg{DelayedArray} machinery.
#' Unary or binary operations will generally create a new DelayedMatrix instance containing a \linkS4class{ResidualMatrixSeed}.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' design <- model.matrix(~gl(5, 50))
#' 
#' library(Matrix)
#' y0 <- rsparsematrix(nrow(design), 200, 0.1)
#' y <- ResidualMatrix(y0, design)
#' y
#' 
#' # For comparison:
#' fit <- lm.fit(x=design, y=as.matrix(y0))
#' DelayedArray(fit$residuals)
#' 
#' # Keeping some of the factors:
#' y2 <- ResidualMatrix(y0, design, keep=1:2)
#' y2
#' DelayedArray(fit$residuals + design[,1:2] %*% fit$coefficients[1:2,])
#' 
#' # Matrix multiplication:
#' crossprod(y)
#' tcrossprod(y)
#' y %*% rnorm(200)
#' 
#' @aliases
#' ResidualMatrix
#' ResidualMatrix-class
#' dimnames<-,ResidualMatrix,ANY-method
#' t,ResidualMatrix-method
#' [,ResidualMatrix,ANY,ANY,ANY-method
#' colSums,ResidualMatrix-method
#' rowSums,ResidualMatrix-method
#' colMeans,ResidualMatrix-method
#' rowMeans,ResidualMatrix-method
#' %*%,ANY,ResidualMatrix-method
#' %*%,ResidualMatrix,ANY-method
#' %*%,ResidualMatrix,ResidualMatrix-method
#' crossprod,ResidualMatrix,missing-method
#' crossprod,ResidualMatrix,ANY-method
#' crossprod,ANY,ResidualMatrix-method
#' crossprod,ResidualMatrix,ResidualMatrix-method
#' tcrossprod,ResidualMatrix,missing-method
#' tcrossprod,ResidualMatrix,ANY-method
#' tcrossprod,ANY,ResidualMatrix-method
#' tcrossprod,ResidualMatrix,ResidualMatrix-method
#' 
#' @name ResidualMatrix-class
#' @docType class
NULL

#' @export
#' @importFrom DelayedArray DelayedArray
ResidualMatrix <- function(x, design=NULL, keep=NULL, restrict=NULL) {
    DelayedArray(ResidualMatrixSeed(x, design, keep, restrict))
}

#' @export
#' @importFrom DelayedArray DelayedArray new_DelayedArray
setMethod("DelayedArray", "ResidualMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="ResidualMatrix")
)

###################################
# Overridden utilities from DelayedArray, for efficiency.

#' @export
#' @importFrom DelayedArray DelayedArray seed
setReplaceMethod("dimnames", "ResidualMatrix", function(x, value) {
    DelayedArray(rename_ResidualMatrixSeed(seed(x), value))
})

#' @export
#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray seed
setMethod("t", "ResidualMatrix", function(x) {
    DelayedArray(transpose_ResidualMatrixSeed(seed(x)))
})

#' @export
#' @importFrom DelayedArray DelayedArray seed
setMethod("[", "ResidualMatrix", function(x, i, j, ..., drop=TRUE) {
    if (missing(i)) i <- NULL
    if (missing(j)) j <- NULL

    rseed <- seed(x)
    if (was_transposed <- is_transposed(rseed)) {
        rseed <- transpose_ResidualMatrixSeed(rseed)
        tmp <- i
        i <- j
        j <- tmp
    }

    rseed <- subset_ResidualMatrixSeed(rseed, i, j)
    if (was_transposed) {
        rseed <- transpose_ResidualMatrixSeed(rseed)
    }

    out <- DelayedArray(rseed)
    if (drop && any(dim(out)==1L)) {
        return(drop(out))
    }
    out
})



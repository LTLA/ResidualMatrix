#' The ResidualMatrixSeed class
#'
#' This is a seed class that powers the \pkg{DelayedArray} machinery underlying the \linkS4class{ResidualMatrix}.
#' 
#' @section Construction:
#' \code{ResidualMatrixSeed(x, design=NULL, keep=NULL)} returns a ResidualMatrixSeed object, given:
#' \itemize{
#' \item \code{x}, a matrix-like object.
#' This can alternatively be a ResidualMatrixSeed, in which case \code{design} is ignored.
#' \item \code{design}, a numeric matrix containing the experimental design,
#' to be used for linear model fitting on each \emph{column} of \code{x}.
#' This defaults to an intercept-only matrix.
#' \item \code{keep}, an integer vector specifying the columns of \code{design} to \emph{not} regress out.
#' By default, all columns of \code{design} are regressed out.
#' }
#' 
#' @section Methods:
#' ResidualMatrixSeed objects are implemented as \linkS4class{DelayedMatrix} backends.
#' They support standard operations like \code{dim}, \code{dimnames} and \code{extract_array}.
#' 
#' Passing a ResidualMatrixSeed object to the \code{\link{DelayedArray}} or \code{\link{ResidualMatrix}} constructors
#' will create a \linkS4class{ResidualMatrix} (which is what most users should be working with, anyway).
#'
#' @examples
#' design <- model.matrix(~gl(5, 50))
#' 
#' library(Matrix)
#' y0 <- rsparsematrix(nrow(design), 200, 0.1)
#' s <- ResidualMatrixSeed(y0, design)
#' s
#' 
#' ResidualMatrix(s)
#'
#' DelayedArray(s)
#' @author Aaron Lun
#'
#' @aliases
#' ResidualMatrixSeed
#' ResidualMatrixSeed-class
#' dim,ResidualMatrixSeed-method
#' dimnames,ResidualMatrixSeed-method
#' extract_array,ResidualMatrixSeed-method
#' DelayedArray,ResidualMatrixSeed-method
#' show,ResidualMatrixSeed-method
#' @docType class
#' @name seed
NULL

#' @export
#' @importFrom methods new is
#' @importFrom Matrix crossprod
ResidualMatrixSeed <- function(x, design=NULL, keep=NULL) {
    if (missing(x)) {
        x <- matrix(0, 0, 0)
    } else if (is(x, "ResidualMatrixSeed")) {
        return(x)
    } 

    if (is.null(design)) {
        design <- matrix(1, nrow(x), 1)
    } else if (nrow(design)!=nrow(x)) {
        stop("'nrow(x)' and 'nrow(design)' should be equal")
    }

    QR <- qr(design)
    if (QR$rank < ncol(design)) {
        stop("'design' does not appear to be of full rank")
    }
    Q <- as.matrix(qr.Q(QR))
    Qty <- as.matrix(crossprod(Q, x))

    if (!is.null(keep)) {
        # We want to add back X %*% coef', where non-kept coefs are set to zero in coef'.
        # As X %*% coef' = QR %*% coef' = Q (R %*% coef'), we can just subtract (R %*% coef') from Qty.
        # We subtract here because Qty is itself subtracted to obtain the residuals,
        # so if we want to add X %*% coef' back to the residuals, we need to subtract here.
        R <- qr.R(QR)
        coefs <- backsolve(R, Qty)
        coefs[!QR$pivot %in% keep,] <- 0 # I've never been able to get a full-rank example that uses non-trivial pivoting.
        Qty <- Qty - R %*% coefs
    }

    new("ResidualMatrixSeed", .matrix=x, Q=Q, Qty=Qty, transposed=FALSE)
}

#' @importFrom S4Vectors setValidity2
#' @importFrom methods is
setValidity2("ResidualMatrixSeed", function(object) {
    msg <- character(0)

    x <- get_matrix2(object)
    Q <- get_Q(object)
    if (nrow(x)!=nrow(Q)) {
        msg <- c(msg, "'nrow(x)' and 'nrow(Q)' are not the same")
    }
    if (!is.numeric(Q)) {
        msg <- c(msg, "'Q' should be a numeric matrix")
    }

    Qty <- get_Qty(object)
    if (ncol(x)!=ncol(Qty)) {
        msg <- c(msg, "'ncol(x)' and 'ncol(Qty)' are not the same")
    }
    if (ncol(Q)!=nrow(Qty)) {
        msg <- c(msg, "'ncol(Q)' and 'nrow(Qty)' are not the same")
    }
    if (!is.numeric(Q)) {
        msg <- c(msg, "'Qty' should be a numeric matrix")
    }

    if (length(is_transposed(object))!=1L) {
        msg <- c(msg, "'transposed' must be a logical scalar")
    } 

    if (length(msg)) {
        return(msg)
    } 
    TRUE
})

#' @export
#' @importFrom methods show
setMethod("show", "ResidualMatrixSeed", function(object) {
    cat(sprintf("%i x %i ResidualMatrixSeed object", nrow(object), ncol(object)),
    sep="\n")
})

###################################
# Internal getters.

get_matrix2 <- function(x) x@.matrix

get_Q <- function(x) x@Q

get_Qty <- function(x) x@Qty

is_transposed <- function(x) x@transposed

###################################
# DelayedArray support utilities. 

#' @export
setMethod("dim", "ResidualMatrixSeed", function(x) {
    d <- dim(get_matrix2(x))
    if (is_transposed(x)) { d <- rev(d) }
    d
})

#' @export
setMethod("dimnames", "ResidualMatrixSeed", function(x) {
    d <- dimnames(get_matrix2(x))
    if (is_transposed(x)) { d <- rev(d) }
    d
})

#' @export
#' @importFrom DelayedArray extract_array
#' @importFrom Matrix t crossprod
setMethod("extract_array", "ResidualMatrixSeed", function(x, index) {
    if (was_transposed <- is_transposed(x)) {
        x <- transpose_ResidualMatrixSeed(x)
        index <- rev(index)
    }
	x2 <- subset_ResidualMatrixSeed(x, index[[1]], index[[2]])
    resid <- get_matrix2(x2) - get_Q(x2) %*% get_Qty(x2)
    if (was_transposed) {
        resid <- t(resid)
    }
    as.matrix(resid)
})

###################################
# Additional utilities for efficiency.

subset_ResidualMatrixSeed <- function(x, i, j) {
    mat <- get_matrix2(x)
    Q <- get_Q(x)
    Qty <- get_Qty(x)

    if (!is.null(i)) {
        if (is.character(i)) {
            i <- match(i, rownames(mat))
        }
        mat <- mat[i,,drop=FALSE]
        Q <- Q[i,,drop=FALSE]
    }

    if (!is.null(j)) {
        if (is.character(j)) {
            j <- match(j, colnames(mat))
        }
        mat <- mat[,j,drop=FALSE]
        Qty <- Qty[,j,drop=FALSE]
    }

    initialize(x, .matrix=mat, Q=Q, Qty=Qty)
}

transpose_ResidualMatrixSeed <- function(x) {
    initialize(x, transposed=!is_transposed(x))
}

rename_ResidualMatrixSeed <- function(x, value) {
    x2 <- get_matrix2(x)
    if (is_transposed(x)) {
        value <- rev(value)
    }
    dimnames(x2) <- value
    initialize(x, .matrix=x2)
}

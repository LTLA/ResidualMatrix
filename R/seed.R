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
#' \item \code{restrict}, an integer or logical vector specifying the rows of \code{x} to use for model fitting.
#' If \code{NULL}, all rows of \code{x} are used.
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
#' @name ResidualMatrixSeed-class
NULL

#' @export
#' @importFrom Matrix crossprod
ResidualMatrixSeed <- function(x, design=NULL, keep=NULL, restrict=NULL) {
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
    if (QR$rank < ncol(design) && nrow(design)!=0) {
        stop("'design' does not appear to be of full rank")
    }

    Q <- as.matrix(qr.Q(QR))

    if (is.null(restrict) && is.null(keep)) {
        Qty <- as.matrix(crossprod(Q, x))
    } else if (!is.null(restrict)) {
        subdesign <- design[restrict,,drop=FALSE]
        subQR <- qr(subdesign)
        if (subQR$rank < ncol(subdesign) && nrow(subdesign)!=0) {
            stop("'design[restrict,]' does not appear to be of full rank")
        }

        subQ <- as.matrix(qr.Q(subQR))
        subQty <- as.matrix(crossprod(subQ, x[restrict,,drop=FALSE]))
        subR <- qr.R(subQR)
        coefs <- backsolve(subR, subQty)

        if (!is.null(keep)) {
            # See the logic in the next clause.
            coefs[subQR$pivot %in% keep,] <- 0 
        }
        coefs <- coefs[match(QR$pivot, subQR$pivot),,drop=FALSE]

        # We want to subtract the fitted values obtained by X %*% coefs,
        # where coefs is computed using only the restricted subset of samples.
        # This leads use to realize that X %*% coefs is just Q (R %*% coefs).
        R <- qr.R(QR)
        Qty <- R %*% coefs
    } else {
        # We want to subtract the fitted values X %*% coef', where kept coefs are set to zero in coef'.
        # As X %*% coef' = QR %*% coef' = Q (R %*% coef'), we can use (R %*% coef') as our Qty.
        R <- qr.R(QR)
        Qty <- as.matrix(crossprod(Q, x))
        coefs <- backsolve(R, Qty)
        coefs[QR$pivot %in% keep,] <- 0 
        Qty <- R %*% coefs
    }

    new("ResidualMatrixSeed", .matrix=x, Q=Q, Qty=Qty, transposed=FALSE)
}

#' @importFrom S4Vectors setValidity2
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
    resid <- as.matrix(get_matrix2(x2)) - as.matrix(get_Q(x2) %*% get_Qty(x2))
    if (was_transposed) {
        resid <- t(resid)
    }
    resid
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

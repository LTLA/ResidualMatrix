###################################
# Matrix multiplication.

# Note that "%*%" methods should NOT call other matrix multiplication methods
# directly on their ResidualMatrix arguments. Rather, any ResidualMatrix should
# be broken down into the seed or the underlying matrix before further multiplication.
# This reduces the risk of infinite S4 recursion when 'y' is also an S4 matrix class. 
# 
# Specifically, the .*_ResidualMatrix functions take a seed object that IGNORES
# any non-FALSE setting of @transposed. It will then break down the seed into 
# its constituents and perform the multiplication, such that any further S4
# dispatch occurs on the lower components of the seed. 
#
# We also coerce each matrix product to a full matrix - which it usually is, anyway 
# - to avoid the unnecessary overhead of multiplying DelayedArray instances.

#' @export
#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray seed
setMethod("%*%", c("ResidualMatrix", "ANY"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- t(.leftmult_ResidualMatrix(t(y), x_seed))
    } else {
        out <- .rightmult_ResidualMatrix(x_seed, y)
    }
    DelayedArray(out)
})

#' @importFrom Matrix crossprod
.rightmult_ResidualMatrix <- function(x_seed, y) {
    # Order of operations chosen to minimize size of intermediates,
    # under the assumption that ncol(y) is very small.
    as.matrix(get_matrix2(x_seed) %*% y) - get_Q(x_seed) %*% as.matrix(get_Qty(x_seed) %*% y)
}

#' @export
#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray seed
setMethod("%*%", c("ANY", "ResidualMatrix"), function(x, y) {
    y_seed <- seed(y)
    if (is_transposed(y_seed)) {
        if (!is.null(dim(x))) {
            # Vectors don't quite behave as 1-column matrices here.
            # so we need to be a bit more careful.
            x <- t(x) 
        }
        out <- t(.rightmult_ResidualMatrix(y_seed, x))
    } else {
        out <- .leftmult_ResidualMatrix(x, y_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix tcrossprod
.leftmult_ResidualMatrix <- function(x, y_seed) {
    # Order of operations chosen to minimize size of intermediates,
    # under the assumption that nrow(x) is very small.
    as.matrix(x %*% get_matrix2(y_seed)) - as.matrix(x %*% get_Q(y_seed)) %*% get_Qty(y_seed)
}

#' @export
#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray seed
setMethod("%*%", c("ResidualMatrix", "ResidualMatrix"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- t(.leftmult_ResidualMatrix(t(y), x_seed))
    } else {
        out <- .rightmult_ResidualMatrix(x_seed, y)
    }
    DelayedArray(out)
})

# Unlike DeferredMatrix, we are happy to delegate to the left/right %*% 
# for dual ResidualMatrix multiplication. This is because there are no
# operations that will collapse a ResidualMatrix instance to a DelayedMatrix
# prior to multiplication. Thus, we do not have to worry about writing
# specialized methods to avoid the overhead of DelayedMatrix multiplication.

###################################
# Crossproduct.

# Technically, this could be defined in terms of '%*%'.
# However, we re-define it in terms of 'crossprod' for efficiency,
# to exploit the potential availability of optimized crossprod for .matrix.

#' @export
#' @importFrom Matrix crossprod
#' @importFrom DelayedArray DelayedArray seed
setMethod("crossprod", c("ResidualMatrix", "missing"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        #  No need for t(), it's symmetric.
        out <- .tcp_ResidualMatrix(x_seed)
    } else {
        out <- .crossprod_ResidualMatrix(x_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix crossprod 
.crossprod_ResidualMatrix <- function(x_seed) {
    mat <- get_matrix2(x_seed)
    Qty <- get_Qty(x_seed)
    Q <- get_Q(x_seed)

    # We assume that nrow(Q) >> ncol(Q) and nrow(mat) >> ncol(mat) 
    # in order for this to be efficient.
    ytQQty <- crossprod(Qty)
    QtQ <- crossprod(Q)

    # Using this addition order to minimize numeric instability.
    (crossprod(mat) - ytQQty) + (crossprod(Qty, QtQ %*% Qty) - ytQQty)
}

#' @export
#' @importFrom Matrix crossprod
#' @importFrom DelayedArray DelayedArray seed
setMethod("crossprod", c("ResidualMatrix", "ANY"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- .rightmult_ResidualMatrix(x_seed, y)
    } else {
        out <- .rightcross_ResidualMatrix(x_seed, y)
    }
    DelayedArray(out)
})

#' @importFrom Matrix crossprod
.rightcross_ResidualMatrix <- function(x_seed, y) {
    # Order of operations chosen to minimize size of intermediates,
    # under the assumption that ncol(y) is very small.
    as.matrix(crossprod(get_matrix2(x_seed), y)) - crossprod(get_Qty(x_seed), as.matrix(crossprod(get_Q(x_seed), y)))
}

#' @export
#' @importFrom Matrix crossprod t
#' @importFrom DelayedArray DelayedArray seed
setMethod("crossprod", c("ANY", "ResidualMatrix"), function(x, y) {
    y_seed <- seed(y)
    if (is_transposed(y_seed)) {
        out <- t(.rightmult_ResidualMatrix(y_seed, x))
    } else {
        out <- .leftcross_ResidualMatrix(x, y_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix crossprod
.leftcross_ResidualMatrix <- function(x, y_seed) {
    # Order of operations chosen to minimize size of intermediates,
    # under the assumption that ncol(x) is very small.
    as.matrix(crossprod(x, get_matrix2(y_seed))) - as.matrix(crossprod(x, get_Q(y_seed))) %*% get_Qty(y_seed)
}

#' @export
#' @importFrom Matrix crossprod
#' @importFrom DelayedArray DelayedArray seed
setMethod("crossprod", c("ResidualMatrix", "ResidualMatrix"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- .rightmult_ResidualMatrix(x_seed, y)
    } else {
        out <- .rightcross_ResidualMatrix(x_seed, y)
    }
    DelayedArray(out)
})

###################################
# Transposed crossproduct.

# Technically, this could be defined in terms of '%*%'.
# However, we re-define it in terms of 'tcrossprod' for efficiency,
# to exploit the potential availability of optimized crossprod for .matrix.

#' @export
#' @importFrom Matrix tcrossprod
#' @importFrom DelayedArray DelayedArray seed
setMethod("tcrossprod", c("ResidualMatrix", "missing"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- .crossprod_ResidualMatrix(x_seed)
    } else {
        out <- .tcp_ResidualMatrix(x_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix tcrossprod 
.tcp_ResidualMatrix <- function(x_seed) {
    mat <- get_matrix2(x_seed)
    Qty <- get_Qty(x_seed)
    Q <- get_Q(x_seed)

    # We assume that ncol(mat) >> nrow(mat) for this to be efficient.
    # We also try to avoid constructing QQt under the assumption that 
    # nrow(Q) >> ncol(Q) for Q derived from a full-rank design matrix.
    QQtyyt<- Q %*% tcrossprod(Qty, mat)
    QQtyytQQt <- tcrossprod(QQtyyt %*% Q, Q)

    # Using this addition order to minimize numeric instability.
    (tcrossprod(mat) - QQtyyt) + (QQtyytQQt - t(QQtyyt))
}

#' @export
#' @importFrom Matrix tcrossprod t
#' @importFrom DelayedArray DelayedArray seed
setMethod("tcrossprod", c("ResidualMatrix", "ANY"), function(x, y) {
    if (is.null(dim(y))) { # for consistency with base::tcrossprod.
        stop("non-conformable arguments")
    }

    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- t(.leftmult_ResidualMatrix(y, x_seed))
    } else {
        out <- .righttcp_ResidualMatrix(x_seed, y)
    }
    DelayedArray(out)
})

#' @importFrom Matrix tcrossprod
.righttcp_ResidualMatrix <- function(x_seed, y) {
    # Order of operations chosen to minimize size of intermediates,
    # assuming that nrow(y) is very small.
    as.matrix(tcrossprod(get_matrix2(x_seed), y)) - get_Q(x_seed) %*% as.matrix(tcrossprod(get_Qty(x_seed), y))
}

#' @export
#' @importFrom Matrix tcrossprod
#' @importFrom DelayedArray DelayedArray seed
setMethod("tcrossprod", c("ANY", "ResidualMatrix"), function(x, y) {
    y_seed <- seed(y)
    if (is_transposed(y_seed)) {
        out <- .leftmult_ResidualMatrix(x, y_seed)
    } else {
        out <- .lefttcp_ResidualMatrix(x, y_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix tcrossprod
.lefttcp_ResidualMatrix <- function(x, y_seed) {
    # Order of operations chosen to minimize size of intermediates.
    # assuming that nrow(x) is very small.
    as.matrix(tcrossprod(x, get_matrix2(y_seed))) - tcrossprod(as.matrix(tcrossprod(x, get_Qty(y_seed))), get_Q(y_seed))
}

#' @export
#' @importFrom Matrix tcrossprod t
#' @importFrom DelayedArray DelayedArray seed
setMethod("tcrossprod", c("ResidualMatrix", "ResidualMatrix"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- t(.leftmult_ResidualMatrix(y, x_seed))
    } else {
        out <- .righttcp_ResidualMatrix(x_seed, y)
    }
    DelayedArray(out)
})

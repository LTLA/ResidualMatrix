###################################
# Basic matrix stats.

#' @export
#' @importFrom Matrix colSums rowSums drop
setMethod("colSums", "ResidualMatrix", function(x, na.rm = FALSE, dims = 1L) {
    if (is_transposed(seed(x))) {
        return(rowSums(t(x)))
    }

    out <- rep(1, nrow(x)) %*% x
    out <- drop(out)
    names(out) <- colnames(x)
    out
})

#' @export
#' @importFrom Matrix colSums rowSums drop
setMethod("rowSums", "ResidualMatrix", function(x, na.rm = FALSE, dims = 1L) {
    if (is_transposed(seed(x))) {
        return(colSums(t(x)))
    }

    out <- x %*% rep(1, ncol(x))
    out <- drop(out)
    names(out) <- rownames(x)
    out
})

#' @export
#' @importFrom Matrix colMeans colSums
setMethod("colMeans", "ResidualMatrix", function(x, na.rm = FALSE, dims = 1L) colSums(x)/nrow(x))

#' @export
#' @importFrom Matrix rowMeans rowSums
setMethod("rowMeans", "ResidualMatrix", function(x, na.rm = FALSE, dims = 1L) rowSums(x)/ncol(x))



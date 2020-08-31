#' @export
#' @import methods
setClass("ResidualMatrixSeed", slots=c(.matrix="ANY", Q="matrix", Qty="matrix", transposed="logical"))

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("ResidualMatrix",
    contains="DelayedMatrix",
    representation(seed="ResidualMatrixSeed")
)

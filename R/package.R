#' The \pkg{ResidualMatrix} package
#'
#' Originally implemented in the \pkg{BiocSingular} package,
#' the \code{ResidualMatrix} class has been placed into its own package to enable greater re-use.
#' This class provides delayed computation of residuals from a linear model fit,
#' allowing us to represent large matrices of residuals without actually calculating them in memory.
#' The idea is to allow us to easily regress out uninteresting factors for big datasets,
#' much like a delayed, scalable version of \pkg{limma}'s venerable \code{removeBatchEffect} function. 
#' 
#' @name ResidualMatrix-package
NULL

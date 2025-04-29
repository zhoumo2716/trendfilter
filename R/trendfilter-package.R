#' @importFrom Rcpp sourceCpp
#' @importFrom checkmate assert_scalar assert_logical assert_numeric
#' @importFrom checkmate assert_character assert_integerish assert_class
#' @importFrom cli cli_abort cli_warn
#' @importFrom rlang arg_match %||%
#' @importFrom dspline dspline_interp
#' @importFrom tvdenoising tvdenoising
## usethis namespace: start
#' @useDynLib trendfilter, .registration = TRUE
## usethis namespace: end
NULL

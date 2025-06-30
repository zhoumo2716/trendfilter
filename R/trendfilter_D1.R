#' Predict for a two-penalty trend-filter fit
#'
#' @param object  a list returned by trendfilter_D1_cpp / trendfilter_D1
#'                Must contain elements  x   and  theta .
#'                theta can be a vector or an n × L matrix.
#' @param x       numeric vector of new locations to predict at.
#'                Default = original grid stored in the object.
#' @param lambdaindex column of theta to use (ignored if theta is a vector).
#' @param ...     unused (keeps R CMD check happy + generic compatibility)
#'
#' @return numeric vector of fitted values at x
#' @export
predict.trendfilter_D1 <- function(object,
                                   x           = object$x,
                                   lambdaindex = NULL,
                                   ...) {

  if (is.null(object$theta))
    stop("object has no 'theta' slot")

  ## ---- pick the right column ------------------------------------------
  th <- object$theta
  if (is.null(dim(th))) {                   # theta is a *vector*
    fitted <- th
  } else {                                  # theta is a matrix
    if (is.null(lambdaindex))
      lambdaindex <- ncol(th)               # default = last column
    if (length(lambdaindex) != 1L ||
        lambdaindex < 1L ||
        lambdaindex > ncol(th))
      stop("'lambdaindex' out of range")
    fitted <- th[, lambdaindex]
  }

  ## ---- interpolate / extrapolate --------------------------------------
  stats::approx(x = object$x, y = fitted, xout = x,
                rule = 2L,  # linear inside range, constant outside
                ties = "ordered")$y
}

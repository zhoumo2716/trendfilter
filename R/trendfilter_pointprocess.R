#' Trend Filtering for Point Process Intensity Estimation
#'
#' This function performs trend filtering
#' on a Poisson point process to estimate the log-intensity function.
#'
#' @param x Numeric vector of event times.
#' @param k Integer, order of the trend filter (e.g., 0 = piecewise constant, 1 = piecewise linear).
#' @param A Numeric, start of the observation window.
#' @param B Numeric, end of the observation window.
#' @param lambda Numeric, regularization parameter (default = 1).
#' @param rho_scale Numeric, scaling factor for ADMM penalty parameter (default = 1e-8).
#' @param max_iter Integer, maximum ADMM iterations (default = 100).
#' @param tol Numeric, tolerance for ADMM convergence (default = 1e-5).
#' @param newton_max_iters Integer, maximum Newton iterations (default = 50).
#' @param newton_tol Numeric, tolerance for Newton solver (default = 1e-5).
#'
#' @return A list containing estimated coefficients, intensity values, diagnostics, etc.
#'
#' When \code{k = 0}, performs total variation denoising (TV denoising);
#' otherwise, calls the Rcpp backend \code{trendfilter_pointProcess()}.
#'
#' @importFrom tvdenoising tvdenoising
#' @export
#'
#' @examples
#' \dontrun{
#'   # Example data
#'   set.seed(1)
#'   x <- sort(runif(50, 0, 1))
#'   result <- trendfilter_pointprocess(x, k = 1, A = 0, B = 1, lambda = 0.1)
#'   print(result$theta)
#' }
trendfilter_pointprocess <- function(data, A, B, k = 0,
    lambda = 1,
    rho_scale = 1e-8,
    max_iter = 100,
    tol = 1e-5,
    newton_max_iters = 50,
    newton_tol = 1e-5
) {
  # Use tvdenosing when k = 0
  if (k==0) {
    ss <- c(data[1] - A, diff(data), B - data[length(data)])
    ss <- pmax(ss, 1e-12)
    y <- 1/ss
    w <- ss/2
    theta <- tvdenoising::tvdenoising(y=y, w=w,lambda=1)
  } else {
    theta <- trendfilter_pointProcess(
      x = data,
      k = k,
      A = A,
      B = B,
      lambda = lambda,
      rho_scale = rho_scale,
      max_iter = max_iter,
      tol = tol,
      newton_max_iters = newton_max_iters,
      newton_tol = newton_tol
    )$theta
  }
}

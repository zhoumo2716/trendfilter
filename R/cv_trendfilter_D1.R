#' Cross-validation for two-penalty trend filtering
#'
#' Performs leave-v-th-out cross-validation over a **2-D grid**:
#' a sequence of \eqn{\lambda_k} values (high-order penalty) and a
#' *scalar* grid that multiplies \eqn{\lambda_k} to obtain
#' \eqn{\lambda_1}.
#'
#' For each combination of \eqn{(\lambda_k,\;\rho)} it
#' \enumerate{
#'   \item fits \code{trendfilter_D1()} on the training folds,
#'   \item predicts on the held-out fold via \code{predict.trendfilter_D1()},
#'   \item stores the mean-squared error.
#' }
#' The output gives the full CV error surface and, for every \eqn{\lambda_k},
#' the scalar that minimises the CV score.
#'
#' @param y Numeric vector of observations.
#' @param x Numeric vector of the same length as \code{y};
#'   defaults to \code{1:length(y)}.
#' @param k Polynomial order of the trend filter (default 3).
#' @param lambda1_scalar Numeric vector of multipliers
#'   \eqn{\rho} so that \eqn{\lambda_1 = \rho\,\lambda_k}.
#' @param nfolds Integer \(\ge 2\).  Number of CV folds (leave-v-th-out).
#' @param nlambda Number of distinct \eqn{\lambda_k} values along the path.
#' @param lambda_max Optional; maximum \eqn{\lambda_k}.
#'   If \code{NULL} it is computed via \code{get_lambda_max()}.
#' @param lambda_min_ratio Smallest \eqn{\lambda_k} is
#'   \code{lambda_max * lambda_min_ratio}.
#' @param ... Additional arguments forwarded to \code{trendfilter_D1()}.
#'
#' @return A list of class \code{"cv_trendfilter_D1"} with elements
#'   \describe{
#'     \item{\code{lambdak}}{numeric vector of \eqn{\lambda_k}.}
#'     \item{\code{lambda1_scalar}}{the scalar grid supplied.}
#'     \item{\code{cv_scores}}{matrix \eqn{nlambda \times length(lambda1_scalar)}
#'       of averaged MSEs.  Row = \eqn{\lambda_k}, column = scalar.}
#'     \item{\code{scalar_min}}{data frame of length \eqn{nlambda} giving,
#'       for each \eqn{\lambda_k}, the scalar that minimises CV error.}
#'     \item{\code{call}}{the matched call.}
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' x <- 1:40
#' y <- c(rep(3, 10), rep(2, 15), seq(2, -2, length = 15))
#' cv2d <- cv_trendfilter_D1(y, x, k = 2)
#' head(cv2d$scalar_min)
#' }
#'
#' @export
cv_trendfilter_D1 <- function(
    y,
    x        = seq_along(y),
    k        = 3L,
    lambda1_scalar = c(0.1, 0.2, 0.5, 0.7, 1, 1.5, 2, 10),  # scalar grid
    nfolds   = 5L,
    nlambda  = 50L,
    lambda_max       = NULL,
    lambda_min_ratio = 1e-5,
    ...
) {
  n <- length(y)

  ## --- 1. build lambdak path ------------------------------------------
  if (is.null(lambda_max)) {
    lambda_max <- get_lambda_max(x, y, rep(1, n), k)
  }

  lambdak <- numeric(nlambda)
  lambdak <- get_lambda_seq_r(lambdak, lambda_max,lambda_min = -1,lambda_min_ratio, nlambda)

  ## --- 2. fold assignment (leave-v-th-out) -----------------------------
  foldid <- c(0, rep_len(1:nfolds, n - 2), 0)   # endpoints never held out
  mse    <- function(a, b) (a - b)^2

  ## cv_score[row = lambdak , col = lambda1_scalar]
  cv_score <- matrix(NA, nrow = nlambda, ncol = length(lambda1_scalar))
  best_scalar <- double(nlambda)

  ## --- 3. double loop: lambdak  ×  scalar ------------------------------
  for (m in seq_along(lambdak)) {
    lamk     <- lambdak[m]
    fold_scores <- matrix(NA, nfolds, length(lambda1_scalar))

    for (s in seq_along(lambda1_scalar)) {
      rho <- lambda1_scalar[s]

      for (f in 1:nfolds) {
        tr <- foldid != f
        te <- !tr

        fit_tr <- trendfilter_D1(
          k        = k,
          lambda1_scalar = rho,
          lambdak  = lamk,
          y        = y[tr],
          x        = x[tr]
        )

        y_hat <- predict.trendfilter_D1(fit_tr, x[te])
        fold_scores[f, s] <- mean(mse(y[te], y_hat))
      }
    }

    ## average along columns (over folds) -> one score per scalar
    cv_score[m, ]       <- colMeans(fold_scores, na.rm = TRUE)
    best_scalar_index <- which.min(cv_score[m, ])
    best_scalar[m] <- lambda1_scalar[best_scalar_index]
  }

  ## --- 4. list of (lambdak, best scalar) pairs -------------------------
  scalar_min <- data.frame(
    lambdak = lambdak,
    best_lambda1_scalar = best_scalar
  )

  ## package output ------------------------------------------------------
  out <- list(
    lambdak         = lambdak,
    lambda1_scalar  = lambda1_scalar,
    cv_scores       = cv_score,
    scalar_min      = scalar_min,
    call            = match.call()
  )
  class(out) <- "cv_trendfilter_D1"
  return(out)
}


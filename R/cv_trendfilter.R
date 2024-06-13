#' Leave-vth-out cross validation for trendfilter
#'
#' @inheritParams trendfilter
#' @param error_measure Metric used to calculate cross validation scores. May
#'   be `mse`, `mae`, or `deviance`.
#' @param nfolds Integer. The number of folds to use. For leave-vth-out cross
#'   validation, every vth `y` value and its corresponding position (and weight)
#'   are placed into the same fold. The first and last observations are not
#'   assigned to any folds. This value must be at least 2. As an example, with
#'   15 data points and `nfolds = 4`, the points are assigned to folds in the
#'   following way:
#'   \deqn{
#'   0 \; 1 \; 2 \; 3 \; 4 \; 1 \; 2 \; 3 \;  4 \; 1 \; 2 \; 3 \; 4 \; 1 \; 0
#'   }{0 1 2 3 4 1 2 3 4 1 2 3 4 1 0} where 0 indicates no assignment.
#'   Therefore, the folds are not random and running `cv_trendfilter()` twice
#'   will give the same result.
#' @param ... additional parameters passet along to `trendfilter()`.
#'
#' @return An object with S3 class `"cv_trendfilter"`. Among the list components:
#' * `full_fit` An object with S3 class `"trendfilter"`, estimated with all
#'   `y` and `lambda`
#' * `cv_scores` leave-vth-out cross validation scores
#' * `cv_se` leave-vth-out cross validation standard error
#' * `lambda_min` lambda which achieved the optimal cross validation score
#' * `lambda_1se` lambda that gives the optimal cross validation score
#' within one standard error.
#' * `lambda` the value of `lambda` used in the algorithm.
#' @export
#'
#' @examples
#' x <- 1:100 / 101 * 2 * pi
#' y <- sin(x) + .2 * rnorm(100)
#' cv <- cv_trendfilter(y, x)
cv_trendfilter <- function(
    y,
    x = seq_along(y),
    weights = rep(1, n),
    k = 2L,
    error_measure = c("deviance", "mse", "mae"),
    nfolds = 5L,
    ...) {
  n <- length(y)
  min_n_train <- (n - 2) %/% nfolds * (nfolds - 1)
  assert_numeric(y, finite = TRUE)
  assert_numeric(x, finite = TRUE, len = n)
  assert_numeric(weights, lower = 0, finite = TRUE, len = n)
  # need k+1 =< n_train forall folds
  assert_integerish(k, lower = 0L, upper = min_n_train - 1L, len = 1L)
  error_measure <- arg_match(error_measure)
  assert_integerish(nfolds, lower = 2L, len = 1L, upper = n - 2L) # 2 for endpoints


  ## Run program one time to create lambda
  full_fit <- trendfilter(y, x, weights, k, ...)
  lambda <- full_fit$lambda
  middle_fold <- rep_len(1:nfolds, n - 2)
  foldid <- c(0, middle_fold, 0)
  cvall <- matrix(NA, nfolds, length(lambda))

  error_measure <- match.arg(error_measure)
  err_fun <- switch(error_measure,
    mse = function(y, m) (y - m)^2,
    mae = function(y, m) abs(y - m),
    deviance = function(y, m) (y - m)^2
  )

  for (f in 1:nfolds) {
    train_idx <- foldid != f
    test_idx <- !train_idx

    mod <- trendfilter(
      y[train_idx],
      x = x[train_idx],
      weights = weights[train_idx],
      k = k,
      lambda = lambda,
      ...
    )

    interp <- predict(mod, x[test_idx])
    score <- colMeans(err_fun(y[test_idx], interp))
    cvall[f, seq_along(score)] <- score
    if (length(lambda) != length(score)) {
      cli::cli_warn(c(
        "Estimation for the full `lambda` sequence did not occur for fold {.val {f}}",
        "because the maximum number of iterations was exhausted.",
        "i" = "You may wish to increase `maxiter` from the current {.val {control$max_iter}}."
      ))
    }
  }
  index <- apply(cvall, 2, function(x) any(is.na(x)))
  cvall <- cvall[, !index]
  lambda <- lambda[!index]


  ### Calculate CV summary
  cv_scores <- colMeans(cvall)
  cv_se <- apply(cvall, FUN = stats::sd, MARGIN = 2) / sqrt(nfolds)
  i0 <- which.min(cv_scores)

  structure(
    enlist(
      full_fit,
      cv_scores,
      cv_se,
      lambda,
      lambda_min = lambda[i0],
      lambda_1se = max(
        lambda[cv_scores <= cv_scores[i0] + cv_se[i0]],
        na.rm = TRUE
      ),
      call = match.call()
    ),
    class = "cv_trendfilter"
  )
}

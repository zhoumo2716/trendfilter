#' @method summary cv_trendfilter
#' @export
summary.cv_trendfilter <- function(object, ...) {
  rlang::check_dots_empty()

  tab <- with(object, data.frame(
    lambda = lambda,
    index = seq_along(lambda),
    cv_scores = cv_scores,
    cv_se = cv_se,
    dof = full_fit$dof[seq_along(lambda)]
  ))
  n <- nrow(tab)
  if (n > 5) {
    l1 <- which(abs(object$lambda - object$lambda_1se) < 1e-10)
    l2 <- which(abs(object$lambda - object$lambda_min) < 1e-10)
    idx <- c(1, l1, l2, n)
    tab <- tab[idx, ]
    rownames(tab) <- c("Max Lambda", "1se Lambda", "CV Minimizer", "Min Lambda")
  }

  out <- structure(
    list(call = object$call, table = tab, korder = object$full_fit$k),
    class = "summary.cv_trendfilter"
  )
  out
}

#' @method print cv_trendfilter
#' @export
print.cv_trendfilter <- function(x, ...) {
  print(summary(x, ...))
}


#' @method print summary.cv_trendfilter
#' @export
print.summary.cv_trendfilter <- function(
    x,
    digits = max(3, getOption("digits") - 3),
    ...) {
  rlang::check_dots_empty()

  lambda_warning <- NULL
  if (x$table$index[2] == 1) lambda_warning <- "smallest"
  if (nrow(x$table) > 3 & x$table$index[2] == x$table$index[4]) {
    lambda_warning <- "largest"
  }

  cat("\nCall:", deparse(x$call), fill = TRUE)
  cat("\nDegree of the estimated piecewise polynomial curve:", x$korder, "\n")
  if (!is.null(lambda_warning)) {
    cat(
      "Warning: the CV minimum occurred at the", lambda_warning,
      "lambda in the path.\n\n"
    )
  }
  cat("\nSummary of cross validation across lambda:\n")
  print(x$tab, digits = digits)
  cat("\n")
}

#' Plot Cross-Validated trendfilter
#'
#' @param x the result of `cv_trendfilter()` of class `cv_trendfilter`
#' @param which_lambda select which solutions to show.
#'   If not provided, the cross validation score will be plotted. If provided a
#'   vector of `lambda` values, the corresponding \eqn{\theta} estimates will be
#'   plotted. If provided a string, it must be either one of `cv_scores`,
#'   `lambda_min`, or `lambda_1se`.
#'
#'    * `cv_scores`: plot the cross validation score curve. Vertical lines
#'      indicate `lambda_min` and `lambda_1se` (from left to right).
#'    * `lambda_min`: plot the estimates corresponding to the lambda
#'      that minimizes the cross validation score.
#'    * `lambda_1se`: plot estimates corresponding to the largest lambda whose
#'      corresponding CV score is within 1 standard error of the
#'      minimal cross validation score.
#'    * If `NULL`, all estimated \eqn{\theta}'s are plotted.
#'
#' @param ... Not used.
#'
#' @return plot of cv scores
#' @exportS3Method
#'
#' @examples
#' x <- 1:100 / 101 * 2 * pi
#' y <- sin(x) + .2 * rnorm(100)
#' cv <- cv_trendfilter(y, x, nlambda = 20L)
#' plot(cv)
#' plot(cv, which_lambda = cv$lambda[1])
#' plot(cv, which_lambda = "lambda_min")
#' plot(cv, which_lambda = "lambda_1se")
#' plot(cv, NULL)
plot.cv_trendfilter <- function(
    x, which_lambda = c("cv_scores", "lambda_min", "lambda_1se"), ...) {
  rlang::check_dots_empty()
  plt_scores <- FALSE
  if (is.character(which_lambda)) {
    which_lambda <- arg_match(which_lambda)
    if (which_lambda == "cv_scores") {
      plt_scores <- TRUE
    } else {
      which_lambda <- x[[which_lambda]]
    }
  } else {
    assert_numeric(which_lambda,
      lower = min(x$full_fit$lambda),
      upper = max(x$full_fit$lambda), null.ok = TRUE
    )
  }

  if (plt_scores) {
    df <- with(x, data.frame(
      cv_scores = cv_scores,
      lambda = lambda,
      cv_se = cv_se,
      upper = cv_scores + cv_se,
      lower = cv_scores - cv_se
    ))
    plt <- ggplot2::ggplot(df) +
      ggplot2::geom_errorbar(ggplot2::aes(
        x = .data$lambda,
        y = .data$cv_scores,
        ymin = .data$lower,
        ymax = .data$upper,
        width = 0.1
      )) +
      ggplot2::geom_point(ggplot2::aes(x = .data$lambda, y = .data$cv_scores),
        color = "darkblue"
      ) +
      ggplot2::geom_vline(xintercept = x$lambda_min, linetype = "dotted") +
      ggplot2::geom_vline(xintercept = x$lambda_1se, linetype = "dotted") +
      ggplot2::theme_bw() +
      ggplot2::labs(x = expression(log(lambda)), y = "CV scores") +
      ggplot2::scale_x_log10()
  } else {
    plt <- plot(x$full_fit, lambda = which_lambda)
  }

  return(plt)
}




#' Predict with trendfilter at new (interior) design points
#'
#' @param object result of `cv_trendfilter()` of class `cv_trendfilter`
#' @param which_lambda select which solutions to show. If provided a
#'   vector of `lambda` values, the corresponding \eqn{\theta} estimates will be
#'   plotted. If a string, it must be either one of `lambda_min`, or `lambda_1se`.
#'
#'    * `lambda_min`: plot the estimates corresponding to the lambda
#'      that minimizes the cross validation score.
#'    * `lambda_1se`: plot estimates corresponding to the largest lambda whose
#'      corresponding CV score is within 1 standard error of the
#'      minimal cross validation score.
#'    * If NULL, all estimated \eqn{\theta}'s are plotted (at `newx`).
#' @inheritParams predict.trendfilter
#' @param ... additional arguments passed to `predict.trendfilter()`
#'
#' @return A vector or matrix of predictions.
#' @export
#' @examples
#' x <- 1:100 / 101 * 2 * pi
#' y <- sin(x) + .2 * rnorm(100)
#' cv <- cv_trendfilter(y, x, nlambda = 20L)
#' p <- predict(cv)
#' p <- predict(cv, which_lambda = cv$lambda[1])
#' p <- predict(cv, which_lambda = "lambda_1se")
#' p <- predict(cv, which_lambda = NULL)
#' plot(y)
#' matlines(p, lty = 2)
predict.cv_trendfilter <- function(object,
                                   newx = NULL,
                                   which_lambda = c("lambda_min", "lambda_1se"),
                                   ...) {
  if (is.character(which_lambda)) {
    which_lambda <- arg_match(which_lambda)
    which_lambda <- object[[which_lambda]]
  } else {
    assert_numeric(which_lambda,
      lower = min(object$full_fit$lambda),
      upper = max(object$full_fit$lambda), null.ok = TRUE
    )
  }
  predict(object$full_fit, newx, which_lambda, ...)
}

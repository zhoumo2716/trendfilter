interpolate_lambda <- function(lambda, new_lambda) {
  if (length(lambda) == 1) {
    nums <- length(new_lambda)
    left <- rep(1, nums)
    right <- left
    sfrac <- rep(1, nums)
  } else {
    k <- length(lambda)
    sfrac <- (lambda[1] - new_lambda) / (lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda) / (lambda[1] - lambda[k])
    sfrac[sfrac < min(lambda)] <- min(lambda)
    sfrac[sfrac > max(lambda)] <- max(lambda)
    coord <- stats::approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
    sfrac[left == right] <- 1
    sfrac[abs(lambda[left] - lambda[right]) < .Machine$double.eps] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}

interpolate_mat <- function(mat, lam_list, lambda) {
  m <- length(lambda)
  out <- mat[, lam_list$left, drop = FALSE] %*% diag(lam_list$frac, m, m) +
    mat[, lam_list$right, drop = FALSE] %*% diag(1 - lam_list$frac, m, m)
  drop(out)
}


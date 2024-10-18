library(Matrix)

test_that("test a single iterate of linear solvers yields same results", {
  rho <- 2

  n <- 30
  x <- 1:n / (n + 1) * 2 * pi
  weights <- rep(1, n)
  # random standard-normally distributed signals
  y0 <- c(0.006, -0.119, 0.262, 0.043, 1.114, -1.020, -0.350, -0.889, 0.476,
          0.865, 0.645, -0.631, 0.277, -0.480, 0.364, 0.599, 0.608, 1.755,
          -0.041, -0.124, -0.031, 0.195, -0.093, -0.207, 2.405, -0.026,
          -1.115, 0.474, 0.129, 0.222)
  y <- sin(x) + .2 * y0
  k <- 2
  # random standard-normally distributed signals
  mn <- c(0.485, -0.047, -0.200, 1.646, -0.594, 1.167, -0.816, 0.295, 0.784,
          1.446, 1.328, 0.852, -0.571, 1.007, -0.741, 1.114, 0.232, -1.886,
          -0.542, 2.111, -0.222, 0.507, -0.111, 0.199, 0.181, -0.363, 0.505,
          -0.868)
  Dkx <- dspline::d_mat(k, x, FALSE)
  theta <- solve(diag(weights) + rho * Matrix::crossprod(Dkx),
                 diag(weights) %*% y + rho * t(as.matrix(Dkx)) %*% mn)[,1]
  theta_sparseQR <- linear_single_solve_test(1, y, weights, x, rho, mn)
  theta_kf <- linear_single_solve_test(2, y, weights, x, rho, mn)
  expect_equal(theta_kf, theta)
  expect_equal(theta_sparseQR, theta)

  n <- 50
  x <- 1:n / (n + 1) * 2 * pi
  weights <- rep(1, n)
  y0 <- c(-0.065, 0.654, 1.397, 0.626, 1.233, 0.818, 0.586, -0.693, -0.345,
          -0.550, 0.618, -0.216, -1.355, 0.131, -0.461, 0.581, -1.709, -0.881,
          0.506, 1.749, -0.449, 1.856, -0.840, -2.980, -0.355, -1.221, 0.450,
          0.006, 0.350, 0.790, 0.405, 1.467, 0.890, 0.282, -1.038, 0.208, -0.039,
          1.166, 1.192, 0.027, 1.078, 2.363, -0.080, -0.074, -0.638, -0.539,
          -1.721, 0.646, -0.899, 0.600)
  y <- sin(x) + .2 * y0
  k <- 3
  mn <- c(0.194, -0.579, 0.240, -1.641, -0.009, 0.865, -0.362, -0.739, 0.463,
          0.854, 0.103, -0.718, -0.110, 0.367, 0.131, 0.368, -1.452, 2.719,
          -0.314, 0.846, 0.017, -0.899, -0.569, 0.116, -0.849, 0.969, -0.614,
          1.181, 1.432, 0.369, 0.427, 0.200, 0.218, -0.678, -0.513, -1.233,
          -1.161, -0.601, 1.064, -3.186, -1.311, 0.603, -0.092, -0.429, 0.997,
          -0.419, -2.362)
  Dkx <- dspline::d_mat(k, x, FALSE)
  theta <- solve(diag(weights) + rho * Matrix::crossprod(Dkx),
                 diag(weights) %*% y + rho * t(as.matrix(Dkx)) %*% mn)[,1]

  theta_sparseQR <- linear_single_solve_test(1, y, weights, x, rho, mn)
  theta_kf <- linear_single_solve_test(2, y, weights, x, rho, mn)
  expect_equal(theta_kf, theta)
  expect_equal(theta_sparseQR, theta)
})

lambda <- 5
n <- 30
x <- 1:n / (n + 1) * 2 * pi
# random Poisson signals with mean 5
weig <- c(5, 5, 1, 5, 8, 3, 2, 6, 4, 5, 4, 3, 4, 6, 5, 2, 9, 4, 1, 2, 5, 4, 3,
          5, 3, 3, 9, 4, 6, 7)
# random standard normal signals
y0 <- c(0.418, 2.954, 0.530, -0.967, 0.446, 1.458, 0.709, -1.424, 0.281, -1.166,
        1.243, -0.039, 1.134, 0.946, 0.033, 0.166, -0.761, -0.885, 1.629, -1.269,
        -0.482, -1.146, 0.468, 1.179, -1.195, 0.416, -0.957, 0.659, -0.606,
        -0.321)
y <- sin(x) + .2 * y0
test_that("test linear solvers run with no errors until convergence for single lambda", {
  k <- 2L
  theta_sparseQR <- trendfilter(
    y, x, weights = weig, lambda = lambda, k = k,
    control = trendfilter_control_list(
      admm_control = admm_control_list(k = k, linear_solver = 1L)))
  theta_kf <- trendfilter(y, x, weights = weig, k = k, lambda = lambda)
  expect_no_error(theta_kf)
  expect_no_error(theta_sparseQR)

  k <- 3L
  theta_sparseQR <- trendfilter(
    y, x, weights = weig, lambda = lambda, k = k,
    control = trendfilter_control_list(
      admm_control = admm_control_list(k = k, linear_solver = 1L)))
  theta_kf <- trendfilter(y, x, weights = weig, k = k, lambda = lambda)
  expect_no_error(theta_kf)
  expect_no_error(theta_sparseQR)
})

test_that("test linear solvers run with no errors for lambda sequences", {
  nlambda <- 10

  k <- 2L
  n <- 30
  x <- 1:n / (n + 1) * 2 * pi
  weig <- rnorm(n, 5, 1)
  set.seed(117)
  y <- sin(x) + .2 * rnorm(n)
  expect_no_error(trendfilter(y, x, weights = weig, k = k, nlambda = nlambda))
  expect_no_error(trendfilter(
    y, x, weights = weig, k = k, nlambda = nlambda,
    control = trendfilter_control_list(
      admm_control = admm_control_list(k = k, linear_solver = 1L))))

  k <- 3L
  n <- 50
  x <- 1:n / (n + 1) * 2 * pi
  weig <- rnorm(n, 5, 1)
  set.seed(121)
  y <- sin(x) + .2 * rnorm(n)
  expect_no_error(trendfilter(y, x, weights = weig, k = k, nlambda = nlambda))
  expect_no_error(trendfilter(
    y, x, weights = weig, k = k, nlambda = nlambda,
    control = trendfilter_control_list(
      admm_control = admm_control_list(k = k, linear_solver = 1L))))
})

test_that("test marginal cases", {
  k <- 2L
  lambda <- 10
  n <- 30
  x <- 1:n / (n + 1) * 2 * pi
  weig <- rep(0, n)
  y <- rep(1, n)
  expect_no_error(trendfilter(y, x, weights = weig, k = k, lambda = lambda))
  expect_no_error(trendfilter(
    y, x, weights = weig, k = k, lambda = lambda,
    control = trendfilter_control_list(
      admm_control = admm_control_list(k = k, linear_solver = 1L))))
})

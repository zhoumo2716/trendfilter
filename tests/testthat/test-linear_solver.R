library(Matrix)
test_that("test a single iterate of linear solvers yields same results", {
  rho <- 2

  n <- 30
  x <- 1:n / (n + 1) * 2 * pi
  weights <- rep(1, n)
  set.seed(839)
  y <- sin(x) + .2 * rnorm(n)

  k <- 2
  c <- rnorm(n - k)
  Dkx <- dspline::d_mat(k, x, FALSE)
  theta <- solve(diag(weights) + rho * Matrix::crossprod(Dkx),
                 diag(weights) %*% y + rho * t(as.matrix(Dkx)) %*% c)[,1]
  theta_sparseQR <- linear_single_solve_test(1, y, weights, x, rho, c)
  theta_kf <- linear_single_solve_test(2, y, weights, x, rho, c)
  expect_equal(theta_kf, theta)
  expect_equal(theta_sparseQR, theta)

  n <- 50
  x <- 1:n / (n + 1) * 2 * pi
  weights <- rep(1, n)
  set.seed(135)
  y <- sin(x) + .2 * rnorm(n)
  k <- 3
  c <- rnorm(n - k)
  Dkx <- dspline::d_mat(k, x, FALSE)
  theta <- solve(diag(weights) + rho * Matrix::crossprod(Dkx),
                 diag(weights) %*% y + rho * t(as.matrix(Dkx)) %*% c)[,1]

  theta_sparseQR <- linear_single_solve_test(1, y, weights, x, rho, c)
  theta_kf <- linear_single_solve_test(2, y, weights, x, rho, c)
  expect_equal(theta_kf, theta)
  expect_equal(theta_sparseQR, theta)
})

test_that("test linear solvers run with no errors until convergence for single lambda", {
  lambda <- 5
  n <- 30
  x <- 1:n / (n + 1) * 2 * pi
  set.seed(910)
  weig <- rpois(n, 5)
  y <- sin(x) + .2 * rnorm(n)

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

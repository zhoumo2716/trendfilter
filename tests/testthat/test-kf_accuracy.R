# common problem designs
n <- 10L
k <- 4L

# V-shape signals
y <- c((n/2):1, 1:(n/2)) + 1

# equal space and weights
x0 <- 1:n
w0 <- rep(1, n)
adj_mean0 <- numeric(n - k)

# arbitrary space and weights
set.seed(400)
x1 <- rpois(n, 5) + 1
x1 <- cumsum(x1)
x1 <- c(x1 / x1[n] * (n - 1)) + 1
w1 <- rpois(n, 5)
w1 <- w1 / sum(w1) * n
adj_mean1 <- rnorm(n - k)

# hyper-parameters
rho <- 1.0

Dkx0 <- dspline::d_mat(k, x0, FALSE)
Dkx1 <- dspline::d_mat(k, x1, FALSE)

# averaged squared error between A*x and b in the linear system A*x=b
f_mse <- function(A, x, b) {
  return(sum((A %*% x - b)^2) / length(x1))
}

test_that("LS: V-shape signals -- Scenario 1", {
  x <- x0
  Dkx <- Dkx0
  w <- w0
  mn <- adj_mean0

  A <- diag(w) + rho * Matrix::crossprod(Dkx)
  b <- diag(w) %*% y + rho * t(as.matrix(Dkx)) %*% mn
  theta <- solve(A, b)[,1]
  theta_kf <- linear_single_solve_test(2, y, w, x, rho, mn)
  expect_equal(theta_kf, theta)

  mse_kf <- f_mse(A, theta_kf, b)
  expect_true(mse_kf < 1e-15)
})
test_that("LS: V-shape signals -- Scenario 2", {
  x <- x1
  Dkx <- Dkx1
  w <- w0
  mn <- adj_mean0

  A <- diag(w) + rho * Matrix::crossprod(Dkx)
  b <- diag(w) %*% y + rho * t(as.matrix(Dkx)) %*% mn
  theta <- solve(A, b)[,1]
  theta_kf <- linear_single_solve_test(2, y, w, x, rho, mn)
  expect_equal(theta_kf, theta)

  mse_kf <- f_mse(A, theta_kf, b)
  expect_true(mse_kf < 1e-15)
})
test_that("LS: V-shape signals -- Scenario 3", {
  x <- x0
  Dkx <- Dkx0
  w <- w0
  mn <- adj_mean1

  A <- diag(w) + rho * Matrix::crossprod(Dkx)
  b <- diag(w) %*% y + rho * t(as.matrix(Dkx)) %*% mn
  theta <- solve(A, b)[,1]
  theta_kf <- linear_single_solve_test(2, y, w, x, rho, mn)
  expect_equal(theta_kf, theta)

  mse_kf <- f_mse(A, theta_kf, b)
  expect_true(mse_kf < 1e-15)
})
test_that("LS: V-shape signals -- Scenario 4", {
  x <- x0
  Dkx <- Dkx0
  w <- w1
  mn <- adj_mean0

  A <- diag(w) + rho * Matrix::crossprod(Dkx)
  b <- diag(w) %*% y + rho * t(as.matrix(Dkx)) %*% mn
  theta <- solve(A, b)[,1]
  theta_kf <- linear_single_solve_test(2, y, w, x, rho, mn)
  expect_equal(theta_kf, theta)

  mse_kf <- f_mse(A, theta_kf, b)
  expect_true(mse_kf < 1e-15)
})
test_that("LS: V-shape signals -- Scenario 5", {
  x <- x0
  Dkx <- Dkx0
  w <- w1
  mn <- adj_mean1

  A <- diag(w) + rho * Matrix::crossprod(Dkx)
  b <- diag(w) %*% y + rho * t(as.matrix(Dkx)) %*% mn
  theta <- solve(A, b)[,1]
  theta_kf <- linear_single_solve_test(2, y, w, x, rho, mn)
  expect_equal(theta_kf, theta)

  mse_kf <- f_mse(A, theta_kf, b)
  expect_true(mse_kf < 1e-15)
})
test_that("LS: V-shape signals -- Scenario 6", {
  x <- x1
  Dkx <- Dkx1
  w <- w0
  mn <- adj_mean1

  A <- diag(w) + rho * Matrix::crossprod(Dkx)
  b <- diag(w) %*% y + rho * t(as.matrix(Dkx)) %*% mn
  theta <- solve(A, b)[,1]
  theta_kf <- linear_single_solve_test(2, y, w, x, rho, mn)
  expect_equal(theta_kf, theta)

  mse_kf <- f_mse(A, theta_kf, b)
  expect_true(mse_kf < 1e-15)
})
test_that("LS: V-shape signals -- Scenario 7", {
  x <- x1
  Dkx <- Dkx1
  w <- w1
  mn <- adj_mean0

  A <- diag(w) + rho * Matrix::crossprod(Dkx)
  b <- diag(w) %*% y + rho * t(as.matrix(Dkx)) %*% mn
  theta <- solve(A, b)[,1]
  theta_kf <- linear_single_solve_test(2, y, w, x, rho, mn)
  expect_equal(theta_kf, theta)

  mse_kf <- f_mse(A, theta_kf, b)
  expect_true(mse_kf < 1e-15)
})
test_that("LS: V-shape signals -- Scenario 8", {
  x <- x1
  Dkx <- Dkx1
  w <- w1
  mn <- adj_mean1

  A <- diag(w) + rho * Matrix::crossprod(Dkx)
  b <- diag(w) %*% y + rho * t(as.matrix(Dkx)) %*% mn
  theta <- solve(A, b)[,1]
  theta_kf <- linear_single_solve_test(2, y, w, x, rho, mn)
  expect_equal(theta_kf, theta)

  mse_kf <- f_mse(A, theta_kf, b)
  expect_true(mse_kf < 1e-15)
})

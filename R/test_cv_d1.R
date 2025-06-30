segment1 <- rep(3, 10)                    # Flat segment
segment2 <- 2 + 0.01 * (1:20)^2 - 0.001*(1:20)^3
segment3 <- rep(8, 10)                     # Flat segment
y <- c(segment1, segment2, segment3)
x <- 1:length(y)

cv2d <- cv_trendfilter_D1(y, x, k = 3)
head(cv2d$scalar_min)

# segment1 <- rep(3, 10)                    # Flat segment
# segment2 <- 2 + 0.01 * (1:20)^2 - 0.001*(1:20)^3
# segment3 <- rep(8, 10)                     # Flat segment
# y <- c(segment1, segment2, segment3)
# x <- 1:length(y)
# k <- 3
# lambda1_scalar <- c(0.1, 0.02, 0.5, 0.7, 1, 1.5, 2, 10)  # scalar grid
# nfolds <- 3
# nlambda <- 10
# lambda_max <- NULL
# lambda_min_ratio <- 1e-5
#
# n <- length(y)
#   if (is.null(lambda_max)) {
#     lambda_max <- get_lambda_max(x, y, rep(1, n), k)
#   }
#
# lambdak <- numeric(nlambda)
# lambdak <- get_lambda_seq_r(lambdak, lambda_max,lambda_min = -1,lambda_min_ratio, nlambda)
#
#
# foldid <- c(0, rep_len(1:nfolds, n - 2), 0)   # endpoints never held out
# mse    <- function(a, b) (a - b)^2
#
# cv_score <- matrix(NA, nrow = nlambda, ncol = length(lambda1_scalar))
# best_scalar <- double(nlambda)
#
#   ## --- 3. double loop: lambdak  ×  scalar ------------------------------
#   for (m in seq_along(lambdak)) {
#     lamk     <- lambdak[m]
#     fold_scores <- matrix(NA, nfolds, length(lambda1_scalar))
#
#     for (s in seq_along(lambda1_scalar)) {
#       rho <- lambda1_scalar[s]
#
#       for (f in 1:nfolds) {
#         tr <- foldid != f
#         te <- !tr
#
#         fit_tr <- trendfilter_D1(
#          k        = k,
#         lambda1_scalar = rho,
#          lambdak  = lamk,
#          y        = y[tr],
#          x        = x[tr]
#         )
#
#         y_hat <- predict.trendfilter_D1(fit_tr, x[te])
#         fold_scores[f, s] <- mean(mse(y[te], y_hat))
#       }
#    }
#
#     ## average along columns (over folds) -> one score per scalar
#     cv_score[m, ]       <- colMeans(fold_scores, na.rm = TRUE)
#     best_scalar_index <- which.min(cv_score[m, ])
#     best_scalar[m] <- lambda1_scalar[best_scalar_index]
#   }
#
#   ## --- 4. list of (lambdak, best scalar) pairs -------------------------
#   scalar_min <- data.frame(
#     lambdak = lambdak,
#     best_lambda1_scalar = best_scalar
#   )
#
# #
# # #print(rho)
# # #print(lamk)
#

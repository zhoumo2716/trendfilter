#devtools::load_all()
segment1 <- rep(3, 10)                    # Flat segment
segment2 <- 2 + 0.01 * (1:20)^2 - 0.001*(1:20)^3
segment3 <- rep(8, 10)                     # Flat segment
y <- c(segment1, segment2, segment3)
x <- 1:length(y)
n <- length(y)

### Plot 1: Using trendfilter_D1 with lambda1 = 1, lambdak = 1, and k = 3
out1 <- trendfilter_D1(k = 3, lambda1_scalar = 0.1, lambdak = 7.599391e-03, y = y, x = x)
plot(x, y, main = "trendfilter_D1: k=3, lambda1_scalar = 0.1",
     xlab = "x", ylab = "y", pch = 19, col = "darkgray")
points(x, out1$theta, col = "blue", lwd = 2)

### Plot 2: Using trendfilter_D1 with lambda1 = 3, lambdak = 1, and k = 3
out2 <- trendfilter_D1(k = 3, lambda1_scalar = 0.5, lambdak = 7.599391e-03, y = y, x = x)
plot(x, y, main = "trendfilter_D1: k=3, lambda1_scalar = 0.5",
     xlab = "x", ylab = "y", pch = 19, col = "darkgray")
points(x, out2$theta, col = "blue", lwd = 2)

### Plot 3: Using trendfilter_D1 with lambda1 = 1, lambdak = 1, and k = 3
out3 <- trendfilter_D1(k = 3, lambda1_scalar = 1, lambdak = 7.599391e-03, y = y, x = x)
plot(x, y, main = "trendfilter_D1: k=3, lambda1_scalar = 1",
     xlab = "x", ylab = "y", pch = 19, col = "darkgray")
points(x, out3$theta, col = "blue", lwd = 2)

### Plot 4: Using trendfilter_D1 with lambda1 = 3, lambdak = 1, and k = 3
out4 <- trendfilter_D1(k = 3, lambda1_scalar = 2, lambdak = 7.599391e-03, y = y, x = x)
plot(x, y, main = "trendfilter_D1: k=3, lambda1_scalar = 2",
     xlab = "x", ylab = "y", pch = 19, col = "darkgray")
points(x, out4$theta, col = "blue", lwd = 2)

### Plot: Regular Trendfilter with k = 3
# Here we assume that trendfilter() is another function (e.g., your baseline trend filter without the extra D1 penalty)
# that returns a list with an element "theta".
# For demonstration, let's assume lambda = 1 for the regular trend filter.
# out_t <- trendfilter(y,x, k=3, lambda = 0.00299044493944087)
out_t <- trendfilter(y,x, k=3, nlambda = 5)
plot(x, y, main = "Regular trendfilter: k = 3",
     xlab = "x", ylab = "y", pch = 19, col = "darkgray")
points(x, out_t$theta, col = "red", lwd = 2)

### Plot: Regular Trendfilter with k = 1
out_t1 <- trendfilter(y,x, k=1, nlambda = 5)
plot(x, y, main = "Regular trendfilter: k = 1",
     xlab = "x", ylab = "y", pch = 19, col = "darkgray")
points(x, out_t1$theta, col = "red", lwd = 2)



#######################################################
## ------------------------------------------------------------
## Make a signal with several waves and sharp jumps between them
## ------------------------------------------------------------

# Segment 1: flat
segment1 <- rep(2.5, 10)

# Segment 2: small wave
t2 <- seq(0, 2*pi, length.out = 18)
segment2 <- 3.0 + 0.8 * sin(t2)

# Sharp jump up
segment3 <- rep(7.5, 8)

# Segment 4: larger wave
t4 <- seq(0, 3*pi, length.out = 20)
segment4 <- 7.2 + 1.0 * sin(t4) + 0.3 * cos(2 * t4)

# Sharp jump down
segment5 <- rep(1.8, 8)

# Segment 6: another wave
t6 <- seq(0, 2.5*pi, length.out = 18)
segment6 <- 2.0 + 0.7 * sin(t6) - 0.2 * cos(1.5 * t6)

# Sharp jump up again
segment7 <- rep(6.8, 8)

# Segment 8: final wave
t8 <- seq(0, 2*pi, length.out = 18)
segment8 <- 6.5 + 0.9 * cos(t8)

# Combine
y <- c(segment1, segment2, segment3, segment4,
       segment5, segment6, segment7, segment8)

x <- seq_along(y)
n <- length(y)

## If you want noisy observations instead, uncomment:
# set.seed(123)
y <- y + rnorm(n, sd = 0.15)


## ------------------------------------------------------------
## Tuning values
## ------------------------------------------------------------
library(ggplot2)
output_dir <- "/Users/xli/Documents/UBC stats PhD RA/proposal-template-main/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

lambda_tf_k3 <- c(0.005, 0.05, 0.01, 0.1, 1)
lambda_tf_k1 <- c(0.005, 0.01, 0.05, 0.1, 0.5)
cols <- c("blue", "red", "forestgreen", "purple", "orange")

y_df <- data.frame(x = x, y = y)

lambda_pairs <- data.frame(
  lambda1 = c(0.001, 0.1, 1, 0.1, 0.1),
  lambdak = c(0.005, 0.005, 0.005, 0.05, 0.1)
)

# Labels for the five penalty combinations
labels <- paste0(
  "\u03bb1 = ", lambda_pairs$lambda1,
  ", \u03bbk = ", lambda_pairs$lambdak
)

# Generate fits
fit_df_mixed <- do.call(rbind, lapply(seq_len(nrow(lambda_pairs)), function(i) {
  fit <- trendfilter_D1(
    k = 3,
    lambda1_scalar = lambda_pairs$lambda1[i],
    lambdak = lambda_pairs$lambdak[i],
    y = y,
    x = x
  )

  data.frame(
    x = x,
    theta = as.numeric(fit$theta),
    tuning = labels[i]
  )
}))

fit_df_mixed$tuning <- factor(fit_df_mixed$tuning, levels = labels)

# Plot with legend outside on the right
p_mixed <- ggplot() +
  geom_point(
    data = data.frame(x = x, y = y),
    aes(x = x, y = y),
    color = "grey70", size = 2
  ) +
  geom_point(
    data = fit_df_mixed,
    aes(x = x, y = theta, color = tuning),
    size = 1,
    alpha = 0.5
  ) +
  scale_color_manual(
    values = setNames(cols, labels),
    breaks = labels
  ) +
  labs(
    title = "Mixed Penalty Trend Filtering (k = 3)",
    x = "x", y = "y",
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

print(p_mixed)

ggsave(
  file.path(output_dir, "mixed_penalty_combinations.pdf"),
  plot = p_mixed,
  width = 11, height = 6
)


## ------------------------------------------------------------
## Regular trend filtering, k = 3 and 1
## ------------------------------------------------------------
plot_tf <- function(k, lambda_values) {

  # Fit models for each lambda
  fit_df <- do.call(rbind, lapply(lambda_values, function(lam) {
    fit <- trendfilter(y = y, x = x, k = k, lambda = lam)

    data.frame(
      x = x,
      theta = as.numeric(fit$theta),
      lambda = paste0("lambda = ", lam)
    )
  }))

  # Observed signal
  y_df <- data.frame(x = x, y = y)

  # Plot
  ggplot() +
    geom_point(
      data = y_df,
      aes(x = x, y = y),
      color = "grey70", size = 2
    ) +
    geom_point(
      data = fit_df,
      aes(x = x, y = theta, color = lambda),
      size = 1,
      alpha = 0.5
    ) +
    scale_color_manual(
      values = setNames(
        cols[seq_along(lambda_values)],
        paste0("lambda = ", lambda_values)
      )
    ) +
    labs(
      title = paste0("Regular Trend Filtering (k = ", k, ")"),
      x = "x",
      y = "y",
      color = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )
}

# Regular trend filtering: vary lambda directly
# k = 3
p_k3 <- plot_tf(k = 3, lambda_values = lambda_tf_k3)

# k = 1
p_k1 <- plot_tf(k = 1, lambda_values = lambda_tf_k1)

# Display
print(p_k3)
print(p_k1)

# Save as PDF
ggsave(
  filename = file.path(output_dir, "regular_trendfilter_k3.pdf"),
  plot = p_k3,
  width = 11, height = 6, units = "in"
)

ggsave(
  filename = file.path(output_dir, "regular_trendfilter_k1.pdf"),
  plot = p_k1,
  width = 11, height = 6, units = "in"
)


###############################################################################################
# start <- Sys.time()       # start timer
# ##--------------------------------------------------------------------------------
# for (i in 1:10000) {
# out_r <- trendfilter(y,x, k=3, lambda = 0.00299044493944087)
# out_d1 <- trendfilter_D1(k = 3, lambda1_scalar = 1, lambdak = 0.2, y = y, x = x)
#
# }
# ##--------------------------------------------------------------------------------
# end   <- Sys.time()       # stop timer
# run   <- end - start
#
# plot(x, y, main = "Trendfilter Comparison",
#      xlab = "x", ylab = "y", pch = 19, col = "darkgray")
# points(x, out_r$theta, col = "red", lwd = 2)
# points(x, out_d1$theta, col = "blue", lwd = 2)
#
#
# mse_global <- mean( (y - out_d1$theta)^2 )

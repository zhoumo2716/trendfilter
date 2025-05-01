segment1 <- rep(3, 10)                    # Flat segment
segment2 <- 2 + 0.01 * (1:20)^2 - 0.001*(1:20)^3
segment3 <- rep(8, 10)                     # Flat segment
y <- c(segment1, segment2, segment3)
x <- 1:length(y)
n <- length(y)

### Plot 1: Using trendfilter_D1 with lambda1 = 1, lambdak = 1, and k = 3
out1 <- trendfilter_D1(n = n, k = 3, lambda1 = 1, lambdak = 1, y = y, x = x)
plot(x, y, main = "trendfilter_D1: k=3, lambda1=1, lambdak=1",
     xlab = "x", ylab = "y", pch = 19, col = "darkgray")
points(x, out1$theta, col = "blue", lwd = 2)

### Plot 2: Using trendfilter_D1 with lambda1 = 3, lambdak = 1, and k = 3
out2 <- trendfilter_D1(n = n, k = 3, lambda1 = 3, lambdak = 1, y = y, x = x)
plot(x, y, main = "trendfilter_D1: k=3, lambda1=3, lambdak=1",
     xlab = "x", ylab = "y", pch = 19, col = "darkgray")
points(x, out2$theta, col = "blue", lwd = 2)

### Plot 3: Regular Trendfilter with k = 3
# Here we assume that trendfilter() is another function (e.g., your baseline trend filter without the extra D1 penalty)
# that returns a list with an element "theta".
# For demonstration, let's assume lambda = 1 for the regular trend filter.
out3 <- trendfilter(y,x, k=3, nlambda = 1)
plot(x, y, main = "Regular trendfilter: k = 3",
     xlab = "x", ylab = "y", pch = 19, col = "darkgray")
points(x, out3$theta, col = "blue", lwd = 2)

### Plot 4: Regular Trendfilter with k = 1
out4 <- trendfilter(y,x, k=1, nlambda = 1)
plot(x, y, main = "Regular trendfilter: k = 1",
     xlab = "x", ylab = "y", pch = 19, col = "darkgray")
points(x, out4$theta, col = "blue", lwd = 2)

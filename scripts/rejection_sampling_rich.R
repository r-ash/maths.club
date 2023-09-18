## 1. Importance sampling
p1 <- function(x) {
  exp(-(x^2)) * 0.5 + exp(-((x - 1)^2)) * 0.3 + exp(-((x + 2)^2)) * 0.2
}

## Compute the mean directly (can't normally do this, note we compute
## the normalizing constant in order to do this).
z <- integrate(function(x) p1(x), -Inf, Inf)$value
integrate(function(x) x * p1(x) / z, -Inf, Inf)

## Use importance sampling, with normal distribution:
f <- identity
q1 <- dnorm
r1 <- rnorm

n <- 10000
r <- r1(n)
w <- p1(r) / q1(r)
plot(cumsum(f(r) * w) / cumsum(w), type = "l")
abline(h = -0.1, lty = 3, col = "red")

## And again with Cauchy - sometimes better, but still really high variance
q2 <- dcauchy
r2 <- rcauchy

n <- 1000
r <- r2(n)
w <- p1(r) / q2(r)
plot(cumsum(f(r) * w) / cumsum(w), type = "l")
abline(h = -0.1, lty = 3, col = "red")

## 2. Rejection sampling, using scaled normal
curve(p1, -5, 5, log = "y")
curve(dnorm(x, 0, 1.5) * 3.5, add = TRUE, col = "red", lty = 2)

q1 <- function(x) dnorm(x, 0, 1.5) * 3.5
r1 <- function(n) rnorm(n, 0, 1.5)

rr <- r1(100000)
uu <- runif(length(rr))
accept <- uu < p1(rr) / q1(rr)
plot(rr, uu * q1(rr), col = ifelse(accept, "blue", "red"), pch = ".")
curve(p1, add = TRUE)
curve(q1, add = TRUE)


plot(density(rr[accept]))

var(rr[accept])
summary(rr[accept])

## 3. MH MCMC, in non-log space (typically we'll want this in log space)
mcmc <- function(f, x, n, proposal) {
  p <- f(x)
  xx <- matrix(NA, n + 1, length(x))
  xx[1, ] <- x
  for (i in seq_len(n)) {
    x_new <- proposal(x)
    p_new <- f(x_new)
    if (runif(1) < p_new / p) {
      x <- x_new
      p <- p_new
    }
    xx[i + 1, ] <- x
  }
  xx
}

x <- mcmc(p1, 0, 1000000, function(x) x + runif(1) - 0.5)
plot(density(x))
curve(p1(x) / z, add = TRUE, col = "red")

plot(drop(mcmc(p1, 0, 10000, function(x) x + runif(1, -2, 2))), pch = ".")

x <- mcmc(p1, 0, 100000, function(x) x + runif(1, -5, 5))
plot(density(x))
curve(p1(x) / z, add = TRUE, col = "red")

## do it properly, in log space

mcmc2 <- function(f, x, n, proposal) {
  p <- f(x)
  xx <- matrix(NA, n + 1, length(x))
  xx[1, ] <- x
  for (i in seq_len(n)) {
    x_new <- proposal(x)
    p_new <- f(x_new)
    if (runif(1) < exp(p_new - p)) {
      x <- x_new
      p <- p_new
    }
    xx[i + 1, ] <- x
  }
  xx
}

p2 <- function(x) log(p1(x))

x <- mcmc2(p2, 0, 10000, function(x) x + runif(1) - 0.5)
plot(density(x))
curve(p1(x) / z, add = TRUE, col = "red")

proposal <- function(w) function(x) x + runif(1, -w, w)
plot(drop(mcmc2(p2, 0, 10000, proposal(100))), pch = ".")

x <- mcmc(p1, 0, 100000, function(x) x + runif(1, -5, 5))
plot(density(x))
curve(p1(x) / z, add = TRUE, col = "red")

## A really simple example, sample from mvn:
mvnorm <- function(mean, sigma) {
  sigma_inv <- solve(sigma)
  function(x, log = TRUE) {
    distval <- mahalanobis(x, center = mean, cov = sigma_inv, TRUE)
    logdet <- as.numeric(determinant.matrix(sigma, TRUE)$modulus)
    ret <- -(length(x) * log(2 * pi) + logdet + distval)/2
    if (log) ret else exp(ret)
  }
}

f2 <- mvnorm(c(0, 0), cbind(c(1, 0.7), c(0.7, 2)))

proposal <- function(sd) {
  function(x) {
    rnorm(length(x), x, sd)
  }
}

x <- mcmc2(f2, c(-1, -2), 100000, proposal(0.5))
plot(x, pch = ".", col = "#00000055")
plot(density(x[, 1]))
plot(density(x[, 2]))

## try something weirder:
g <- function(a, b) {
  function(theta) {
    x <- theta[[1]]
    y <- theta[[2]]
    -log(a + (1 - x)^2 + b * (y - x^2)^2)
  }
}

x <- mcmc2(g(1, 5), c(-1, -2), 100000, proposal(0.1))
plot(x, pch = ".", col = "#00000055")

x <- mcmc2(g(1, 100), c(-1, -2), 100000, proposal(0.1))
plot(x, pch = ".", col = "#00000055")

library(ggplot2)
data <- data.frame(x = x[, 1], y = x[, 2])
ggplot(data, aes(x=x, y=y) ) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

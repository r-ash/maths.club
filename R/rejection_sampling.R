f <- function(x) {
  exp(-(x^2)) * 0.5 + exp(-((x - 1)^2)) * 0.3 + exp(-((x + 2)^2)) * 0.2
}

covering <- function(x, scale = 1) {
  scale * dnorm(x, sd = 1.5)
}

# curve(f, -10, 10, n = 1000, log = "y")
# curve(dnorm(x, sd = 4, log = TRUE) * 2, add = TRUE, col = "red")

plot_f <- function(f, x_min, x_max) {
  x <- seq(x_min, x_max, 0.1)
  plot(x, f(x), type = "l")
}

plot_f_g <- function(f, g, x_min, x_max) {
  x <- seq(x_min, x_max, 0.1)
  data <- data.frame(x = x)
  data$f <- f(x)
  data$g <- g(x) * 3
  colours <- wesanderson::wes_palette("Zissou1", 3)
  ggplot2::ggplot(data) +
    ggplot2::geom_line(ggplot2::aes(x = x, y = f), color = colours[1]) +
    ggplot2::geom_line(ggplot2::aes(x = x, y = g), color = colours[3])
}

rejection_sampling <- function(p, q, c, x_min, x_max, n_samples) {
  res <- rep(NA_real_, n_samples)
  for (i in seq(1, n_samples)) {
    x <- rnorm(1, sd = 1.5)
    sample_point <- runif(1, min = 0, max = c * q(x))
    p_val <- p(x)
    if (sample_point < p_val) {
      res[i] <- x
    }
  }
  data.frame(sample = res[!is.na(res)])
}

plot_rejection_sampling <- function(samples, p, x_min, x_max) {
  x <- seq(x_min, x_max, 0.1)
  p_val <- p(x)
  normalise <- integrate(p, x_min, x_max)
  p_val <- p(x)/normalise$value
  actual <- data.frame(x = x, p = p_val)
  ggplot2::ggplot(samples, ggplot2::aes(sample, ggplot2::after_stat(density))) +
    ggplot2::geom_histogram(bins = 200) +
    ggplot2::geom_line(data = actual, ggplot2::aes(x = x, y = p) )
}

# x <- rejection_sampling(f, covering, 3, -5, 5, 10000)
# plot_rejection_sampling(x, f, -5, 5)


rejection_sampling_all_points <- function(p, q, c, x_min, x_max, n_samples) {
  x <- rnorm(n_samples, sd = 1.5)
  sample_point <- runif(n_samples, min = 0, max = c * q(x))
  p_val <- p(x)
  data.frame(sample = x, value = sample_point, accepted = sample_point < p_val)
}


plot_animated <- function(samples) {
  # x <- seq(x_min, x_max, 0.1)
  # p_val <- p(x)
  # actual <- data.frame(x = x, p = p_val)

  samples$states <- seq(1, nrow(samples))
  ggplot2::ggplot(samples, ggplot2::aes(x = sample, y = value, color = ifelse(accepted, "green", "red"))) +
    ggplot2::geom_point(size = 0.5)

  # anim <- ggplot2::ggplot(x, aes(x, y)) +
  #   geom_point() +
  #   transition_states(z, state_length = 0) +
  #   shadow_mark() +
  #   labs(title = "Point num: {closest_state}"
  #
  # gganimate::animate(anim)
}

# x <- rejection_sampling_all_points(f, covering, 3, -5, 5, 1000)
# plot_animated(x)

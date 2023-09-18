f <- function(x) {
  exp(-(x^2)) * 0.5 + exp(-((x - 1)^2)) * 0.3 + exp(-((x + 2)^2)) * 0.2
}

# plot_f(f, -5, 5)

metropolis <- function(f, initial, steps, propsal_sd) {
  current_x <- initial
  history <- rep(NA_real_, steps + 1)
  history[1] <- initial
  for (i in seq(2, steps + 1)) {
    proposal <- rnorm(1, mean = current_x, sd = propsal_sd)
    acceptance_ratio <- f(proposal) / f(current_x)
    if (runif(1) <= acceptance_ratio) {
      current_x <- proposal
    }
    history[i] <- current_x
  }
  data.frame(sample = history)
}

plot_samples <- function(samples, p, x_min, x_max) {
  x <- seq(x_min, x_max, 0.1)
  normalise <- integrate(p, x_min, x_max)
  p_val <- p(x)/normalise$value
  actual <- data.frame(x = x, p = p_val)
  ggplot2::ggplot(samples, ggplot2::aes(sample, ggplot2::after_stat(density))) +
    ggplot2::geom_histogram(bins = 200) +
    ggplot2::geom_line(data = actual, ggplot2::aes(x = x, y = p, color = "#3B9AB2")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

plot_line <- function(samples) {
  ggplot2::ggplot(samples, ggplot2::aes(y = sample, x = as.numeric(row.names(samples)))) +
    ggplot2::geom_line()
}

plot_both <- function(samples, p, x_min, x_max) {
  x <- seq(x_min, x_max, 0.1)
  normalise <- integrate(p, x_min, x_max)
  p_val <- p(x)/normalise$value
  actual <- data.frame(x = x, p = p_val)

  s <- ggplot2::ggplot(samples, ggplot2::aes(sample, ggplot2::after_stat(density))) +
    ggplot2::geom_histogram(bins = 200) +
    ggplot2::geom_line(data = actual, ggplot2::aes(x = x, y = p, color = "#3B9AB2")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none",
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    ggplot2::coord_flip()

  line <- ggplot2::ggplot(samples, ggplot2::aes(y = sample, x = as.numeric(row.names(samples)))) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ylim(x_min, x_max) +
    ggplot2::xlab("Sample number")

  gridExtra::grid.arrange(line, s, nrow = 1)
}

# points <- metropolis(f, 0, 1000, 0.3)
# plot_samples(points, f, -5, 5)
# plot_line(points)
# plot_both(points, f, -5, 5)

f_3d <- function(x, y) {
  sin(10 * (x^2 + y^2)) / 10
}

plot_3d <- function(f, x_min, x_max, y_min, y_max) {
  x <- seq(x_min, x_max, 0.05)
  y <- seq(y_min, y_max, 0.05)
  z <- outer(x, y, f)
  fig <- plotly::plot_ly(x = x, y = y, z = z) %>% plotly::add_surface()
  fig
}

# plot_3d(f_3d, -1, 1, -1, 1)

fig <- plotly::plot_ly(z = ~volcano)
fig <- fig %>% plotly::add_surface()

plot_3d <- function(f, x_min, x_max, y_min, y_max) {
  x <- seq(x_min, x_max, 0.05)
  y <- seq(y_min, y_max, 0.05)
  z <- outer(x, y, f)
  fig <- plotly::plot_ly(x = x, y = y, z = z) %>% plotly::add_surface()
  fig
}

f_3d <- function(x, y) {
  ## Himmelblau's function
  z <- -(x^2 + y - 11)^2 - (x + y^2 -7)^2
  update <- x < -5 | x > 5 | y < -5 | y > 5
  z[update] <- NA_real_
  z
}

shekel <- function(x, y) {

}

plot_3d(f_3d, -6, 5, -5, 5)

metropolis_3d <- function(f, initial_x, initial_y, steps, propsal_sd) {
  current_x <- initial_x
  current_y <- initial_y
  history_x <- rep(NA_real_, steps + 1)
  history_x[1] <- initial_x
  history_y <- rep(NA_real_, steps + 1)
  history_y[1] <- initial_y
  for (i in seq(2, steps + 1)) {
    proposal_x <- rnorm(1, mean = current_x, sd = propsal_sd)
    proposal_y <- rnorm(1, mean = current_y, sd = propsal_sd)
    message(f(proposal_x, proposal_y))
    acceptance_ratio <- f(proposal_x, proposal_y) / f(current_x, current_y)
    if (!is.na(acceptance_ratio) && runif(1) <= acceptance_ratio) {
      current_x <- proposal_x
      current_y <- proposal_y
    }
    history_x[i] <- current_x
    history_y[i] <- current_y
  }
  data.frame(sample_x = history_x, sample_y = history_y)
}

data <- metropolis_3d(f_3d, 0, 0, 100000, 0.2)

plot_3d_histogram <- function(data) {
  ggplot2::ggplot(data = data, ggplot2::aes(x = sample_x, y = sample_y)) +
    ggplot2::geom_bin2d(bins = 100) +
    ggplot2::scale_fill_continuous(type = "viridis") +
    ggplot2::theme_minimal()
}

plot_3d_histogram(data)

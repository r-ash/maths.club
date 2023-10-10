## This example is adapted from Lilith's continuous time vignette
## https://mrc-ide.github.io/mcstate/articles/continuous.html
##
## To run it, you will need any version of odin, and a C compiler.
##
## We'll avoid the high level support in odin.dust and mcstate and
## instead write some wrappers ourselves.

## The data to be fitted to:
data <- data.frame(
  day = c(42, 49, 56, 63, 70, 77, 84, 91, 98),
  n_tested = c(107, 105, 96, 92, 93, 97, 99, 99, 109),
  n_positive = c(6, 14, 25, 31, 14, 15, 11, 6, 3))


par(bty = "n", mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0))
test_cols <- c(n_tested = "darkorange2", n_positive = "green4")
plot(data$day, data$n_tested, type = "b",
     xlab = "Day", ylab = "Number of individuals",
     ylim = c(0, max(data$n_tested)), col = test_cols[1], pch = 20)
lines(data$day, data$n_positive, type = "b", col = test_cols[2], pch = 20)
legend("left", names(test_cols), fill = test_cols, bty = "n")

## Here is the generating model; this defines a set of ODEs for S, I
## and R that run forward in time.  We can write this out without odin
## pretty easily too to remove some of the mystery, though it will be
## a bit slower.
sir <- odin::odin({
  initial(S) <- N - I_init
  initial(I) <- I_init
  initial(R) <- 0

  deriv(S) <- -beta * S * I / N
  deriv(I) <- (beta * S * I / N) - gamma * I
  deriv(R) <- gamma * I

  beta <- user(0.4)
  gamma <- user(0.2)
  N <- user(1e5)
  I_init <- user(1)

  output(prevalence) <- I / N
})

## Remember Bayes theorem
## P(theta|x) proportional to P(theta)*P(x|theta)
## Where theta are our parameters and x is the data
## We want to know the probability of theta given the observed data
## Probability of theta is given by our prior belief
## Probability of x|theta data given the parameters is given by the log
## likelihood from our model.


## We need a likelihood for the above model.
## Recall the binomial distribution, this is for parameters n and p
## it gives the probability of the number of successes in a sequence of n
## independent experiments. Where each experiment has a successor failure
## outcome with probability of success p
## So for our example, each experiment is a person tested and p is the modelled
## prevalence

## So the likelihood in Lilith's model is defined per time point by the data
## for number of positive people (`data$n_positive`) being *biniomially
## distributed* with parameters n as the number of people tested
## (`data$n_tested`) and the modelled prevalance (`I / N`).  We then
## take the product of these probabilities over all time points. This
## is much easier to do by getting the log-probability from the
## binomial distribution and taking sum, using the fact that
##
##   log(a * b) = log(a) + log(b)
##
## So for a single value of theta (which is beta and gamma in our example)
## we can work out the log likelihood
pars <- list(beta = 0.4, gamma = 0.2)
mod <- sir$new(user = pars)
y <- mod$run(c(0, data$day))[-1, ]
sum(dbinom(data$n_positive, data$n_tested, y[, "prevalence"], log = TRUE))

## We can then wrap that up as a little function. In order to work
## with mcmc we need an unstructured vector of parameters to consider,
## so let's write a function that takes some parameter vector theta
## which contains beta as the first element and gamma as the second.
##
## The conditional statement at the top of the function means that if
## we propose negative values of parmeters to unconditionally reject
## them.
loglikelihood <- function(theta) {
  if (any(theta <= 0)) {
    return(-Inf)
  }
  ## Note that this is finding 'mod' and 'data' from the global
  ## environment, which is gross but keeps things going here. It's
  ## easy enough to tidy that away if you want to...
  mod$set_user(beta = theta[[1]], gamma = theta[[2]])
  y <- mod$run(c(0, data$day))[-1, ]
  if (any(y[, "prevalence"] < 0)) {
    return(-Inf)
  }
  sum <- sum(dbinom(data$n_positive, data$n_tested, y[, "prevalence"], log = TRUE))
  if (is.nan(sum)) {
    sum = -Inf
  }
  sum
}

loglikelihood(c(0.4, 0.2))

## We have loglikelhood, going back to Bayes
## ## P(theta|x) proportional to P(theta)*P(x|theta) we need the prior P(theta)
## The modellers use intuition to define some priors for these
## parameters (if you are really
## stuck, consider a uniform prior on 0 to some upper bound, say 0.9,
## or an exponential prior with a mean of about 0.3). Your prior
## function should take the same theta vector and return the log prior
## density for the combination (i.e., the sum of log prior
## probabilities of each parameter).
##
## For example, a prior with a mean of 1/3 for *both* parameters might
## be defined as:
logprior <- function(theta) {
  sum(dexp(theta, 3, log = TRUE))
}

plot(dexp, 0, 5)

## Recall when we move in MCMC we move based on a proposal distribution
## You can do a multivariate
## normal proposal using the mvtnorm package, but independent normal
## proposals for each parameter is a reasonable starting point. Write
## a function that takes theta and returns a pair of parameters. For
## example, a really bad (but valid) proposal function which updates
## the parmeters by independently moving them by +/- 0.1 would look
## like:
proposal <- function(theta) {
  theta + runif(length(theta), -0.1, 0.1)
}


## The log posterior is then logprior(theta) + loglikelihood(theta) and we
## have all the ingredients to write an MCMC

log_posterior <- function(theta) {
  logprior(theta) + loglikelihood(theta)
}

output <- metropolis2(log_posterior, c(beta = 0.4, gamma = 0.2), 10000, proposal)
plot_3d_histogram(output)

## In the vignette, Lilith defines the model using the data comparison
## DSL which is part of odin.dust; this gives the same answer but with
## more magic.
sir2 <- odin.dust::odin_dust({
  initial(S) <- N - I_init
  initial(I) <- I_init
  initial(R) <- 0

  deriv(S) <- -beta * S * I / N
  deriv(I) <- beta * S * I / N - gamma * I
  deriv(R) <- gamma * I

  beta <- user(0.4)
  gamma <- user(0.2)
  N <- user(1e5)
  I_init <- user(1)

  ## likelihood
  n_positive <- data()
  n_tested <- data()
  compare(n_positive) ~ binomial(n_tested, I / N)
})
mod2 <- sir2$new(list(), 0, 1)
mod2$set_data(dust::dust_data(data, "day"))
mod2$filter()$log_likelihood

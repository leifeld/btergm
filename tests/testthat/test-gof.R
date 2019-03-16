context("test goodness-of-fit functions")

library("btergm")
library("network")

set.seed(12345)
networks <- list()
for(i in 1:10){            # create 10 random networks with 10 actors
  mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
  diag(mat) <- 0           # loops are excluded
  nw <- network(mat)       # create network object
  networks[[i]] <- nw      # add network to the list
}

covariates <- list()
for (i in 1:10) {          # create 10 matrices as covariate
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  covariates[[i]] <- mat   # add matrix to the list
}

fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), R = 100, usefastglm = TRUE, verbose = FALSE)

test_that("basic GOF functionality works", {
  g <- gof(fit, nsim = 2, MCMC.burnin = 1000, MCMC.interval = 500, verbose = FALSE)
  expect_equal(length(g), 7)
  expect_equal(class(g), "gof")
})
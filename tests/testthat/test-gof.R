context("test goodness-of-fit functions")

library("btergm")
library("network")

set.seed(12345)
networks <- list()
for(i in 1:10) {            # create 10 random networks with 10 actors
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

test_that("edgeprob works with ergm, btergm, and mtergm object with curved terms", {
  ep <- edgeprob(fit)
  expect_s3_class(ep, "data.frame")
  expect_equivalent(dim(ep), c(900, 10))
  fit <- btergm(networks ~ edges + gwidegree(1.2, fixed = FALSE) + edgecov(covariates), R = 100, usefastglm = TRUE, verbose = FALSE)
  ep <- edgeprob(fit)
  expect_s3_class(ep, "data.frame")
  expect_equivalent(dim(ep), c(900, 10))
  
  skip_on_cran()
  sink("/dev/null")
  fit2 <- mtergm(networks ~ edges + gwidegree(1.2, fixed = FALSE) + edgecov(covariates))
  sink()
  ep2 <- edgeprob(fit2)
  expect_s3_class(ep2, "data.frame")
  expect_equivalent(dim(ep2), c(900, 10))
  
  skip_if_not_installed("ergm", minimum_version = "3.10.4")
  sink("/dev/null")
  require("ergm")
  nw1 <- networks[[1]]
  cov1 <- covariates[[1]]
  fit3 <- ergm(nw1 ~ edges + gwidegree(1.2, fixed = FALSE) + edgecov(cov1))
  sink()
  ep3 <- edgeprob(fit3)
  expect_s3_class(ep3, "data.frame")
  expect_equivalent(dim(ep3), c(90, 10))
})
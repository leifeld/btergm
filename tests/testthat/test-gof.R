context("test goodness-of-fit functions")


test_that("basic GOF functionality works", {
  set.seed(12345)
  networks <- list()
  for(i in 1:10) {                # create 10 random networks with 10 actors
    mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
    diag(mat) <- 0               # loops are excluded
    nw <- network::network(mat)  # create network object
    networks[[i]] <- nw          # add network to the list
  }
  
  covariates <- list()
  for (i in 1:10) {              # create 10 matrices as covariate
    mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
    covariates[[i]] <- mat       # add matrix to the list
  }
  
  fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), R = 100, verbose = FALSE)
  
  sink(nullfile())
  g <- gof(fit, nsim = 2, MCMC.burnin = 1000, MCMC.interval = 500, verbose = FALSE)
  sink()
  expect_equal(length(g), 7)
  expect_equal(class(g), "gof")
  path <- tempfile(fileext = ".png")
  png(path)
  plot(g)
  dev.off()
  fs <- file.size(path)
  fs <- fs < 100000 && fs > 1000
  expect_true(fs, 5)
})

test_that("edgeprob works with ergm, btergm, and mtergm objects, with and without curved terms", {
  skip_on_cran()
  
  # simulate networks with fixed and changing covariate and gwidegree
  set.seed(12345)
  nnodes <- 30
  T <- 15
  sim <- list()
  fixed_covariate <- matrix(rnorm(nnodes^2), nrow = nnodes, ncol = nnodes)
  changing_covariate <- list()
  for (t in 1:T) {
    changing <- matrix(rnorm(nnodes^2), nrow = nnodes, ncol = nnodes)
    changing_covariate[[t]] <- changing
    sim[[t]] <- ergm::simulate_formula(network::network(nnodes) ~ edges + edgecov(fixed_covariate) + edgecov(changing) + gwidegree(0.5, fixed = TRUE),
                                       nsim = 1,
                                       coef = c(-2, 0.3, 0.6, 1.4))
  }
  
  # btergm with fixed GW decay
  expect_silent({
    fit1 <- btergm(sim ~ edges + edgecov(fixed_covariate) + edgecov(changing_covariate) + gwidegree(1.0, fixed = TRUE),
                   R = 100, verbose = FALSE)
  })
  expect_length(coef(fit1), 4)
  expect_silent(ep1 <- edgeprob(fit1))
  expect_s3_class(ep1, "data.frame")
  expect_equivalent(dim(ep1), c(13050, 11))
  
  # btergm with variable GW decay: currently unsupported in edgeprob
  expect_error(fit2 <- btergm(sim ~ edges + edgecov(fixed_covariate) + edgecov(changing_covariate) + gwidegree(fixed = FALSE),
                   R = 20, verbose = FALSE), "NAs generated during bootstrap")
  
  # mtergm with fixed GW decay
  expect_silent({
    fit3 <- mtergm(sim ~ edges + edgecov(fixed_covariate) + edgecov(changing_covariate) + gwidegree(1.0, fixed = TRUE),
                   verbose = FALSE)
  })
  expect_length(coef(fit3), 4)
  expect_silent(ep3 <- edgeprob(fit3))
  expect_s3_class(ep3, "data.frame")
  expect_equivalent(dim(ep3), c(13050, 11))
  
  # mtergm with variable GW decay: currently unsupported in edgeprob
  expect_silent({
    fit4 <- mtergm(sim ~ edges + edgecov(fixed_covariate) + edgecov(changing_covariate) + gwidegree(fixed = FALSE),
                   verbose = FALSE)
  })
  expect_length(coef(fit4), 5)
  expect_error(ep4 <- edgeprob(fit4), "MPLE-based \\(T\\)ERGMs with variable GW\\* decay are currently not supported")
  
  # ergm with fixed GW decay
  nnodes <- 50
  set.seed(12345)
  cov1 <- matrix(rnorm((nnodes)^2), nrow = nnodes, ncol = nnodes)
  cov2 <- matrix(rnorm((nnodes)^2), nrow = nnodes, ncol = nnodes)
  sim <- ergm::simulate_formula(network::network(nnodes) ~ edges + edgecov(cov1) + edgecov(cov2) + gwidegree(1.0, fixed = TRUE),
                                     nsim = 1,
                                     coef = c(-3, 0.3, 0.6, 0.8))
  expect_silent({
    suppressMessages(fit5 <- ergm::ergm(sim ~ edges + edgecov(cov1) + edgecov(cov2) + gwidegree(1.0, fixed = TRUE),
                                        verbose = FALSE))
  })
  expect_length(coef(fit5), 4)
  expect_silent(ep5 <- edgeprob(fit5))
  expect_s3_class(ep5, "data.frame")
  expect_equivalent(dim(ep5), c(2450, 11))
  
  # test validity of coefficients
  expect_equivalent(coef(fit1) - coef(fit3), rep(0, 4), tolerance = 0.1)
})

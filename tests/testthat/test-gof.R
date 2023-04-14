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

  fit <- suppressWarnings(btergm(networks ~ edges + istar(2) + edgecov(covariates), R = 100, verbose = FALSE))

  sink(nullfile())
  g <- btergm::gof(fit, nsim = 2, MCMC.burnin = 1000, MCMC.interval = 500, verbose = FALSE)
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

test_that("checkdegeneracy works with btergm and mtergm", {
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

  fit <- suppressWarnings(btergm(networks ~ edges + istar(2) + edgecov(covariates), R = 100, verbose = FALSE))
  d <- checkdegeneracy(fit)
  expect_equal(class(d), "degeneracy")
  expect_length(d$target.stats, 10)
  expect_length(d$sim, 10)
  expect_length(d$target.stats[[1]], 3)
  expect_equal(dim(d$sim[[1]]), c(1000, 3))

  fit2 <- mtergm(networks ~ edges + istar(2) + edgecov(covariates), verbose = FALSE)
  expect_output(checkdegeneracy(fit2), "Sample statistics summary")
})
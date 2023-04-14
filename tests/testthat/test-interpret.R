context("test micro-level interpretation functions")

test_that("edgeprob works with one-mode ERGM with network object", {
  skip_on_cran()
  set.seed(12345)
  data("chemnet")
  com <- committee %*% t(committee)
  nw <- network(pol, directed = TRUE)
  model1 <- ergm::ergm(nw ~ edges + edgecov(com) + istar(2))
  expect_error(ep <- edgeprob(model1), NA)
  expect_s3_class(ep, "data.frame")
  expect_equal(dim(ep), c(870, 10))
  expect_equal(colnames(ep), c("tie", "edges", "edgecov.com[[i]]", "istar2", "i", "j", "t", "i.name", "j.name", "probability"))
  expect_equal(ep$i[1], 1)
  expect_equal(ep$j[1], 2)
  expect_equal(ep$i.name[1], "BMA")
  expect_equal(ep$j.name[1], "BML")
  expect_false(any(is.na(ep$i.name)))
  expect_false(any(is.na(ep$j.name)))
  expect_equal(class(ep$probability), "numeric")
  expect_equal(class(ep$i), "integer")
  expect_equal(class(ep$j), "integer")
  expect_equal(class(ep$t), "integer")
  expect_equal(class(ep$i.name), "character")
  expect_equal(class(ep$j.name), "character")
  expect_equal(class(ep$tie), "integer")
  expect_lte(max(ep$probability), 1)
  expect_gte(max(ep$probability), 0)
})

test_that("edgeprob works with one-mode ERGM with matrix", {
  skip_on_cran()
  set.seed(12345)
  data("chemnet")
  com <- committee %*% t(committee)
  model1 <- ergm::ergm(pol ~ edges + edgecov(com) + istar(2))
  expect_error(ep <- edgeprob(model1), NA)
  expect_s3_class(ep, "data.frame")
  expect_equal(dim(ep), c(870, 10))
  expect_equal(colnames(ep), c("tie", "edges", "edgecov.com[[i]]", "istar2", "i", "j", "t", "i.name", "j.name", "probability"))
  expect_equal(ep$i[1], 1)
  expect_equal(ep$j[1], 2)
  expect_equal(ep$i.name[1], "BMA")
  expect_equal(ep$j.name[1], "BML")
  expect_false(any(is.na(ep$i.name)))
  expect_false(any(is.na(ep$j.name)))
  expect_equal(class(ep$probability), "numeric")
  expect_equal(class(ep$i), "integer")
  expect_equal(class(ep$j), "integer")
  expect_equal(class(ep$t), "integer")
  expect_equal(class(ep$i.name), "character")
  expect_equal(class(ep$j.name), "character")
  expect_equal(class(ep$tie), "integer")
  expect_lte(max(ep$probability), 1)
  expect_gte(max(ep$probability), 0)
})

test_that("edgeprob works with bipartite ERGM", {
  skip_on_cran()
  set.seed(12345)
  data("chemnet")
  cm <- network(committee, bipartite = TRUE, directed = FALSE)
  set.vertex.attribute(cm, "type", types[, 1])
  suppressMessages(model1 <- ergm::ergm(cm ~ edges + nodefactor("type", levels = 1) + b1star(2)))
  expect_error(ep <- edgeprob(model1), NA)
  expect_s3_class(ep, "data.frame")
  expect_equal(dim(ep), c(600, 10))
  expect_equal(colnames(ep), c("tie", "edges", "nodefactor.type.gov", "b1star2", "i", "j", "t", "i.name", "j.name", "probability"))
  expect_equal(ep$i[1], 1)
  expect_equal(ep$j[1], 31)
  expect_equal(ep$i.name[1], "BMA")
  expect_equal(ep$j.name[1], "IPU")
  expect_false(any(is.na(ep$i.name)))
  expect_false(any(is.na(ep$j.name)))
  expect_equal(class(ep$probability), "numeric")
  expect_equal(class(ep$i), "integer")
  expect_equal(class(ep$j), "integer")
  expect_equal(class(ep$t), "integer")
  expect_equal(class(ep$i.name), "character")
  expect_equal(class(ep$j.name), "character")
  expect_equal(class(ep$tie), "integer")
  expect_lte(max(ep$probability), 1)
  expect_gte(max(ep$probability), 0)
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
    fit1 <- suppressWarnings(btergm(sim ~ edges + edgecov(fixed_covariate) + edgecov(changing_covariate) + gwidegree(1.0, fixed = TRUE),
                                    R = 100, verbose = FALSE))
  })
  expect_length(coef(fit1), 4)
  expect_silent(ep1 <- edgeprob(fit1))
  expect_s3_class(ep1, "data.frame")
  expect_equivalent(dim(ep1), c(13050, 11))

  # btergm with variable GW decay: currently unsupported in edgeprob
  expect_warning(fit2 <- btergm(sim ~ edges + edgecov(fixed_covariate) + edgecov(changing_covariate) + gwidegree(fixed = FALSE),
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
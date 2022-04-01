context("test btergm estimation")

set.seed(12345)
networks <- list()
for(i in 1:10) {               # create 10 random networks with 10 actors
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

test_that("btergm estimation works", {
  set.seed(12345)
  fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), R = 100, verbose = FALSE)
  expect_equal(round(unname(coef(fit)), 4), c(-1.1707, 0.0543, 0.0045))
  expect_equal(names(coef(fit)), c("edges", "istar2", "edgecov.covariates[[i]]"))
  expect_equal(class(fit@boot), "boot")
  expect_equal(fit@boot$R, 100)
  expect_equal(fit@R, 100)
  expect_equal(fit@nobs, 900)
  expect_equal(fit@time.steps, 10)
  expect_equal(class(fit@formula), "formula")
  expect_equal(class(fit@formula2), "character")
  expect_equal(fit@formula, as.formula("networks ~ edges + istar(2) + edgecov(covariates)"))
  expect_equal(fit@formula2, "networks[[i]] ~ edges + istar(2) + edgecov(covariates[[i]])")
  expect_equal(length(fit@response), 900)
  expect_equal(is.numeric(fit@response), TRUE)
  expect_equal(class(fit@effects), "data.frame")
  expect_equal(dim(fit@effects), c(900, 3))
  expect_equal(unique(fit@effects$edges), 1)
  expect_equal(median(fit@effects$istar2), 2)
  expect_equal(round(mean(fit@effects$`edgecov.covariates[[i]]`), 4), -0.0144)
  expect_equal(unique(fit@weights), 1)
  expect_equal(fit@auto.adjust, FALSE)
  expect_equal(fit@offset, FALSE)
  expect_equal(fit@directed, TRUE)
  expect_equal(fit@bipartite, FALSE)
  expect_equal(unname(rowSums(fit@nvertices)), c(100, 100))
  expect_equal(round(confint(fit)[1, 3], 1), -1.4)
  expect_equal(round(confint(fit)[1, 4], 1), -0.8)
  expect_equal(round(confint(fit)[2, 3], 0), 0)
  expect_equal(round(confint(fit)[2, 4], 1), 0.1)
  expect_equal(round(confint(fit)[3, 3], 1), -0.1)
  expect_equal(round(confint(fit)[3, 4], 1), 0.1)
  expect_true(all(round(confint(fit)[, 1] - confint(fit)[, 2], 1) == 0))
})

test_that("fastglm works like speedglm", {
  skip_if_not_installed("fastglm", minimum_version = "0.0.1")

  set.seed(12345)
  fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), R = 100, verbose = FALSE)
  set.seed(12345)
  fit2 <- btergm(networks ~ edges + istar(2) + edgecov(covariates), R = 100, usefastglm = TRUE, verbose = FALSE)
  expect_equal(all(round(confint(fit), 4) == round(confint(fit2), 4)), TRUE)
})

test_that("offset argument in btergm works without composition change", {
  set.seed(12345)
  fit1 <- btergm(networks ~ edges + istar(2) + edgecov(covariates),
                 R = 100,
                 offset = FALSE,
                 usefastglm = TRUE,
                 verbose = FALSE)
  set.seed(12345)
  fit2 <- btergm(networks ~ edges + istar(2) + edgecov(covariates),
                 R = 100,
                 offset = TRUE,
                 usefastglm = TRUE,
                 verbose = FALSE)
  expect_equal(confint(fit1), confint(fit2))
})

test_that("offset argument in btergm works with composition change", {
  skip_on_cran()

  # example taken from 2018 JSS article
  require("sna")
  data("knecht")

  # step 1: make sure the network matrices have node labels
  for (i in 1:length(friendship)) {
    rownames(friendship[[i]]) <- 1:nrow(friendship[[i]])
    colnames(friendship[[i]]) <- 1:ncol(friendship[[i]])
  }
  rownames(primary) <- rownames(friendship[[1]])
  colnames(primary) <- colnames(friendship[[1]])
  sex <- demographics$sex
  names(sex) <- 1:length(sex)

  # step 2: imputation of NAs and removal of absent nodes:
  suppressMessages(friendship <- handleMissings(friendship, na = 10, method = "remove"))
  suppressMessages(friendship <- handleMissings(friendship, na = NA, method = "fillmode"))

  # step 3: add nodal covariates to the networks
  for (i in 1:length(friendship)) {
    s <- adjust(sex, friendship[[i]])
    friendship[[i]] <- network(friendship[[i]])
    friendship[[i]] <- set.vertex.attribute(friendship[[i]], "sex", s)
    idegsqrt <- sqrt(degree(friendship[[i]], cmode = "indegree"))
    friendship[[i]] <- set.vertex.attribute(friendship[[i]],
                                            "idegsqrt", idegsqrt)
    odegsqrt <- sqrt(degree(friendship[[i]], cmode = "outdegree"))
    friendship[[i]] <- set.vertex.attribute(friendship[[i]],
                                            "odegsqrt", odegsqrt)
  }
  expect_equal(unname(sapply(friendship, network.size)), c(26, 26, 25, 25))

  # step 4: estimate models
  set.seed(12345)
  m1 <- btergm(friendship ~ edges + mutual + ttriple +
                 transitiveties + ctriple + nodeicov("idegsqrt") +
                 nodeicov("odegsqrt") + nodeocov("odegsqrt") +
                 nodeofactor("sex") + nodeifactor("sex") + nodematch("sex") +
                 edgecov(primary) + delrecip + memory(type = "stability"),
               R = 100, usefastglm = TRUE, offset = TRUE, verbose = FALSE)
  m2 <- btergm(friendship ~ edges + mutual + ttriple +
                 transitiveties + ctriple + nodeicov("idegsqrt") +
                 nodeicov("odegsqrt") + nodeocov("odegsqrt") +
                 nodeofactor("sex") + nodeifactor("sex") + nodematch("sex") +
                 edgecov(primary) + delrecip + memory(type = "stability"),
               R = 100, usefastglm = TRUE, offset = FALSE, verbose = FALSE)

  # test results
  expect_equal(dim(confint(m1)), c(14, 4))
  expect_equal(dim(confint(m2)), c(14, 4))
  expect_equal(all(confint(m1)[, 4] - confint(m1)[, 3] > 0), TRUE)
  expect_equal(all(confint(m2)[, 4] - confint(m2)[, 3] > 0), TRUE)
  expect_equal(m1@offset, TRUE)
  expect_equal(m2@offset, FALSE)
  expect_equal(sapply(m1@data$offsmat, sum), c(0, 51, 51))
  expect_equal(sapply(m2@data$offsmat, sum), c(0, 0, 0))
  expect_equal(unname(nobs(m1)), c(3, 1850, 100))
  expect_equal(nobs(m1), nobs(m2))
})

test_that("mtergm estimation works", {
  skip_if_not_installed("fastglm", minimum_version = "0.0.1")

  set.seed(12345)
  fit1 <- btergm(networks ~ edges + istar(2) + edgecov(covariates), usefastglm = TRUE, verbose = FALSE)
  set.seed(12345)
  fit2 <- mtergm(networks ~ edges + istar(2) + edgecov(covariates), verbose = FALSE)
  expect_equivalent(coef(fit1) - coef(fit2), rep(0, 3), tolerance = 0.05)
  expect_equivalent(coef(fit2), c(-1.17, 0.06, 0.00), tolerance = 0.05)
  expect_equivalent(fit2@se, c(0.193, 0.079, 0.074), tolerance = 0.05)
  expect_equal(class(fit2)[1], "mtergm")
  expect_equal(class(fit2@ergm), "ergm")
})

test_that("simulation of new networks works", {
  skip_if_not_installed("fastglm", minimum_version = "0.0.1")

  # for btergm
  fit1 <- btergm(networks ~ edges + istar(2) + edgecov(covariates), usefastglm = TRUE, verbose = FALSE)
  sim1 <- simulate(fit1, 5, verbose = FALSE)
  expect_equal(length(sim1), 5)
  expect_equal(class(sim1[[1]]), "network")

  # with an index
  sim1b <- simulate(fit1, index = 3, verbose = FALSE)
  expect_equal(class(sim1b), "network")

  # for mtergm
  fit2 <- mtergm(networks ~ edges + istar(2) + edgecov(covariates), verbose = FALSE)
  sim2 <- simulate(fit2, 5, verbose = FALSE)
  expect_equal(length(sim2), 5)
  expect_equal(class(sim2[[1]]), "network")
})

test_that("tbergm estimation works", {
  skip_on_cran()
  skip_if_not_installed("Bergm", minimum_version = "5.0.2")
  set.seed(12345)
  fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), R = 100, verbose = FALSE)
  suppressMessages(suppressWarnings(fit_b <- tbergm(networks ~ edges + istar(2) + edgecov(covariates), verbose = FALSE)))
  expect_s3_class(fit_b@bergm, "bergm")
  expect_length(fit_b@bergm$ess, 4)
  expect_length(fit_b@data, 3)
  expect_equal(dim(as.matrix(fit_b@bergm$Theta)), c(8000, 4))
  expect_silent(g <- btergm::gof(fit_b, verbose = FALSE))
  expect_s3_class(g, "gof")
  expect_length(g, 7)
  expect_equivalent(nobs(fit_b), c(10, 900))
  expect_equivalent(coef(fit) - apply(fit_b@bergm$Theta, 2, mean)[-4], rep(0, 3), tolerance = 0.2)
})

test_that("bipartite network objects work, including composition change, with and without offset", {
  skip_on_cran()
  set.seed(12345)
  nw_matrix <- lapply(1:20, function(x) {
    mat <- matrix(rbinom(200, 1, 0.2), nrow = 20)
    rownames(mat) <- letters[1:20]
    colnames(mat) <- letters[1:10]
    mat <- mat[sample(1:20, 18), sample(1:10, 9)]
  })
  nw_network <- lapply(nw_matrix, function(x) {
    network::network(x, bipartite = TRUE, directed = FALSE)
  })
  expect_warning(model1 <- btergm(nw_matrix ~ edges + threetrail + kstar(2) + memory("autoregression"), R = 50, offset = FALSE, verbose = FALSE), regexp = NA)
  expect_warning(model2 <- btergm(nw_network ~ edges + threetrail + kstar(2) + memory("autoregression"), R = 50, offset = FALSE, verbose = FALSE), regexp = NA)
  expect_warning(model3 <- btergm(nw_matrix ~ edges + threetrail + kstar(2) + memory("autoregression"), R = 50, offset = TRUE, verbose = FALSE), regexp = NA)
  expect_warning(model4 <- btergm(nw_network ~ edges + threetrail + kstar(2) + memory("autoregression"), R = 50, offset = TRUE, verbose = FALSE), regexp = NA)
  expect_length(model1@coef, 4)
  expect_length(model2@coef, 4)
  expect_length(model3@coef, 4)
  expect_length(model4@coef, 4)
  expect_equal(dim(confint(model1)), c(4, 4))
  expect_equal(dim(confint(model2)), c(4, 4))
  expect_equal(dim(confint(model3)), c(4, 4))
  expect_equal(dim(confint(model4)), c(4, 4))
})
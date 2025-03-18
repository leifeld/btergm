context("test goodness-of-fit functions")

test_that("basic GOF functionality works", {
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

test_that("gof statistics for vectors work", {
  set.seed(12345)
  mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
  diag(mat) <- 0               # loops are excluded
  nw <- network::network(mat)  # create network object
  
  # dsp
  expect_no_error(s <- dsp(mat))
  expect_equal(attributes(s)$label, "Dyad-wise shared partners")
  expect_equal(names(s), as.character(0:8))
  expect_length(s, 9)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 90)
  
  # esp
  expect_no_error(s <- esp(mat))
  expect_equal(attributes(s)$label, "Edge-wise shared partners")
  expect_equal(names(s), as.character(0:8))
  expect_length(s, 9)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 23)
  
  # nsp
  expect_no_error(s <- nsp(mat))
  expect_equal(attributes(s)$label, "Non-edge-wise shared partners")
  expect_equal(names(s), as.character(0:8))
  expect_length(s, 9)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 67)
  
  # deg
  expect_no_error(s <- deg(mat))
  expect_equal(attributes(s)$label, "Degree")
  expect_equal(names(s), as.character(0:9))
  expect_length(s, 10)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 10)
  
  # b1deg
  expect_no_error(s <- b1deg(mat))
  expect_equal(attributes(s)$label, "Degree (first mode)")
  expect_equal(names(s), as.character(0:10))
  expect_length(s, 11)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 10)
  
  # b2deg
  expect_no_error(s <- b2deg(mat))
  expect_equal(attributes(s)$label, "Degree (second mode)")
  expect_equal(names(s), as.character(0:10))
  expect_length(s, 11)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 10)
  
  # odeg
  expect_no_error(s <- odeg(mat))
  expect_equal(attributes(s)$label, "Outdegree")
  expect_equal(names(s), as.character(0:9))
  expect_length(s, 10)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 10)
  
  # ideg
  expect_no_error(s <- ideg(mat))
  expect_equal(attributes(s)$label, "Indegree")
  expect_equal(names(s), as.character(0:9))
  expect_length(s, 10)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 10)
  
  # kstar
  expect_no_error(s <- kstar(mat))
  expect_equal(attributes(s)$label, "k-star")
  expect_equal(names(s), as.character(0:9))
  expect_length(s, 10)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 442)
  
  # b1star
  expect_no_error(s <- b1star(mat))
  expect_equal(attributes(s)$label, "k-star (first mode)")
  expect_equal(names(s), as.character(0:10))
  expect_length(s, 11)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 100)
  
  # b2star
  expect_no_error(s <- b2star(mat))
  expect_equal(attributes(s)$label, "k-star (second mode)")
  expect_equal(names(s), as.character(0:10))
  expect_length(s, 11)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 71)
  
  # ostar
  expect_no_error(s <- ostar(mat))
  expect_equal(attributes(s)$label, "Outgoing k-star")
  expect_equal(names(s), as.character(0:9))
  expect_length(s, 10)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 100)
  
  # istar
  expect_no_error(s <- istar(mat))
  expect_equal(attributes(s)$label, "Incoming k-star")
  expect_equal(names(s), as.character(0:9))
  expect_length(s, 10)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 71)
  
  # kcycle
  expect_warning(s <- kcycle(mat), "cycles of length less than 2 cannot exist")
  expect_equal(attributes(s)$label, "Cycle")
  expect_equal(names(s), as.character(0:7))
  expect_length(s, 8)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 14)
  
  # geodesic
  expect_no_error(s <- geodesic(mat))
  expect_equal(attributes(s)$label, "Geodesic distances")
  expect_equal(names(s), c(as.character(1:9), "Inf"))
  expect_length(s, 10)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 90)
  
  # triad.directed
  expect_no_error(s <- triad.directed(mat))
  expect_equal(attributes(s)$label, "Triad census")
  expect_equal(names(s), c("003", "012", "102", "021D", "021U", "021C", "111D", "111U", "030T", "030C", "201", "120D", "120U", "120C", "210", "300"))
  expect_length(s, 16)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 120)
  
  # triad.undirected
  expect_no_error(s <- triad.undirected(mat))
  expect_equal(attributes(s)$label, "Triad census")
  expect_equal(names(s), as.character(0:3))
  expect_length(s, 4)
  expect_contains(class(s), "numeric")
  expect_equal(sum(s), 120)
})

test_that("gof statistics for single values work", {
  set.seed(12345)
  mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
  diag(mat) <- 0               # loops are excluded
  nw <- network::network(mat)  # create network object
  
  # walktrap.modularity
  expect_no_error(s <- walktrap.modularity(mat))
  expect_equal(attributes(s)$label, "Modularity (walktrap)")
  expect_length(s, 1)
  expect_contains(class(s), "numeric")
  expect_gt(s, 0)
  expect_lt(s, 1)
  
  # fastgreedy.modularity
  expect_no_error(s <- fastgreedy.modularity(mat))
  expect_equal(attributes(s)$label, "Modularity (fast & greedy)")
  expect_length(s, 1)
  expect_contains(class(s), "numeric")
  expect_gt(s, 0)
  expect_lt(s, 1)
  
  # louvain.modularity
  expect_no_error(s <- louvain.modularity(mat))
  expect_equal(attributes(s)$label, "Modularity (Louvain)")
  expect_length(s, 1)
  expect_contains(class(s), "numeric")
  expect_gt(s, 0)
  expect_lt(s, 1)
  
  # maxmod.modularity
  expect_no_error(s <- maxmod.modularity(mat))
  expect_equal(attributes(s)$label, "Maximum modularity")
  expect_length(s, 1)
  expect_contains(class(s), "numeric")
  expect_gt(s, 0)
  expect_lt(s, 1)
  
  # edgebetweenness.modularity
  expect_no_error(s <- edgebetweenness.modularity(mat))
  expect_equal(attributes(s)$label, "Modularity (edge betweenness)")
  expect_length(s, 1)
  expect_contains(class(s), "numeric")
  expect_gt(s, 0)
  expect_lt(s, 1)
  
  # spinglass.modularity
  expect_no_error(s <- spinglass.modularity(mat))
  expect_equal(attributes(s)$label, "Modularity (spinglass)")
  expect_length(s, 1)
  expect_contains(class(s), "numeric")
  expect_gt(s, 0)
  expect_lt(s, 1)
})

test_that("tie prediction gof statistics work", {
  set.seed(16) # note that a specific random seed is used because the spinglass algorithm requires connected networks in both the observed and target networks
  networks <- list()
  for(i in 1:7) {                # create 7 random networks with 10 actors
    mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
    diag(mat) <- 0               # loops are excluded
    nw <- network::network(mat)  # create network object
    networks[[i]] <- nw          # add network to the list
  }
  
  covariates <- list()
  for (i in 1:7) {              # create 7 matrices as covariate
    mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
    covariates[[i]] <- mat       # add matrix to the list
  }
  
  fit <- suppressWarnings(btergm(networks ~ edges + istar(2) + edgecov(covariates), R = 100, verbose = FALSE))
  
  # spinglass.roc
  expect_no_error(g <- btergm::gof(fit, statistics = spinglass.roc, nsim = 10, verbose = FALSE))
  expect_contains(class(g), "gof")
  expect_length(g, 1)
  expect_length(g[[1]], 6)
  expect_equal(names(g), "Tie prediction")
  expect_equal(names(g[[1]]), c("type", "auc.roc", "auc.roc.rgraph", "roc", "roc.rgraph", "label"))
  expect_equal(g[[1]]$label, "Spinglass community comembership prediction")
  expect_gt(g[[1]]$auc.roc, 0)
  expect_lt(g[[1]]$auc.roc, 1)
  
  # spinglass.pr
  set.seed(3) # note that a specific random seed is used because the spinglass algorithm requires connected networks in both the observed and target networks
  expect_no_error(g <- btergm::gof(fit, statistics = spinglass.pr, nsim = 10, verbose = FALSE))
  expect_contains(class(g), "gof")
  expect_length(g, 1)
  expect_length(g[[1]], 6)
  expect_equal(names(g), "Tie prediction")
  expect_equal(names(g[[1]]), c("type", "auc.pr", "auc.pr.rgraph", "pr", "pr.rgraph", "label"))
  expect_equal(g[[1]]$label, "Spinglass community comembership prediction")
  expect_gt(g[[1]]$auc.pr, 0)
  expect_lt(g[[1]]$auc.pr, 1)
  
  # walktrap.roc
  expect_no_error(g <- btergm::gof(fit, statistics = walktrap.roc, nsim = 10, verbose = FALSE))
  expect_contains(class(g), "gof")
  expect_length(g, 1)
  expect_length(g[[1]], 6)
  expect_equal(names(g), "Tie prediction")
  expect_equal(names(g[[1]]), c("type", "auc.roc", "auc.roc.rgraph", "roc", "roc.rgraph", "label"))
  expect_equal(g[[1]]$label, "Walktrap community comembership prediction")
  expect_gt(g[[1]]$auc.roc, 0)
  expect_lt(g[[1]]$auc.roc, 1)
  
  # walktrap.pr
  expect_no_error(g <- btergm::gof(fit, statistics = walktrap.pr, nsim = 10, verbose = FALSE))
  expect_contains(class(g), "gof")
  expect_length(g, 1)
  expect_length(g[[1]], 6)
  expect_equal(names(g), "Tie prediction")
  expect_equal(names(g[[1]]), c("type", "auc.pr", "auc.pr.rgraph", "pr", "pr.rgraph", "label"))
  expect_equal(g[[1]]$label, "Walktrap community comembership prediction")
  expect_gt(g[[1]]$auc.pr, 0)
  expect_lt(g[[1]]$auc.pr, 1)
  
  # fastgreedy.roc
  expect_no_error(g <- btergm::gof(fit, statistics = fastgreedy.roc, nsim = 10, verbose = FALSE))
  expect_contains(class(g), "gof")
  expect_length(g, 1)
  expect_length(g[[1]], 6)
  expect_equal(names(g), "Tie prediction")
  expect_equal(names(g[[1]]), c("type", "auc.roc", "auc.roc.rgraph", "roc", "roc.rgraph", "label"))
  expect_equal(g[[1]]$label, "Fast & greedy community comembership prediction")
  expect_gt(g[[1]]$auc.roc, 0)
  expect_lt(g[[1]]$auc.roc, 1)
  
  # fastgreedy.pr
  expect_no_error(g <- btergm::gof(fit, statistics = fastgreedy.pr, nsim = 10, verbose = FALSE))
  expect_contains(class(g), "gof")
  expect_length(g, 1)
  expect_length(g[[1]], 6)
  expect_equal(names(g), "Tie prediction")
  expect_equal(names(g[[1]]), c("type", "auc.pr", "auc.pr.rgraph", "pr", "pr.rgraph", "label"))
  expect_equal(g[[1]]$label, "Fast & greedy community comembership prediction")
  expect_gt(g[[1]]$auc.pr, 0)
  expect_lt(g[[1]]$auc.pr, 1)
  
  # louvain.roc
  expect_no_error(g <- btergm::gof(fit, statistics = louvain.roc, nsim = 10, verbose = FALSE))
  expect_contains(class(g), "gof")
  expect_length(g, 1)
  expect_length(g[[1]], 6)
  expect_equal(names(g), "Tie prediction")
  expect_equal(names(g[[1]]), c("type", "auc.roc", "auc.roc.rgraph", "roc", "roc.rgraph", "label"))
  expect_equal(g[[1]]$label, "Louvain community comembership prediction")
  expect_gt(g[[1]]$auc.roc, 0)
  expect_lt(g[[1]]$auc.roc, 1)
  
  # louvain.pr
  expect_no_error(g <- btergm::gof(fit, statistics = louvain.pr, nsim = 10, verbose = FALSE))
  expect_contains(class(g), "gof")
  expect_length(g, 1)
  expect_length(g[[1]], 6)
  expect_equal(names(g), "Tie prediction")
  expect_equal(names(g[[1]]), c("type", "auc.pr", "auc.pr.rgraph", "pr", "pr.rgraph", "label"))
  expect_equal(g[[1]]$label, "Louvain community comembership prediction")
  expect_gt(g[[1]]$auc.pr, 0)
  expect_lt(g[[1]]$auc.pr, 1)
  
  # maxmod.roc
  expect_no_error(g <- btergm::gof(fit, statistics = maxmod.roc, nsim = 10, verbose = FALSE))
  expect_contains(class(g), "gof")
  expect_length(g, 1)
  expect_length(g[[1]], 6)
  expect_equal(names(g), "Tie prediction")
  expect_equal(names(g[[1]]), c("type", "auc.roc", "auc.roc.rgraph", "roc", "roc.rgraph", "label"))
  expect_equal(g[[1]]$label, "Maximum modularity community comembership prediction")
  expect_gt(g[[1]]$auc.roc, 0)
  expect_lt(g[[1]]$auc.roc, 1)
  
  # maxmod.pr
  expect_no_error(g <- btergm::gof(fit, statistics = maxmod.pr, nsim = 10, verbose = FALSE))
  expect_contains(class(g), "gof")
  expect_length(g, 1)
  expect_length(g[[1]], 6)
  expect_equal(names(g), "Tie prediction")
  expect_equal(names(g[[1]]), c("type", "auc.pr", "auc.pr.rgraph", "pr", "pr.rgraph", "label"))
  expect_equal(g[[1]]$label, "Maximum modularity community comembership prediction")
  expect_gt(g[[1]]$auc.pr, 0)
  expect_lt(g[[1]]$auc.pr, 1)
  
  # edgebetweenness.roc
  expect_no_error(g <- btergm::gof(fit, statistics = edgebetweenness.roc, nsim = 10, verbose = FALSE))
  expect_contains(class(g), "gof")
  expect_length(g, 1)
  expect_length(g[[1]], 6)
  expect_equal(names(g), "Tie prediction")
  expect_equal(names(g[[1]]), c("type", "auc.roc", "auc.roc.rgraph", "roc", "roc.rgraph", "label"))
  expect_equal(g[[1]]$label, "Edge betweenness community comembership prediction")
  expect_gt(g[[1]]$auc.roc, 0)
  expect_lt(g[[1]]$auc.roc, 1)
  
  # edgebetweenness.pr
  expect_no_error(g <- btergm::gof(fit, statistics = edgebetweenness.pr, nsim = 10, verbose = FALSE))
  expect_contains(class(g), "gof")
  expect_length(g, 1)
  expect_length(g[[1]], 6)
  expect_equal(names(g), "Tie prediction")
  expect_equal(names(g[[1]]), c("type", "auc.pr", "auc.pr.rgraph", "pr", "pr.rgraph", "label"))
  expect_equal(g[[1]]$label, "Edge betweenness community comembership prediction")
  expect_gt(g[[1]]$auc.pr, 0)
  expect_lt(g[[1]]$auc.pr, 1)
})

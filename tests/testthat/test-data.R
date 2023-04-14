context("test btergm estimation")

test_that("chemnet example can be replicated", {
  skip_on_cran()

  # preparatory steps
  library("network")
  library("sna")
  library("ergm")
  seed <- 12345
  set.seed(seed)
  data("chemnet")

  # create confirmed network relation
  sci <- scito * t(scifrom)  # equation 1 in the AJPS paper
  prefsim <- dist(intpos, method = "euclidean")  # equation 2
  prefsim <- max(prefsim) - prefsim  # equation 3
  prefsim <- as.matrix(prefsim)
  committee <- committee %*% t(committee)  # equation 4
  diag(committee) <- 0 # the diagonal has no meaning
  types <- types[, 1]  # convert to vector

  # create network objects and store attributes
  nw.pol <- network(pol) # political/stratgic information exchange
  set.vertex.attribute(nw.pol, "orgtype", types)
  set.vertex.attribute(nw.pol, "betweenness",
      betweenness(nw.pol)) # centrality

  nw.sci <- network::network(sci) # technical/scientific information exchange
  set.vertex.attribute(nw.sci, "orgtype", types)
  set.vertex.attribute(nw.sci, "betweenness",
      betweenness(nw.sci)) # centrality

  # ERGM: model 1 in the AJPS paper; only preference similarity
  suppressMessages({
    model1 <- ergm(nw.pol ~ edges + edgecov(prefsim),
                   control = control.ergm(seed = seed))
  })
  expect_equal(class(model1), "ergm")
  expect_length(coef(model1), 2)

  # ERGM: model 2 in the AJPS paper; complete model
  suppressMessages({
    model2 <- ergm(nw.pol ~
                     edges +
                     edgecov(prefsim) +
                     mutual +
                     nodemix("orgtype", base = -7) +
                     nodeifactor("orgtype", base = -1) +
                     nodeofactor("orgtype", base = -5) +
                     edgecov(committee) +
                     edgecov(nw.sci) +
                     edgecov(infrep) +
                     gwesp(0.1, fixed = TRUE) +
                     gwdsp(0.1, fixed = TRUE),
                   control = control.ergm(seed = seed)
    )
  })
  expect_equal(class(model2), "ergm")
  expect_length(coef(model2), 11)

  # ERGM: model 3 in the AJPS paper; only preference similarity
  suppressMessages({
    model3 <- ergm(nw.sci ~ edges + edgecov(prefsim),
                   control = control.ergm(seed = seed))
  })
  expect_equal(class(model3), "ergm")
  expect_length(coef(model3), 2)

  # ERGM: model 4 in the AJPS paper; complete model
  suppressMessages({
    model4 <- ergm(nw.sci ~
                     edges +
                     edgecov(prefsim) +
                     mutual +
                     nodemix("orgtype", base = -7) +
                     nodeifactor("orgtype", base = -1) +
                     nodeofactor("orgtype", base = -5) +
                     edgecov(committee) +
                     edgecov(nw.pol) +
                     edgecov(infrep) +
                     gwesp(0.1, fixed = TRUE) +
                     gwdsp(0.1, fixed = TRUE),
                   control = control.ergm(seed = seed)
    )
  })
  expect_equal(class(model4), "ergm")
  expect_length(coef(model4), 11)

  # goodness of fit using the btergm package
  gof4 <- btergm::gof(model4, verbose = FALSE)
  expect_equal(length(gof4), 7)
  plot(gof4)
})
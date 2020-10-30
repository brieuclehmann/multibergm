test_that("Default priors for single group return expected results", {
  n_nets <- 10
  n_nodes <- 10
  n_iters <- 5
  
  nets <- lapply(seq_len(n_nets), function(x) network(n_nodes))
  
  ergm_formula <- nets ~ edges

  prior <- set_priors(ergm_formula)
  
  expect_equal(prior$muPop$mean, rep(0, 1))
  expect_equal(prior$muPop$cov, diag(100, 1))
  expect_equal(prior$covMuGroup, NULL)
  expect_equal(prior$covTheta$df, 2)
  expect_equal(prior$covTheta$scale, diag(1))
})

test_that("Default priors for multiple groups return expected results", {
  n_nets <- 10
  n_nodes <- 10
  n_iters <- 5
  
  nets <- lapply(seq_len(n_nets), function(x) network(n_nodes))
  
  ergm_formula <- nets ~ edges
  groups <- rep(c(1, 2), 5)
  
  prior <- set_priors(ergm_formula, groups)
  
  expect_equal(prior$muPop$mean, rep(0, 1))
  expect_equal(prior$muPop$cov, diag(100, 1))
  expect_equal(prior$covMuGroup$df, 2)
  expect_equal(prior$covMuGroup$scale, diag(1))
  expect_equal(prior$covTheta$df, 2)
  expect_equal(prior$covTheta$scale, diag(1))
})

test_that("User-specified priors fail when invalid", {
  n_nets <- 10
  n_nodes <- 10
  n_iters <- 5
  
  nets <- lapply(seq_len(n_nets), function(x) network(n_nodes))
  
  ergm_formula <- nets ~ edges
  
  prior <- list(muPop = list(cov = diag(2)))
  expect_error(set_priors(ergm_formula, groups, prior))
  
  prior <- list(muPop = list(mean = c(0, 0)))
  expect_error(set_priors(ergm_formula, groups, prior))
  
  prior <- list(covTheta = list(df = -1))
  expect_error(set_priors(ergm_formula, groups, prior))
  
  prior <- list(covTheta = list(scale = c(0, 0)))
  expect_error(set_priors(ergm_formula, groups, prior))
})
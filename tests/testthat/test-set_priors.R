test_that("Default priors for single group return expected results", {
  n_nets <- 10
  n_nodes <- 10
  n_iters <- 5

  nets <- lapply(seq_len(n_nets), function(x) network(n_nodes))

  ergm_formula <- nets ~ edges
  control <- control_multibergm(ergm_formula)
  prior <- set_priors(ergm_formula, control)

  expect_equal(prior$mu_pop$mean, rep(0, 1))
  expect_equal(prior$mu_pop$cov, diag(100, 1))
  expect_equal(prior$cov_mu_group, NULL)
  expect_equal(prior$cov_theta$df, 2)
  expect_equal(prior$cov_theta$scale, diag(1))
})

test_that("Default priors for multiple groups return expected results", {
  n_nets <- 10
  n_nodes <- 10
  n_iters <- 5

  nets <- lapply(seq_len(n_nets), function(x) network(n_nodes))

  ergm_formula <- nets ~ edges
  groups <- rep(c(1, 2), 5)
  control = control_multibergm(ergm_formula, groups = groups)

  prior <- set_priors(ergm_formula, control, groups)

  expect_equal(prior$mu_pop$mean, rep(0, 1))
  expect_equal(prior$mu_pop$cov, diag(100, 1))
  expect_equal(prior$cov_mu_group$df, 2)
  expect_equal(prior$cov_mu_group$scale, diag(1))
  expect_equal(prior$cov_theta$df, 2)
  expect_equal(prior$cov_theta$scale, diag(1))
})

test_that("User-specified priors fail when invalid", {
  n_nets <- 10
  n_nodes <- 10
  n_iters <- 5

  nets <- lapply(seq_len(n_nets), function(x) network(n_nodes))

  ergm_formula <- nets ~ edges
  control <- control_multibergm(ergm_formula)

  prior <- list(mu_pop = list(cov = diag(2)))
  expect_error(set_priors(ergm_formula, control, groups, prior))

  prior <- list(mu_pop = list(mean = c(0, 0)))
  expect_error(set_priors(ergm_formula, control, groups, prior))

  prior <- list(cov_theta = list(df = -1))
  expect_error(set_priors(ergm_formula, control, groups, prior))

  prior <- list(cov_theta = list(scale = c(0, 0)))
  expect_error(set_priors(ergm_formula, control, groups, prior))
})

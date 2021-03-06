### ERGM wrapper

test_that("ergm wrapper returns correct dimension", {
  n_nets <- 10L
  n_nodes <- 10L
  n_iters <- 5L

  nets <- lapply(seq_len(n_nets),
                 function(x) network(n_nodes, directed = FALSE))

  model <- nets ~ edges
  control <- control_multibergm(model)
  coefs <- matrix(1, n_nets, 1)

  expect_equal(dim(ergm_wrapper(coefs, control)), dim(coefs))
})

### Exchange update

test_that("exchange_update() returns correct dimension", {
  n <- 3
  p <- 2
  curr <- rmvnorm(n, rep(0, p))
  prop <- rmvnorm(n, rep(0, p))
  prior_cov <- diag(1, p)
  delta <- rmvnorm(n, rep(0, p))
  etamap <- NULL

  expect_equal(dim(exchange_update(curr, prop, delta, prior_cov, etamap)),
               c(n, p))
  expect_equal(dim(exchange_update(curr[1, ], prop[1, ],
                                   delta[1, , drop = FALSE], prior_cov, etamap)),
               NULL)
})

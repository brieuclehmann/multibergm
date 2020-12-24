
n_nets <- 10L
n_nodes <- 10L
n_iters <- 5L

nets <- lapply(seq_len(n_nets), function(x) network(n_nodes, directed = FALSE))

ergm_formula_single <- nets ~ edges
fit_single <- multibergm(ergm_formula_single, main_iters = n_iters)

ergm_formula <- nets ~ edges + triangles
fit <- multibergm(ergm_formula, main_iters = n_iters)

test_that("Basic single group multibergm runs without errors", {

  expect_equal(dim(fit_single$params$theta), c(1, n_iters, n_nets, 1))
  expect_equal(length(fit_single$networks), n_nets)
  
  expect_equal(dim(fit$params$theta), c(1, n_iters, n_nets, 2))
  expect_equal(length(fit$networks), n_nets)

})

test_that("print.multibergm() returns correct output", {

  verify_output(test_path("print_one_group.txt"), fit_single)
})


test_that("summary.multibergm() returns correct output", {

  verify_output(test_path("summary_one_group.txt"), summary(fit_single))
})


test_that("plot.multibergm() runs as expected", {

  vdiffr::expect_doppelganger("multibergm MCMC plots", plot(fit_single))
})

test_that("gof.multibergm() runs as expected", {

  vdiffr::expect_doppelganger("multibergm MCMC gof",
                              gof(fit_single, sample_size = 5))
})

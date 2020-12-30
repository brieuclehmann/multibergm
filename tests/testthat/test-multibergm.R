set.seed(1)
n_nets <- 10L
n_nodes <- 10L
n_iters <- 40L

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

  verify_output(test_path("print_two_group.txt"), fit_single)
})


test_that("summary.multibergm() returns correct output", {

  verify_output(test_path("summary_two_group.txt"), summary(fit_single))
})


group_ind <- rep(c(1,2), each = n_nodes / 2)
fit_edges_two <- multibergm(ergm_formula_single, 
                            groups = group_ind, 
                            main_iters = n_iters)

ergm_formula <- nets ~ edges + triangles
fit_edges_triangle_two <- multibergm(ergm_formula, main_iters = n_iters)

test_that("Two group multibergm runs without errors", {
  
  expect_equal(dim(fit_single$params$theta), c(1, n_iters, n_nets, 1))
  expect_equal(length(fit_single$networks), n_nets)
  
  expect_equal(dim(fit$params$theta), c(1, n_iters, n_nets, 2))
  expect_equal(length(fit$networks), n_nets)
  
})
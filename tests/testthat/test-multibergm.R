data(nets)

formula <- nets ~ edges

fit <- multibergm(formula, mainIters = 10)

test_that("Basic single group multibergm runs without errors", {
  
  expect_equal(dim(fit$params$theta), c(10, 10, 1))
  expect_equal(length(fit$networks), length(nets))
  
})

test_that("print.multibergm() returns correct output", {
  
  verify_output(test_path("print_one_group.txt"), fit)
})
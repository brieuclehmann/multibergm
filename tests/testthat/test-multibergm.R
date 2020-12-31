set.seed(1)

### GENERATE NETWORKS

n_nets <- 10L
n_nodes <- 10L
n_iters <- 30L

nets <- lapply(seq_len(n_nets), function(x) network(n_nodes, directed = FALSE))
group_ind <- rep(c(1, 2), each = n_nodes / 2)

### SET FORMULAS

ergm_formula_one_stat <- nets ~ edges
ergm_formula_two_stat <- nets ~ edges + triangles
ergm_formula_curved   <- nets ~ edges + gwesp()

formulas <- c(ergm_formula_one_stat, ergm_formula_two_stat, ergm_formula_curved)

### RUN TESTS FOR SINGLE GROUP

for (f in formulas) {
  fit <- multibergm(f, main_iters = n_iters)
  
  # Main function
  test_name <- paste("Single group multibergm with", 
                     format(fit$formula), "runs without errors")
  test_that(test_name, {
    n_terms <- nparam(fit$control$model)
    expect_equal(dim(fit$params$theta), c(1, n_iters, n_nets, n_terms))
    expect_equal(length(fit$networks), n_nets)
  })
  
  # Print
  test_name <- paste("print.multibergm() returns correct output for",
                      format(fit$formula), "with a single group")
  params <- paste(param_names(fit$control$model), collapse = "_")
  txt_path <- paste("print", params, "single.txt", sep = "_")
  test_that(test_name, {
    verify_output(test_path(txt_path), fit)
  })
  
  # Summary
  test_name <- paste("summary.multibergm() returns correct output for",
                     format(fit$formula), "with a single group")
  #TODO: Make this test more robust
  set.seed(1)
  fit$params$mu_pop[] <- rnorm(length(fit$params$mu_pop[]))
  fit$accepts$theta[] <- rbinom(length(fit$accepts$theta[]), 1, 0.5)
  fit$accepts$mu[]    <- rbinom(length(fit$accepts$mu[]), 1, 0.5)
  txt_path <- paste("summary", params, "single.txt", sep = "_")
  test_that(test_name, {
    verify_output(test_path(txt_path), summary(fit))
  })
}


### RUN TESTS FOR TwO GROUP

for (f in formulas) {
  set.seed(1)
  fit <- multibergm(f, main_iters = n_iters, groups = group_ind)
  
  # Main function
  test_name <- paste("Two group multibergm with", 
                     format(fit$formula), "runs without errors")
  test_that(test_name, {
    n_terms <- nparam(fit$control$model)
    expect_equal(dim(fit$params$theta), c(1, n_iters, n_nets, n_terms))
    expect_equal(length(fit$networks), n_nets)
  })
  
  # Print
  test_name <- paste("print.multibergm() returns correct output for",
                     format(fit$formula), "with two groups")
  params <- paste(param_names(fit$control$model), collapse = "_")
  txt_path <- paste("print", params, "twogrp.txt", sep = "_")
  test_that(test_name, {
    verify_output(test_path(txt_path), fit)
  })
  
  # Summary
  test_name <- paste("summary.multibergm() returns correct format for",
                     format(fit$formula), "with two groups")
  #TODO: Make this test more robust
  set.seed(1)
  fit$params$mu_group[] <- rnorm(length(fit$params$mu_group[]))
  fit$accepts$theta[]   <- rbinom(length(fit$accepts$theta[]), 1, 0.5)
  fit$accepts$mu[]      <- rbinom(length(fit$accepts$mu[]), 1, 0.5)
  txt_path <- paste("summary", params, "twogrp.txt", sep = "_")
  test_that(test_name, {
    verify_output(test_path(txt_path), summary(fit))
  })
}
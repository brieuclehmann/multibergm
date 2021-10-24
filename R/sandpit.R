ergm_formula <- nets ~ edges + triangles + gwesp(0.5, fixed = TRUE)
model_formula <- ~ 1 + age
n_nets <- 10L
n_nodes <- 10L
n_iters <- 5L
nets <- lapply(seq_len(n_nets), function(x) network(n_nodes, directed = FALSE))
nets <- lapply(nets, function(x) set.network.attribute(x, "age", rnorm(1)))

model_matrix <- get_model_matrix(ergm_formula, model_formula)
control <- control_multibergm(ergm_formula, model_formula,
                             ~ .)
prior <- set_priors(ergm_formula, 
                    model_matrix, 
                    control)



nets <- list()

for (n in 1:20){
  nets[[n]] <- simulate(y ~ edges, coef = -2,
                        basis = network(n_nodes, directed=FALSE))
}

  
fit <- multibergm(ergm_formula, model_formula)

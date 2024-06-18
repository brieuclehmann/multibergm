
devtools::load_all()

ergm_formula <- nets ~ edges + triangles + gwesp(0.5, fixed = TRUE)
ergm_formula <- nets ~ edges
model_formula <- ~ 1 + group
n_nets <- 10L
n_nodes <- 10L
n_iters <- 5L
ages <- rnorm(n_nets)
groups <- rep(c(0, 1), each = n_nets / 2)

nets <- lapply(seq_len(n_nets), function(x) network(n_nodes, directed = FALSE))
nets <- lapply(1:n_nets, function(x) set.network.attribute(nets[[x]], "group", groups[x]))
nets <- lapply(1:n_nets, function(x) set.network.attribute(nets[[x]], "age", ages[x]))

model_matrix <- get_model_matrix(ergm_formula, model_formula)
model_matrix[ ,1] <- model_matrix[ ,1] - model_matrix[ ,2]
true_mu <- c(-3, -2)
true_theta <- model_matrix %*% true_mu

for (n in 1:n_nets) {
  nets[[n]] <- simulate(y ~ edges, coef = true_theta[n],
                        basis = network(n_nodes, directed=FALSE))
}

nets <- lapply(1:n_nets, function(x) set.network.attribute(nets[[x]], "group", groups[x]))
nets <- lapply(1:n_nets, function(x) set.network.attribute(nets[[x]], "age", ages[x]))
control <- control_multibergm(ergm_formula, model_matrix, ~.)


prior <- set_priors(ergm_formula,
                    model_matrix,
                    control)

#init <- set_init(ergm_formula, prior, model_matrix)

fit <- multibergm(ergm_formula, model_matrix = model_matrix, main_iters = 1000)

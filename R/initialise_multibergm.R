# Set initial values ==========================================================
#' @importFrom mvtnorm rmvnorm
#' @importFrom MCMCpack riwish

set_init <- function(formula, prior, model_matrix, init = list()){
  
  networks <- statnet.common::eval_lhs.formula(formula)
  if (is.null(groups))
    groups <- rep(1, length(networks))
  
  n_networks <- length(networks)
  n_terms    <- length(attr(terms(formula), "term.labels"))
  
  # Use default if not specified explicitly
  if (is.null(init$mu)) {
    init$mu <- prior$mu$mean
    for (q in 1:ncol(prior$mu$mean))
      init$mu[,q] <- init$mu[,q] + rmvnorm(1, sigma = as.matrix(prior$mu$cov))
  }
  
  if (is.null(init$cov_theta))
    init$cov_theta   <- riwish(prior$cov_theta$df, prior$cov_theta$scale)
  
  if (is.null(init$theta))
    init$theta      <- rmvnorm(n_networks, sigma = init$cov_theta)
  
  if (!is.null(prior$cov_mu$df)) {
    if (is.null(init$cov_mu))
      init$cov_mu <- riwish(prior$cov_mu$df, prior$cov_mu$scale)
    
    if (is.null(init$muGroup))
      init$muGroup    <- rmvnorm(n_groups, init$mu, init$cov_mu)
  }

  mcmcr::as.mcmcr(init)
}

#' Initialise a multibergm fit
#' 
#' Initialise a multibergm fit by specifying the ERGM formula, group
#' memberships, prior hyperparameters, initial values and tuning parameters for 
#' the MCMC.
#' 
#' @export

initialise_multibergm <- function(formula,
                                  constraints = ~.,
                                  groups = NULL,
                                  prior = set_priors(formula, groups),
                                  init = set_init(formula, prior, groups),
                                  control = control_multibergm(formula,
                                                               constraints)) {
  
  networks <- statnet.common::eval_lhs.formula(formula)
  
  if (is.null(groups))
    groups <- rep(1, length(networks))
  
  out = list(formula = formula,
             constraints = constraints,
             groups = groups,
             prior = prior,
             params = init,
             main_iters = 1, 
             control = control,
             accepts = list(),
             time = as.difftime(0, units = "secs"))
  
  structure(out,
            class = "multibergm")
  
}

preallocate <- function(object, ...) {
  UseMethod("preallocate")
}

preallocate.mcmcarray <- function(x, iters) {
  
  empty <- array(NA, c(nchains(x), iters, pdims(x)))
  class(empty) <- "mcmcarray"
  
  bind_iterations(x, empty)
}

preallocate.mcmcr <- function(x, iters) {
  
  x <- lapply(x, preallocate, iters = iters)
  class(x) <- "mcmcr"
  
  x
}

set_class <- function(x, class) 
{
  class(x) <- class
  x
}
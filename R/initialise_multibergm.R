# Set initial values ==========================================================
#' @importFrom mvtnorm rmvnorm
#' @importFrom MCMCpack riwish
#'
#' @inheritParams multibergm

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
  
  # if (!is.null(prior$cov_mu$df)) {
  #   if (is.null(init$cov_mu))
  #     init$cov_mu <- riwish(prior$cov_mu$df, prior$cov_mu$scale)
  #   
  #   if (is.null(init$muGroup))
  #     init$muGroup    <- rmvnorm(n_groups, init$mu, init$cov_mu)
  # }

  mcmcr::as.mcmcr(init)
}

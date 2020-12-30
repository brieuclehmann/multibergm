# Set initial values ==========================================================
#' @importFrom mvtnorm rmvnorm
#' @importFrom MCMCpack riwish
#'
#' @inheritParams multibergm

set_init <- function(formula, prior, groups = NULL, init = list()) {

  networks <- statnet.common::eval_lhs.formula(formula)
  if (is.null(groups))
    groups <- rep(1, length(networks))

  n_networks <- length(networks)
  n_groups   <- length(unique(groups))

  # Use default if not specified explicitly
  if (is.null(init$mu_pop))
    init$mu_pop      <- rmvnorm(1, prior$mu_pop$mean, prior$mu_pop$cov)[1, ]

  if (is.null(init$cov_theta))
    init$cov_theta   <- riwish(prior$cov_theta$df,
                               prior$cov_theta$scale)

  if (is.null(init$theta))
    init$theta       <- rmvnorm(n_networks, sigma = init$cov_theta)

  if (n_groups > 1) {
    if (is.null(init$cov_mu_group))
      init$cov_mu_group <- riwish(prior$cov_mu_group$df,
                                  prior$cov_mu_group$scale)

    if (is.null(init$mu_group))
      init$mu_group    <- rmvnorm(n_groups, init$mu_pop, init$cov_mu_group)
  }

  mcmcr::as.mcmcr(init)
}

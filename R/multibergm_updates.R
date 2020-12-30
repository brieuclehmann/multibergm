#' MCMC update for single-group multi-BERGM
#'
#' Function used to perform one exchange-within-Gibbs MCMC update for
#' all the parameters in a single-group multi-BERGM. This applies a standard
#' Gibbs update for the network-level covariance parameter before using the
#' exchange algorithm within an Ancillarity-Sufficiency Interweaving Strategy
#' (ASIS) to update the remaining parameters.
#'
#' @inheritParams multibergm
#' @param curr A list of the current values of the model parameters
#' @param proposals A list of the current RW proposal parameters
#' @param control A list of parameters set by \code{\link{control_multibergm}}
#'  specifying priors, proposal variances, and group labels.
#'
#' @return A list of the updated values of the model parameters and the
#'   acceptance counts for the exchange updates.
#' @importFrom mvtnorm rmvnorm
singlegroup_update <- function(curr, prior, groups, proposals, control) {

  # Preallocate parameter values for next iteration
  nxt     <- curr
  accepts <- list()

  # First, update network-level covariance parameter
  nxt$cov_theta <- cov_update(curr$theta, prior$cov_theta$df,
                              prior$cov_theta$scale,
                              rep(0, dim(curr$theta)[2]))

  # Update network-level mean parameters in centered parameterisation
  theta_prop <- array(NA, dim(curr$theta))
  for (i in 1:nrow(curr$theta)) {
    this_proposal <- as.matrix(proposals$theta[i,,])
    theta_prop[i,] <- rmvnorm(1, curr$theta[i,], sigma = this_proposal)
  }
  coefs       <- sweep(theta_prop, 2, curr$mu_pop, "+")
  delta_theta <- ergm_wrapper(coefs, control)
  theta_mid   <- exchange_update(curr$theta, theta_prop,
                                 delta_theta, nxt$cov_theta)

  # Update population-level mean parameter in centered parameterisation
  mu_pop_curr <- sweep(theta_mid, 2, curr$mu_pop, "+")
  mu_pop_mid  <- mean_update(mu_pop_curr, nxt$cov_theta,
                             prior$mu_pop$mean, prior$mu_pop$cov)

  # Switch to non-centered parameterisation
  nxt$theta <- sweep(theta_mid, 2, curr$mu_pop - mu_pop_mid, "+")

  # Update population-level parameters in non-centered parameterisation
  this_proposal <- as.matrix(proposals$mu[1, , ])
  mu_pop_prop  <- rmvnorm(1, sigma = this_proposal) + mu_pop_mid
  coefs        <- sweep(nxt$theta, 2, mu_pop_prop, "+")
  delta_mu     <- ergm_wrapper(coefs, control)
  nxt$mu_pop   <- exchange_update(mu_pop_mid, mu_pop_prop, delta_mu,
                                  prior$mu_pop$cov,
                                  prior_mean = prior$mu_pop$mean,
                                  labels = groups)

  # Track acceptance counts
  accepts$theta <- apply(theta_mid != curr$theta, 1, any)
  accepts$mu    <- any(nxt$mu_pop != mu_pop_mid)

  list(params = nxt, accepts = accepts)
}


#' MCMC update for multiple-group multi-BERGM
#'
#' Function used to perform one exchange-within-Gibbs MCMC update for
#' all the parameters in a multiple-group multi-BERGM. This applies a standard
#' Gibbs update for the population-level parameters and the network-level
#' covariance parameter before using the exchange algorithm within an
#' Ancillarity-Sufficiency Interweaving Strategy (ASIS) to update the
#' remaining parameters.
#'
#' @inheritParams multibergm
#' @inheritParams singlegroup_update
#'
#' @return A list of the updated values of the model parameters and the
#'   acceptance counts for the exchange updates.
#' @importFrom mvtnorm rmvnorm
multigroup_update <- function(curr, prior, groups, proposals, control) {

  n_groups    <- dim(curr$mu_group)[1]

  # Preallocate parameter values for next iteration
  nxt     <- lapply(curr, function(x) array(NA, dim(x)))
  accepts <- list()

  # First, update covariance parameters and population-level mean parameter
  nxt$cov_mu_group     <- cov_update(curr$mu_group, prior$cov_mu_group$df,
                                     prior$cov_mu_group$scale,
                                     curr$mu_pop)

  nxt$mu_pop           <- mean_update(curr$mu_group, nxt$cov_mu_group,
                                      prior$mu_pop$mean, prior$mu_pop$cov)

  nxt$cov_theta        <- cov_update(curr$theta, prior$cov_theta$df,
                                     prior$cov_theta$scale,
                                     rep(0, dim(curr$theta)[2]))

  # Update network-level mean parameters in centered parameterisation
  theta_prop <- array(NA, dim(curr$theta))
  for (i in 1:nrow(curr$theta)) {
    this_proposal <- proposals$theta[i,,]
    theta_prop[i,] <- rmvnorm(1, curr$theta[i,], sigma = this_proposal)
  }
  coefs       <- array(NA, dim(theta_prop))
  for (n in 1:dim(coefs)[1])
    coefs[n, ] <- theta_prop[n, ] + curr$mu_group[groups[n], ]

  delta_theta <- ergm_wrapper(coefs, control)
  theta_mid   <- exchange_update(curr$theta, theta_prop,
                                 delta_theta, nxt$cov_theta)

  # Update group-level mean parameters in centered parameterisation
  mu_group_mid <- array(NA, dim(curr$mu_group))
  for (g in 1:n_groups) {
    nws                <- which(groups == g)
    mu_group_curr      <- sweep(theta_mid[nws, ], 2, curr$mu_group[g, ], "+")
    mu_group_mid[g, ]  <- mean_update(mu_group_curr, nxt$cov_theta,
                                      nxt$mu_pop, nxt$cov_mu_group)

    # Switch to non-centered parameterisation
    nxt$theta[nws, ] <- sweep(theta_mid[nws, ], 2,
                              curr$mu_group[g, ] - mu_group_mid[g, ], "+")
  }

  # Update group-level parameters in non-centered parameterisation
  mu_group_prop <- array(NA, dim(mu_group_mid))
  for (i in 1:n_groups) {
    this_proposal <- as.matrix(proposals$mu[i,,])
    mu_group_prop[i,] <- rmvnorm(1, mu_group_mid[i,], sigma = this_proposal)
  }
  coefs         <- array(NA, dim(nxt$theta))
  for (n in 1:dim(coefs)[1])
    coefs[n, ] <- nxt$theta[n, ] + mu_group_prop[groups[n], ]

  delta_mu     <- ergm_wrapper(coefs, control)
  nxt$mu_group <- exchange_update(mu_group_mid, mu_group_prop, delta_mu,
                                  nxt$cov_mu_group, prior_mean = nxt$mu_pop,
                                  labels = groups)

  # Track acceptance counts
  accepts$theta <- apply(theta_mid != curr$theta, 1, any)
  accepts$mu    <- apply(nxt$mu_group != mu_group_mid, 1, any)

  list(params = nxt, accepts = accepts)
}

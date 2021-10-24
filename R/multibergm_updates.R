#' MCMC update for single-group multi-BERGM
#'
#' Function used to perform one exchange-within-Gibbs MCMC update for
#' all the parameters in a single-group multi-BERGM. This applies a standard
#' Gibbs update for the network-level covariance parameter before using the
#' exchange algorithm within an Ancillarity-Sufficiency Interweaving Strategy
#' (ASIS) to update the remaining parameters.
#'
#' @param curr A list of the current values of the model parameters
#' @param control A list of parameters set by \code{\link{control_multibergm}}
#'  specifying priors, proposal variances, and group labels.
#'
#' @return A list of the updated values of the model parameters and the
#'   acceptance counts for the exchange updates.
#' @importFrom mvtnorm rmvnorm
single_var_update <- function(curr, prior, groups, control) {

  proposal   <- control$proposal

  # Preallocate parameter values for next iteration
  nxt     <- curr
  accepts <- list()

  # First, update network-level covariance parameter
  nxt$cov_theta <- cov_update(curr$theta,
                              prior$cov_theta$df,
                              prior$cov_theta$scale,
                              rep(0, dim(curr$theta)[2]))

  # Update network-level mean parameters in centered parameterisation
  theta_prop  <- rmvnorm(dim(curr$theta)[1], sigma=proposal$theta) + curr$theta
  coefs       <- sweep(theta_prop, 2, curr$mu, "+")
  delta_theta <- ergm_wrapper(coefs, control)
  theta_mid   <- exchange_update(curr$theta,
                                 theta_prop,
                                 delta_theta,
                                 nxt$cov_theta)

  # Update population-level mean parameter in centered parameterisation
  mu_curr <- sweep(theta_mid, 2, curr$mu, "+")
  mu_mid  <- mean_update(mu_curr, nxt$cov_theta, prior$mu$mean, prior$mu$cov)

  # Switch to non-centered parameterisation
  nxt$theta <- sweep(theta_mid, 2, curr$mu - mu_mid, "+")

  # Update population-level parameters in non-centered parameterisation
  mu_prop <- rmvnorm(1, sigma = proposal$mu) + mu_mid
  coefs     <- sweep(nxt$theta, 2, mu_prop, "+")
  delta_mu   <- ergm_wrapper(coefs, control)
  nxt$mu <- exchange_update(mu_mid,
                            mu_prop,
                            delta_mu,
                            prior$mu$cov,
                            prior_mean = prior$mu$mean,
                            labels = groups)
  
  # Track acceptance counts
  accepts$theta <- apply(theta_mid != curr$theta, 1, any)
  accepts$mu    <- any(nxt$mu != mu_mid)

  list(params = nxt, accepts = accepts)
}

#' MCMC update for multi-BERGM
#'
#' Function used to perform one exchange-within-Gibbs MCMC update for
#' all the parameters in a multi-BERGM. This applies a standard
#' Gibbs update for the network-level covariance parameter before using the
#' exchange algorithm within an Ancillarity-Sufficiency Interweaving Strategy
#' (ASIS) to update the remaining parameters.
#'
#' @param curr A list of the current values of the model parameters
#' @param control A list of parameters set by \code{\link{control_multibergm}}
#'  specifying priors, proposal variances, and group labels.
#'
#' @return A list of the updated values of the model parameters and the
#'   acceptance counts for the exchange updates.
#' @importFrom mvtnorm rmvnorm
asis_update <- function(curr, prior, model_matrix, control) {
  
  proposal   <- control$proposal
  
  # Preallocate parameter values for next iteration
  nxt     <- curr
  accepts <- list()
  
  # First, update network-level covariance parameter
  nxt$cov_theta <- cov_update(curr$theta,
                              array(0, dim(curr$theta)),
                              prior$cov_theta$df,
                              prior$cov_theta$scale)
  
  # Update network-level mean parameters in centered parameterisation
  theta_prop  <- rmvnorm(dim(curr$theta)[1], sigma=proposal$theta) + curr$theta
  coefs       <- theta_prop + (model_matrix %*% curr$mu)
  delta_theta <- ergm_wrapper(coefs, control)
  theta_mean_curr <- model_matrix %*% curr$mu
  theta_mid   <- exchange_update_cp(curr$theta,
                                    theta_prop,
                                    array(0, dim(curr$theta)),
                                    delta_theta,
                                    nxt$cov_theta)
  
  # Update population-level mean parameter in centered parameterisation
  #mu_curr <- theta_mid + theta_mean_curr
  mu_mid  <- mean_update(theta_mid,
                         nxt$cov_theta,
                         prior$mu$mean,
                         prior$mu$cov,
                         model_matrix)
  
  # Switch to non-centered parameterisation
  #nxt$theta <- sweep(theta_mid, 2, curr$mu - mu_mid, "+")
  nxt$theta <- theta_mid + (model_matrix %*% (curr$mu - mu_mid))
  
  # Update population-level parameters in non-centered parameterisation
  mu_prop <- mu_mid
  for (j in seq(nrow(curr$mu))) {
    mu_prop[j, ] <- mu_prop[j, ] + rmvnorm(1, sigma = as.matrix(proposal$mu[ , ,j]))
  }
  coefs <- nxt$theta + (model_matrix %*% mu_prop)
  delta_mu   <- ergm_wrapper(coefs, control)
  nxt$mu <- exchange_update_ncp(mu_mid,
                                mu_prop,
                                delta_mu,
                                nxt$cov_theta,
                                prior$mu$cov,
                                prior$mu$mean,
                                model_matrix)

  #mu_prop <- rmvnorm(nrow(mu_mid), sigma = proposal$mu) + mu_mid
  #coefs     <- sweep(nxt$theta, 2, mu_prop, "+")
  
  
  # Track acceptance counts
  accepts$theta <- apply(theta_mid != curr$theta, 1, any)
  accepts$mu    <- any(nxt$mu != mu_mid)
  
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
#' @param curr A list of the current values of the model parameters
#' @param control A list of parameters set by \code{\link{control_multibergm}}
#'  specifying priors, proposal variances, and group labels.
#'
#' @return A list of the updated values of the model parameters and the
#'   acceptance counts for the exchange updates.

#' @importFrom mvtnorm rmvnorm
multi_var_update <- function(curr, prior, groups, control){

  proposal   <- control$proposal
  n_groups    <- dim(curr$muGroup)[1]

  # Preallocate parameter values for next iteration
  nxt     <- lapply(curr, function(x) array(NA, dim(x)))
  accepts <- list()

  # First, update covariance parameters and population-level mean parameter
  nxt$covMuGroup      <- cov_update(curr$muGroup, prior$covMuGroup$df,
                                    prior$covMuGroup$scale,
                                    curr$mu)

  nxt$mu           <- mean_update(curr$muGroup, nxt$covMuGroup,
                                     prior$mu$mean, prior$mu$cov)

  nxt$cov_theta        <- cov_update(curr$theta, prior$cov_theta$df,
                                    prior$cov_theta$scale,
                                    rep(0, dim(curr$theta)[2]))

  # Update network-level mean parameters in centered parameterisation
  theta_prop  <- rmvnorm(dim(curr$theta)[1], sigma=proposal$theta) + curr$theta
  coefs      <- array(NA, dim(theta_prop))
  for (n in 1:dim(coefs)[1])
    coefs[n, ] <- theta_prop[n, ] + curr$muGroup[groups[n], ]

  delta_theta <- ergm_wrapper(coefs, control)
  theta_mid   <- exchange_update(curr$theta, theta_prop, delta_theta, nxt$cov_theta)

  # Update group-level mean parameters in centered parameterisation
  muGroupMid <- array(NA, dim(curr$muGroup))
  for (g in 1:n_groups){
    nws              <- which(groups == g)
    muGroupCurr      <- sweep(theta_mid[nws, ], 2, curr$muGroup[g, ], "+")
    muGroupMid[g, ]  <- mean_update(muGroupCurr, nxt$cov_theta,
                                    nxt$mu, nxt$covMuGroup)

    # Switch to non-centered parameterisation
    nxt$theta[nws, ] <- sweep(theta_mid[nws, ], 2,
                              curr$muGroup[g, ] - muGroupMid[g, ], "+")
  }

  # Update group-level parameters in non-centered parameterisation
  muGroupProp <- rmvnorm(n_groups, sigma = proposal$mu) + muGroupMid
  coefs       <- array(NA, dim(nxt$theta))
  for (n in 1:dim(coefs)[1])
    coefs[n, ] <- nxt$theta[n, ] + muGroupProp[groups[n], ]

  delta_mu     <- ergm_wrapper(coefs, control)
  nxt$muGroup <- exchange_update(muGroupMid, muGroupProp, delta_mu,
                                 nxt$covMuGroup, prior_mean = nxt$mu,
                                 labels = groups)

  # Track acceptance counts
  accepts$theta <- apply(theta_mid != curr$theta, 1, any)
  accepts$mu    <- apply(nxt$muGroup != muGroupMid, 1, any)

  list(params = nxt, accepts = accepts)
}

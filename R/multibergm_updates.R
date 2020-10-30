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
singleGroup_update <- function(curr, prior, groups, control){

  proposal   <- control$proposal

  # Preallocate parameter values for next iteration
  nxt     <- curr
  accepts <- list()

  # First, update network-level covariance parameter
  nxt$covTheta <- cov_update(curr$theta, prior$covTheta$df,
                             prior$covTheta$scale,
                             rep(0, dim(curr$theta)[2]))

  # Update network-level mean parameters in centered parameterisation
  thetaProp  <- rmvnorm(dim(curr$theta)[1], sigma=proposal$theta) + curr$theta
  coefs      <- sweep(thetaProp, 2, curr$muPop, "+")
  deltaTheta <- ergm_wrapper(coefs, control)
  thetaMid   <- exchange_update(curr$theta, thetaProp, deltaTheta, nxt$covTheta)

  # Update population-level mean parameter in centered parameterisation
  muPopCurr <- sweep(thetaMid, 2, curr$muPop, "+")
  muPopMid  <- mean_update(muPopCurr, nxt$covTheta,
                           prior$muPop$mean, prior$muPop$cov)

  # Switch to non-centered parameterisation
  nxt$theta <- sweep(thetaMid, 2, curr$muPop - muPopMid, "+")

  # Update population-level parameters in non-centered parameterisation
  muPopProp <- rmvnorm(1, sigma = proposal$mu) + muPopMid
  coefs     <- sweep(nxt$theta, 2, muPopProp, "+")
  deltaMu   <- ergm_wrapper(coefs, control)
  nxt$muPop <- exchange_update(muPopMid, muPopProp, deltaMu,
                               prior$muPop$cov,
                               prior_mean = prior$muPop$mean,
                               labels = groups)

  # Track acceptance counts
  accepts$theta <- apply(thetaMid != curr$theta, 1, any)
  accepts$mu    <- any(nxt$muPop != muPopMid)

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
multiGroup_update <- function(curr, prior, groups, control){

  proposal   <- control$proposal
  n_groups    <- dim(curr$muGroup)[1]

  # Preallocate parameter values for next iteration
  nxt     <- lapply(curr, function(x) array(NA, dim(x)))
  accepts <- list()

  # First, update covariance parameters and population-level mean parameter
  nxt$covMuGroup      <- cov_update(curr$muGroup, prior$covMuGroup$df,
                                    prior$covMuGroup$scale,
                                    curr$muPop)

  nxt$muPop           <- mean_update(curr$muGroup, nxt$covMuGroup,
                                     prior$muPop$mean, prior$muPop$cov)

  nxt$covTheta        <- cov_update(curr$theta, prior$covTheta$df,
                                    prior$covTheta$scale,
                                    rep(0, dim(curr$theta)[2]))

  # Update network-level mean parameters in centered parameterisation
  thetaProp  <- rmvnorm(dim(curr$theta)[1], sigma=proposal$theta) + curr$theta
  coefs      <- array(NA, dim(thetaProp))
  for (n in 1:dim(coefs)[1])
    coefs[n, ] <- thetaProp[n, ] + curr$muGroup[groups[n], ]

  deltaTheta <- ergm_wrapper(coefs, control)
  thetaMid   <- exchange_update(curr$theta, thetaProp, deltaTheta, nxt$covTheta)

  # Update group-level mean parameters in centered parameterisation
  muGroupMid <- array(NA, dim(curr$muGroup))
  for (g in 1:n_groups){
    nws              <- which(groups == g)
    muGroupCurr      <- sweep(thetaMid[nws, ], 2, curr$muGroup[g, ], "+")
    muGroupMid[g, ]  <- mean_update(muGroupCurr, nxt$covTheta,
                                    nxt$muPop, nxt$covMuGroup)

    # Switch to non-centered parameterisation
    nxt$theta[nws, ] <- sweep(thetaMid[nws, ], 2,
                              curr$muGroup[g, ] - muGroupMid[g, ], "+")
  }

  # Update group-level parameters in non-centered parameterisation
  muGroupProp <- rmvnorm(n_groups, sigma = proposal$mu) + muGroupMid
  coefs       <- array(NA, dim(nxt$theta))
  for (n in 1:dim(coefs)[1])
    coefs[n, ] <- nxt$theta[n, ] + muGroupProp[groups[n], ]

  deltaMu     <- ergm_wrapper(coefs, control)
  nxt$muGroup <- exchange_update(muGroupMid, muGroupProp, deltaMu,
                                 nxt$covMuGroup, prior_mean = nxt$muPop,
                                 labels = groups)

  # Track acceptance counts
  accepts$theta <- apply(thetaMid != curr$theta, 1, any)
  accepts$mu    <- apply(nxt$muGroup != muGroupMid, 1, any)

  list(params = nxt, accepts = accepts)
}

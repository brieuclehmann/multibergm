#' Internal functions for MCMC updates
#'
#' Functions used internally by to perform Gibbs or exchange updates.
#'
#' @name mcmc_updates
#'
#' @param curr A vector or matrix of current mean parameter values.
#' @param prop A vector or matrix of proposed mean parameter values
#' @param delta Change in summary statistics as produced by
#'   \code{\link{ergm_wrapper}}
#' @param priorMean Prior mean of parameter
#' @param priorCov Prior covariance of parameter
#' @param netLabels Labels specifying which networks to associate with
#'   each row of the the parameter matrices
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats runif

exchange_update <- function(curr, 
                            prop, 
                            delta, 
                            prior_cov,
                            prior_mean = double(ncol(curr)), 
                            labels = seq_len(nrow(curr))) {

  if (all(curr == prop)) {
    return(curr)
  }

  if (is.vector(curr)) curr <- matrix(curr, nrow=1)
  if (is.vector(prop)) prop <- matrix(prop, nrow=1)

  n   <- nrow(curr)
  new <- curr
  
  for (i in seq_len(n)) {
    pr_prop <- dmvnorm(prop[i, ], prior_mean, prior_cov, log = TRUE)
    pr_curr <- dmvnorm(curr[i, ], prior_mean, prior_cov, log = TRUE)

    this_delta <- colSums(delta[which(labels == i), , drop = FALSE])
    beta      <- sum((curr[i, ] - prop[i, ]) * this_delta) + pr_prop - pr_curr

    if (beta >= log(runif(1)))
      new[i, ] <- prop[i, ]
  }

  new
}

#' @param obsData Matrix of observed data
#' @param obsCov Fixed (known) covariance of observed data
#'
#' @describeIn mcmc_updates Gibbs update of a mean parameter
#' @importFrom mvtnorm rmvnorm

mean_update <- function(obs_data, obs_cov, prior_mean, prior_cov) {

  if (all(prior_cov == 0)) {
    
    post_cov <- prior_cov
    post_mean <- prior_mean
    
  } else {
    
    n <- dim(obs_data)[1]
    
    post_cov  <- solve(solve(prior_cov) + n*solve(obs_cov))
    post_mean <- post_cov %*% ((solve(prior_cov) %*% prior_mean) +
                               (n*(solve(obs_cov) %*% colMeans(obs_data))))
    
  }

  rmvnorm(1, post_mean, post_cov)[1, ]
}

#' @param priorDf Prior degrees of freedom in Inverse-Wishart
#' @param priorScale Prior scale matrix in Inverse-Wishart
#' @param obsMean Fixed (known) mean of observed data
#' @param dataLabels Labels to associate each observation with a grouping
#' @param currCov Current value of covariance parameter
#'
#' @describeIn mcmc_updates Gibbs update of a covariance parameter
#'
#' @importFrom MCMCpack riwish

cov_update <- function(obsData, priorDf, priorScale, obsMean,
                       dataLabels = NULL, currCov = NULL) {

  # Return current covariance if this level of hierarchy is fixed
  if (all(priorScale == 0))
    return(currCov)

  # Else perform Gibbs update from conditional posterior (inverse-Wishart)
  nObs <- dim(obsData)[1]
  if (is.null(dataLabels))
    dataLabels <- rep(1, nObs)

  if (!is.matrix(obsMean))
    obsMean <- matrix(obsMean, nrow = 1)

  nGroups <- dim(obsMean)[1]
  pairwiseDev <- 0
  for (g in 1:nGroups) {
    ind <- which(dataLabels == g)
    pairwiseDev <- pairwiseDev + crossprod(sweep(obsData[ind, , drop = F],
                                                 2, obsMean[g, ]))
  }

  postDf    <- priorDf + nObs
  postScale <- priorScale + pairwiseDev

  riwish(postDf, postScale)
}

#' @param coefs Matrix of coefficients to simulate ERGs from
#' @param control A list of parameters set by \code{\link{control_multibergm}}
#'   containing settings for the ERGM simulations
#' @describeIn mcmc_updates Gibbs update of a mean parameter
#' @import ergm
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
ergm_wrapper <- function(coefs, control) {

  n_nets <- dim(coefs)[1]
  seeds  <- rngtools::RNGseq(n_nets)

  # Parallel call to ergm_MCMC_slave
  delta <- foreach(n = seq_len(n_nets), 
                   r = seeds,
                   .combine = rbind, 
                   .packages = "ergm") %dopar% {

                     rngtools::setRNG(r)
                     ergm_MCMC_slave(control$Clists[[n]],
                                     control$MHproposals,
                                     coefs[n, ],
                                     control.simulate.formula(),
                                     burnin = control$auxIters,
                                     samplesize = 1,
                                     interval = 1,
                                     verbose = FALSE)$s
                   }

  return(delta)
}


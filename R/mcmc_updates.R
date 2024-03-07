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
#' @param prior_mean Prior mean of parameter
#' @param prior_cov Prior covariance of parameter
#' @param etamap The list of values that constitutes the theta-> eta mapping and
#'   is returned by ergm.etamap
#' @param labels Labels specifying which networks to associate with
#'   each row of the the parameter matrices
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom matrixNormal dmatnorm
#' @importFrom stats runif

exchange_update_cp <- function(curr, 
                               prop,
                               prior_mean,
                               delta, 
                               prior_cov,
                               etamap) {
  
  if (all(curr == prop)) {
    return(curr)
  }

  d <- dim(curr)

  if (is.vector(curr)) curr <- matrix(curr, nrow = 1)
  if (is.vector(prop)) prop <- matrix(prop, nrow = 1)

  n   <- nrow(curr)
  new <- curr

  # eta-theta mapping for curved ERGMs
  if (length(etamap$curved) > 0) {
    curr_eta <- t(apply(curr, 1, ergm.eta, etamap))
    prop_eta <- t(apply(prop, 1, ergm.eta, etamap))
  } else {
    curr_eta <- curr
    prop_eta <- prop
  }

  
  for (i in seq_len(n)) {
    #pr_prop <- dmvnorm(prop[i, ], prior_mean, prior_cov, log = TRUE)
    #pr_curr <- dmvnorm(curr[i, ], prior_mean, prior_cov, log = TRUE)
    pr_prop <- dmvnorm(prop[i, ], prior_mean[i, ], prior_cov, log = TRUE)
    pr_curr <- dmvnorm(curr[i, ], prior_mean[i, ], prior_cov, log = TRUE)

    #this_delta <- colSums(delta[which(labels == i), , drop = FALSE])
    beta      <- sum((curr[i, ] - prop[i, ]) * delta[i, , drop = FALSE]) + pr_prop - pr_curr

    if (beta >= log(runif(1)))
      new[i, ] <- prop[i, ]
  }
  dim(new) <- d
  new
}

exchange_update_ncp <- function(curr, 
                                prop, 
                                delta,
                                obs_cov,
                                prior_cov,
                                prior_mean, 
                                model_matrix) {
  
  if (all(curr == prop)) {
    return(curr)
  }
  
  d <- dim(curr)
  
  if (is.vector(curr)) curr <- matrix(curr, nrow = 1)
  if (is.vector(prop)) prop <- matrix(prop, nrow = 1)
  
  n   <- nrow(model_matrix)
  new <- curr
  
  pr_prop <- dmatnorm(prop, prior_mean, prior_cov, obs_cov, log = TRUE)
  pr_curr <- dmatnorm(curr, prior_mean, prior_cov, obs_cov, log = TRUE)
  
  beta <- pr_prop - pr_curr
  theta_diff <- model_matrix %*% (curr - prop)
  for (i in seq_len(n)) {
    #pr_prop <- dmvnorm(prop[i, ], prior_mean, prior_cov, log = TRUE)
    #pr_curr <- dmvnorm(curr[i, ], prior_mean, prior_cov, log = TRUE)
    #beta      <- sum((curr[i, ] - prop[i, ]) * this_delta[i, ]) + pr_prop - pr_curr
    beta <- beta + sum(theta_diff[i, ] * delta[i, ])
  }
  if (beta >= log(runif(1)))
    new <- prop
  dim(new) <- d
  new
}

#' @param obs_data Matrix of observed data
#' @param obs_cov Fixed (known) covariance of observed data
#'
#' @describeIn mcmc_updates Gibbs update of a mean parameter
#' @importFrom mvtnorm rmvnorm
#' @importFrom matrixNormal rmatnorm

mean_update <- function(obs_data, obs_cov, prior_mean, prior_cov, model_matrix) {

  if (all(prior_cov == 0)) {

    post_cov  <- prior_cov
    post_mean <- prior_mean

  } else {

    n <- dim(obs_data)[1]
    
    prior_prec <- solve(prior_cov)
    U_inv <- prior_prec + (t(model_matrix) %*% model_matrix)
    U <- solve(U_inv)
    V <- obs_cov
    
    M <- U %*% ((t(model_matrix) %*% obs_data) + (prior_prec %*% prior_mean))
    
    return(rmatnorm(1, M, U, V))
    
    data_prec <- apply(model_matrix, 1,  
                       function(x) t(x) %*% solve(obs_cov) %*% x) 
    
    post_cov  <- solve(solve(prior_cov) + n * solve(obs_cov))
    post_mean <- post_cov %*% ((solve(prior_cov) %*% prior_mean) +
                               (n * (solve(obs_cov) %*% colMeans(obs_data))))
    
  }

  rmvnorm(1, post_mean, post_cov)[1, ]
}

#' @param prior_df Prior degrees of freedom in Inverse-Wishart
#' @param prior_scale Prior scale matrix in Inverse-Wishart
#' @param var_mean Fixed (known) mean of parameters
#' @param curr_cov Current value of covariance parameter
#'
#' @describeIn mcmc_updates Gibbs update of a covariance parameter
#'
#' @importFrom MCMCpack riwish

cov_update <- function(obs_data, # Y - XB
                       obs_mean, 
                       prior_df, # nu_0
                       prior_scale,
                       curr_cov = NULL) {

  #if (all(prior_obs_scale == 0)) {
  #  return(curr_cov)
  #}

  n <- nrow(obs_data)

  if (!is.matrix(obs_mean))
    obs_mean <- matrix(obs_mean, nrow = 1)
  
  data_centered <- obs_data - obs_mean
  
  pairwise_dev <- t(data_centered) %*% data_centered
  
  #pairwise_dev <- 0
  #for (g in seq_len(nrow(obs_mean))) {
  #  ind <- which(labels == g)
  #  pairwise_dev <- pairwise_dev + crossprod(sweep(obs_data[ind, , drop = F],
  #                                               2, obs_mean[g, ]))
  #}

  post_df    <- prior_df + n
  post_scale <- prior_scale + pairwise_dev

  riwish(post_df, post_scale)
}

#' @param prior_df Prior degrees of freedom in Inverse-Wishart
#' @param prior_scale Prior scale matrix in Inverse-Wishart
#' @param var_mean Fixed (known) mean of parameters
#' @param curr_cov Current value of covariance parameter
#'
#' @describeIn mcmc_updates Gibbs update of a covariance parameter
#'
#' @importFrom MCMCpack riwish

cov_mv_update <- function(obs_data, # Y - XB
                          obs_mean, 
                          prior_df, # nu_0
                          prior_scale,
                          curr_mu,
                          mu_prior_mean,
                          mu_prior_scale) {
  
  #if (all(prior_obs_scale == 0)) {
  #  return(curr_cov)
  #}
  
  n <- nrow(obs_data)
  
  if (!is.matrix(obs_mean))
    obs_mean <- matrix(obs_mean, nrow = 1)
  
  data_centered <- obs_data - obs_mean
  mu_centered <- curr_mu - mu_prior_mean
  
  pairwise_dev <- t(data_centered) %*% data_centered
  mu_dev <- t(mu_centered) %*% mu_prior_scale %*% mu_centered
  
  #pairwise_dev <- 0
  #for (g in seq_len(nrow(obs_mean))) {
  #  ind <- which(labels == g)
  #  pairwise_dev <- pairwise_dev + crossprod(sweep(obs_data[ind, , drop = F],
  #                                               2, obs_mean[g, ]))
  #}
  
  post_df    <- prior_df + n
  post_scale <- prior_scale + pairwise_dev + mu_dev
  
  riwish(post_df, post_scale)
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
  r <- n <- NULL
  delta <- foreach(n = seq_len(n_nets),
                   r = seeds,
                   .combine = rbind,
                   .packages = "ergm") %dopar% {

                     rngtools::setRNG(r)
                     ergm_MCMC_slave(control$clists[[n]],
                                     control$mh_proposals,
                                     coefs[n, ],
                                     control.simulate.formula(),
                                     burnin = control$aux_iters,
                                     samplesize = 1,
                                     interval = 1,
                                     verbose = FALSE)$s
                   }

  return(delta)
}

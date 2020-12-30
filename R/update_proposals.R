#' Update Metropolis random-walk covariance parameters
#' 
#' This is an internal function used to update the covariance parameters for the
#' mu_pop and theta proposals. 
#' 
#' @param proposals Current set of proposals
#' @param accept_rate Acceptance rates since the last proposal adaptation
#' @param params Posterior samples
#' @param control A list of parameters set by \code{\link{control_multibergm}}
#' 
#' @importFrom stats cov

update_proposals <- function(proposals, accept_rate, params, control) {
  
  beta <- 0.05
  delta_scale <- 100 * min(0.01, 1 / sqrt(mcmcr::niters(params$theta)))
  # Update theta proposals
  n_networks <- dim(proposals$theta)[1]
  n_stats    <- dim(proposals$theta)[2]
  for (i in 1:n_networks) {
    if (accept_rate$theta[i] > 0.3) {
      proposals$theta_scale[i] <- proposals$theta_scale[i] + delta_scale
    } else {
      proposals$theta_scale[i] <- proposals$theta_scale[i] - delta_scale
    }
    this_sample <- as.matrix(params$theta[1, ,i, ])
    posterior_cov <- (((1 - beta) * 2.38) ^ 2) * cov(this_sample) / n_stats
    regular_cov <- ((beta * 0.1) ^ 2) * diag(n_stats) / n_stats
    proposals$theta[i,,] <- exp(proposals$theta_scale[i]) * (posterior_cov + regular_cov)
  }
  
  # Update mu proposals
  n_groups <- dim(proposals$mu)[1]
  if (n_groups == 1) {
    if (accept_rate$mu[1] > 0.3) {
      proposals$mu_scale[1] <- proposals$mu_scale[1] + delta_scale
    } else {
      proposals$mu_scale[1] <- proposals$mu_scale[1] - delta_scale
    }
    this_sample <- as.matrix(params$mu[1, , ])
    posterior_cov <- (((1 - beta) * 2.38) ^ 2) * cov(this_sample) / n_stats
    regular_cov <- ((beta * 0.1) ^ 2) * diag(n_stats) / n_stats
    proposals$mu[1,,] <- exp(proposals$mu_scale[1]) * (posterior_cov + regular_cov)
  } else {
    for (i in 1:n_groups) {
      if (accept_rate$mu[i] > 0.3) {
        proposals$mu_scale[i] <- proposals$mu_scale[i] - delta_scale
      } else {
        proposals$mu_scale[i] <- proposals$mu_scale[i] + delta_scale
      }
      this_sample <- as.matrix(params$mu[1, ,i, ])
      posterior_cov <- (((1 - beta) * 2.38) ^ 2) * cov(this_sample) / n_stats
      regular_cov <- ((beta * 0.1) ^ 2) * diag(n_stats) / n_stats
      proposals$mu[i,,] <- exp(proposals$mu_scale[i]) * (posterior_cov + regular_cov)
    }
  }

  proposals
}
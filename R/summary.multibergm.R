#' Summarising multibergm model fits
#'
#' \link{summary} method for \link{multibergm} fits
#'
#' @param object A \link{multibergm} object
#' @param param Multibergm parameter to be summarised
#' @param thin Amount of thinning to apply to the posterior samples
#' @param burn_in Amount of burn-in to remove from the start of posterior
#'   samples (pre-thinning)
#' @param ... Additional parameters to be passed on to lower-level functions.
#'
#' @return The function computes and prints posterior means and quantiles for
#'   the specified parameter, as well as the acceptance rates for parameters
#'   updated via the exchange algorithm.
#'
#' @export
summary.multibergm <- function(object,
                               param = "mu",
                               thin = 1L,
                               burn_in = 0L,
                               ...) {
  n_groups <- length(unique(object$groups))
  thin    <- as.integer(thin)
  burn_in <- as.integer(burn_in)

  # Remove burn_in and apply thinning
  post_iters      <- seq(burn_in + 1L, object$main_iters, thin)
  output  <- get(param, object$params)
  output  <- subset(output, iterations = post_iters)
  
  model_vars <- colnames(control$mod_mat)
  n_vars <- length(model_vars)
  ergm_terms <- param_names(object$control$model)
  n_terms     <- nparam(object$control$model)
  n_groups   <- length(unique(object$groups))
  
  cat("\n", "Posterior Density Estimate for Model: \ny ~",
      paste(object$formula[3]), "\n", "\n")

  if (param == "mu") {
    FFmu   <- mcmcr::as.mcmc(output, start = burn_in + 1, thin = thin)

    table1 <- summary(FFmu)$statistics
    #rnames <- paste0("mu", seq_len(n_terms), " (", ergm_terms, ")")
    rnames <- paste("mu", 
                    rep(model_vars, n_terms), 
                    rep(ergm_terms, each = n_vars), sep = "_")
    table1 <- matrix(table1, n_terms * n_vars,
                     dimnames = list(rnames, colnames(table1)))

    table2 <- summary(FFmu)$quantiles
    table2 <- matrix(table2, n_terms * n_vars,
                     dimnames = list(rnames, colnames(table2)))

    print(table1, digits = 4)
    cat("\n")
    print(table2, digits = 4)

  } else {
    
    ff_mu   <- mcmcr::as.mcmc(output, start = burn_in + 1, thin = thin)
    
    table1 <- summary(ff_mu)$statistics
    rnames <- paste0("mu", rep(seq_len(n_terms), each = n_groups),
                     " (", rep(model_terms, each = n_groups), 
                     ", G", rep(1:n_groups, n_terms), ")")
    table1 <- matrix(table1, n_terms * n_groups,
                     dimnames = list(rnames, colnames(table1)))
    
    table2 <- summary(ff_mu)$quantiles
    table2 <- matrix(table2, n_terms * n_groups,
                     dimnames = list(rnames, colnames(table2)))
    
    print(table1, digits = 4)
    cat("\n")
    print(table2, digits = 4)
  }

  # Print acceptance rates
  cat("\n Theta acceptance rate: \n")
  theta_ar <- colSums(object$accepts$theta) / (length(post_iters) - 1)
  cat(formatC(mean(theta_ar), digits = 3, format = "f"),
      " (", formatC(min(theta_ar), digits = 3, format = "f"), ", ",
      formatC(max(theta_ar), digits = 3, format = "f"), ")", sep = "")

  cat("\n Mu acceptance rate: \n")
  mu_ar <- colSums(object$accepts$mu) / (length(post_iters) - 1)
  cat(formatC(mean(mu_ar), digits = 3, format = "f"),
      " (", formatC(min(mu_ar), digits = 3, format = "f"), ", ",
      formatC(max(mu_ar), digits = 3, format = "f"), ") \n", sep = "")

}

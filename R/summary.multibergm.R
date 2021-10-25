#' Summarising multibergm model fits
#'
#' \link{summary} method for \link{multibergm} fits
#'
#' @param object A \link{multibergm} object
#' @param param Multibergm parameter to be summarised
#' @param thin Amount of thinning to apply to the posterior samples
#' @param burn_in Amount of burn-in to remove from the start of posterior
#'   samples (pre-thinning)
#'
#' @return The function computes and prints posterior means and quantiles for
#'   the specified parameter, as well as the acceptance rates for parameters
#'   updated via the exchange algorithm.
#'
#' @export
summary.multibergm <- function(object,
                               param = "mu",
                               thin = 1L,
                               burn_in = 0L) {
  thin    <- as.integer(thin)
  burn_in <- as.integer(burn_in)
  
  # Remove burn_in and apply thinning
  post_iters      <- seq(burn_in + 1L, object$main_iters, thin)
  output  <- get(param, object$params)
  output  <- subset(output, iterations = post_iters)
  
  model_vars <- colnames(control$mod_mat)
  n_vars <- length(model_vars)
  ergm_terms <- object$control$model$coef.names
  n_terms     <- length(attr(terms(object$formula), "term.labels"))
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
    for (g in 1:object$init$n_groups) {
      FFmu   <- coda::mcmc(output, start = burn_in + 1, thin = thin)

      table1 <- summary(FFmu)$statistics
      rnames <- paste0("mu", seq_len(n_terms), " (", model_terms, ", G",
                       rep(g, each = n_terms), ")")
      table1 <- matrix(table1, n_terms, 
                       dimnames = list(rnames, colnames(table1)))

      table2 <- summary(FFmu)$quantiles
      table2 <- matrix(table2, n_terms,
                       dimnames = list(rnames, colnames(table2)))

      print(table1, digits = 4)
      print("\n")
      print(table2, digits = 4)
    }
  }

  # Print acceptance rates
  cat("\n Theta acceptance rate: \n")
  thetaAR <- colSums(object$accepts$theta)/(object$main_iters - 1)
  cat(formatC(mean(thetaAR), digits = 3, format = "f"),
      " (", formatC(min(thetaAR), digits = 3, format = "f"), ", ",
      formatC(max(thetaAR), digits = 3, format = "f"), ")", sep="")

  cat("\n Mu acceptance rate: \n")
  muAR <- colSums(object$accepts$mu)/(object$main_iters - 1)
  cat(formatC(mean(muAR), digits = 3, format = "f"),
      " (", formatC(min(muAR), digits = 3, format = "f"), ", ",
      formatC(max(muAR), digits = 3, format = "f"), ") \n", sep="")

}

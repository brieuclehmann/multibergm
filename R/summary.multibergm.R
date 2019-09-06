#' Summarising multibergm model fits
#'
#' \link{summary} method for \link{multibergm} fits
#'
#' @param object A \link{multibergm} object
#' @param param Multibergm parameter to be summarised
#' @param thin Amount of thinning to apply to the posterior samples
#' @param burnIn Amount of burn-in to remove from the start of posterior
#'   samples (pre-thinning)
#'
#' @return The function computes and prints posterior means and quantiles for
#'   the specified parameter, as well as the acceptance rates for parameters
#'   updated via the exchange algorithm.
#'
#' @export
summary.multibergm <- function(object,
                               param = "muPop",
                               thin = 1,
                               burnIn = 0) {

  # Remove burnIn iterations and apply thinning (default: no thinning)
  postIters      <- seq(burnIn + 1, object$mainIters, thin)
  object$params  <- lapply(object$params,
                          function(x) abind::asub(x, postIters, 1))

  modelTerms <- object$control$model$coef.names
  nTerms     <- length(attr(terms(object$formula), "term.labels"))
  output     <- get(param, object$params)

  cat("\n", "Posterior Density Estimate for Model: \ny ~", paste(object$formula[3]),
      "\n", "\n")

  if (param == "muPop") {
    FFmu   <- coda::mcmc(output, start=burnIn+1, thin = thin)

    table1 <- summary(FFmu)$statistics
    rnames <- paste0("mu", 1:nTerms, " (", modelTerms, ")")
    table1 <- matrix(table1, nTerms,
                     dimnames = list(rnames, colnames(table1)))

    table2 <- summary(FFmu)$quantiles
    table2 <- matrix(table2, nTerms,
                     dimnames = list(rnames, colnames(table2)))

    print(table1, digits = 4)
    cat("\n")
    print(table2, digits = 4)

  } else {
    for (g in 1:object$init$nGroups) {
      FFmu   <- coda::mcmc(output, start=burnIn+1, thin = thin)

      table1 <- summary(FFmu)$statistics
      rnames <- paste0("mu", 1:nTerms, " (", modelTerms, ", G",
                       rep(g, each = nTerms), ")")
      table1 <- matrix(table1, object$init$nStats,
                       dimnames = list(rnames, colnames(table1)))

      table2 <- summary(FFmu)$quantiles
      table2 <- matrix(table2, object$init$nStats,
                       dimnames = list(rnames, colnames(table2)))

      print(table1, digits = 4)
      print("\n")
      print(table2, digits = 4)
    }
  }

  # Print acceptance rates
  cat("\n Theta acceptance rate: \n")
  thetaAR <- colSums(object$accepts$theta)/(object$mainIters - 1)
  cat(formatC(mean(thetaAR), digits = 3, format = "f"),
      " (", formatC(min(thetaAR), digits = 3, format = "f"), ", ",
      formatC(max(thetaAR), digits = 3, format = "f"), ")", sep="")

  cat("\n Mu acceptance rate: \n")
  muAR <- colSums(object$accepts$mu)/(object$mainIters - 1)
  cat(formatC(mean(muAR), digits = 3, format = "f"),
      " (", formatC(min(muAR), digits = 3, format = "f"), ", ",
      formatC(max(muAR), digits = 3, format = "f"), ") \n", sep="")

}

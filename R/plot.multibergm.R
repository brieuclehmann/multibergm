#' Plot BERGM posterior output
#'
#' This function creates MCMC diagnostic plots for \code{multibergm} objects.
#'
#' @inheritParams summary.multibergm
#' @importFrom ggplot2 ggplot
#'
#' @return Outputs a density plot, trace plot, and autocorrelation plot for
#'   each of the model parameters.
#'
#' @export
plot.multibergm <- function(object,
                            param = "muPop",
                            burnIn = 0,
                            thin = 1){

  # Remove burnIn iterations and apply thinning (default: no thinning)
  postIters      <- seq(burnIn + 1, object$mainIters, thin)
  object$params  <- lapply(object$params,
                           function(x) abind::asub(x, postIters, 1))

  output <- get(param, object$params)

  densityFig <- densityplot(output)
  traceFig <- traceplot(output)
  autocorrFig <- autocorrplot(output)

  cowplot::plot_grid(densityFig, traceFig, autocorrFig, nrow=1)
}


#' @import ggplot2
#' @importFrom reshape2 melt
traceplot <- function(output) {

  varNames <- c("Iteration", "Stat")

  trace_df <- melt(output, varnames = varNames, value.name = "Posterior")
  trace_df$Iteration <- as.integer(trace_df$Iteration)

  ggplot(trace_df, aes(x = Iteration, y = Posterior)) +
    geom_line(colour = "black") +
    facet_wrap("Stat", ncol = 1, scales = "free_y")
}


#' @import ggplot2
#' @importFrom reshape2 melt
densityplot <- function(output) {

  varNames <- c("Iteration", "Stat")
  out_df <- melt(output, varnames = varNames, value.name = "Posterior")

  ggplot(out_df, aes(x = Posterior)) +
    geom_density(fill = "black", alpha = 0.5, colour = NA) +
    facet_wrap("Stat", ncol = 1, scales = "free") +
    ylab("Density")
}

#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats acf
groupDensityplot <- function(output) {

  varNames <- c("Iteration", "Group", "Stat")
  out_df <- melt(output, varnames = varNames, value.name = "Posterior")
  out_df$Group <- factor(out_df$Group)

  ggplot(out_df, aes(x = Posterior)) +
    geom_density(aes(fill = Group), alpha = 0.5, colour = NA) +
    facet_wrap("Stat", ncol = 1, scales = "free") +
    ylab("Density")
}


#' @import ggplot2
#' @importFrom reshape2 melt
autocorrplot <- function(output, lag.max = 40) {

  # Get autocorrelation for each model term
  autocorr <- apply(output, 2, function(x) acf(x, lag.max=lag.max, plot=F)$acf)

  varNames <- c("Lag", "Stat")
  ac_df <- melt(autocorr, varnames = varNames, value.name = "Autocorrelation")
  ac_df$Lag <- ac_df$Lag - 1

  ggplot(ac_df, aes(x = Lag, y = Autocorrelation)) +
    geom_bar(stat = "identity", width = 0.3) +
    facet_wrap("Stat", ncol = 1, scales = "free_y")
}

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
                            param = "mu",
                            burn_in = 0L,
                            thin = 1L){
  thin    <- as.integer(thin)
  burn_in <- as.integer(burn_in)
  
  # Remove burnIn iterations and apply thinning (default: no thinning)
  post_iters      <- seq(burn_in + 1L, object$main_iters, thin)
  object$params  <- subset(object$params, iterations = post_iters)
  output <- get(param, object$params)
  output <- abind::adrop(unclass(output), 1)
  
  ergm_terms <- object$control$model$coef.names
  model_vars <- colnames(control$mod_mat)
  n_vars <- length(model_vars)
  n_dim <- length(dim(output))
  dim_names <- vector("list", n_dim)
  dim_names[[n_dim]] <- ergm_terms
  dim_names[[n_dim - 1]] <- model_vars
  dimnames(output) <- dim_names

  density_fig <- densityplot(output)
  trace_fig <- traceplot(output)
  autocorr_fig <- autocorrplot(output)

  cowplot::plot_grid(density_fig, trace_fig, autocorr_fig, nrow=1)
}


#' @import ggplot2
#' @importFrom reshape2 melt
traceplot <- function(output) {

  varnames <- c("Iteration", "Var", "Stat")

  trace_df <- melt(output, varnames = varnames, value.name = "Estimate")
  trace_df$Iteration <- as.integer(trace_df$Iteration)

  ggplot(trace_df, aes(x = Iteration, y = Estimate)) +
    geom_line(colour = "black") +
    facet_wrap(c("Var", "Stat"), ncol = 1, scales = "free_y")
}


#' @import ggplot2
#' @importFrom reshape2 melt
densityplot <- function(output) {

  varnames <- c("Iteration", "Var", "Stat")
  out_df <- melt(output, varnames = varnames, value.name = "Estimate")

  ggplot(out_df, aes(x = Estimate)) +
    geom_density(fill = "black", alpha = 0.5, colour = NA) +
    facet_wrap(c("Var", "Stat"), ncol = 1, scales = "free") +
    ylab("Density")
}

#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats acf
groupDensityplot <- function(output) {

  varnames <- c("Iteration", "Group", "Stat")
  out_df <- melt(output, varnames = varnames, value.name = "Estimate")
  out_df$Group <- factor(out_df$Group)

  ggplot(out_df, aes(x = Estimate)) +
    geom_density(aes(fill = Group), alpha = 0.5, colour = NA) +
    facet_wrap("Stat", ncol = 1, scales = "free") +
    ylab("Density")
}


#' @import ggplot2
#' @importFrom reshape2 melt
autocorrplot <- function(output, lag.max = 40) {

  # Get autocorrelation for each model term
  autocorr <- apply(output, c(2, 3), function(x) acf(x, lag.max=lag.max, plot=F)$acf)

  varnames <- c("Lag", "Var", "Stat")
  ac_df <- melt(autocorr, varnames = varnames, value.name = "Autocorrelation")
  ac_df$Lag <- ac_df$Lag - 1

  ggplot(ac_df, aes(x = Lag, y = Autocorrelation)) +
    geom_bar(stat = "identity", width = 0.3) +
    facet_wrap(c("Var", "Stat"), ncol = 1, scales = "free_y")
}

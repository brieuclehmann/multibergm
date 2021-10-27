#' Plot BERGM posterior output
#'
#' This function creates MCMC diagnostic plots for \code{multibergm} objects.
#'
#' @inheritParams summary.multibergm
#' @param x A \link{multibergm} object
#' @param ind Indices to subset when param is theta or mu_pop
#' @param ... Additional parameters to be passed on to lower-level functions.
#' @importFrom ggplot2 ggplot
#'
#' @return Outputs a density plot, trace plot, and autocorrelation plot for
#'   each of the model parameters.
#'
#' @export

plot.multibergm <- function(x,
                            param = "mu",
                            burn_in = 0L,
                            thin = 1L,
                            ind  = 1L,
                            ...) {
  thin    <- as.integer(thin)
  burn_in <- as.integer(burn_in)

  # Remove burn_in iterations and apply thinning (default: no thinning)
  post_iters      <- seq(burn_in + 1L, x$main_iters, thin)
  x$params  <- subset(x$params, iters = post_iters)
  output <- get(param, x$params)
  output <- abind::adrop(unclass(output), 1)

  ergm_terms <- param_names(x$control$model)
  model_vars <- colnames(control$mod_mat)
  n_vars <- length(model_vars)
  n_dim <- length(dim(output))
  dim_names <- vector("list", n_dim)
  dim_names[[n_dim]] <- ergm_terms
  dim_names[[n_dim - 1]] <- model_vars
  dimnames(output) <- dim_names
  
  if (param %in% c("theta", "mu_group")) {
    output <- output[ ,ind, ,drop=FALSE]
    if (length(ind) == 1)
      output <- abind::adrop(output, drop = 2)
  }
  
  density_fig <- densityplot(output)
  trace_fig <- traceplot(output)
  autocorr_fig <- autocorrplot(output)

  cowplot::plot_grid(density_fig, trace_fig, autocorr_fig, nrow=1)
}


#' @import ggplot2
#' @importFrom reshape2 melt
traceplot <- function(output) {
  
  var_names <- c("iteration", "var", "stat")
  trace_df <- melt(output, varnames = var_names, value.name = "estimate")
  trace_df$iteration <- as.integer(trace_df$iteration)
  p <- ggplot(trace_df, aes(x = .data$iteration, y = .data$estimate)) +
    geom_line(colour = "black") +
    facet_wrap(c("var", "stat"), ncol = 1, scales = "free_y")
  
  p
}


#' @import ggplot2
#' @importFrom reshape2 melt
densityplot <- function(output, param) {
  
  if (length(dim(output)) == 3) {
    p <- group_densityplot(output)
  } else {
    var_names <- c("iteration", "var", "stat")
    out_df <- melt(output, varnames = var_names, value.name = "estimate")
    p <- ggplot(out_df, aes(x = .data$estimate)) +
      geom_density(fill = "black", alpha = 0.5, colour = NA) +
      facet_wrap(c("var", "stat"), ncol = 1, scales = "free") +
      ylab("density")
  }

  p 
}

#' @import ggplot2
#' @importFrom reshape2 melt
group_densityplot <- function(output) {

  var_names <- c("iteration", "group", "stat")
  out_df <- melt(output, varnames = var_names, value.name = "estimate")
  out_df$group <- factor(out_df$group)

  ggplot(out_df, aes(x = .data$estimate)) +
    geom_density(aes(fill = .data$group), alpha = 0.5, colour = NA) +
    facet_wrap("stat", ncol = 1, scales = "free") +
    ylab("density")
}


#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats acf
autocorrplot <- function(output, lagmax = 40) {

  if (length(dim(output)) == 3) {
    p <- group_autocorrplot(output, lagmax)
  } else {
    
    # Get autocorrelation for each model term
    autocorr <- apply(output, 2,
                      function(x) acf(x, lag.max = lagmax, plot = F)$acf)
    
    var_names <- c("lag", "stat")
    ac_df <- melt(autocorr, varnames = var_names, value.name = "autocorrelation")
    ac_df$lag <- ac_df$lag - 1
    
    ggplot(ac_df, aes(x = .data$lag, y = .data$autocorrelation)) +
      geom_bar(stat = "identity", width = 0.3) +
      facet_wrap("stat", ncol = 1, scales = "free_y")
  }
}

#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats acf
group_autocorrplot <- function(output, lagmax) {
  var_names <- c("lag", "var", "stat")
  
  # Get autocorrelation for each model term
  autocorr <- apply(output, c(2, 3),
                    function(x) acf(x, lag.max=lagmax, plot=F)$acf)

  ac_df <- melt(autocorr, varnames = var_names, value.name = "autocorrelation")
  ac_df$lag <- ac_df$lag - 1

  ggplot(ac_df, aes(x = lag, y = autocorrelation)) +
    geom_bar(stat = "identity", width = 0.3) +
    facet_wrap(c("var", "stat"), ncol = 1, scales = "free_y")
}

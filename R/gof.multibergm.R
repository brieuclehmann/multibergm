#' Goodness-of-fit diagnostics for multibergm
#'
#' Function to calculate summaries for degree, minimum geodesic distances,
#' and edge-wise shared partner distributions to diagnose the Bayesian
#' goodness-of-fit of exponential random graph models fit to multiple
#' networks.
#'
#' @inheritParams summary.multibergm
#' @param sample_size Number of networks to be simulated and compared to the
#'   observed networks.
#' @param aux_iters Number of iterations used for network simulation.
#' @param coefs Optional set of model coefficients to use for network
#'   simulation. The default is to use the posterior samples of the
#'   population-level parameter.
#' @param ... Additional parameters to be passed on to lower-level functions.
#'
#' @return Outputs bar plots (for observed networks) overlayed with ribbons
#'   (for simulated network) of degree, minimum geodesic distances,
#'   and edge-wise shared partner distributions.
#'
#' @export
#'
#' @importFrom statnet.common nonsimp_update.formula
#' @importFrom stats quantile
#' @import dplyr
#' @import ggplot2
#' @import ergm

gof.multibergm <- function(object,
                           coefs = NULL,
                           sample_size = 100L,
                           aux_iters = 1.5 * object$control$aux_iters,
                           burn_in = 0L,
                           thin = 1L,
                           ...) {

  thin    <- as.integer(thin)
  burn_in <- as.integer(burn_in)

  # Remove burn_in iterations and apply thinning (default: no thinning)
  post_iters     <- seq(burn_in + 1L, object$main_iters, thin)
  object$params  <- subset(object$params, iters = post_iters)

  # Get statistics for observed networks
  obs_df <- get_net_stats(object$networks, object$formula, "gof")

  n_iters <- mcmcr::niters(object$params$theta)
  if (is.null(coefs)) {
    coefs  <- subset(object$params$mu_pop, iters = sample(n_iters, sample_size))
  }

  sim_df <- obs_df[0, ]
  for (i in seq_len(sample_size)) {
    y         <- object$networks[[sample(length(object$networks), 1)]]
    myformula <- nonsimp_update.formula(object$formula, y ~.,
                                        from.new = "y")

    net_sim <- ergm::simulate_formula(myformula,
                                      coef = coefs[1, i, ],
                                      nsim = 1,
                                      constraints = object$constraints,
                                      control = control.simulate.formula(
                                        MCMC.burnin = aux_iters)
    )

    this_df <- get_net_stats(net_sim, myformula, "gof")
    sim_df  <- rbind(sim_df, this_df)
  }

  sim_df <- sim_df %>%
    group_by(.data$Stat, .data$n) %>%
    summarise(Group = mean(.data$Value),
              Lower = quantile(.data$Value, 0.025),
              Upper = quantile(.data$Value, 0.975)) %>%
    mutate(nMax = max(.data$n[.data$Upper > 0 & is.finite(.data$n)])) %>%
    filter(is.na(.data$n) | .data$n <= .data$nMax + 1)

  obs_df <- obs_df %>%
    group_by(.data$Stat) %>%
    mutate(nMax = max(.data$n[.data$Value > 0 & is.finite(.data$n)])) %>%
    filter(is.na(.data$n) | .data$n <= .data$nMax + 1)

  ggplot(obs_df, aes(x = .data$n, y = .data$Value, group = .data$n)) +
    geom_boxplot() +
    geom_line(data = sim_df, aes(y = .data$Group, x = .data$n),
              inherit.aes = FALSE) +
    geom_ribbon(data = sim_df,
                aes(ymin = .data$Lower, ymax = .data$Upper, x = .data$n),
                alpha = 0.4, inherit.aes = FALSE) +
    facet_wrap("Stat", scales = "free", ncol = 1) +
    scale_linetype_discrete(name = NULL, labels = "Observed data") +
    xlab("k") +
    ylab("Probability")
}

#' Compute network statistics
#'
#' GetNetStats is an internal function used to compute network statistics.
#'
#' @param object A network or list of networks
#' @param formula The ERGM formula containing the summary statistics to be
#'   computed
#' @param which_stats A string specifying which network statistics to be
#'   computed ("all", "model", "gof", or "other").
#' @param ... Arguments to be passed to methods.
#'
#' @export
get_net_stats <- function(object, ...)
  UseMethod("get_net_stats")


#' @export
#'
#' @importFrom plyr ldply
#' @describeIn get_net_stats Network statistics for a list of networks.
get_net_stats.list <- function(object, formula, which_stats, ...)
  ldply(object, get_net_stats, formula, which_stats)

###############################################################################

#' @export
#'
#' @importFrom statnet.common nonsimp_update.formula
#' @importFrom tidyr separate
#' @import network
#'
#' @describeIn get_net_stats Network statistics for a single network
get_net_stats.network <- function(object, formula, which_stats, ...) {

  model_stats <- summary(nonsimp_update.formula(formula, object ~ .,
                                                from.new = "object"))

  # Edgewise shared partners distribution
  n_nodes  <- network.edgecount(object)
  esp      <- summary(object ~ esp(0:(n_nodes - 1))) / n_nodes
  esp_df   <- data.frame(Stat = names(esp), Value = unname(esp))
  esp_df   <- separate(esp_df, Stat, c("Stat", "n"), sep = 3)

  # Geodesic distance distribution
  n_dyads      <- network.dyadcount(object)
  geo_dist     <- ergm.geodistdist(object) / n_dyads
  geo_dist_df  <- data.frame(Stat  = "distance",
                             n     = names(geo_dist),
                             Value = unname(geo_dist), stringsAsFactors = FALSE)
  geo_dist_df$n[geo_dist_df$n == Inf] <- NA_integer_

  # Degree distribution
  if (is.directed(object)) {
    n_nodes  <- network.size(object)
    ideg_dist <- summary(object ~ idegree(0:(n_nodes - 1))) / n_nodes
    ideg_df   <- data.frame(Stat = names(ideg_dist), Value = unname(ideg_dist))
    ideg_df   <- separate(ideg_df, Stat, c("Stat", "n"), sep = 7)

    odeg_dist <- summary(object ~ odegree(0:(n_nodes - 1))) / n_nodes
    odeg_df   <- data.frame(Stat = names(odeg_dist), Value = unname(odeg_dist))
    odeg_df   <- separate(odeg_df, Stat, c("Stat", "n"), sep = 7)

    gof_stats   <- bind_rows(ideg_df, odeg_df, esp_df, geo_dist_df)
  } else {
    n_nodes  <- network.size(object)
    deg_dist <- summary(object ~ degree(0:(n_nodes - 1))) / n_nodes
    deg_df   <- data.frame(Stat = names(deg_dist), Value = unname(deg_dist))
    deg_df   <- separate(deg_df, Stat, c("Stat", "n"), sep = 6)

    gof_stats   <- bind_rows(deg_df, esp_df, geo_dist_df)
  }

  gof_stats$n <- as.integer(gof_stats$n)

  switch(which_stats,
         model = data.frame(t(model_stats)),
         gof   = gof_stats)
}

###############################################################################

#' @export
#'
#' @describeIn get_net_stats Network statistics for a network.list
get_net_stats.network.list <- function(object, formula, ...)
  plyr::ldply(object, get_net_stats, formula)

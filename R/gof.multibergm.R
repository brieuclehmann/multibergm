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
                           param = NULL,
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
  object$networks <- mapply(function(x,y) set.network.attribute(x, "group", y),
                            object$networks, object$groups, SIMPLIFY = FALSE)
  obs_df <- get_net_stats(object$networks, object$formula, "gof")

  # Get posterior samples
  if (is.null(coefs)) {
    if (is.null(param)) {
      n_groups <- length(unique(object$groups))
      param <- ifelse(n_groups > 1, "mu_group", "mu_pop")
    }
    output <- get(param, object$params)
    
    n_iters <- mcmcr::niters(object$params$theta)
    coefs   <- subset(output, iters = sample(n_iters, sample_size))
    coefs   <- abind::adrop(unclass(coefs), 1)
  }

  sim_df <- obs_df[0, ]
  for (g in unique(object$groups)) {
    group_nets <- object$networks[object$groups == g]
    
    for (i in seq_len(sample_size)) {
      y         <- group_nets[[sample(length(group_nets), 1)]]
      myformula <- nonsimp_update.formula(object$formula, y ~.,
                                          from.new = "y")
      
      if (length(dim(coefs)) == 2) {
        this_coef <- coefs[i,]
      } else {
        this_coef <- coefs[i,g,]
      }
      
      net_sim <- ergm::simulate_formula(myformula,
                                        coef = this_coef,
                                        nsim = 1,
                                        constraints = object$constraints,
                                        control = control.simulate.formula(
                                          MCMC.burnin = aux_iters
                                        )
      )
      
      this_df <- get_net_stats(net_sim, myformula, "gof")
      this_df$group <- g
      sim_df  <- rbind(sim_df, this_df)
    }
  }

  sim_df <- sim_df %>%
    group_by(.data$Stat, .data$n, .data$group) %>%
    summarise(Group = quantile(.data$Value, 0.5),
              Lower = quantile(.data$Value, 0.05),
              Upper = quantile(.data$Value, 0.95))
  
  # Degree distribution
  if (is.directed(object$networks[[1]])) {
    ideg_obs <- filter(obs_df, Stat == "idegree")
    ideg_sim <- filter(sim_df, Stat == "idegree")
    p0 <- plot_dist_gof(ideg_obs, ideg_sim, "In-degree")
    
    odeg_obs <- filter(obs_df, Stat == "odegree")
    odeg_sim <- filter(sim_df, Stat == "odegree")
    p1 <- plot_dist_gof(odeg_obs, odeg_sim, "Out-degree")
  } else {
    deg_obs <- filter(obs_df, Stat == "degree")
    deg_sim <- filter(sim_df, Stat == "degree")
    p1 <- plot_dist_gof(deg_obs, deg_sim, "Degree")
  }

  
  # Geodesic distance distrobution
  geodist_obs <- filter(obs_df, Stat == "distance")
  geodist_sim <- filter(sim_df, Stat == "distance")
  p2 <- plot_dist_gof(geodist_obs, geodist_sim, "Geodesic distance")
  
  # Esp distrobution
  espdist_obs <- filter(obs_df, Stat == "esp")
  espdist_sim <- filter(sim_df, Stat == "esp")
  p3 <- plot_dist_gof(espdist_obs, espdist_sim, "Edgewise shared partners")

  if (is.directed(object$networks[[1]])) {
    cowplot::plot_grid(p0, p1, p2, p3, ncol = 1)
  } else {
    cowplot::plot_grid(p1, p2, p3, ncol = 1)
  }
  
}

geodist_gof <- function(obs_df, sim_df) {
  
  n_max_sim <- max(1L, sim_df$n[sim_df$Upper > 0], na.rm = TRUE)
  n_max_obs <- max(1L, obs_df$n[obs_df$Value > 0], na.rm = TRUE)
  n_max <- max(n_max_sim, n_max_obs)
  
  sim_df <- sim_df %>%
    filter(n <= n_max | is.na(n)) %>%
    mutate(n = if_else(is.na(n), n_max + 1L, n))
  obs_df <- obs_df %>%
    filter(n <= n_max | is.na(n)) %>%
    mutate(n = if_else(is.na(n), n_max + 1L, n))
  
  finite_breaks <- pretty(seq(n_max))
  finite_breaks <- finite_breaks[finite_breaks < n_max]
  breaks <- c(finite_breaks, n_max + 1)
  labels <- c(finite_breaks, "Inf")
  p1 <- ggplot(obs_df, aes(x = .data$n, y = .data$Value, group = .data$n)) +
    geom_boxplot() +
    geom_line(data = sim_df, aes(y = .data$Group, x = .data$n),
              inherit.aes = FALSE) +
    geom_ribbon(data = sim_df,
                aes(ymin = .data$Lower, ymax = .data$Upper, x = .data$n),
                alpha = 0.4, inherit.aes = FALSE) +
    scale_linetype_discrete(name = NULL, labels = "Observed data") +
    scale_x_continuous(breaks = breaks, minor_breaks = seq(n_max + 1),
                       labels = labels) +
    xlab("Geodesic distance") +
    ylab("Probability")
  
  p1
}

plot_dist_gof <- function(obs_df, sim_df, name) {
  
  n_max_sim <- max(1L, sim_df$n[sim_df$Upper > 0], na.rm = TRUE)
  n_max_obs <- max(1L, obs_df$n[obs_df$Value > 0], na.rm = TRUE)
  n_max <- max(n_max_sim, n_max_obs)
  
  sim_df <- sim_df %>%
    filter(n <= n_max | is.na(n)) %>%
    mutate(n = if_else(is.na(n), n_max + 1L, n))
  obs_df <- obs_df %>%
    filter(n <= n_max | is.na(n)) %>%
    mutate(n = if_else(is.na(n), n_max + 1L, n))
  
  if (name == "Geodesic distance") {
    finite_breaks <- pretty(seq(n_max))
    finite_breaks <- finite_breaks[finite_breaks <= n_max]
    breaks <- c(finite_breaks, n_max + 1)
    labels <- c(finite_breaks, "Inf")
  } else {
    finite_breaks <- pretty(seq(0, n_max))
    finite_breaks <- round(finite_breaks[finite_breaks <= n_max])
    breaks <- finite_breaks
    labels <- finite_breaks
  }

  sim_df$group <- factor(sim_df$group)
  obs_df$group <- factor(obs_df$group)
  
  n_groups <- length(unique(sim_df$group))
  if (n_groups == 1) {
    ggplot(obs_df, aes(x = .data$n, y = .data$Value, group = .data$n)) +
      geom_boxplot() +
      geom_line(data = sim_df, aes(y = .data$Group, x = .data$n),
                inherit.aes = FALSE) +
      geom_ribbon(data = sim_df,
                  aes(ymin = .data$Lower, ymax = .data$Upper, x = .data$n),
                  alpha = 0.4, inherit.aes = FALSE) +
      scale_linetype_discrete(name = NULL, labels = "Observed data") +
      scale_x_continuous(breaks = breaks, minor_breaks = seq(n_max + 1),
                         labels = labels) +
      xlab(name) +
      ylab("Probability")
  } else {
    ggplot(obs_df, aes(x = .data$n, y = .data$Value, group = .data$n)) +
      geom_boxplot() +
      geom_line(data = sim_df, aes(y = .data$Group, x = .data$n, color = .data$group),
                inherit.aes = FALSE) +
      geom_ribbon(data = sim_df,
                  aes(ymin = .data$Lower, ymax = .data$Upper, x = .data$n,
                      fill = .data$group),
                  alpha = 0.4, inherit.aes = FALSE) +
      scale_linetype_discrete(name = NULL, labels = "Observed data") +
      scale_x_continuous(breaks = breaks, minor_breaks = seq(n_max + 1),
                         labels = labels) +
      xlab(name) +
      ylab("Probability") +
      facet_wrap("group", ncol = n_groups)
  }
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
get_net_stats.list <- function(object, formula, which_stats, ...) {
  out <- ldply(object, get_net_stats, formula, which_stats)
  out$group <- sapply(object, get.network.attribute, "group")
  
  out
}

###############################################################################

#' @export
#'
#' @importFrom statnet.common nonsimp_update.formula
#' @import network
#'
#' @describeIn get_net_stats Network statistics for a single network
get_net_stats.network <- function(object, formula, which_stats, ...) {

  switch(which_stats,
         model    = get_model_stats(object, formula),
         deg      = get_deg_dist(object),
         esp      = get_esp_dist(object),
         geodist  = get_geodist_dist(object),
         gof      = bind_rows(get_deg_dist(object),
                              get_esp_dist(object),
                              get_geodist_dist(object)))
}

###############################################################################

#' @export
#'
#' @describeIn get_net_stats Network statistics for a network.list
get_net_stats.network.list <- function(object, formula, which_stats, ...) {
  out <- plyr::ldply(object, get_net_stats, formula, which_stats)
  out$group <- sapply(object, get.network.attribute, "group")
  out
}
  

###############################################################################

#' Get model summary statistics
#' 
#' Internal function to get the model summary statistics for a single network
#' @export
#' 
#' @importFrom statnet.common nonsimp_update.formula
#' @import network
#' 
#' @param network A network object
#' 
#' @rdname get_net_stats
get_model_stats <- function(network, formula) {
  formula <- nonsimp_update.formula(formula, network ~ ., from.new = "network")
  model_stats <- summary(formula)
  
  data.frame(t(model_stats))
}

#' Get degree distribution
#' 
#' Internal function to get the degree distribution for a single network
#' @export
#'
#' @importFrom tidyr separate
#' @import network
#'
#' @rdname get_net_stats
get_deg_dist <- function(network) {
  
  if (is.directed(network)) {
    n_nodes  <- network.size(network)
    ideg_dist <- summary(network ~ idegree(0:(n_nodes - 1))) / n_nodes
    ideg_df   <- data.frame(Stat = names(ideg_dist), Value = unname(ideg_dist))
    ideg_df   <- separate(ideg_df, Stat, c("Stat", "n"), sep = 7)
    
    odeg_dist <- summary(network ~ odegree(0:(n_nodes - 1))) / n_nodes
    odeg_df   <- data.frame(Stat = names(odeg_dist), Value = unname(odeg_dist))
    odeg_df   <- separate(odeg_df, Stat, c("Stat", "n"), sep = 7)
    
    out  <- bind_rows(ideg_df, odeg_df)
  } else {
    n_nodes  <- network.size(network)
    deg_dist <- summary(network ~ degree(0:(n_nodes - 1))) / n_nodes
    deg_df   <- data.frame(Stat = names(deg_dist), Value = unname(deg_dist))
    out      <- separate(deg_df, Stat, c("Stat", "n"), sep = 6)
  }
  out$n <- as.integer(out$n)
  
  out
}

#' Get edgewise-shared partner distribution
#' 
#' Internal function to get the edgewise-shared partner distribution for a single network
#' @export
#'
#' @importFrom tidyr separate
#' @import network
#'
#' @rdname get_net_stats
get_esp_dist <- function(network) {
  
  # Edgewise shared partners distribution
  n_nodes  <- network.size(network)
  n_edges  <- max(1, network.edgecount(network))
  esp      <- summary(network ~ esp(0:(n_nodes - 2))) / n_edges
  esp_df   <- data.frame(Stat = names(esp), Value = unname(esp))
  
  esp_df <- separate(esp_df, Stat, c("Stat", "n"), sep = 3)
  esp_df$n <- as.integer(esp_df$n)
  
  esp_df
}

#' Get geodesic distance distribution
#' 
#' Internal function to get the geodesic distance distribution for a single network
#' 
#' @export
#'
#' @importFrom tidyr separate
#' @import network
#'
#' @rdname get_net_stats 
get_geodist_dist <- function(network) {
  
  # Geodesic distance distribution
  n_nodes      <- network.size(network)
  n_dyads      <- network.dyadcount(network)
  geo_dist     <- ergm.geodistdist(network) / n_dyads
  geo_dist_df  <- data.frame(Stat  = "distance",
                             n     = c(1:(n_nodes - 1), NA_integer_),
                             Value = unname(geo_dist), stringsAsFactors = FALSE)
  
  geo_dist_df$Value <- geo_dist_df$Value / sum(geo_dist_df$Value)
  geo_dist_df
}

#' Goodness-of-fit diagnostics for multibergm
#'
#' Function to calculate summaries for degree, minimum geodesic distances,
#' and edge-wise shared partner distributions to diagnose the Bayesian
#' goodness-of-fit of exponential random graph models fit to multiple
#' networks.
#'
#' @inheritParams summary.multibergm
#' @param sampleSize Number of networks to be simulated and compared to the
#'   observed networks.
#' @param auxIters Number of iterations used for network simulation.
#' @param coefs Optional set of model coefficients to use for network
#'   simulation. The default is to use the posterior samples of the
#'   population-level parameter.
#'
#' @return Outputs bar plots (for observed networks) overlayed with ribbons
#'   (for simulated network) of degree, minimum geodesic distances,
#'   and edge-wise shared partner distributions.
#'
#' @export gof.multibergm
#'
#' @importFrom statnet.common nonsimp_update.formula
#' @importFrom stats quantile
#' @import dplyr
#' @import ggplot2
#' @import ergm

gof.multibergm <- function(object,
                           coefs = NULL,
                           sampleSize = 100,
                           auxIters = 1.5*object$control$auxIters,
                           burnIn = 0,
                           thin = 1){

  # Remove burnIn iterations and apply thinning (default: no thinning)
  postIters      <- seq(burnIn + 1, object$mainIters, thin)
  object$params  <- lapply(object$params,
                           function(x) abind::asub(x, postIters, 1))

  # Get statistics for observed networks
  obs_df <- GetNetStats(object$networks, object$formula, "gof")

  nIters <- dim(object$params$theta)[1]
  if (is.null(coefs))
    coefs  <- as.matrix(object$params$muPop[sample(nIters, sampleSize), ])

  sim_df <- obs_df[0, ]
  for (i in 1:sampleSize) {
    y         <- object$networks[[sample(length(object$networks),1)]]
    myformula <- nonsimp_update.formula(object$formula, y ~.,
                                        from.new = "y")

    netSim <- simulate(myformula, coef = coefs[i, ],
                       nsim = 1, constraints = object$constraints,
                       control = control.simulate.formula(MCMC.burnin=auxIters))

    this_df <- GetNetStats(netSim, myformula, "gof")
    sim_df  <- rbind(sim_df,this_df)
  }

  sim_df <- sim_df %>%
    group_by(Stat, n) %>%
    summarise(Group = mean(Value),
              Lower = quantile(Value, 0.025),
              Upper = quantile(Value, 0.975)) %>%
    mutate(nMax = max(n[Upper>0 & is.finite(n)])) %>%
    filter(is.infinite(n) | n <= nMax + 1)

  obs_df <- obs_df %>%
    group_by(Stat) %>%
    mutate(nMax = max(n[Value>0 & is.finite(n)])) %>%
    filter(is.infinite(n) | n <= nMax + 1)

  ggplot(obs_df, aes(x=n, y = Value, group = n)) +
    geom_boxplot() +
    geom_line(data = sim_df, aes(y=Group, x=n),
              inherit.aes = FALSE) +
    geom_ribbon(data = sim_df,
                aes(ymin=Lower, ymax=Upper, x=n), alpha=0.4,
                inherit.aes=F) +
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
#' @param whichStats A string specifying which network statistics to be
#'   computed ("all", "model", "gof", or "other").
#'
#' @export
GetNetStats <- function(object, formula, whichStats)
  UseMethod("GetNetStats")


#' @export
#'
#' @importFrom plyr ldply
GetNetStats.list <- function(object, formula, whichStats)
  ldply(object, GetNetStats, formula, whichStats)

##########################################################################################

#' @export
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph assortativity.degree
#' @importFrom igraph average.path.length
#' @importFrom igraph transitivity
#' @importFrom statnet.common nonsimp_update.formula
#' @importFrom tidyr separate
#' @import network
GetNetStats.network <- function(object, formula, whichStats) {

  modelStats <- summary(nonsimp_update.formula(formula, object ~ .,
                                               from.new = "object"))

  graph      <- graph_from_adjacency_matrix(as.matrix(object))
  graphStats <- data.frame(assortativity     = assortativity.degree(graph),
                           transitivity      = transitivity(graph),
                           averagePathLength = average.path.length(graph,
                                                                   directed=F))

  # Degree distribution
  nNodes  <- network.size(object)
  degDist <- summary(object ~ degree(0:(nNodes - 1)))/nNodes
  deg_df  <- data.frame(Stat=names(degDist), Value=unname(degDist))
  deg_df  <- separate(deg_df, Stat, c("Stat", "n"), sep = 6)

  # Edgewise shared partners distribution
  nEdges  <- network.edgecount(object)
  esp     <- summary(object ~ esp(0:(nNodes - 1)))/nEdges
  esp_df  <- data.frame(Stat=names(esp), Value=unname(esp))
  esp_df  <- separate(esp_df, Stat, c("Stat", "n"), sep = 3)

  # Geodesic distance distribution
  nDyads     <- network.dyadcount(object)
  geoDist    <- ergm.geodistdist(object)/nDyads
  geoDist_df <- data.frame(Stat  = "distance",
                           n     = names(geoDist),
                           Value = unname(geoDist), stringsAsFactors = FALSE)

  gofStats   <- bind_rows(deg_df, esp_df, geoDist_df)
  gofStats$n <- as.numeric(gofStats$n)

  switch(whichStats,
         all   = data.frame((t(modelStats)), graphStats),
         model = data.frame(t(modelStats)),
         gof   = gofStats,
         other = data.frame(graphStats))

}

###############################################################################
GetNetStats.network.list <- function(object, ergmFormula)
  plyr::ldply(object, GetNetStats, ergmFormula)

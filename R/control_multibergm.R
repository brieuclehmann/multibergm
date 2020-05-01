#' Set multibergm options
#'
#' Auxiliary function used to specify settings controlling multibergm fitting.
#' This allows the user to control various aspects of the model and fitting
#' process, including prior specification, parameter initialisation, MCMC
#' proposal covariance, auxiliary iterations for ERGM simulation, and number of
#' cores to be used in parallel.
#'
#' @param formula An \R \code{\link{formula}} object, of the form
#'   \code{y ~ <model terms>}, where \code{y} is a
#'   \code{\link[network]{network}} object or a
#'   \code{\link[ergm]{network.list}} object.
#' @param constraints A one-sided formula specifying one or more constraints
#'   on the support of the distribution of the networks being simulated.
#' @param priors A list of prior settings.
#' @param init A list of initial values.
#' @param auxIters Number of internal (auxiliary) MCMC iterations used to
#'   simulate a network from the model.
#' @param nBatches Number of computing cores available to simulate networks
#'   from the model in parallel.
#' @param proposal A list of covariance matrices used to generate the MCMC
#'   proposals.
#'
#' @return A list containing the following control parameters:
#'   \itemize{
#'      \item \code{auxIters}: Number of internal (auxiliary) MCMC iterations
#'      used to simulate a network from the model
#'      \item \code{prior}: A list of prior settings
#'      \item \code{init}: A list of initial values
#'      \item \code{proposal}: A list of covariance matrices used to
#'      generate the MCMC proposals.
#'      \item \code{groupLabel}: A vector containing  group
#'      labels for each network.
#'      \item \code{batches}: Used to batch parallel runs of ergm simulation
#'      \item \code{Clists}: ERGM parameters passed to the ergm simulation
#'      function
#'      \item \code{MHproposals}: ERGM parameters passed to the ergm simulation
#'      function
#'   }
#'
#' @import ergm
#' @importFrom stats terms
#' @export
control_multibergm <- function(formula, 
                               constraints = ~. ,
                               priors      = NULL,
                               init        = NULL, 
                               proposal    = NULL,
                               auxIters    = 10000, 
                               nBatches    = 1) {

  networks <- statnet.common::eval_lhs.formula(formula)

  # For compatibility with a single network
  if (is.network(networks))  networks <- list(networks)

  if (!is.network(networks[[1]]))
    stop("Input must be a network or list of networks")

  groupLabel    <- get_labels(networks)
  nNetworks <- length(networks)
  nTerms    <- length(attr(terms(formula), "term.labels"))

  # Set default settings for priors unless explicitly specified
  prior <- set_priors(formula, priors)

  # Set default settings for proposals unless explicitly specified
  if (is.null(proposal$theta)) proposal$theta <- diag(0.1^2, nTerms)
  if (is.null(proposal$mu))    proposal$mu    <- diag(0.1^2, nTerms)

  # Set default settings for initial values unless explicitly specified
  init <- set_init(formula, init, prior)

  # Specify network batching for parallelisation
  if (nBatches == 1) {
    batches <- list(1:nNetworks)
  } else {
    batches <- split(1:nNetworks, cut(seq_along(1:nNetworks),
                                      nBatches, labels = FALSE))
  }

  # Set up ergm parameters
  model       <- ergm_model(formula, networks[[1]])
  Clists      <- lapply(networks, function(x) ergm.Cprepare(x, model))
  MHproposals <- ergm_proposal(constraints, control.ergm()$MCMC.prop.args,
                               networks[[1]])

  list(auxIters    = auxIters,
       prior       = prior,
       proposal    = proposal,
       groupLabel  = groupLabel,
       init        = init,
       batches     = batches,
       model       = model,
       Clists      = Clists,
       MHproposals = MHproposals)

}


# Set priors ==================================================================

set_priors <- function(formula, prior = list()){

  nTerms <- length(attr(terms(formula), "term.labels"))

  # Use defaults if not specified explicitly

  # Set priors for population mean parameter
  if (is.null(prior$muPop$mean))
    prior$muPop$mean <- rep(0, nTerms)
  if (is.null(prior$muPop$cov))
    prior$muPop$cov <- diag(100, nTerms)

  # Set priors for network-level covariance parameter
  if (is.null(prior$covTheta$df))
    prior$covTheta$df <- nTerms + 1
  if (is.null(prior$covTheta$scale))
    prior$covTheta$scale <- diag(nTerms)

  # Set priors for group-level covariance parameter
  if (is.null(prior$covMuGroup$df))
    prior$covMuGroup$df <- nTerms + 1
  if (is.null(prior$covMuGroup$scale))
    prior$covMuGroup$scale <- diag(nTerms)

  prior
}


# Set initial values ==========================================================
#' @importFrom mvtnorm rmvnorm
#' @importFrom MCMCpack riwish

set_init <- function(formula, init, prior){

  networks <- statnet.common::eval_lhs.formula(formula)

  nNetworks <- length(networks)
  nTerms    <- length(attr(terms(formula), "term.labels"))

  # Get group assignments for each network
  groupLabel  <- get_labels(networks)
  nGroups     <- length(unique(groupLabel))

  # Use default if not specified explicitly
  if (is.null(init$muPop))
    init$muPop      <- rmvnorm(1, prior$muPop$mean, prior$muPop$cov)[1, ]

  if (is.null(init$covMuGroup))
    init$covMuGroup <- riwish(prior$covMuGroup$df, prior$covMuGroup$scale)

  if (is.null(init$covTheta))
    init$covTheta   <- riwish(prior$covTheta$df, prior$covTheta$scale)

  if (is.null(init$muGroup))
    init$muGroup    <- rmvnorm(nGroups, init$muPop, init$covMuGroup)

  if (is.null(init$theta))
    init$theta      <- rmvnorm(nNetworks, sigma=init$covTheta)

  # Add dimension to initial array
  # DIM ensures compatibility with vectors
  DIM <- function(x) {
    if (is.null(dim(x))) return(length(x))
    dim(x)
  }

  lapply(init, function(x) array(x, c(1, DIM(x))))
}



# Function to set/extract labels for/from a group of networks
get_labels <- function(networks){

  nNetworks <- length(networks)

  # Get group assignments for each network
  group <- unlist(sapply(networks, get.network.attribute, "group"))

  if (length(group) == 0){
#    warning(strwrap("Group network assignments have not been specified -
#            assuming each network belongs to the same group."))
    group <- rep(1, nNetworks)
  }

  if (length(group) != nNetworks)
    stop("Not all networks have been assigned a group")

  group
}

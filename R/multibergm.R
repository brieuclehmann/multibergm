# Main function ===============================================================
#
#' Fit multibergm to multiple networks
#'
#' Main function to fit a hierarchical framework of Bayesian exponential random
#' graph models to a population of networks. Each network \eqn{i} is modelled
#' as an ERGM with network-level parameters \eqn{\theta_i}, each of which come
#' from a common group- or population-level distribution. An
#' exchange-within-Gibbs algorithm is used to generate samples from the joint
#' posterior. The group memberships, priors, and other model fitting
#' hyperparameters are set through the \code{\link{control_multibergm}}
#' function.
#'
#' @param object A multibergm object, or a \R \code{\link{formula}} object,
#'   of the form \code{y ~ <model terms>}, where \code{y} is a
#'   \code{\link[network]{network}} object or a
#'   \code{\link[ergm]{network.list}} object.
#' @param main_iters Number of (outer) MCMC iterations
#' @param control A list of parameters set by \code{\link{control_multibergm}}
#' @param groups A vector of group memberships
#' @param prior A list of explicit prior specifications.
#' @param init A list of initial values.
#' @param ... Arguments to be passed to methods.
#'
#' @return A list containing the following elements:
#'   \itemize{
#'     \item \code{networks} - to which the model was fit
#'     \item \code{formula} - specifiying the ERGM
#'     \item \code{constraints} - used to fix any summary statistics
#'     \item \code{main_iters} - the number of MCMC iterations used
#'     \item \code{control} parameters used to fit the model
#'     \item \code{params} - list containing MCMC output for each variable
#'     \item \code{accepts} - list containing acceptance counts at each
#'     iteration
#'     }
#'
#' @export
multibergm <- function(object, ...) {
  UseMethod("multibergm")
}


#' @param constraints A one-sided formula specifying one or more constraints
#'   on the support of the distribution of the networks being simulated.
#'
#' @importFrom abind abind
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom mcmcr bind_iterations is.mcmcr as.mcmcr
#' @describeIn multibergm S3 method for class 'formula'
#' @export
multibergm.formula <- function(object,
                               constraints = ~.,
                               groups = NULL,
                               main_iters = 1000L,
                               prior = set_priors(object, groups),
                               init = set_init(object, prior, groups),
                               control = control_multibergm(object,
                                                            constraints),
                               ...) {

  start_time <- Sys.time()
  main_iters <- as.integer(main_iters)
  if (!mcmcr::is.mcmcr(init)) init <- mcmcr::as.mcmcr(init)

  # TODO: add checks for inputs
  networks <- statnet.common::eval_lhs.formula(object)
  n_networks <- length(networks)
  #prior <- set_priors(object, groups, prior)

  if (is.null(groups))
    groups <- rep(1, n_networks)

  # Set up parallel options
  n_batches <- length(control$batches)
  if (n_batches == 1) {
    foreach::registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(n_batches)
    parallel::clusterEvalQ(cl, {
      library(ergm)
    })
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  }

  # Preallocate and initialise variables for MCMC
  first_iter <- dim(init$theta)[1] + 1
  # TODO: change this to abind::afill

  add_iters <- main_iters - first_iter + 1
  params <- init

  # Preallocate acceptance counts for MCMC
  n_groups <- length(unique(groups))
  accepts <- list(theta = array(0, c(main_iters, n_networks)),
                  mu    = array(0, c(main_iters, n_groups)))

  pb <- txtProgressBar(min = 0, max = main_iters, style = 3)

  # Run the MCMC
  for (k in seq(first_iter, main_iters)) {
    setTxtProgressBar(pb, k)

    curr <- subset(params, iters = k - 1L)
    curr <- lapply(curr, function(x) abind::adrop(unclass(x), c(1, 2)))

    # Perform Gibbs update for all variables
    if (n_groups == 1) {
      mcmc_update <- singlegroup_update(curr, prior, groups, control)
    } else {
      mcmc_update <- multigroup_update(curr, prior, groups, control)
    }

    params <- mcmcr::bind_iterations(params,
                                     mcmcr::as.mcmcr(mcmc_update$params))

    accepts$theta[k, ] <- mcmc_update$accepts$theta
    accepts$mu[k, ]    <- mcmc_update$accepts$mu

  }

  close(pb)

  end_time <- Sys.time()

  out <- list(networks    = networks,
              formula     = object,
              groups      = groups,
              main_iters  = main_iters,
              control     = control,
              prior       = prior,
              accepts     = accepts,
              params      = params,
              constraints = constraints,
              time        = end_time - start_time)

  class(out) <- "multibergm"

  return(out)
}


#' @describeIn multibergm S3 method for class 'multibergm', used to continue
#'   generating posterior samples from a previous fit
#' @export
multibergm.multibergm <- function(object, main_iters = 1000, ...) {

    old_iters <- dim(object$params$theta)[1]

    # Pass on control parameters from previous run
    control <- control_multibergm(object$formula,
                                  object$constraints,
                                  object$control$proposal,
                                  object$control$aux_iters,
                                  length(object$control$batches))

    new_iters <- main_iters + old_iters

    multibergm.formula(object$formula, object$constraints,
                       new_iters, control)
}

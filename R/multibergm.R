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
#' hyperparameters are set through the \code{\link{control_multibergm}} function.
#'
#' @param object A multibergm object, or a \R \code{\link{formula}} object,
#'   of the form \code{y ~ <model terms>}, where \code{y} is a
#'   \code{\link[network]{network}} object or a
#'   \code{\link[ergm]{network.list}} object.
#' @param mainIters Number of (outer) MCMC iterations
#' @param control A list of parameters set by \code{\link{control_multibergm}}
#'
#' @return A list containing the following elements:
#'   \itemize{
#'     \item \code{networks} - to which the model was fit
#'     \item \code{formula} - specifiying the ERGM
#'     \item \code{constraints} - used to fix any summary statistics
#'     \item \code{mainIters} - the number of MCMC iterations used
#'     \item \code{control} parameters used to fit the model
#'     \item \code{params} - list containing MCMC output for each variable
#'     \item \code{accepts} - list containing acceptance counts at each iteration
#'     }
#'
#' @export
multibergm <- function(object,
                       mainIters = 1000,
                       control = control_multibergm(formula))
  UseMethod("multibergm")

#' @param constraints A one-sided formula specifying one or more constraints
#'   on the support of the distribution of the networks being simulated.
#'
#' @importFrom abind abind
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @describeIn multibergm S3 method for class 'formula'
#' @export
multibergm.formula <- function(object,
                               constraints = ~.,
                               mainIters = 1000,
                               control = control_multibergm(formula, constraints)) {

  startTime <- Sys.time()

  # TODO: add check for constraints
  networks <- statnet.common::eval_lhs.formula(formula)

  # Set up parallel options
  nBatches <- length(control$batches)
  if (nBatches == 1) {
    foreach::registerDoSEQ()
  } else {
    if (!is.null(nodeList)) {
      cl <- parallel::makeCluster(nodeList)
    } else cl <- parallel::makeCluster(nBatches)

    parallel::clusterEvalQ(cl, {
      library(ergm)
    })
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  }

  # Preallocate and initialise variables for MCMC
  firstIter <- dim(control$init$theta)[1] + 1
  # TODO: change this to abind::afill
  params <- lapply(control$init,
                   function(x) abind(x,
                                     array(NA,
                                           c(mainIters-firstIter+1, dim(x)[-1])),
                                     along = 1))

  # Preallocate acceptance counts for MCMC
  nGroups <- length(unique((control$groupLabel)))
  accepts <- list(theta = array(0, dim(params$theta)[c(1,2)]),
                  mu = array(0, c(mainIters, nGroups)))

  pb <- txtProgressBar(min = 0, max = mainIters, style = 3)

  # Run the MCMC
  for (k in firstIter:mainIters) {
    setTxtProgressBar(pb, k)

    # Get the variables from the last iteration
    curr <- lapply(params,
                   function(x) array(x[slice.index(x,1) == k-1], dim(x)[-1]))

    # Perform Gibbs update for all variables
    if (nGroups == 1) {
      mcmcUpdate <- singleGroup_update(curr, control)
    } else {
      mcmcUpdate <- multiGroup_update(curr, control)
    }

    # Assign values and acceptance counts for this iteration
    for (x in names(params))
      params[[x]][slice.index(params[[x]],1) == k] <- mcmcUpdate$params[[x]]

    accepts$theta[k, ] <- mcmcUpdate$accepts$theta
    accepts$mu[k, ]    <- mcmcUpdate$accepts$mu

  }

  close(pb)

  endTime <- Sys.time()

  out = list(networks = networks, formula = formula,
             mainIters = mainIters, control = control,
             accepts = accepts, params = params,
             constraints = constraints, time = endTime - startTime)

  class(out) <- "multibergm"

  return(out)
}


#' @describeIn multibergm S3 method for class 'multibergm', used to continue
#'   generating posterior samples from a previous fit
#' @export
multibergm.multibergm <- function(object, mainIters = 1000) {

    oldIters <- dim(object$params$theta)[1]

    # Pass on control parameters from previous run
    control <- control_multibergm(object$formula,
                                  priors = object$control$prior,
                                  proposal = object$control$proposal,
                                  nBatches = length(object$control$batches),
                                  auxIters = object$control$auxIters,
                                  init = object$params)

    newIters <- mainIters + oldIters

    multibergm.formula(object$formula, object$constraints,
                       newIters, control)
}

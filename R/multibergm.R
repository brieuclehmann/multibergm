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
#'     \item \code{ergm_formula} - specifying the ERGM
#'     \item \code{model_formula} - specifying the linear model relating network covariates to ERGM summary statistics
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
                               model_formula = ~ 1,
                               constraints = ~ .,
                               main_iters = 1000L,
                               model_matrix = get_model_matrix(object, 
                                                               model_formula),
                               control = control_multibergm(object,
                                                            model_matrix,
                                                            constraints),
                               prior = set_priors(object, 
                                                  model_matrix, 
                                                  control),
                               init = set_init(object, prior, model_matrix),
                               ...) {

  start_time <- Sys.time()
  main_iters <- as.integer(main_iters)
  if (!mcmcr::is.mcmcr(init)) init <- mcmcr::as.mcmcr(init)

  # TODO: add checks for inputs
  networks <- statnet.common::eval_lhs.formula(object)
  n_networks <- length(networks)

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

  # Preallocate acceptance counts for MCMC
  n_vars <- nrow(model_matrix)
  accepts <- list(theta = array(0, c(main_iters, n_networks)),
                  mu = array(0, c(main_iters, n_vars)))
  
  # Preallocate and initialise variables for MCMC
  first_iter <- dim(init$theta)[1] + 1

  add_iters <- main_iters - first_iter + 1
  # TODO: change this to abind::afill and preallocate
  params <- init
  proposals <- list(theta = control$init_proposals$theta,
                    theta_scale = rep(0, n_networks),
                    mu       = control$init_proposals$mu,
                    mu_scale = rep(0, n_vars))

  pb <- txtProgressBar(min = 0, max = main_iters, style = 3)

  # Run the MCMC
  for (k in seq(first_iter, main_iters)) {
    setTxtProgressBar(pb, k)

    curr <- subset(params, iters = k - 1L)
    curr <- lapply(curr, function(x) abind::adrop(unclass(x), c(1, 2)))

    # Perform Gibbs update for all variables
    mcmc_update <- asis_update(curr, prior, model_matrix, proposals, control)
    #if (n_vars == 1) {
    #  mcmc_update <- single_var_update(curr, prior, model_matrix, control)
    #} else {
    #  mcmc_update <- multi_var_update(curr, prior, model_matrix, control)
    #}
    
    # Assign values and acceptance counts for this iteration
    # for (x in names(params))
    #   params[[x]][slice.index(params[[x]],1) == k] <- mcmc_update$params[[x]]
    
    params <- mcmcr::bind_iterations(params,
                                     mcmcr::as.mcmcr(mcmc_update$params))

    accepts$theta[k, ] <- mcmc_update$accepts$theta
    accepts$mu[k, ]    <- mcmc_update$accepts$mu
    
    # Update RW proposals (adaptive MCMC step)
    if (k %% control$proposal_update_freq == 0 & k < control$proposal_update_max) {
      last_iters <- k:(k - control$proposal_update_freq + 1)
      accept_rate <- list(theta = colMeans(accepts$theta[last_iters, ,drop=FALSE]),
                          mu    = colMeans(accepts$mu[last_iters, ,drop=FALSE]))
      proposals <- update_proposals(proposals, accept_rate, params, control)
    }

  }

  close(pb)

  end_time <- Sys.time()

  out = list(networks = networks,
             formula = object,
             main_iters = main_iters,
             control = control,
             prior = prior,
             accepts = accepts,
             params = params,
             model_matrix = model_matrix,
             constraints = constraints,
             time = end_time - start_time)

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

    multibergm.formula(object$formula,
                       object$model_formula,
                       object$constraints,
                       new_iters, 
                       control)
}



get_model_matrix <- function(ergm_formula, model_formula, ...) {
  networks <- statnet.common::eval_lhs.formula(ergm_formula)
  net_attributes <- unique(c(sapply(networks, list.network.attributes)))
  
  net_df <- lapply(networks, function(x) {
    sapply(net_attributes, function(y) x %n% y)
  })
  net_df <- as.data.frame(do.call(rbind, net_df))
  
  model.matrix(model_formula, net_df, ...)
}

#' Set multibergm MCMC options
#'
#' Auxiliary function used to specify settings controlling multibergm fitting.
#' This allows the user to control various aspects of the fitting algorithm,
#' including MCMC proposal covariance, auxiliary iterations for ERGM simulation,
#' and number of cores to be used in parallel.
#'
#' @param ergm_formula An \R \code{\link{formula}} object, of the form
#'   \code{y ~ <model terms>}, where \code{y} is a
#'   \code{\link[network]{network}} object or a
#'   \code{\link[ergm]{network.list}} object.
#' @param constraints A one-sided formula specifying one or more constraints
#'   on the support of the distribution of the networks being simulated.
#' @param groups A vector of group memberships
#' @param aux_iters Number of internal (auxiliary) MCMC iterations used to
#'   simulate a network from the model.
#' @param n_batches Number of computing cores available to simulate networks
#'   from the model in parallel.
#' @param init_proposals A list of covariance matrices used to initialise the 
#'   MCMC proposals.
#' @param proposal_update_freq How often to update the MCMC proposal 
#'   covariances.
#' @param proposal_update_max When to stop adaptive the MCMC proposals.
#' @param proposal_rescale When to reset adaptation following an initial burn_in
#'
#' @return A list containing the following control parameters:
#'   \itemize{
#'      \item \code{aux_iters}: Number of internal (auxiliary) MCMC iterations
#'      used to simulate a network from the model
#'      \item \code{proposal}: A list of covariance matrices used to
#'      generate the MCMC proposals.
#'      \item \code{batches}: Used to batch parallel runs of ergm simulation
#'      \item \code{model}: Internal representation of the ergm network model
#'      \item \code{clists}: ERGM parameters passed to the ergm simulation
#'      function
#'      \item \code{mh_proposals}: ERGM parameters passed to the ergm simulation
#'      function
#'   }
#'
#' @import ergm
#' @importFrom stats terms
#' @export

control_multibergm <- function(ergm_formula,
                               mod_mat,
                               constraints   = ~ . ,
                               proposal_update_freq = 20,
                               proposal_update_max = 1000,
                               proposal_rescale  = 500,
                               init_proposals   = NULL,
                               aux_iters     = 1000, 
                               n_batches     = 1) {

  networks <- statnet.common::eval_lhs.formula(ergm_formula)

  # For compatibility with a single network
  if (is.network(networks))  networks <- list(networks)

  if (!is.network(networks[[1]]))
    stop("Input must be a network or list of networks")
  
  n_nets <- length(networks)
  
  # Set up ergm parameters
  model       <- ergm_model(ergm_formula, networks[[1]])
  clists      <- lapply(networks, function(x) ergm.Cprepare(x, model))
  mh_proposals <- ergm_proposal(constraints, control.ergm()$MCMC.prop.args,
                                networks[[1]])
  n_terms    <- nparam(model)
  n_vars     <- ncol(mod_mat)
  etamap     <- ergm.etamap(model)
  
  # Set default settings for proposals unless explicitly specified
  if (is.null(init_proposals$theta)) {
    init_proposals$theta <- array(NA, c(n_nets, n_terms, n_terms))
    for (i in 1:n_nets) {
      init_proposals$theta[i, , ] <- diag(0.1^2, n_terms)
    }
  }
  
  #TODO: FIX THIS
  if (is.null(init_proposals$mu)) {
    init_proposals$mu    <- array(NA, c(n_vars, n_terms, n_terms))
      for (i in 1:n_vars) {
        init_proposals$mu[i, , ] <- diag(0.1^2, n_terms)
      }
  }

  # Specify network batching for parallelisation
  if (n_batches == 1) {
    batches <- list(seq_len(n_nets))
  } else {
    batches <- split(seq_len(n_nets), cut(seq_len(n_nets),
                                      n_batches, labels = FALSE))
  }

  # Set up ergm parameters
  model       <- ergm_model(ergm_formula, networks[[1]])
  Clists      <- lapply(networks, function(x) ergm.Cprepare(x, model))
  MHproposals <- ergm_proposal(constraints, control.ergm()$MCMC.prop.args,
                               networks[[1]])

  list(aux_iters   = aux_iters,
       init_proposals       = init_proposals,
       proposal_update_freq = proposal_update_freq,
       proposal_update_max  = proposal_update_max,
       proposal_rescale     = proposal_rescale,
       batches     = batches,
       model       = model,
       clists      = clists,
       mh_proposals = mh_proposals,
       mod_mat     = mod_mat,
       etamap               = etamap)

}

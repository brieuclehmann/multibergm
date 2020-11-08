#' Set multibergm MCMC options
#'
#' Auxiliary function used to specify settings controlling multibergm fitting.
#' This allows the user to control various aspects of the fitting algorithm,
#' including MCMC proposal covariance, auxiliary iterations for ERGM simulation,
#' and number of cores to be used in parallel.
#'
#' @param formula An \R \code{\link{formula}} object, of the form
#'   \code{y ~ <model terms>}, where \code{y} is a
#'   \code{\link[network]{network}} object or a
#'   \code{\link[ergm]{network.list}} object.
#' @param constraints A one-sided formula specifying one or more constraints
#'   on the support of the distribution of the networks being simulated.
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
                               constraints = ~.,
                               proposal    = NULL,
                               aux_iters   = 1000,
                               n_batches   = 1) {

  networks <- statnet.common::eval_lhs.formula(formula)

  # For compatibility with a single network
  if (is.network(networks))  networks <- list(networks)

  if (!is.network(networks[[1]]))
    stop("Input must be a network or list of networks")

  n_nets <- length(networks)
  n_terms    <- length(attr(terms(formula), "term.labels"))

  # Set default settings for proposals unless explicitly specified
  if (is.null(proposal$theta)) proposal$theta <- diag(0.1^2, n_terms)
  if (is.null(proposal$mu))    proposal$mu    <- diag(0.1^2, n_terms)

  # Specify network batching for parallelisation
  if (n_batches == 1) {
    batches <- list(seq_len(n_nets))
  } else {
    batches <- split(seq_len(n_nets), cut(seq_len(n_nets),
                                      n_batches, labels = FALSE))
  }

  # Set up ergm parameters
  model       <- ergm_model(formula, networks[[1]])
  clists      <- lapply(networks, function(x) ergm.Cprepare(x, model))
  mh_proposals <- ergm_proposal(constraints, control.ergm()$MCMC.prop.args,
                               networks[[1]])

  list(aux_iters    = aux_iters,
       proposal     = proposal,
       batches      = batches,
       model        = model,
       clists       = clists,
       mh_proposals = mh_proposals)

}

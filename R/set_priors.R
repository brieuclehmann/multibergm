# Set priors ==================================================================
#' Set priors for multibergm
#'
#' This is a utility function to set priors for a multibergm. When not
#' explicitly specified, set_priors() uses default (vague) priors. A
#' validity check is also performed to ensure the compatibility of the priors
#' with the model.
#' 
#' @inheritParams multibergm
#' 
#' @param ergm_formula An \R \code{\link{formula}} object, of the form
#'   \code{y ~ <model terms>}, where \code{y} is a
#'   \code{\link[network]{network}} object or a
#'   \code{\link[ergm]{network.list}} object.
#' @param prior A list of explicit prior specifications.
#' 
#' @export
#' @importFrom statnet.common eval_lhs.formula

set_priors <- function(ergm_formula,
                       model_matrix,
                       control,
                       prior = list()) {

  n_terms <- length(attr(terms(ergm_formula), "term.labels"))
  n_vars  <- ncol(model_matrix)
  
  # Model: \theta = X\mu + \epsilon
  
  # Set priors for \mu parameter
  if (is.null(prior$mu$mean))
    prior$mu$mean <- matrix(0, n_vars, n_terms)

  if (is.null(prior$mu$cov)) {
    #prior$mu$cov <- array(diag(100, n_terms), c(n_terms, n_terms, n_vars))
    prior$mu$cov <- diag(100, n_vars)
  }

  # Set priors for network-level covariance parameter
  if (is.null(prior$cov_theta$df))
    prior$cov_theta$df <- n_terms + 1

  if (is.null(prior$cov_theta$scale))
    prior$cov_theta$scale <- diag(n_terms)

  # if (!is.null(hyper_model_matrix)) {
  #   # Set priors for group-level covariance parameter
  #   #TODO: figure this out...
  #   n_hyper <- ncol(hyper_model_matrix)
  #   if (is.null(prior$cov_mu$df))
  #     prior$cov_mu$df <- rep(n_terms + 1, n_hyper)
  # 
  #   if (is.null(prior$cov_mu$scale))
  #     prior$cov_mu$scale <- array(diag(n_terms), c(n_terms, n_terms, n_hyper))
  # }

  #check_prior(prior, n_terms, n_groups)

  prior
}

#' Check validity of multibergm prior
#'
#' Internal functions to check compatibility of the prior with the model.
#'
#' @inheritParams set_priors
#'
#' @param n_terms Number of terms (summary statistics) in the exponential random
#'     graph model
#' @param n_groups Number of distinct groups
#' @param x Prior mean or covariance to be checked

check_prior <- function(prior, n_terms, n_groups) {

  check_prior_mean(prior$mu$mean, n_terms)

  check_prior_cov(prior$mu$cov, n_terms)

  check_prior_df(prior$cov_theta$df, n_terms)

  check_prior_scale(prior$cov_theta$scale, n_terms)

  if (n_groups > 1) {
    check_prior_df(prior$cov_mu$df, n_terms)

    check_prior_scale(prior$cov_mu$scale, n_terms)
  }

  invisible(NULL)
}

#' @rdname check_prior
check_prior_mean <- function(x, n_terms) {

  obj_name <- deparse(substitute(x))

  if (is.null(x))
    stop(paste(obj_name, "must be specified."))

  if (length(x) != n_terms | !is.numeric(x)) {
    stop(paste(obj_name, "must be a vector of length n_terms =", n_terms))
  }


}

#' @rdname check_prior
is_covmat <- function(x) {

  out <- TRUE
  if (!is.matrix(x) | !isSymmetric(x) | !is.numeric(x)) {
    return(FALSE)
  }

  eigenvalues <- eigen(x, symmetric = TRUE, only.values = TRUE)$values

  if (any(eigenvalues <= 0)) {
    out <- FALSE
  }

  out
}

#' @rdname check_prior
check_prior_cov <- function(x, n_terms) {

  obj_name <- deparse(substitute(x))

  if (is.null(x)) {
    stop(paste(obj_name, "must be specified."))
  }

  if (any(dim(x) != c(n_terms, n_terms)) || !is_covmat(x)) {
    stop(paste(obj_name, "must be a positive semidefinite matrix of",
               "dimension (n_terms, n_terms) = (", n_terms, n_terms, ")"))
  }

  invisible(NULL)
}

#' @rdname check_prior
check_prior_scale <- function(x, n_terms) {

  obj_name <- deparse(substitute(x))

  if (is.null(x)) {
    stop(paste(obj_name, "must be specified."))
  }

  if (!is_covmat(x) | any(dim(x) != c(n_terms, n_terms))) {
    stop(paste(obj_name, "must be a positive definite matrix of",
               "dimension (n_terms, n_terms) = (", n_terms, n_terms, ")"))
  }

  invisible(NULL)
}

#' @rdname check_prior
check_prior_df <- function(x, n_terms) {

  is_scalar <- function(x) is.atomic(x) && length(x) == 1L && is.numeric(x)

  obj_name <- deparse(substitute(x))

  if (is.null(x)) {
    stop(paste(obj_name, "must be specified."))
  }

  if (x <= n_terms - 1 | !is_scalar(x)) {
    stop(paste(obj_name, "must be a scalar greater than n_terms - 1 =",
               n_terms - 1))
  }

  invisible(NULL)
}

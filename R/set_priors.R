# Set priors ==================================================================
#' Set priors for multibergm
#' 
#' This is a utility function to set priors for a multibergm. When not 
#' explicitly specified, set_priors() uses default (vague) priors. A
#' validity check is also performed to ensure the compatibility of the priors
#' with the model.
#' 
#' @param formula An \R \code{\link{formula}} object, of the form
#'   \code{y ~ <model terms>}, where \code{y} is a
#'   \code{\link[network]{network}} object or a
#'   \code{\link[ergm]{network.list}} object.
#' @param groups A vector of group memberships
#' @param prior A list of explicit prior specifications.
#' 
#' @export
#' @importFrom statnet.common eval_lhs.formula

set_priors <- function(formula, groups = NULL, prior = list()) {

  n_terms  <- length(attr(terms(formula), "term.labels"))
  if (is.null(groups)) {
    networks <- statnet.common::eval_lhs.formula(formula)
    groups <- rep(1, length(networks))
  }
  n_groups <- length(unique(groups))

  # Set priors for population mean parameter
  if (is.null(prior$muPop$mean))
    prior$muPop$mean <- rep(0, n_terms)

  if (is.null(prior$muPop$cov))
    prior$muPop$cov <- diag(100, n_terms)

  # Set priors for network-level covariance parameter
  if (is.null(prior$covTheta$df))
    prior$covTheta$df <- n_terms + 1

  if (is.null(prior$covTheta$scale))
    prior$covTheta$scale <- diag(n_terms)

  if (n_groups > 1) {
    # Set priors for group-level covariance parameter

    if (is.null(prior$covMuGroup$df))
      prior$covMuGroup$df <- n_terms + 1

    if (is.null(prior$covMuGroup$scale))
      prior$covMuGroup$scale <- diag(n_terms)
  }

  check_prior(prior, n_terms, n_groups)

  prior
}

#' Check validity of multibergm prior
#' 
#' Internal functions to check compatibility of the prior with the model.

check_prior <- function(prior, n_terms, n_groups) {

  check_prior_mean(prior$muPop$mean, n_terms)

  check_prior_cov(prior$muPop$cov, n_terms)

  check_prior_df(prior$covTheta$df, n_terms)

  check_prior_scale(prior$covTheta$scale, n_terms)

  if (n_groups > 1) {
    check_prior_df(prior$covMuGroup$df, n_terms)

    check_prior_scale(prior$covMuGroup$scale, n_terms)
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
is_covmat <- function(x, semi = FALSE) {

  if (!is.matrix(x) | !isSymmetric(x) | !is.numeric(x)) {
    return(FALSE)
  }

  eigenvalues <- eigen(x, symmetric = TRUE, only.values = TRUE)$values

  if (any(eigenvalues < 0)) {
    return(FALSE)
  }

  if (semi & any(eigenvalues == 0)) {
    return(FALSE)
  }

  TRUE
}

#' @rdname check_prior
check_prior_cov <- function(x, n_terms) {

  obj_name <- deparse(substitute(x))

  if (is.null(x)) {
    stop(paste(obj_name, "must be specified."))
  }

  if (any(dim(x) != c(n_terms, n_terms)) || !is_covmat(x, semi = TRUE)) {
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

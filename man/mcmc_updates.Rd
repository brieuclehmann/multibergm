% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_updates.R
\name{mcmc_updates}
\alias{mcmc_updates}
\alias{exchange_update_cp}
\alias{mean_update}
\alias{cov_update}
\alias{ergm_wrapper}
\title{Internal functions for MCMC updates}
\usage{
exchange_update_cp(curr, prop, prior_mean, delta, prior_cov, etamap)

mean_update(obs_data, obs_cov, prior_mean, prior_cov, model_matrix)

cov_update(
  obs_data,
  obs_mean,
  prior_df,
  prior_scale,
  labels = rep(1, nrow(obs_data)),
  curr_cov = NULL
)

ergm_wrapper(coefs, control)
}
\arguments{
\item{curr}{A vector or matrix of current mean parameter values.}

\item{prop}{A vector or matrix of proposed mean parameter values}

\item{prior_mean}{Prior mean of parameter}

\item{delta}{Change in summary statistics as produced by
\code{\link{ergm_wrapper}}}

\item{prior_cov}{Prior covariance of parameter}

\item{etamap}{The list of values that constitutes the theta-> eta mapping and
is returned by ergm.etamap}

\item{obs_data}{Matrix of observed data}

\item{obs_cov}{Fixed (known) covariance of observed data}

\item{prior_df}{Prior degrees of freedom in Inverse-Wishart}

\item{labels}{Labels to associate each observation with a grouping}

\item{curr_cov}{Current value of covariance parameter}

\item{coefs}{Matrix of coefficients to simulate ERGs from}

\item{control}{A list of parameters set by \code{\link{control_multibergm}}
containing settings for the ERGM simulations}

\item{prior_obs_scale}{Prior scale matrix in Inverse-Wishart}

\item{var_mean}{Fixed (known) mean of parameters}
}
\description{
Functions used internally by to perform Gibbs or exchange updates.
}
\section{Functions}{
\itemize{
\item \code{mean_update}: Gibbs update of a mean parameter

\item \code{cov_update}: Gibbs update of a covariance parameter

\item \code{ergm_wrapper}: Gibbs update of a mean parameter
}}


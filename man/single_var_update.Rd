% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multibergm_updates.R
\name{single_var_update}
\alias{single_var_update}
\title{MCMC update for single-group multi-BERGM}
\usage{
single_var_update(curr, prior, groups, proposals, control)
}
\arguments{
\item{curr}{A list of the current values of the model parameters}

\item{prior}{A list of explicit prior specifications.}

\item{groups}{A vector of group memberships}

\item{proposals}{A list of the current RW proposal parameters}

\item{control}{A list of parameters set by \code{\link{control_multibergm}}
specifying priors, proposal variances, and group labels.}
}
\value{
A list of the updated values of the model parameters and the
  acceptance counts for the exchange updates.
}
\description{
Function used to perform one exchange-within-Gibbs MCMC update for
all the parameters in a single-group multi-BERGM. This applies a standard
Gibbs update for the network-level covariance parameter before using the
exchange algorithm within an Ancillarity-Sufficiency Interweaving Strategy
(ASIS) to update the remaining parameters.
}

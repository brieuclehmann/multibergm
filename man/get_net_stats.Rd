% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gof.multibergm.R
\name{get_net_stats}
\alias{get_net_stats}
\alias{get_net_stats.list}
\alias{get_net_stats.network}
\alias{get_net_stats.network.list}
\alias{get_model_stats}
\alias{get_deg_dist}
\alias{get_esp_dist}
\alias{get_geodist_dist}
\title{Compute network statistics}
\usage{
get_net_stats(object, ...)

\method{get_net_stats}{list}(object, formula, which_stats, ...)

\method{get_net_stats}{network}(object, formula, which_stats, ...)

\method{get_net_stats}{network.list}(object, formula, which_stats, ...)

get_model_stats(network, formula)

get_deg_dist(network)

get_esp_dist(network)

get_geodist_dist(network)
}
\arguments{
\item{object}{A network or list of networks}

\item{...}{Arguments to be passed to methods.}

\item{formula}{The ERGM formula containing the summary statistics to be
computed}

\item{which_stats}{A string specifying which network statistics to be
computed ("all", "model", "gof", or "other").}

\item{network}{A network object}
}
\description{
GetNetStats is an internal function used to compute network statistics.

Internal function to get the model summary statistics for a single network

Internal function to get the degree distribution for a single network

Internal function to get the edgewise-shared partner distribution for a single network

Internal function to get the geodesic distance distribution for a single network
}
\section{Methods (by class)}{
\itemize{
\item \code{list}: Network statistics for a list of networks.

\item \code{network}: Network statistics for a single network

\item \code{network.list}: Network statistics for a network.list
}}


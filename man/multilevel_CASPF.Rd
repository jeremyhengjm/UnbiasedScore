% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multilevel_CASPF.R
\name{multilevel_CASPF}
\alias{multilevel_CASPF}
\title{Multilevel Conditional Ancestor Sampling Particle Filter}
\usage{
multilevel_CASPF(
  model,
  theta,
  discretization,
  observations,
  nparticles,
  resampling_threshold,
  coupled_resampling,
  ref_trajectory_coarse = NULL,
  ref_trajectory_fine = NULL,
  treestorage = FALSE
)
}
\arguments{
\item{model}{a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}}

\item{theta}{a vector of parameters as input to model functions}

\item{discretization}{lists containing stepsize, nsteps, statelength, obstimes for fine and coarse levels, 
and coarsetimes of length statelength_fine indexing time steps of coarse level}

\item{observations}{a matrix of observations, of size nobservations x ydimension}

\item{nparticles}{number of particles}

\item{resampling_threshold}{ESS proportion below which resampling is triggered (always resample at observation times by default)}

\item{coupled_resampling}{a 2-way coupled resampling scheme, such as \code{\link{coupled2_maximal_independent_residuals}}}

\item{ref_trajectory_coarse}{a matrix of reference trajectory for coarser discretization level, of size xdimension x statelength_coarse}

\item{ref_trajectory_fine}{a matrix of reference trajectory for finer discretization level, of size xdimension x statelength_fine}

\item{treestorage}{logical specifying tree storage of Jacob, Murray and Rubenthaler (2013);
if missing, this function store all states and ancestors}
}
\value{
two new trajectories stored as matrices of size xdimension x statelength_coarse/fine.
}
\description{
Runs two coupled conditional particle filters (one at each discretization level) with ancestor sampling (Lindsten, Jordan and Schon, 2014).
}

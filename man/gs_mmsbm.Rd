% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{gs_mmsbm}
\alias{gs_mmsbm}
\title{Wrapper of MCMC algorithm for mixed membership stochastic block model (MMSBM)}
\usage{
gs_mmsbm(g, alpha_init, s_alpha_init, a, b, K, n_swap, N, thin, burnin,
  print_freq, seed, f)
}
\arguments{
\item{g}{a graph object, assumed to be a DAG}

\item{alpha_init}{initial value of alpha}

\item{s_alpha_init}{initial value of proposal standard deviation for alpha's Metropolis step}

\item{a, b}{hyperparameters for alpha's prior}

\item{K}{number of clusters, fixed}

\item{n_swap}{number of swaps in order per iteration}

\item{N, thin, burnin, print_freq}{MCMC quantities}

\item{seed}{random seed to initiate MCMC}

\item{f}{MCMC algo in C++; currently supports rgs_mmsbm()}
}
\value{
a list of data frames and matrices representing the MCMC output
}
\description{
Wrapper of MCMC algorithm for mixed membership stochastic block model (MMSBM)
}

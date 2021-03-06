% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baptrapezoid.R
\name{laplaceBAPmontecarlo.trapezoid}
\alias{laplaceBAPmontecarlo.trapezoid}
\title{Laplace transform for BAP process by simulation with trapezoid rule.}
\usage{
laplaceBAPmontecarlo.trapezoid(X0, horizon, mu, kappa, sigma, lambda = 0,
  zeta = 1, w = -1, trials = 10000, ticks = 100)
}
\arguments{
\item{X0}{Initial condition.}

\item{horizon}{Scalar horizon time.}

\item{mu}{Drift of process.}

\item{kappa}{Mean reversion of process.}

\item{sigma}{Volatility of process.}

\item{lambda}{Jump arrival rate of process.}

\item{zeta}{Mean of jump size distribution.}

\item{w=-1}{Auxiliary parameter in transform.}

\item{trials=10000}{Number of trials in simulation}

\item{ticks=100}{Number of ticks in horizon}
}
\value{
psi vector of equal length to tt.
}
\description{
Laplace transform of Y[t]
}

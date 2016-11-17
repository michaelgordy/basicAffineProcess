#' basicAffineProcess: A package of tools for the Basic Affine Process.
#'
#' The basicAffineProcess package provides Laplace transforms and
#' random number generators for the basic affine process and its important
#' special case of the CIR process.
#'
#' A basic affine process (BAP) follows stochastic differential equation
#' \deqn{dX[t] = (\mu-\kappa X[t]) dt + \sigma X[t]^{1/2} dW[t] + dJ[t]}
#' where W[t] is a Brownian motion and J[t] is a compound Poisson process
#' with arrival rate \eqn{\lambda} and jump sizes distributed exponential with
#' mean \eqn{\zeta}.
#'
#' We have separate routines for two special cases:
#' \itemize{
#'  \item{\eqn{\lambda=0}: }{X[t] is a CIR process.}
#'  \item{\eqn{\sigma=0}: }{X[t] is a mean-reverting compound Poisson (MRCP) process.}
#' }
#'
#' Let Y[t] be the time-integral of X[s] over \eqn{s=(0,t)}.  Notation X[t] and Y[t]
#' and parameters \eqn{(\mu,\kappa,\sigma,\lambda,\zeta)} are used throughout
#' the documentation of this package.
#'
#' @section Laplace transform functions:
#' The extended transform is
#' \deqn{\psi(t;u,w,X0) = E[\exp(wY[t]+uX[t])| X0]}
#' Given a term-structure vector of times tt, the Laplace transform returns
#' a list (tt, A0, B0, A1, B1), all vectors of the same length. The returns
#' are defined by
#' \eqn{\psi(t;0,w,X0)=exp(A0+B0*X0)} and
#' \deqn{\psi_u(t;0,w,X0)=exp(A0+B0*X0)(A1+B1*X0)}
#' where \eqn{\psi_u} is the partial derivative with respect to u.
#' Equivalently, A1 is the time-derivative of A0, and B1 the time-derivative
#' of B0.
#'
#' @section Simulation functions:
#' The simulation functions ...
#'
#' @section Cumulants and multi-dimensional functions:
#' From my paper ...
#'
#' @docType package
#' @name basicAffineProcess
NULL

#' Generic CIR function for inheritance of parameter documentation
#'
#' @param mu Drift of process.
#' @param kappa Mean reversion of process.
#' @param sigma Volatility of process.
#' @return NULL
.cirGenericDummyFunction <- function(mu,kappa,sigma)
  NULL

#' Generic BAP function for inheritance of parameter documentation
#'
#' @inheritParams .cirGenericDummyFunction
#' @param lambda Jump arrival rate of process.
#' @param zeta Mean of jump size distribution.
#' @return NULL
.bapGenericDummyFunction <- function(mu,kappa,sigma,lambda,zeta)
  NULL

#' Generic CIR function for inheritance of parameter documentation
#'
#' @param mu Drift of process.
#' @param kappa Mean reversion of process.
#' @param lambda Jump arrival rate of process.
#' @param zeta Mean of jump size distribution.
#' @return NULL
.mrcpGenericDummyFunction <- function(mu,kappa,lambda,zeta)
  NULL

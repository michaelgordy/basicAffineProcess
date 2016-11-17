# These routines use trapezoid rule to integrate X[t] to get Y[t]

#' Draw (X[t],Y[t]) jointly from BAP transition density.
#'
#' Integral Y[t] obtained via simple trapezoid rule.
#'
#' @param X0 Initial condition X=X0
#' @inheritParams .bapGenericDummyFunction
#' @param horizon time t
#' @param ticks=100 is approximate number of knot points in [0,horizon]
#' @return List (X=Xt,Y=Yt) of process at time t.
rBAPtransition2.trapezoid <- function(X0,mu,kappa,sigma,lambda=0,zeta=1,horizon,ticks=100) {
  dtguide <- horizon/ticks  # maximum dt
  # Draw jump times and sizes.  Final jump size is zero because it is beyond horizon.
  jumptimes <- rPoissonTimes(lambda,horizon)
  jumpsizes <- c(rexp(length(jumptimes)-1,1/zeta),0)
  # Between jump times, we have CIR process
  Yt <- jumplag <- 0
  Xlag <- X0
  for (j in seq_along(jumptimes)) {
    horizonj <- min(jumptimes[j],horizon)-jumplag
    ticksj <- ceiling(horizonj/dtguide)
    XXtt <- rCIRpath(Xlag,mu,kappa,sigma,horizonj,ticksj)
    Xt <- tail(XXtt$Xt,1)
    Yt <- Yt + (sum(XXtt$Xt) - Xt/2 - Xlag/2)*(horizonj/ticksj)
    # Now include the jump at the end of the segment
    Xt <- Xlag <- Xt + jumpsizes[j]
    jumplag <- jumptimes[j]
  }
  return(list(X=Xt,Y=Yt))
}

#' Laplace transform for BAP process by simulation with trapezoid rule.
#'
#' Laplace transform of Y[t]
#'
#' @param X0 Initial condition.
#' @param horizon Scalar horizon time.
#' @inheritParams .bapGenericDummyFunction
#' @param w=-1 Auxiliary parameter in transform.
#' @param trials=10000 Number of trials in simulation
#' @param ticks=100 Number of ticks in horizon
#' @return psi vector of equal length to tt.
#' @export
laplaceBAPmontecarlo.trapezoid <- function(X0,horizon,mu,kappa,sigma,lambda=0,zeta=1,w=-1,trials=10000,ticks=100) {
  psi <- replicate(trials,
                   exp(w*rBAPtransition2.trapezoid(X0,mu,kappa,sigma,lambda,zeta,horizon,ticks=100)$Y))
  return(mean(psi))
}

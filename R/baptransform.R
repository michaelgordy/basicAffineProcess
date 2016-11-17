#' Laplace transform for CIR process
#'
#' Laplace transform of Y[t] and time-derivative of transform.
#'
#' @param tt Vector of times.
#' @inheritParams .cirGenericDummyFunction
#' @param w=-1 Auxiliary parameter in transform.
#' @return List with (tt, A0, B0, A1, B1), all vectors of equal length.
#' @export
laplaceCIR <- function(tt,mu,kappa,sigma,w=-1)
   laplaceBAP(tt,mu,kappa,sigma,lambda=0,zeta=1,w)

#' Laplace transform for MRCP process
#'
#' Laplace transform of Y[t] and time-derivative of transform.
#'
#' @param tt Vector of times.
#' @inheritParams .mrcpGenericDummyFunction
#' @param w=-1 Auxiliary parameter in transform.
#' @return List with (tt, A0, B0, A1, B1), all vectors of equal length.
#' @export
laplaceMRCP <- function(tt,mu,kappa,lambda,zeta,w=-1)
  laplaceBAP(tt,mu,kappa,sigma=0,lambda,zeta,w)


#' Laplace transform for BAP
#'
#' Laplace transform of Y[t] and time-derivative of transform.
#'
#' @param tt Vector of times.
#' @inheritParams .bapGenericDummyFunction
#' @param w=-1 Auxiliary parameter in transform.
#' @return List with (tt, A0, B0, A1, B1), all vectors of equal length.
#' @export
laplaceBAP <- function(tt,mu,kappa,sigma,lambda=0,zeta=1,w=-1) {
  # define constants
  gma <- sqrt(kappa^2-2*w*sigma^2)
  c1 <- 0.5*(gma+kappa)/w
  d1 <- 0.5*(gma-kappa)/w
  if (gma==0) {
    # Then sigma==kappa==0. Take limits.
    g1 <- rep(1,length(tt))
    ftt <- tt
    gtt <- tt^2/2
  } else {
    g1 <- w*(c1+d1*exp(-gma*tt))/gma
    ftt <- (1-exp(-gma*tt))/gma
    gtt <- (tt-ftt/g1)/kappa
  }
  # First handle diffusive component
  if (sigma==0) {
    A0 <- w*mu*gtt
  } else {
    h1 <- -2/sigma^2
    A0 <- mu*(h1*log(g1)+tt/c1)
  }
  A1 <- w*mu*ftt/g1  # because g1dot <- ftt/h1
  B0 <- w*ftt/g1
  B1 <- (1+B0*d1)*w*exp(-gma*tt)/g1
  if (lambda>0) {
    h2overjs <- 1/(w*(c1-zeta)*(d1+zeta))  # h2/zeta
    g2 <- g1-w*zeta*ftt
    A0 <- A0 + lambda*zeta*(h2overjs*log(g2)+tt/(c1-zeta))
    A1 <- A1 + w*lambda*zeta*ftt/g2 # because g2dot <- ftt/h2overjs
  }
  return(list(tt=tt, A0=A0, B0=B0, A1=A1, B1=B1))
}

#' Laplace transform for CIR process
#'
#' Laplace transform of Y[t] and time-derivative of transform.
#'
#' @param tt Vector of times.
#' @inheritParams .cirGenericDummyFunction
#' @param w=-1 Auxiliary parameter in transform.
#' @return List with (tt, A0, B0, A1, B1), all vectors of equal length.
#' @export
laplaceCIR <- function(tt,mu,kappa,sigma,w=-1) {
  if (sigma==0) {
     # handle limiting cases for ftt, gtt
     ftt <- ifelse(kappa==0,tt,(1-exp(-kappa*tt))/kappa)
     gtt <- ifelse(kappa==0,tt^2/2,(tt-ftt)/kappa)
     A0 <- w*mu*gtt
     B0 <- w*ftt
     A1 <- w*mu*ftt
     B1 <- w*exp(-kappa*tt)
  } else {
    gma <- sqrt(kappa^2-2*w*sigma^2)
    c1 <- -sigma^2/(gma-kappa)
    d1 <- -sigma^2/(gma+kappa)

    h1 <- -2/sigma^2
    p1 <- c1+d1*exp(-gma*tt)
    g1 <- w*p1/gma

    A0 <- mu*(h1*log(g1)+tt/c1)
    B0 <- (1-exp(-gma*tt))/p1

    q1 <- gma*exp(-gma*tt)/p1
    A1 <- mu*(-h1*d1*q1 + 1/c1)
    B1 <- q1*(1+B0*d1)
  }

  return(list(tt=tt, A0=A0, B0=B0, A1=A1, B1=B1))
}

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
  c1 <- ifelse(sigma>0,-sigma^2/(gma-kappa),kappa/w)
  d1 <- -sigma^2/(gma+kappa)
  if (sigma==0) {
    # handle limiting cases for ftt, gtt
    ftt <- ifelse(kappa==0,tt,(1-exp(-kappa*tt))/kappa)
    gtt <- ifelse(kappa==0,tt^2/2,(tt-ftt)/kappa)
    A0 <- w*mu*gtt
    B0 <- w*ftt
    A1 <- w*mu*ftt
    B1 <- w*exp(-kappa*tt)
  } else {
    h1 <- -2/sigma^2
    p1 <- c1+d1*exp(-gma*tt)
    g1 <- w*p1/gma
    A0 <- mu*(h1*log(g1)+tt/c1)
    B0 <- (1-exp(-gma*tt))/p1
    q1 <- gma*exp(-gma*tt)/p1
    A1 <- mu*(-h1*d1*q1 + 1/c1)
    B1 <- q1*(1+B0*d1)
  }
  if (lambda>0) {
    a2 <- d1/c1
    c2 <- 1-zeta/c1
    d2 <- (d1+zeta)/c1
    h2overjs <- 1/(w*(c1-zeta)*(d1+zeta))  # h2/zeta
    p2 <- c2+d2*exp(-gma*tt)
    g2 <- p2/(c2+d2)
    A0 <- A0 + lambda*zeta*(h2overjs*log(g2)+tt/(c1-zeta))
    g2dot <- (1-exp(-gma*tt))/(h2overjs*gma)
    A1 <- A1 + w*lambda*zeta*h2overjs*g2dot/g2
  }
  return(list(tt=tt, A0=A0, B0=B0, A1=A1, B1=B1))
}

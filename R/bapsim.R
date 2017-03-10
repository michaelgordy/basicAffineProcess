#' Draw from CIR transition density
#'
#' Draw X[t] given X[0].
#'
#' If \code{X0} is a vector, return \code{Xt} will be a vector of the
#' same length. The vector elements are draw independently. This allows
#' us to simulate many independent paths simulataneously.
#'
#' @param X0 Initial condition.
#' @inheritParams .cirGenericDummyFunction
#' @param horizon time t
#' @return Value Xt of process at time t.
#' @export
rCIRtransition <- function(X0,mu,kappa,sigma,horizon) {
  ea1 <- exp(-kappa*horizon)
  if (sigma==0) {
    if (kappa==0) {
      ftt <- tt
    } else {
      ftt <- (1-ea1)/kappa
    }
     Xt <- mu*ftt + X0*ea1
  } else {
    # This fails to handle kappa==0
     cinv <- 0.5*sigma^2*(1-ea1)/kappa
     nu <- 4*mu/sigma^2
     Xt <- 0.5*cinv*rchisq(n=length(X0),df=nu,ncp=2*X0*ea1/cinv)
  }
  return(Xt)
}

#' Draw path from CIR transition density on [0,horizon]
#'
#' Evenly spaced times points 0,dt,...,horizon, where dt=horizon/ticks.
#'
#' @param X0 Initial condition.
#' @inheritParams .cirGenericDummyFunction
#' @param tt horizon Final time.
#' @param ticks Number of slices in [0,horizon]
#' @return List of (tt,Xt), each of length ticks+1.
#' @export
rCIRpath <- function(X0,mu,kappa,sigma,horizon,ticks) {
  dt <- horizon/ticks
  Xt <- rep(X0,ticks+1)
  for (k in 1:ticks) {
    Xt[k+1] <- rCIRtransition(Xt[k],mu,kappa,sigma,dt)
  }
  return(list(tt=seq(0,horizon,length.out=ticks+1),Xt=Xt))
}

#' Draw ([X[t],Y[t]) jointly from CIR transition density
#'
#' Draw (X[t],Y[t]) given (X[0],Y[0]).
#'
#' @param XY0 Initial condition list(X=X0,Y=Y0)
#' @inheritParams .cirGenericDummyFunction
#' @param horizon time t
#' @return List (X=Xt,Y=Yt) of process at time t.
# #' @export
rCIRtransition2 <- function(XY0,mu,kappa,sigma,horizon) {
  Xt <- rCIRtransition(XY0$X,mu,kappa,sigma,horizon)
  Yt <- XY0$Y + rCIRintegrated(XY0$X,Xt,mu,kappa,sigma,horizon)
  return(list(X=Xt,Y=Yt))
}


#' Draw from CIR stationary distribution.
#'
#' If \code{n}>1, return will be a vector of length \code{n},
#' with vector elements draw independently. This is convenient
#' for simulation of many independent paths simulataneously.
#'
#' @inheritParams .cirGenericDummyFunction
#' @param n=1 number of observations.
#' @return X Value of process.
#' @export
rCIRstationary <- function(mu,kappa,sigma,n=1) {
  rgamma(n,shape=2*mu/sigma^2,scale=sigma^2/(2*kappa))
}

#' Draw Y[t] conditional on (X[0],X[t]).
#'
#' Glasserman and Kim FS 2011, Theorem 2.2
#'
#' @param X0 Initial condition.
#' @param Xt Final value at horizon.
#' @inheritParams .cirGenericDummyFunction
#' @param horizon time t
#' @param K number of terms in G-K expansion
#' @return Value Yt of process at time t.
#' @export
rCIRintegrated <- function(X0,Xt,mu,kappa,sigma,horizon,K=10) {
  if (sigma==0) {
    # handle special case of kappa==0
    if (kappa==0) {
      ftt <- horizon
      gtt <- horizon^2/2
    } else {
      ftt <- (1-exp(-kappa*horizon))/kappa
      gtt <- (tt-ftt)/kappa
    }
    Yt <- mu*gtt + X0*ftt
  } else {
    # Calculate delta, lambda_n and gamma_n series
    ktpinsq <- (kappa^2*horizon^2 + 4*pi^2*(1:K)^2)
    gmma <- ktpinsq/(2*sigma^2*horizon^2)
    lmbda <- (4*pi*(1:K))^2/(sigma^2*horizon*ktpinsq)
    dlta <- 4*mu/sigma^2
    nu <- dlta/2-1
    z <- sqrt(X0*Xt)*(2*kappa/sigma^2)/sinh(kappa*horizon/2)
    eta <- rBessel(nu,z)
    gkX1 <- 0
    for (n in (1:K)) {
      Nn <- rpois(1,(X0+Xt)*lmbda[n])
      gkX1 <- gkX1 + sum(rexp(Nn))/gmma[n]
    }
    gkX2 <- sum(rgamma(K,dlta/2)/gmma)
    # gkZ is vector of length eta, each element of which is sum of eta
    # Ga(2,1)/gmma[n] random variables
    gkZ <- as.vector((1/gmma) %*% matrix(rgamma(eta*K,2),nrow = K))
    # Now we need to sample the remainder terms as in Lemma 3.1
    m1 <- 2*(X0+Xt)*horizon/(pi^2*K)
    scale1 <- (1/3)*(sigma*horizon/(pi*K))^2
    gkX1 <- gkX1 + rgamma(1,shape=m1/scale1,scale=scale1)
    m2 <- (dlta/K)*(sigma*horizon/(2*pi))^2
    scaleZ <- scale2 <- 0.5*scale1
    gkX2 <- gkX2 + rgamma(1,shape=m2/scale2,scale=scale2)
    mZ <- (1/K)*(sigma*horizon/pi)^2
    gkZ <- gkZ + rgamma(eta,shape=mZ/scaleZ,scale=scaleZ)
    Yt <- gkX1+gkX2+sum(gkZ)
  }
  return(Yt)
}

#' Draw a Bessel variate.
#'
#' Discrete Bessel distribution of Yuan and Kalbfleish (2000).
#' Method on page 285 of G-K.
#'
#' @param nu First parameter of distribution
#' @param z  Second parameter of distribution
#' @return BES(nu,z) random variate
#' @export
rBessel <- function(nu,z) {
  m <- (sqrt(z^2+nu^2)-nu)/2  # mode
  if (m<10) {
    return(rBessel.exact(nu,z))
  } else {
    return(rBessel.gaussian(nu,z))
  }
}

# For use in rBessel
rBessel.exact <- function(nu,z) {
  U <- runif(1)
  Fn <- pn <- exp(nu*log(z/2)-lgamma(nu+1)-z)/besselI(z,nu,expon.scaled=T)
  # Fn <- pn <- (z/2)^nu/(gamma(nu+1)*besselI(z,nu))
  n <- 0
  while (U>Fn) {
      n <- n+1
      pn <- pn*(z/2)^2/(n*(n+nu))
      Fn <- Fn + pn
  }
  return(n)
}

# Approximation for use in rBessel
rBessel.gaussian <- function(nu,z) {
  U <- runif(1)
  # Get Bessel quotients
  I0 <- besselI(z,nu,expon.scaled=T)
  I1 <- besselI(z,nu+1,expon.scaled=T)
  I2 <- besselI(z,nu+2,expon.scaled=T)
  R0 <- I1/I0
  R1 <- I2/I1
  mu <- R0*(z/2)
  sigma <- sqrt(mu - mu^2*(1-R1/R0))
  X <- mu+sigma*qnorm(U+(1-U)*pnorm(-mu/sigma))
  return(round(X))
}

#' Laplace transform for CIR process by simulation
#'
#' Laplace transform of Y[t]
#'
#' @param X0 Initial condition.
#' @param tt Vector of times.
#' @inheritParams .cirGenericDummyFunction
#' @param w=-1 Auxiliary parameter in transform.
#' @return psi vector of equal length to tt.
# #' @export
laplaceCIRmontecarlo <- function(X0,tt,mu,kappa,sigma,w=-1,trials=10000) {
   psi <- rep(0,length(tt))
   dt <- diff(c(0,tt))
   for (j in 1:trials) {
     Yt <- rep(0,length(tt))
     XY <- list(X=X0,Y=0)
     for (k in 1:length(tt)) {
       XY <- rCIRtransition2(XY,mu,kappa,sigma,dt[k])
       Yt[k] <- XY$Y
     }
     psi <- psi+exp(w*Yt)
   }
   return(psi/trials)
}

#' Laplace transform for BAP process by simulation
#'
#' Laplace transform of Y[t]
#'
#' @param X0 Initial condition.
#' @param tt Vector of times.
#' @inheritParams .bapGenericDummyFunction
#' @param w=-1 Auxiliary parameter in transform.
#' @param trials=10000 Number of trials in simulation
#' @return psi vector of equal length to tt.
laplaceBAPmontecarlo <- function(X0,tt,mu,kappa,sigma,lambda=0,zeta=1,w=-1,trials=10000) {
  psi <- rep(0,length(tt))
  dt <- diff(c(0,tt))
  for (j in 1:trials) {
    Yt <- rep(0,length(tt))
    XY <- list(X=X0,Y=0)
    for (k in 1:length(tt)) {
      XY <- rBAPtransition2(XY,mu,kappa,sigma,lambda,zeta,dt[k])
      Yt[k] <- XY$Y
    }
    psi <- psi+exp(w*Yt)
  }
  return(psi/trials)
}

#' Helper function to draw Poisson jump times out to horizon.
#'
#' Final jump time will be beyond horizon.
#'
#' @param lambda Jump arrival rate of process.
#' @param horizon time t
#' @return Jump times
rPoissonTimes <- function(lambda,horizon) {
  if (lambda>0) {
    jumptimes <- NULL
    tt <- 0
    while (tt<horizon) {
      tt <- tt+rexp(1,lambda)
      jumptimes <- c(jumptimes,tt)
    }
  } else {
    jumptimes <- horizon+1
  }
  return(jumptimes)
}

#' Helper function to draw from BAP transition density, i.e., X[t] given X[0]
#'
#' @param X0 Initial conditions (scalar).
#' @inheritParams .bapGenericDummyFunction
#' @param horizon time t
#' @return Xlist list containing X[horizon], X[jumptimes-], jump times and sizes.
.rBAPdraw <- function(X0,mu,kappa,sigma,lambda,zeta,horizon) {
  # Last element in jumpsizes corresponds to the jump beyond horizon
  jumptimes <- rPoissonTimes(lambda,horizon)
  jumpsizes <- c(rexp(length(jumptimes)-1,1/zeta),0)
  Xminus <- rep(0,length(jumptimes))
  Xlag <- X0
  jumplag <- 0
  for (j in seq_along(Xminus)) {
      dt <- min(jumptimes[j],horizon)-jumplag
      Xminus[j] <- rCIRtransition(Xlag,mu,kappa,sigma,dt)
      jumplag <- jumptimes[j]
      Xlag <- Xminus[j]+jumpsizes[j]
  }
  return(list(Xt=Xlag,Xminus=Xminus,
                      jumpsizes=jumpsizes,jumptimes=jumptimes))
}

#' Draw from MRCP transition density, i.e., X[t] given X[0]
#'
#' @param X0 Initial conditions (scalar).
#' @inheritParams .mrcpGenericDummyFunction
#' @param horizon time t
#' @return Xlist list containing X[horizon], X[jumptimes-], jump times and sizes.
#' @return Xt Value of processes at time t.
#' @export
rMRCPtransition <- function(X0,mu,kappa,lambda=0,zeta=1,horizon) {
  .rBAPdraw(X0,mu,kappa,sigma=0,lambda,zeta,horizon)$Xt
}

#' Draw (X[t],Y[t]) jointly from MRCP transition density.
#'
#' @param XY0 Initial condition list(X=X0,Y=Y0)
#' @inheritParams .mrcpGenericDummyFunction
#' @param horizon time t
#' @return List (X=Xt,Y=Yt) of process at time t.
#' @export
rMCRPtransition2 <- function(XY0,mu,kappa,lambda=0,zeta=1,horizon)
  rBAPtransition2(XY0,mu,kappa,sigma=0,lambda,zeta,horizon)

#' Draw from BAP transition density, i.e., X[t] given X[0]
#'
#' @inheritParams .rBAPdraw
#' @return Xt Value of processes at time t.
#' @export
rBAPtransition <- function(X0,mu,kappa,sigma,lambda=0,zeta=1,horizon) {
  .rBAPdraw(X0,mu,kappa,sigma,lambda,zeta,horizon)$Xt
}

#' Draw (X[t],Y[t]) jointly from BAP transition density.
#'
#' @param XY0 Initial condition list(X=X0,Y=Y0)
#' @inheritParams .bapGenericDummyFunction
#' @param horizon time t
#' @return List (X=Xt,Y=Yt) of process at time t.
#' @export
rBAPtransition2 <- function(XY0,mu,kappa,sigma,lambda=0,zeta=1,horizon) {
  bap1 <- .rBAPdraw(XY0$X,mu,kappa,sigma,lambda,zeta,horizon)
  Xt <- bap1$Xt
  Xlag <- XY0$X
  Yt <- XY0$Y   # initial condition
  jumplag <- 0
  for (j in seq_along(bap1$Xminus)) {
    dt <- min(bap1$jumptimes[j],horizon)-jumplag
    Yt <- XY0$Y + rCIRintegrated(Xlag,bap1$Xminus[j],mu,kappa,sigma,dt)
    jumplag <- bap1$jumptimes[j]
    Xlag <- bap1$Xminus[j]+bap1$jumpsizes[j]
  }
  return(list(X=Xt,Y=Yt))
}



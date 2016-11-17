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
#' @param timestep time t
#' @return Value Xt of process at time t.
#' @export
rCIRtransition <- function(X0,mu,kappa,sigma,timestep) {
  ea1 <- exp(-kappa*timestep)
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

#' Draw ([X[t],Y[t]) jointly from CIR transition density
#'
#' Draw (X[t],Y[t]) given (X[0],Y[0]).
#'
#' @param XY0 Initial condition list(X=X0,Y=Y0)
#' @inheritParams .cirGenericDummyFunction
#' @param timestep time t
#' @return List (X=Xt,Y=Yt) of process at time t.
# #' @export
rCIRtransition2 <- function(XY0,mu,kappa,sigma,timestep) {
  Xt <- rCIRtransition(XY0$X,mu,kappa,sigma,timestep)
  Yt <- XY0$Y + rCIRintegrated(XY0$X,Xt,mu,kappa,sigma,timestep)
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
#' @param Xt Final value at timestep.
#' @inheritParams .cirGenericDummyFunction
#' @param timestep time t
#' @param K number of terms in G-K expansion
#' @return Value Yt of process at time t.
#' @export
rCIRintegrated <- function(X0,Xt,mu,kappa,sigma,timestep,K=10) {
  if (sigma==0) {
    # handle special case of kappa==0
    if (kappa==0) {
      ftt <- timestep
      gtt <- timestep^2/2
    } else {
      ftt <- (1-exp(-kappa*timestep))/kappa
      gtt <- (tt-ftt)/kappa
    }
    Yt <- mu*gtt + X0*ftt
  } else {
    # Calculate delta, lambda_n and gamma_n series
    ktpinsq <- (kappa^2*timestep^2 + 4*pi^2*(1:K)^2)
    gmma <- ktpinsq/(2*sigma^2*timestep^2)
    lmbda <- (4*pi*(1:K))^2/(sigma^2*timestep*ktpinsq)
    dlta <- 4*mu/sigma^2
    nu <- dlta/2-1
    z <- sqrt(X0*Xt)*(2*kappa/sigma^2)/sinh(kappa*timestep/2)
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
    m1 <- 2*(X0+Xt)*timestep/(pi^2*K)
    scale1 <- (1/3)*(sigma*timestep/(pi*K))^2
    gkX1 <- gkX1 + rgamma(1,shape=m1/scale1,scale=scale1)
    m2 <- (dlta/K)*(sigma*timestep/(2*pi))^2
    scaleZ <- scale2 <- 0.5*scale1
    gkX2 <- gkX2 + rgamma(1,shape=m2/scale2,scale=scale2)
    mZ <- (1/K)*(sigma*timestep/pi)^2
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
# #' @export
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

#' Laplace transform for BAP process by simulation with trapezoid rule.
#'
#' Laplace transform of Y[t]
#'
#' @param X0 Initial condition.
#' @param timestep Scalar horizon time.
#' @inheritParams .bapGenericDummyFunction
#' @param w=-1 Auxiliary parameter in transform.
#' @param trials=10000 Number of trials in simulation
#' @param ticks=100 Number of ticks in timestep
#' @return psi vector of equal length to tt.
# #' @export
laplaceBAPmontecarlo.trapezoid <- function(X0,timestep,mu,kappa,sigma,lambda=0,zeta=1,w=-1,trials=10000,ticks=100) {
   psi <- replicate(trials,
            exp(w*rBAPtransition2.trapezoid(X0,mu,kappa,sigma,lambda,zeta,timestep,ticks=100)$Y))
   return(mean(psi))
}

#' Helper function to draw from BAP transition density, i.e., X[t] given X[0]
#'
#' @param X0 Initial conditions (scalar).
#' @inheritParams .bapGenericDummyFunction
#' @param timestep time t
#' @return Xlist list containing X[timestep], X[jumptimes-], jump times and sizes.
.rBAPdraw <- function(X0,mu,kappa,sigma,lambda,zeta,timestep) {
  if(lambda>0) {
    jumptimes <- NULL
    tt <- 0
    while (tt<timestep) {
      tt <- tt+rexp(1,lambda)
      jumptimes <- c(jumptimes,tt)
    }
  } else {
    jumptimes <- timestep+1
  }
  # Last element in jumpsizes corresponds to the jump beyond timestep
  jumpsizes <- c(rexp(length(jumptimes),1/zeta),0)
  Xminus <- rep(0,length(jumptimes))
  Xlag <- X0
  jumplag <- 0
  for (j in seq_along(Xminus)) {
      dt <- min(jumptimes[j],timestep)-jumplag
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
#' @param timestep time t
#' @return Xlist list containing X[timestep], X[jumptimes-], jump times and sizes.
#' @return Xt Value of processes at time t.
#' @export
rMRCPtransition <- function(X0,mu,kappa,lambda=0,zeta=1,timestep) {
  .rBAPdraw(X0,mu,kappa,sigma=0,lambda,zeta,timestep)$Xt
}

#' Draw (X[t],Y[t]) jointly from MRCP transition density.
#'
#' @param XY0 Initial condition list(X=X0,Y=Y0)
#' @inheritParams .mrcpGenericDummyFunction
#' @param timestep time t
#' @return List (X=Xt,Y=Yt) of process at time t.
#' @export
rMCRPtransition2 <- function(XY0,mu,kappa,lambda=0,zeta=1,timestep)
  rBAPtransition2(XY0,mu,kappa,sigma=0,lambda,zeta,timestep)

#' Draw from BAP transition density, i.e., X[t] given X[0]
#'
#' @inheritParams .rBAPdraw
#' @return Xt Value of processes at time t.
#' @export
rBAPtransition <- function(X0,mu,kappa,sigma,lambda=0,zeta=1,timestep) {
  .rBAPdraw(X0,mu,kappa,sigma,lambda,zeta,timestep)$Xt
}

#' Draw (X[t],Y[t]) jointly from BAP transition density.
#'
#' @param XY0 Initial condition list(X=X0,Y=Y0)
#' @inheritParams .bapGenericDummyFunction
#' @param timestep time t
#' @return List (X=Xt,Y=Yt) of process at time t.
#' @export
rBAPtransition2 <- function(XY0,mu,kappa,sigma,lambda=0,zeta=1,timestep) {
  bap1 <- .rBAPdraw(XY0$X,mu,kappa,sigma,lambda,zeta,timestep)
  Xt <- bap1$Xt
  Xlag <- XY0$X
  Yt <- XY0$Y   # initial condition
  jumplag <- 0
  for (j in seq_along(bap1$Xminus)) {
    dt <- min(bap1$jumptimes[j],timestep)-jumplag
    Yt <- XY0$Y + rCIRintegrated(Xlag,bap1$Xminus[j],mu,kappa,sigma,dt)
    jumplag <- bap1$jumptimes[j]
    Xlag <- bap1$Xminus[j]+bap1$jumpsizes[j]
  }
  return(list(X=Xt,Y=Yt))
}

#' Draw (X[t],Y[t]) jointly from BAP transition density.
#'
#' Integral Y[t] obtained via simple trapezoid rule.
#'
#' @param X0 Initial condition X=X0
#' @inheritParams .bapGenericDummyFunction
#' @param timestep time t
#' @param ticks=100 is approximate number of knot points in [0,timestep]
#' @return List (X=Xt,Y=Yt) of process at time t.
# #' @export
rBAPtransition2.trapezoid <- function(X0,mu,kappa,sigma,lambda=0,zeta=1,timestep,ticks=100) {
  dtguide <- timestep/ticks  # maximum dt
  # draw jump times
  if(lambda>0) {
    jumptimes <- NULL
    ts <- 0
    while (ts<timestep) {
      ts <- ts+rexp(1,lambda)
      jumptimes <- c(jumptimes,ts)
    }
  } else {
    jumptimes <- timestep+1
  }
  # Between jump times, we have CIR process
  Yt <- jumplag <- 0
  Xlag <- X0
  nu <- 4*mu/sigma^2
  for (j in seq_along(jumptimes)) {
     timestepj <- min(jumptimes[j],timestep)-jumplag
     ticksj <- ceiling(timestepj/dtguide)
     dtj <- timestepj/ticksj
     ea1 <- exp(-kappa*dtj)
     cinv <- 0.5*sigma^2*(1-ea1)/kappa
     for (k in 1:ticksj) {
        Xt <- 0.5*cinv*rchisq(n=1,df=nu,ncp=2*Xlag*ea1/cinv)
        Yt <- Yt + dtj*(Xt+Xlag)/2
        Xlag <- Xt
     }
     # Now include the jump
     if (jumptimes[j]<=timestep) {
         Xt <- Xlag <- Xlag + rexp(1,1/zeta)
         jumplag <- jumptimes[j]
     }
  }
  return(list(X=Xt,Y=Yt))
}


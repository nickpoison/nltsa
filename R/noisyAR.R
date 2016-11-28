#' Create a Noisy Auto-Regressive(1) NLSS model
#' 
#' @param phi State auto-regressive parameter
#' @param sigmaW Dynamic noise standard deviation
#' @param sigmaV Observation error standard deviation
#' @return S3 object of class NoisyAR
#' @export
NoisyAR <- function(phi, sigmaW, sigmaV) {  
  self <- list(phi = phi, sigmaW = sigmaW, sigmaV = sigmaV)
  attr(self,"spaces") <- list(dimx=1,dimy=1)
  class(self) <- c("list","NoisyAR")
  self
}

#' Checks if parameters are valid for a Noisy AR model
#' 
#' @param nlss Instance of the NoisyAR class
#' @return TRUE if the parameters are valid, FALSE otherwise
#' @export
in.support.NoisyAR <- function(nlss) {
  !(is.nan(nlss$phi)  | is.nan(nlss$sigmaW) | is.nan(nlss$sigmaV)) &
    !(is.na(nlss$phi)  | is.na(nlss$sigmaW) | is.na(nlss$sigmaV)) &
    abs(nlss$phi) < 1 & nlss$sigmaW > 0 & nlss$sigmaV > 0
}

#' Compute the local likelihood of a NoisyAR model
#' 
#' @param self Instance of the NoisyAR class
#' @param x Current state
#' @param y Current observation
#' 
#' @return Logarithm of the local likelihood p(y|x)
#' @export
loclike.logpdf.NoisyAR <- function (self, x, y, ...) { dnorm(y, mean=x, sd=self$sigmaV, log = TRUE) }

loclike.rnd.NoisyAR <- function (self, x, ...) { rnorm(length(x), mean=x, sd=self$sigmaV) }

prior.logpdf.NoisyAR <- function (self, xpast, xnew, ...) { dnorm(xnew, mean=self$phi*xpast, sd=self$sigmaW, log = TRUE) }

prior.rnd.NoisyAR <- function (self, xpast, ...) { rnorm(length(xpast), mean=self$phi*xpast, sd=self$sigmaW) }

initial.rnd.NoisyAR <- function (self, N=1, ...) { rnorm(N, mean=0, sd=self$sigmaW/sqrt(1-self$phi^2)) }

initial.logpdf.NoisyAR <- function (self, x, ...) { dnorm(x, mean=0, sd=self$sigmaW/sqrt(1-self$phi^2), log=TRUE) }

NoisyAR.optimal.params <- function(nlss, xpast, y) {  
  sd <- nlss$sigmaW * nlss$sigmaV / sqrt(nlss$sigmaW^2 + nlss$sigmaV^2)
  mean <- sd^2 * (nlss$phi*xpast / nlss$sigmaW^2 +  y/nlss$sigmaV^2)
  data.frame(mean=mean, sd=sd) 
}

NoisyAR.optimal.rnd <- function(nlss, xpast, y, ...) {
  params <- NoisyAR.optimal.params(nlss, xpast, y)  
  xnew <- rnorm(length(xpast),mean=params$mean,sd=params$sd)
}

NoisyAR.optimal.logpdf <- function(nlss, xpast, xnew, y, ...) {  
  params <- NoisyAR.optimal.params(nlss, xpast, y)
  dnorm(xnew, mean=params$mean, sd=params$sd, log = TRUE)
}

NoisyAR.optimal.logadjweights <- function(nlss, xpast, y, ...) {
  .5*log(nlss$sigmaW^2+nlss$sigmaV^2) - .5*(y - nlss$phi*xpast)^2/(nlss$sigmaW^2+nlss$sigmaV^2)
}

#' @import astsa
Kfilter0.NoisyAR <- function(nlss,y) {  
  var.init <- nlss$sigmaW^2/(1-nlss$phi^2)
  Kfilter0(num = length(y),
           y = y,
           A = 1,
           mu0 = 0,
           Sigma0 =var.init,
           Phi = nlss$phi,
           cQ = nlss$sigmaW,
           cR = nlss$sigmaV)
}

#' @import astsa
Ksmooth0.NoisyAR <- function(nlss,y) {  
  var.init <- nlss$sigmaW^2/(1-nlss$phi^2)
  Ksmooth0(num = length(y),
           y = y,
           A = 1,
           mu0 = 0,
           Sigma0 =var.init,
           Phi = nlss$phi,
           cQ = nlss$sigmaW,
           cR = nlss$sigmaV)
}


KalmanLike.NoisyAR <- function(nlss,y) {
  var.init <- nlss$sigmaW^2/(1-nlss$phi^2)
  KalmanLike(
    y = y,
    mod=list(
      T=nlss$phi, 
      Z=1, 
      h=nlss$sigmaV^2, 
      V=nlss$sigmaW^2, 
      a=0, 
      P=var.init, 
      Pn=var.init
      )
  )
}

#' @import FKF 
fkf.NoisyAR <- function(nlss,y) {
  var.init <- nlss$sigmaW^2/(1-nlss$phi^2)
  fkf(a0=0, 
      P0=matrix(var.init), 
      dt=matrix(0), 
      ct=matrix(0), 
      Tt=matrix(nlss$phi), 
      Zt=matrix(1), 
      HHt=matrix(nlss$sigmaW^2), 
      GGt=matrix(nlss$sigmaV^2), 
      yt=matrix(y,nrow=1))
}

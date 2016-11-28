#' Parameters for ARCH optimal kernel 
#' 
#' Compute the mean and standard deviation of p(X_k|x_{k-1}, y_k) for the ARCH model
#' 
#' @param arch Instance of the ARCH class
#' @param xpast Value of the state at the former time-step
#' @param y Observation at the current time-step
#' @return List with fields
#' \item{mean}{Mean of the Gaussian kernel}
#' \item{sd}{Standard deviation of the Gaussian kernel}
#' @export
ARCH.optimal.params <- function(arch, xpast, y) {
  sigmaW <- sqrt(arch$alpha0+arch$alpha1*xpast^2)
  mean <- (sigmaW^2*y) / (sigmaW^2 + arch$sigmaV^2)
  sd <- sigmaW * arch$sigmaV / sqrt(sigmaW^2 + arch$sigmaV^2)
  list(mean=mean, sd=sd) 
}

#' Sample from ARCH optimal kernel 
#' 
#' Sample from p(X_k|x_{k-1}, y_k) for the ARCH model
#' 
#' @param arch Instance of the ARCH class
#' @param xpast Value of the state at the former time-step
#' @param y Observation at the current time-step
#' @return Sample from p(X_k|x_{k-1} = xpast, y_k = y) 
#' @export
ARCH.optimal.rnd <- function(arch, xpast, y, ...) {
  params <- ARCH.optimal.params(arch, xpast, y)  
  xnew <- rnorm(length(xpast), mean=params$mean, sd=params$sd)
}

#' Log-likelihood of  ARCH optimal kernel 
#' 
#' Evaluate  p(X_|x_{k-1}, y_k) for the ARCH model
#' 
#' @param arch Instance of the ARCH class
#' @param xpast Value of the state at the former time-step
#' @param xnew Value of the state at the current time-step
#' @param y Observation at the current time-step
#' @return Logarithm of the optimal kernel density p(x_k=xnew|x_{k-1}=xpast, y_k=y) 
#' @export
ARCH.optimal.logpdf <- function(arch, xpast, xnew, y, ...) {  
  params <- ARCH.optimal.params(arch, xpast, y)
  dnorm(xnew, mean=params$mean, sd=params$sd, log = TRUE)
}

ARCH.optimal.logadjustment <- function(arch, xpast, y, ...) {
  sigmaW <- sqrt(arch$alpha0+arch$alpha1*xpast^2)
  .5*log(sigmaW^2 + arch$sigmaV^2) -.5 * y^2/(sigmaW^2 + arch$sigmaV^2)
}

ARCH.optimal.initial.params <- function(arch, y) {
  sigmaW <- sqrt(arch$alpha0/(1-arch$alpha1))
  mean <- (sigmaW^2*y) / (sigmaW^2 + arch$sigmaV^2)
  sd <- sigmaW * arch$sigmaV / sqrt(sigmaW^2 + arch$sigmaV^2)
  list(mean=mean, sd=sd) 
}

ARCH.optimal.initial.rnd <- function(arch, N=1, y, ...) {
  params <- ARCH.optimal.initial.params(arch, y)
  rnorm(N, mean=params$mean, sd=params$sd)
}
 
ARCH.optimal.initial.logpdf <- function(arch, xnew, y, ...) {  
  params <- ARCH.optimal.initial.params(arch, y)  
  dnorm(xnew, mean=params$mean, sd=params$sd, log=TRUE)
}


ARCH.laplace.initial.rnd <- function(arch, N=1, y, ...) {  
  require(VGAM, warn.conflicts=FALSE)
  params <- ARCH.optimal.initial.params(arch, y)
  rlaplace(N, location=params$mean, scale=params$sd/sqrt(2))
}

ARCH.laplace.initial.logpdf <- function(arch, xnew, y, ...) {    
  require(VGAM, warn.conflicts=FALSE)
  params <- ARCH.optimal.initial.params(arch, y)  
  dlaplace(xnew, location=params$mean, scale=params$sd/sqrt(2), log=TRUE)
}
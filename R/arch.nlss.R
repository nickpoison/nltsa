#' Create an  ARCH(1) 
#' 
#' Example NLSS model
#' 
#' Xt = sigmaW(Xt-1) Wt 
#' 
#' Yt = Xt + sigmaV Vt, 
#' 
#' Wt and Vt are iid N(0,1) 
#' and sigmaW(x) = sqrt(alpha0 + alpha1*x^2) where alpha0 and alpha1 are positive. 
#' 
#' @param alpha0 Intercept of the variance as affine function of the square of the state 
#' @param alpha1 Linear coefficient of the variance as affine function of the square of the state
#' @param sigmaV Observation error standard deviation
#' @return S3 object of class ARCH
#' @export
ARCH <- function(alpha0, alpha1, sigmaV) {  
  self <- list(alpha0=alpha0, alpha1=alpha1, sigmaV=sigmaV)
  attr(self,"spaces") <- list(dimx=1,dimy=1)
  class(self) <- c("list","ARCH")
  self
}
loclike.logpdf.ARCH <- function (self, x, y, ...) { 
  dnorm(y, mean=x, sd=self$sigmaV, log = TRUE) 
}
loclike.rnd.ARCH <- function (self, x, ...) { 
  rnorm(length(x), mean=x, sd=self$sigmaV) 
}
prior.logpdf.ARCH <- function (self, xpast, xnew, ...) { 
  dnorm(xnew, mean=0, sd=sqrt(self$alpha0+self$alpha1*xpast^2), log = TRUE) 
}
prior.rnd.ARCH <- function (self, xpast, ...) { 
  rnorm(length(xpast), mean=0, sd=sqrt(self$alpha0+self$alpha1*xpast^2)) 
}
initial.rnd.ARCH <- function (self, N=1, ...) { 
  rnorm(N, mean=0, sd=sqrt(self$alpha0/(1-self$alpha1)))
}
initial.logpdf.ARCH <- function (self, x, N=1, ...) { 
  dnorm(x, mean=0, sd=sqrt(self$alpha0/(1-self$alpha1)), log=TRUE)
}
print.ARCH <- function (self, ...) {
  cat("Gaussian ARCH nlss:\n")
  cat("Alpha0:", self$alpha0, "\n");
  cat("Alpha1:", self$alpha1, "\n");
  cat("SigmaV:", self$sigmaV, "\n");
}
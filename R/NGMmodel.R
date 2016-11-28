NGMmodel <- function(alpha, beta, gamma, sig2w, sig2v) {  
   sigmaW = sqrt(sig2w) 
   sigmaV = sqrt(sig2v)
  self <- list(a0=alpha, a1=beta, a2=gamma, sigmaW=sigmaW, sigmaV=sigmaV)
  attr(self, "spaces") <- list(dimx=1, dimy=1, dimalpha=3)
  class(self) <- c("list", "NGM")
  self
}
loclike.h.NGM <- function(self, x, t, ...) {
  .05*x^2
}
prior.f.NGM <- function(self, x, t) {
  f1 <- x
  f2 <- x/(1+x^2)
  f3 <- rep(cos(1.2*(t-1)), length(x))
  as.array(rbind(f1,f2,f3))
}
loclike.logpdf.NGM <- function (self, x, y, ...) { 
  dnorm(y, mean=.05*x^2, sd=self$sigmaV, log = TRUE) 
}
loclike.rnd.NGM <- function (self, x, t, ...) { 
  rnorm(length(x), mean=.05*x^2, sd=self$sigmaV) 
}
meanPrior <- function(self, x, t) {
  self$a0*x + self$a1*x/(1+x^2) + self$a2*cos(1.2*(t-1))
}
prior.logpdf.NGM <- function (self, xpast, xnew, t, ...) { 
  dnorm(xnew, mean=meanPrior(self, xpast, t), sd=self$sigmaW, log = TRUE) 
}
prior.rnd.NGM <- function (self, xpast, t, ...) { 
  rnorm(length(xpast),  mean=meanPrior(self, xpast, t), sd=self$sigmaW) 
}
initial.rnd.NGM <- function (self, N=1, ...) { 
  rep(0.1, N)
}
initial.logpdf.NGM <- function (self, x, ...) { 
  rep(0, length(x))
}
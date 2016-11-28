#' Simulation of a dataset for a NLSS
#' 
#' Generates a sequence of hidden states and of observations for a given NLSS.
#' This NLSS can be any class that supports the methods loclike.rnd, prior.rnd, and 
#' has attributes spaces with fields dimx and dimy.
#' 
#' @param nlss An instance of a NLSS, e.g. ARCH or NoisyAR
#' @param x1 Initial hidden state
#' @param T Number of time-steps including the initial one
#' @return 
#' \item{x}{Array (T, dimx) of hidden states}
#' \item{y}{Array (T, dimy) of observations}
#' \item{y}{Array (T, dimy) of observations}
#' @export
simulate.data <- function(nlss, x1, T) {
  x <- array(dim=c(T,attr(nlss,"spaces")$dimx))
  y <- array(dim=c(T,attr(nlss,"spaces")$dimy))
  x[1,] <- x1
  y[1,] <- loclike.rnd(nlss,x[1,])  
  if (T > 1) {
    for(t in seq(2,T,by=1)) {
      x[t,] <- prior.rnd(nlss,x[t-1,],t)
      y[t,] <- loclike.rnd(nlss,x[t,],t)
    }
  }
  t<-1:T
  list(x=x,y=y,t=t)
}

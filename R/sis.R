#' Sequential Importance Sampling with Prior kernel for 1-D NLSS
#'
#' Runs sequential importance sampling  \strong{without} resamplingon a given Non-Linear State Space models with prior kernel as proposal.
#'
#' This algorithm is not recommended and included for illustrative purposes only. This version
#'  is therefore  a minimalistic and only supports NLLS with univariate states. Use \code{\link{sisr}} instead.
#' 
#' The variant \code{\link{siskernel}} allows for arbitrary proposal kernel.
#'
#' @param nlss Non-linear state space model
#' @param y Sequence of observations. Its length T is the number of timesteps.
#' @param N Number of particles
#' @param resampling Resampling scheme to use
#' @return A list with the following components:
#'   \item{particles}{Array (T, N) of the sampled particles}
#'   \item{logweights}{Array (T, N) of the \strong{logarithm} of the \strong{non-normalized} importance weights of the particles}
#'   \item{weights}{Array (T, N) of the \strong{normalized} importance weights of the particles}
#'    \item{t}{Indices 1 to T, included for ease of plotting}
#'
#' @export
#' @seealso \code{\link{siskernel}} \code{\link{sisr}}    
#'    
sis <- function (nlss, y, N) {
  T <- length(y); 
  p <- array(dim=c(T,N)); 
  lw <- array(dim=c(T,N)); 
  w <- array(dim=c(T,N))

  
  particles <- initial.rnd(nlss, N)
  if (!is.na(y[1]))    
    logweights <- loclike.logpdf(nlss, particles, y=y[1], t=1)
  else    
    logweights <- array(0, dim=c(1,N))
  weights <- normalized.exponential(logweights)
  # Store the history
  p[1,] <- particles; lw[1,] <- logweights; w[1,] <- weights
  
  if (T > 1)
    for (t in 2:T) {
      particles <- prior.rnd(nlss, particles, t=t)
      if (!is.na(y[t]))
          logweights <- logweights + loclike.logpdf(nlss, particles, y[t], t=t)      
      weights <- normalized.exponential(logweights)
      p[t,] <- particles; lw[t,] <- logweights; w[t,] <- weights
    }
  
  # Return the weighted sample
  list(particles=p, logweights=lw, weights=w, t=1:T)  
}
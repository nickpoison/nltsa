
#' Sequential Importance Sampling with arbitrary kernel for 1-D NLSS
#'
#' Runs sequential importance sampling \strong{without} resampling on a given Non-Linear State Space models with user-specified kernel as proposal.
#'
#' This algorithm is a slightly more generic version of \code{\link{sis}}. It is not recommended and included for illustrative purposes only. This version
#'  is therefore  a minimalistic and only supports NLLS with univariate states. Use \code{\link{sisr}} instead.
#'
#'
#' @param nlss Non-linear state space model
#' @param y Sequence of observations. Its length T is the number of timesteps.
#' @param N Number of particles
#' @param proposal.rnd Function sampling from the proposal kernel to use
#' @param proposal.logpdf Function computing the log-pdf of the proposal kernel
#' @param initial.proposal.rnd Function sampling from the proposal kernel to use at initial timestep
#' @param initial.proposal.logpdf Function computing the log-pdf of the proposal kernel at initial timestep
#' @return A list with the following components:
#'   \item{particles}{Array (T, N, D) of the sampled particles}
#'   \item{logweights}{Array (T, N) of the \strong{logarithm} of the \strong{non-normalized} importance weights of the particles}
#'   \item{weights}{Array (T, N) of the \strong{normalized} importance weights of the particles}
#'    \item{t}{Indices 1 to T, included for ease of plotting}
#'    
#' @export
#' @seealso \code{\link{sisr}}
siskernel <- function (nlss, y, N,
  proposal.rnd=prior.rnd,
  proposal.logpdf=prior.logpdf,
  initial.proposal.rnd=initial.rnd,                       
  initial.proposal.logpdf=initial.logpdf) {
  
  T <- length(y); 
  p <- array(dim=c(T,N)); 
  lw <- array(dim=c(T,N)); 
  w <- array(dim=c(T,N))
  
  particles <- initial.proposal.rnd(nlss, N=N, y[1])
  if(!is.na(y[1]))
    logweights <- initial.logpdf(nlss, particles) - initial.proposal.logpdf(nlss, particles, y[1], t=1) + loclike.logpdf(nlss, particles, y[1]) 
  else
    logweights <- array(0,dim=c(1,N))
  weights <- normalized.exponential(logweights)
  # Store the history
  p[1,] <- particles; lw[1,] <- logweights; w[1,] <- weights
  
  if (T > 1)
    for (t in 2:T) {
      xpast <- particles
      particles <- proposal.rnd(nlss, xpast, y[t], t=t)
      if (!is.na(y[t]))      
        logweights <- logweights + loclike.logpdf(nlss, particles, y=y[t], t=t)
      logweights <- logweights + 
        prior.logpdf(nlss, xpast, particles, y=y[t], t=t) - 
        proposal.logpdf(nlss, xpast, particles, y=y[t], t=t)
      weights <- normalized.exponential(logweights)
      # Store the history
      p[t,] <- particles; lw[t,] <- logweights; w[t,] <- weights
    }
    
  # Return the weighted sample
  list(particles=p, logweights=lw, weights=w, t=1:T)
}
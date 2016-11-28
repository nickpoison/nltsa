#' Particle filter via Sequential Importance Sampling with Resampling
#'
#' Runs a particle filter on a given Non-Linear State Space models
#'
#' @param nlss Non-linear state space model
#' @param y Sequence of observations. Its length T is the number of timesteps.
#' @param N Number of particles
#' @param proposal.rnd Function sampling from the proposal kernel to use
#' @param proposal.logpdf Function computing the log-pdf of the proposal kernel
#' @param initial.proposal.rnd Function sampling from the proposal kernel to use at initial timestep
#' @param initial.proposal.logpdf Function computing the log-pdf of the proposal kernel at initial timestep
#' @param resampling Resampling scheme to use
#' @return A list with the following components:
#'   \item{particles}{Array (T, N, D) of the sampled particles}
#'   \item{logweights}{Array (T, N) of the \strong{logarithm} of the \strong{non-normalized} importance weights of the particles}
#'   \item{weights}{Array (T, N) of the \strong{normalized} importance weights of the particles}
#'   \item{ancestors}{Array (T, N) of the index of the ancestor of each particle in the previous generation}
#'    \item{t}{Indices 1 to T, included for ease of plotting}
#'    \item{success}{TRUE if the filtering succeeded}
sisr  <- function (nlss,y,N,
                    proposal.rnd=prior.rnd,
                    proposal.logpdf=prior.logpdf, 
                    resampling=ResidualMultinomialR,
                   initial.proposal.rnd=initial.rnd,                       
                   initial.proposal.logpdf=initial.logpdf
                   ) {
  
  T <- length(y)
  p <- array(dim=c(T,N,attr(nlss,"spaces")$dimx))
  lw <- array(dim=c(T,N)); w <- array(dim=c(T,N)) 
  a <- array(dim=c(T,N))
  particles <- array(dim=c(N,attr(nlss,"spaces")$dimx))
  
  ancestors <- 1:N
  particles[,] <- initial.proposal.rnd(nlss, N=N, y[1,])
  if(!any(is.na(y[1,])))
    logweights <- initial.logpdf(nlss, particles) - initial.proposal.logpdf(nlss, particles, y[1,], t=1) + loclike.logpdf(nlss, particles, y[1,]) 
  else
    logweights <- array(0,dim=c(1,N))
  weights <- normalized.exponential(logweights)
  p[1,,] <- particles; lw[1,] <- logweights; w[1,] <- weights; a[1,] <- ancestors
  
  if (T > 1)
    for (t in 2:T) {  
      # Resampling is optional and can be triggered by, e.g., a test for a low ESS
      ancestors <- resampling(weights, N)
      logweights <- array(0,dim=c(1,N))  
      
      xpast <- particles[ancestors,,drop=FALSE]
      particles[,] <- proposal.rnd(nlss, xpast, y=y[t,], t=t)      
      if(!any(is.na(y[t])))
        logweights <- loclike.logpdf(nlss, particles, y=y[t,], t=t)
      logweights <- logweights +
        prior.logpdf(nlss, xpast, particles, t=t) - 
        proposal.logpdf(nlss, xpast, particles, y=y[t,], t=t)
      weights <- normalized.exponential(logweights)
      p[t,,] <- particles; lw[t,] <- logweights; w[t,] <- weights; a[t,] <- ancestors
    }
  
  # Return the weighted sample
  list(particles=p, logweights=lw, weights=w, ancestors=a, t=1:T, success=TRUE)
}

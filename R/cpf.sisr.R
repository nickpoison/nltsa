#' Conditional Particle filter via Sequential Importance Sampling with Resampling
#'
#' Runs a conditional particle filter on a given Non-Linear State Space models
#'
#' Warning: the resampling is Multinomial. Residual will bias the sampler, due to the 
#' way we overried the Nth particle by the conditioned trajectory.
#'
#' @param x Trajectory on which to condition
#' @param ... same as classical sisr 
#' @return A list with the same components as SISR, plus the  components:
#'   \item{x}{Array (T, D) of the new selected trajectory}
cpf.sisr  <- function (x,nlss,y,N,
                    proposal.rnd=prior.rnd,
                    proposal.logpdf=prior.logpdf, 
                    resampling=MultinomialR,
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
  
  # Insert the conditioned particle at position N
  particles[N,] <- x[1,]
  
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
      
      # Insert the conditioned particle at position N
      particles[N,] <- x[t,]
      # Sample its ancestor
      inserted.ancestors.logweights <- lw[t-1,] + prior.logpdf(nlss, p[t-1,,], x[t,], t=t) 
      inserted.ancestors.weights <- normalized.exponential(inserted.ancestors.logweights)
      ancestors[N] <- resampling(inserted.ancestors.weights, 1)      
      xpast[N,] <- p[t-1,ancestors[N],]

      if(!any(is.na(y[t])))
        logweights <- loclike.logpdf(nlss, particles, y=y[t,], t=t)
      logweights <- logweights +
        prior.logpdf(nlss, xpast, particles, t=t) - 
        proposal.logpdf(nlss, xpast, particles, y=y[t,], t=t)
      weights <- normalized.exponential(logweights)
      p[t,,] <- particles; lw[t,] <- logweights; w[t,] <- weights; a[t,] <- ancestors
    }
  
  # Sample one of the trajectories according to the final weights
  i <- resampling(weights, 1)
  for(t in T:1) {
    x[t,] <- p[t,i,]
    i <- a[t,i]
  }
    
  # Return the weighted sample
  list(particles=p, logweights=lw, weights=w, ancestors=a, t=1:T, x=x, success=TRUE)
}

FFBSm <- function(nlss,sisr.output,pairwise=FALSE) {
  T <- dim(sisr.output$ancestors)[1]
  N <- dim(sisr.output$ancestors)[2]
  
  particles <- sisr.output$particles
  weights <- array(0,dim=dim(sisr.output$weights))
  # Pairwise weights only if required
  if (pairwise) {
    pairweights <- array(0,dim=c(T,N,N))    
    pairlogweights <- array(0,dim=c(T,N,N))
  }  else {
    pairweights <- NA    
    pairlogweights <- NA
  }

  weights[T,] <- sisr.output$weights[T,]
  for(k in (T-1):1) { 
    for (j in 1:N) {
      prior.pdf.j <- exp(prior.logpdf(nlss,particles[k,,],particles[k+1,j,],k+1))
      # Normalization factor of backward kernel (loop over N is hidden within sum())
      normalize.B.j <- sum(sisr.output$weights[k,]*prior.pdf.j)
      # Update of filtering weights at time k into smoothing weights (loop over N is vectorized)
      weights[k,] <- weights[k,] + weights[k+1,j]*prior.pdf.j/normalize.B.j      
    }
    weights[k,] <- sisr.output$weights[k,]*weights[k,]
  }
  
  logweights <- log(weights)
  results <- list(particles=particles, logweights=logweights, weights=weights, t=1:T)
  
  if (pairwise) {
    for (k in 2:T) {
      for (j in 1:N) {
            pairlogweights[k,,j] <- sisr.output$logweights[k-1,] +
              prior.logpdf(nlss, particles[k-1,,], particles[k,j,], k)
        pairweights[k, , j] <- normalized.exponential(pairlogweights[k,,j]) * weights[k,j] 
        pairlogweights[k, , j] <- log(pairweights[k, , j])
      }
    }
    results <- c(results, list(pairweights=pairweights, pairlogweights=pairlogweights))
  }
  
  results
}

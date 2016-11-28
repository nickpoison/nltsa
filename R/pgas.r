PGAS <- function(nlss, y,
                          n.mcmc, n.particles,
                          .progress=progress_text, ...) {
  
  T <- length(y)
  d.x <- attr(nlss, "spaces")$dimx
  
  # Progress bar, EMMCMC run can be long
  progress.bar <- .progress()
  progress.bar$init(n.mcmc)
  trace <- array(0, dim=c(T, n.mcmc, d.x))
  
  # Initial trajectory for the CPF is taken from a plain SISR
  p <- sisr(nlss,y=y,N=n.particles,...)
  i <- MultinomialR(p$weights[T,], 1)
  x <- array(0, dim=c(T,d.x))
  for(t in T:1) {
    x[t,] <- p$particles[T,i,]
    i <- p$ancestors[t,i]
  }  
  progress.bar$step()
  
  # Then just roll on CPF
  for (k in 1:n.mcmc) {
    x <- cpf.sisr(x, nlss, y, N=n.particles, ...)$x
    trace[,k,] <- x
    progress.bar$step()
  }

 list(particles=trace, logweights=matrix(0, ncol=n.mcmc), weights=matrix(1, ncol=n.mcmc), t=1:T)                 
}
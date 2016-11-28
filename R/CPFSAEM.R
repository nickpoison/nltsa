#' Run a Conditional Particle Stochastic Approximation EM chain
#' 
#' This function implements Particle Stochastic Approximation EM
#' 
#' @param initial.nlss NLSS whose parameters value will serve as starting point
#' @param y Observations
#' @param n.particles Number of particles
#' @param n.em Number of SAEM iterations
#' @param sa.rate Learning rate as a function of the iteration
#' @param .progress Progress bar to use, from plyr package, default to progress_text
#' @param ... Extra parameters to pass to the SISR filter, see function \code{\link{sisr}}
#' 
#' @return A list with the same components as \code{\link{MH}}
#' @seealso \code{\link{MH}}, \code{\link{sisr}}, \code{\link{random.walk}}
#' 
#' @import plyr
#' @export
CPFSAEM <- function(initial.nlss, 
                              y, 
                              n.particles=function(i) { 30 },                         
                              sa.rate=function(i) { 1/(i+1)^(2/3) },
                              n.em,
                              .progress=progress_text, ...) {
  
  skeleton <- as.relistable(initial.nlss)    
  nlss <- initial.nlss  
  T <- length(y)
  d.alpha <- attr(nlss,"spaces")$dimalpha
  d.x <- attr(nlss, "spaces")$dimx
  
  # Progress bar, EM run can be long
  progress.bar <- .progress()
  progress.bar$init(n.em)
  trace <- array(0, dim=c(n.em+1, length(nlss)))
  colnames(trace) <- names(nlss)
  trace[1,] <- unlist(nlss)
  
  if (!is.function(n.particles)) {
    fixed.n <- n.particles
    n.particles <- function(it) { fixed.n }
  }
  
  # Sufficient statistics
  EXX <- 0
  EXF <- array(0,dim=c(d.alpha,1))
  EFF <- array(0, dim=c(d.alpha,d.alpha))  
  
  sigmaV2 <- 0
  sigmaW2 <- 0
  
  # Initial trajectory for the CPF is taken from a plain SISR
  p <- sisr(nlss,y=y,N=n.particles(it.em),...)
  i <- MultinomialR(p$weights[T,], 1)
  x <- array(0, dim=c(T,d.x))
  for(t in T:1) {
    x[t,] <- p$particles[T,i,]
    i <- p$ancestors[t,i]
  }
  
  for (it.em in 1:n.em) {
    N <- n.particles(it.em)
  #  print(c('Iteration', it.em))
    smoothed<-cpf.sisr(x,nlss,y=y,N=n.particles(it.em),...)    
    x <- smoothed$x
    
    # Untwist the CPF
    trajectories <- array(0, dim=c(T,N,d.x))
    a <- 1:N
    for(t in T:1) {
      trajectories[t,,] <- smoothed$particles[t,a,]
      a <- smoothed$a[t,a]
    }
        
    EXX.new <- 0
    EXF.new <- array(0,dim=c(d.alpha,d.x))
    EFF.new <- array(0, dim=c(d.alpha,d.alpha))      
    for (k in 2:T) {
      EXX.new <- EXX.new + sum(smoothed$weights[T,] * (trajectories[k,,])^2)
      f <- prior.f(nlss, trajectories[k-1,,], k) 
            
      wf <- (array(1, dim=c(d.alpha,1)) %*% smoothed$weights[T,]) * f # this is d.alpha-by-N

      tk <-  trajectories[k,,,drop=FALSE]
      dim(tk) <- c(N, d.x)
      EXF.new <- EXF.new + wf %*% tk      
      EFF.new <- EFF.new + wf %*% t(f)
    }

    if (it.em == 1) {
      lr <- 1
    } else {
      lr <- sa.rate(it.em)
    }
    EXX <- (1-lr)*EXX + lr*EXX.new
    EXF <- (1-lr)*EXF + lr*EXF.new
    EFF <- (1-lr)*EFF + lr*EFF.new
        
    
    alpha <- solve(EFF, EXF)
    sigmaW <- sqrt((EXX - t(EXF) %*% alpha)/T)
    
    sigmaV2.new <- 0
    for (k in 1:T)  {
      sigmaV2.new <- sigmaV2.new + 
        sum(smoothed$weights[T,] * (y[k] - loclike.h(nlss, trajectories[k,,]))^2)
    }
    sigmaV2.new <- sigmaV2.new / T
    sigmaV2 <- (1-lr)*sigmaV2 + lr*sigmaV2.new    
    sigmaV <- sqrt(sigmaV2)
    
    nlss <- relist(flesh=c(alpha, sigmaW, sigmaV), skeleton=skeleton)  
    trace[it.em+1,] <- unlist(nlss)
      
    progress.bar$step()
  }
  list(nlss=nlss, trace=trace)
}
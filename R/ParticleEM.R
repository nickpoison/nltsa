#' Run a Particle EM chain
#' 
#' This function implements Particle EM
#' 
#' @param initial.nlss NLSS whose parameters value will serve as starting point
#' @param y Observations
#' @param n.particles Number of particles
#' @param n.em Number of EM iterations
#' @param .progress Progress bar to use, from plyr package, default to progress_text
#' @param ... Extra parameters to pass to the SISR filter, see function \code{\link{sisr}}
#' 
#' @return A list with the same components as \code{\link{MH}}
#' @seealso \code{\link{MH}}, \code{\link{sisr}}, \code{\link{random.walk}}
#' 
#' @import plyr
#' @export
ParticleEM <- function(initial.nlss, 
                              y, 
                              n.particles=function(i) { 100 + i^2 }, 
                              n.em,
                              .progress=progress_text, ...) {
  
  skeleton <- as.relistable(initial.nlss)    
  nlss <- initial.nlss  
  T <- length(y)
  d.alpha <- attr(nlss,"spaces")$dimalpha
    
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
  
  for (it.em in 1:n.em) {
     
    p<-sisr(nlss,y=y,N=n.particles(it.em),...)
    smoothed <- FFBSm(nlss, p, pairwise=TRUE)
    
    EXX <- 0
    EXF <- array(0,dim=c(d.alpha,1))
    EFF <- array(0, dim=c(d.alpha,d.alpha))
    for (k in 2:T) {
      EXX <- EXX + sum(smoothed$weights[k,] * (smoothed$particles[k,,])^2)
      f <- prior.f(nlss, smoothed$particles[k-1,,], k) 
      for (d in 1:d.alpha) {
        EXF[d] <- EXF[d] + 
          sum(smoothed$pairweights[k,,] * (t(f[d,,drop=FALSE]) %*% smoothed$particles[k,,]))
      }
      wf <- (array(1, dim=c(d.alpha,1)) %*% smoothed$weights[k-1,]) * f
      EFF <- EFF + wf %*% t(f)
    }
    
    alpha.new <- solve(EFF, EXF)
    sigmaW.new <- sqrt((EXX - t(EXF) %*% alpha.new)/T)
    
    sigmaV2.new <- 0
    for (k in 1:T)  {
      sigmaV2.new <- sigmaV2.new + 
        sum(smoothed$weights[k,] * (y[k] - loclike.h(nlss, smoothed$particles[k,,]))^2)
    }
    sigmaV2.new <- sigmaV2.new / T
    sigmaV.new <- sqrt(sigmaV2.new)
    
    
    nlss <- relist(flesh=c(alpha.new, sigmaW.new, sigmaV.new), skeleton=skeleton)  
    trace[it.em+1,] <- unlist(nlss)
      
    progress.bar$step()
  }
  list(nlss=nlss, trace=trace)
}
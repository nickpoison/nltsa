FFBSi <- function(nlss,
                  sisr.output,
                  resampling=ResidualMultinomialR,
                  n.trajectories=dim(sisr.output$particles)[2]) {
  T <- dim(sisr.output$particles)[1]
  N <- dim(sisr.output$particles)[2]
  d <- dim(sisr.output$particles)[3]  
  trajectories <- array(NA_real_, dim=c(T, n.trajectories, d))
  
  indices <- resampling(sisr.output$weights[T, ], n.trajectories)
  trajectories[T, , ] <- sisr.output$particles[T, indices, ]
  for(k in (T-1):1) { 
    for (j in 1:n.trajectories) {
      back.logweight.j <- sisr.output$logweights[k, ]+
        prior.logpdf(nlss, sisr.output$particles[k, , ], trajectories[k+1, j, ], k+1)
      son.j <- resampling(exp(back.logweight.j), 1)
      trajectories[k, j, ] <- sisr.output$particles[k, son.j, ]
    }
  }  
  list(particles=trajectories, logweights=matrix(0, ncol=n.trajectories), weights=matrix(1, ncol=n.trajectories), t=1:T)
}

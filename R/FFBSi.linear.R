FFBSi.linear <- function(nlss,
                  sisr.output,
                  sigma.plus,
                  resampling=ResidualMultinomialR,
                  n.trajectories=dim(sisr.output$particles)[2]) {
  T <- dim(sisr.output$particles)[1]
  N <- dim(sisr.output$particles)[2]  
  d <- dim(sisr.output$particles)[3]    
  trajectories <- array(NA_real_, dim=c(T, n.trajectories, d))
  J <- matrix(NA_real_, nrow=T, ncol=n.trajectories)
  
  J[T, ] <- resampling(sisr.output$weights[T, ], n.trajectories)
  trajectories[T, , ] <- sisr.output$particles[T, J[T, ], ]
  
  for(s in (T-1):1) { 
     L <- 1:n.trajectories
     while (length(L) > 0) {
       K <- length(L)
       I <- resampling(sisr.output$weights[s, ], K)
       U <- runif(K)
       nL <- NULL
       for (k in 1:K) {
         from <- sisr.output$particles[s, I[k], ]
         to <- sisr.output$particles[s+1, J[s+1, L[k]], ]
         if (U[k] < exp(prior.logpdf(nlss, from, to, s))/sigma.plus)
           J[s, L[k]] <- I[k]
         else
           nL <- c(nL, L[k])           
       }
       L <- nL
       trajectories[s, , ] <- sisr.output$particles[s, J[s, ], ]
    }     
  }
  list(particles=trajectories, logweights=matrix(0, ncol=n.trajectories), weights=matrix(1, ncol=n.trajectories), t=1:T)
}

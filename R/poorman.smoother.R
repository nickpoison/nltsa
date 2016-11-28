poorman.smoother <- function(sisr.output) {
  # Flatten the genealogic tree of particles and the weights according to the ancestors
  T <- dim(sisr.output$ancestors)[1]
  N <- dim(sisr.output$ancestors)[2]
  weights <- sisr.output$weights[T,]
  logweights <- sisr.output$logweights[T,]
  
  trajectories <- sisr.output$particles
  ancestors <- 1:N
  for(k in T:1) { 
    trajectories[k,,] <- sisr.output$particles[k,ancestors,]
    ancestors <- sisr.output$ancestors[k,ancestors]
  }
  list(particles=trajectories, logweights=logweights, weights=weights, t=1:T)
}
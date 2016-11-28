ResidualMultinomialR <- function (p, n=length(p)) {  
  if(sum(p) == 0) stop("too few positive probabilities")
  p <- p/sum(p)
  determinist <- floor(n*p)  
  residual.weights <- n*p - determinist
  nleft <- n-sum(determinist)
  if (nleft>0) 
    residual <- rmultinom(1,nleft,residual.weights)
  else 
    residual <- rep(0,length(p))
  rep(1:length(p),times=determinist+residual)
}
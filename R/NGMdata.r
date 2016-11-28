NGMdata <-
function(n.obs=100, alpha=.5, beta=25, gamma=8, sig2w=10, sig2v=1, seed=90210)
{
set.seed(seed)  
g = function(arg){c(arg[1],arg[1]/(1+arg[1]^2),cos(1.2*arg[2]))}
n         = n.obs
sig2      = sig2v
tau2      = sig2w
sig       = sqrt(sig2)
tau       = sqrt(tau2)
theta     = c(alpha,beta,gamma)
ptr       = c(theta,tau2,sig2)
y         = rep(0,n)
x         = rep(0,n)
x0        = 0
Z         = g(c(x0,0))
Fx        = rep(0,n)
Fx[1]     = sum(Z*theta)
x[1]      = rnorm(1,Fx[1],tau)
y[1]      = rnorm(1,x[1]^2/20,sig)
 for (t in 2:n){
   Z    = g(c(x[t-1],t-1))
   Fx[t] = sum(Z*theta)
   x[t] = rnorm(1,Fx[t],tau)
   y[t] = rnorm(1,x[t]^2/20,sig)
 }
 list(x=x, y=y)
}
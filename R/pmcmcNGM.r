##############################################################
pmcmcNGM <-
function(nmcmc, burnin, step, y, alpha00, beta00, gamma00, sig2w00, sig2v00, npart, mcmseed){
###############################################################
###############################################################
# nmcmc = number of iterations for mcmc
# burnin = # of iterations used for burnin
# step =  # of iterations to skip  
# y = data 
# sig2w00 = hyperpameters for sigma_w^2
# sig2v00 =  ditto for sigma_v^2
# npart =  number of particles
##################################################################

niter = nmcmc*step + burnin         #total number of iterations

pr <- progress_text()               # displays progress (from plyr)
pr$init(niter) 

##### setup #########
set.seed(mcmseed)
nobs    = length(y)
alpha   = rep(0, niter)
beta    = rep(0, niter)
gamma   = rep(0, niter)
sig2w   = rep(0,niter)
sig2v   = rep(0,niter)
X       = matrix(0,niter, nobs)
loglik  = rep(0,niter)
accrate = rep(0, niter) # accept rate
#########

# Initialize the parameters
alpha[1] =  alpha00[1]                    # initial value of regr parameters
beta[1]  =  beta00[1]
gamma[1] =  gamma00[1]
sig2w[1] =  sig2w00[2]/(sig2w00[1]-1)      # initial value of sig2w (mean of IG prior)
sig2v[1] =  sig2v00[2]/(sig2v00[1]-1)      # initial value of sig2v  ( ditto )

# tuning parameters
tuninga = alpha00[2]/5              # for gauss rw (each are .2 of their prior stdevs)
tuningb = beta00[2]/5
tuningg = gamma00[2]/5
tuningw = sqrt(sig2w[1]^2/(sig2w00[1]-2))/5     # for sig2w                
tuningv = sqrt(sig2v[1]^2/(sig2v00[1]-2))/5     # for sig2v
#############

   

# Run pfilter to initialize state and likelihood 
u = pfilter(y, alpha[1], beta[1], gamma[1], sig2w[1], sig2v[1], npart)
 particles = u$particles   # returned particles
 w = u$w                   # returned weights
 loglik.p = u$loglik.p     # returned log-likelihood
# Draw J
 J = which( (runif(1)-cumsum(w[,nobs])) < 0 )[1]
 X[1,] = particles[J,]

 # Run MCMC loop
 for(k in 2:niter){
   pr$step()  
   alpha.p = alpha[k-1] + tuninga*rnorm(1)
   beta.p  = beta[k-1] + tuningb*rnorm(1)
   gamma.p = gamma[k-1] + tuningg*rnorm(1)
   sig2w.p = abs(sig2w[k-1] + tuningw*rnorm(1))    # abs to insure >= 0 - should rarely happen
   sig2v.p = abs(sig2v[k-1] + tuningv*rnorm(1))
   # Run a pfilter to evaluate the likelihood
    u = pfilter(y, alpha.p, beta.p, gamma.p, sig2w.p, sig2v.p, npart)
    particles = u$particles
    w = u$w
    loglik.p = u$loglik.p         
					  
    acceptprob = exp(loglik.p - loglik[k-1])   # Likelihood contribution
	numer = dnorm(alpha.p, alpha00[1], alpha00[1])*
	        dnorm(beta.p, beta00[1], beta00[2])*
			dnorm(gamma.p, gamma00[1], gamma00[2])*
	        ipgamma(sig2w.p, sig2w00[1], sig2w00[2])*
			ipgamma(sig2v.p, sig2v00[1], sig2v00[2]) 
	denom = dnorm(alpha[k-1], alpha00[1], alpha00[1])*
	        dnorm(beta[k-1], beta00[1], beta00[2])*
			dnorm(gamma[k-1], gamma00[1], gamma00[2])*
	        ipgamma(sig2w[k-1], sig2w00[1], sig2w00[2])*
			ipgamma(sig2v[k-1], sig2v00[1], sig2v00[2])
	
    acceptprob = acceptprob * numer/denom                 # Prior contribution
    accept = runif(1) < acceptprob
	         
    if(accept){  
	    accrate[k]=1    # count successes 
        alpha[k] = alpha.p
        beta[k]  = beta.p
        gamma[k] = gamma.p		
        sig2w[k] = sig2w.p
        sig2v[k] = sig2v.p
        # Draw J (extract a particle trajectory)
		J = which( (runif(1)-cumsum(w)) < 0 )[1]          # w's are the weights
        X[k,] = particles[J,]
        loglik[k] = loglik.p
     }else{
	    alpha[k] = alpha[k-1]
        beta[k]  = beta[k-1]
        gamma[k] = gamma[k-1]
        sig2w[k] = sig2w[k-1]
        sig2v[k] = sig2v[k-1]
        X[k,] = X[k-1,]
        loglik[k] = loglik[k-1]        
    }
  }
ind = seq(burnin+1, niter, by=step)  
list(alpha=alpha[ind], beta=beta[ind], gamma=gamma[ind], sig2w=sig2w[ind], sig2v=sig2v[ind], accrate=sum(accrate[ind])/nmcmc, X=X[ind,])
}
#########################################################################################



#######################################################################################
pfilter = function(y, alpha, beta, gamma, sig2w, sig2v, npart){   #####################
#######################################################################################
# particle filter
# y = observations
# sig2w = state noise variance
# sig2v = obs noise variance
# npart = number of particles
#########################################################
#########################################################

nobs = length(y);
x = matrix(0, npart, nobs) # particles
a = matrix(0, npart, nobs) # ancestor indices
w = matrix(0, npart, nobs) # weights
x[,1] = 0;                 # initial condition
loglik = 0

 for(t in 1:nobs){
    if(t != 1){
        ind = resamplew(w[,t-1])           
		xt1 = x[,t-1]
        xpred = alpha*xt1 + beta*xt1/(1+xt1^2) + gamma*cos(1.2*(t-1))  
        x[,t] = xpred + sqrt(sig2w)*rnorm(npart)        
        a[,t] = ind        # ancestor indices
    }
  # importance weights
    ypred = x[,t]^2/20
    logweights = -1/2*log(2*pi*sig2v) - (1/(2*sig2v))*(y[t] - ypred)^2
    const = max(logweights)    # subtract the maximum value for numerical stability
    weights = exp(logweights-const)
    # Compute loglikelihood
    loglik = loglik + const + log(sum(weights))        
    w[,t] = weights/sum(weights)    # normalized weights
 }

 # Generate trajectories from ancestor indices
  ind = a[,nobs]
  for(t in (nobs-1):1){
    x[,t] = x[ind,t]
    ind = a[ind,t]
   }
list(particles=x, w=w, loglik.p=loglik)
}


######################################################
ipgamma = function(x,a,b){
######################################################
# inverse gamma pdf at x with parameters shape = a and scale = b 
px = exp(a*log(b) - log(gamma(a)) - (a+1)*log(x) - b/x)
return(px)
}
#######################################################
#######################################################

#########################################################
resamplew = function(w){
##########################################################
# multinomial resampling
N = length(w)
u = runif(N)
cw = cumsum(w)
cw = cw/cw[N]
 ucw = c(u,cw)
 ind1 = sort(ucw, index.return=TRUE)$ix
 ind2 = which(ind1<=N)
i = ind2-(0:(N-1))
return(i)
}
#############################################################




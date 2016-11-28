#function [q,r,X] = pgas(numMCMC, y, prior, N, qinit, rinit, q0, r0, plotOn)

pgasNGM = function(nmcmc, burnin, y, alpha00, beta00, gamma00, sig2w00, sig2v00, npart, mcmseed){

if(sig2w00[1] < 2 || sig2v00[1] < 2 ) stop("priors on variances must have at least a second moment")

numMCMC = nmcmc+burnin
N = npart
set.seed(mcmseed)

pr <- progress_text()               # displays progress (from plyr)
pr$init(numMCMC) 

T = length(y)
q = rep(0,numMCMC)             # q = sig2w
r = rep(0,numMCMC)             # r = sig2v
X = matrix(0,numMCMC,T)
theta = matrix(0,numMCMC, 3)   #  cols are alpha, beta, gamma

# Initialize the parameters
priorq = sig2w00
priorr = sig2v00
q[1] = sig2w00[2]/(sig2w00[1]-1)   # mean of igamma(a,b) = b/(a-1) ... require a>2 for a variance 
r[1] = sig2v00[2]/(sig2v00[1]-1)
theta[1,] = c(alpha00[1], beta00[1], gamma00[1])

# Initialize the state by running a PF
#[particles, w] = cpf_as(y, q(1), r(1), N, X);
u = cpf_as(y, alpha00[1], beta00[1], gamma00[1], q[1], r[1], N, X);
   particles = u$x           # returned particles
   w = u$w                   # returned weights
# Draw J
 J = which( (runif(1)-cumsum(w[,T])) < 0 )[1]
 X[1,] = particles[J,]


# Run MCMC loop
 for(k in 2:numMCMC){
    pr$step()                                                     
    # Sample the parameters (inverse gamma posteriors)
	   t1 = 1:(T-1)
	   alpha = theta[k-1,1]; beta = theta[k-1,2]; gamma = theta[k-1,3]
	   xpred = alpha*X[k-1, t1] + beta*X[k-1, t1]/(1+X[k-1, t1]^2) + gamma*cos(1.2*(t1))
      err_q = X[k-1,2:T] - xpred                   
   # q(k) = igamrnd(prior[1] + (T-1)/2, prior[2] + sum(err_q)^2/2);
     q[k] = 1/rgamma(1, priorq[1] + (T-1)/2, priorq[2] + sum(err_q^2)/2)
      err_r = y - X[k-1,]^2/20 
    #r(k) = igamrnd(prior.a + T/2, prior.b + err_r*err_r'/2);
     r[k] =  1/rgamma(1, priorr[1] + T/2, priorr[2] + sum(err_r^2)/2)
	# update theta
	  Z = cbind(X[k-1,2:T], X[k-1, t1], (X[k-1, t1]/(1+X[k-1, t1]^2)), cos(1.2*(t1)))
	    B = 1/(1/alpha00[2] + (1/q[k])*sum(Z[,2]^2))
        b = alpha00[1]/alpha00[2] + (1/q[k])*t(Z[,2])%*%(Z[,1]-beta*Z[,3]-gamma*Z[,4])	   
     theta[k,1] = rnorm(1, B*b, B)
	    B = 1/(1/beta00[2] + (1/q[k])*sum(Z[,3]^2))
	    b = beta00[1]/beta00[2] + (1/q[k])*t(Z[,3])%*%(Z[,1]-alpha*Z[,2]-gamma*Z[,4])	
	 theta[k,2] = rnorm(1, B*b, B)
	    B = 1/(1/beta00[2] + (1/q[k])*sum(Z[,4]^2))
	    b = gamma00[1]/gamma00[2] + (1/q[k])*t(Z[,4])%*%(Z[,1]-alpha*Z[,2]-beta*Z[,3])
     theta[k,3] = rnorm(1, B*b, B)	
	 
	# Run CPF-AS
    # [particles, w] = cpf_as(y, q(k), r(k), N, X(k-1,:));
      u = cpf_as(y, theta[k,1], theta[k,2], theta[k,3], q[k], r[k], N, X[k-1,])
	  particles = u$x   # returned particles
      w = u$w                   # returned weight
	# Draw J (extract a particle trajectory)
	 # Draw J
        J = which( (runif(1)-cumsum(w[,T])) < 0 )[1]
        X[k,] = particles[J,]    
 }#end
bi = 1:burnin 
list(theta=theta[-bi,], sig2w=q[-bi], sig2v=r[-bi], X=X[-bi,])
}#end

###############################################################################
#--------------------------------------------------------------------------
# function [x,w] = cpf_as(y, q, r, N, X)
  cpf_as = function(y, alpha, beta, gamma, q, r, N, X){
# Conditional particle filter with ancestor sampling
# Input:
#   y - measurements
#   q - process noise variance
#   r - measurement noise variance
#   N - number of particles
#   X - conditioned particles - if not provided... 
unconditional = is.null(X)   # ... an unconditional PF is run

T = length(y);
x = matrix(0,N, T); # Particles
a = matrix(0,N, T); # Ancestor indices
w = matrix(0, N, T); # Weights
x[,1] = 0; # Deterministic initial condition
x[N,1] = X[1] 

for (t in 1:T){
    if(t != 1){
        ind = resamplew(w[,t-1]);             
        ind = ind[sample.int(N)];
           t1 = t-1
	       xpred = alpha*x[, t1] + beta*x[, t1]/(1+x[, t1]^2) + gamma*cos(1.2*(t1))
        x[,t] = xpred[ind] + sqrt(q)*rnorm(N);
	
            x[N,t] = X[t]; 
            # Ancestor sampling                                              
            m = exp(-1/(2*q)*(X[t]-xpred)^2);
            w_as = w[,t-1]*m
            w_as = w_as/sum(w_as);
            #ind(N) = find(runif(1) < cumsum(w_as),1,'first');
			ind[N] =  which( (runif(1)-cumsum(w_as)) < 0 )[1]
        
		# Store the ancestor indices
        a[,t] = ind;
    }#end
    # Compute importance weights
    ypred = x[,t]^2/20
    logweights = -1/(2*r)*(y[t] - ypred)^2; # (up to an additive constant)
    const = max(logweights); # Subtract the maximum value for numerical stability
    weights = exp(logweights-const);
    w[,t] = weights/sum(weights); # Save the normalized weights
}#end

# Generate the trajectories from ancestor indices
ind = a[,T];
for(t in (T-1):1){
    x[,t] = x[ind,t];
    ind = a[ind,t];
}#end
list(x=x, w=w)
}#end
#-------------------------------------------------------------------
#-------------------------------------------------------------------

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

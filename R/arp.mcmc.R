arp.mcmc <- 
function(xdata, porder, n.iter=1000, n.warmup=100){
#
if (porder > 9) {stop("WARNING: The AR order must be less than 10")}
#
yy=xdata
nobs=length(yy)
lagp=porder
lagp1=lagp+1
nwarmup = n.warmup
niter = n.iter + nwarmup
y=yy[lagp1:nobs];
x=matrix(1,nobs-lagp,lagp1)
for (j in lagp1:2) {x[,j]=yy[(lagp-j+2):(nobs-j+1)]}
tt=(lagp1:nobs)/(nobs-lagp)
phi=matrix(0, lagp1, niter+1)
sigma=rep(1,niter+1)
sigphi=50
fit= matrix(0,nobs-lagp,niter)
pred_den_hat=matrix(0,nobs-lagp,1);
for (p in 1:niter){
   # Drawing the phis
    var_phi=solve((1/sigma[p])*(t(x)%*%x+(1/sigphi)*diag(1,lagp1)))
   # mean 
    mu_phi=(1/sigma[p])*var_phi%*%t(x)%*%y
	  Z = as.matrix(rnorm(lagp1,0,1))
      eV=eigen(var_phi, symmetric=TRUE)
     phi[,p+1]= mu_phi + eV$vectors%*%diag(sqrt(pmax(eV$values,0)),lagp1)%*%Z
   ### phi[,p+1]=mvrnorm(mu=mu_phi,Sigma=var_phi) # replace need for mvrnorm above 3 lines
   # Drawing sigma
    siga=(nobs-lagp)/2+0.01
    sigb=1/(sum((y-x%*%phi[,p+1])^2)/2+.01+ (1/(2*sigphi))*(t(phi[,p+1])%*%phi[,p+1]))
    sigma[p+1]=1/rgamma(1,shape=siga,scale=sigb)
   } 
# plot results	
indx = (nwarmup+1):niter
phit=t(phi[-1,indx])
sigma=sigma[indx]
u = ts.union(phi=ts(phit), tauinv=ts(sigma))
#u = window(u, start=(nwarmup+1), end=niter)
plot(u, main="", xlab="Iteration")
print(summary(u))
u = list(phi=phit, tauinv=sigma)
return(invisible(u))
}





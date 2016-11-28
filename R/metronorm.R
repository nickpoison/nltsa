metronorm <-
function(CovProp,CovPi,nbiter,Xinit){
#
accept= 0
CovProp= as.matrix(CovProp)
CovPi= as.matrix(CovPi)
Xinit= as.matrix(Xinit)
d= ncol(CovPi)
X= matrix(0,d,nbiter)
X[,1]= Xinit
CholFactorCovProp= t(chol(CovProp))
CholFactorInvCovPi=t(chol(solve(CovPi)))
for (iter in 2:nbiter){
    Xcur= X[,iter-1]
    Y= Xcur+CholFactorCovProp%*%as.matrix(rnorm(d),d)
    alpha= 0.5*(-(vecnorm(CholFactorInvCovPi%*%Y))^2+(vecnorm(CholFactorInvCovPi%*%Xcur))^2)
    if (alpha < -10) {alpha=0;  X[,iter]= Xcur}
     else if (alpha > 0) {accept= accept+1; X[,iter]= Y}
     else {alpha=exp(alpha); U= (runif(1) < alpha)
           if (U == TRUE){accept= accept+1; X[,iter]= Y}
           else{X[,iter]= Xcur}} 
}       
accept= accept/nbiter
X = ts(t(X))
x.out = list(X=X)
cat("Acceptance Rate = ", accept, "\n")
return(invisible(x.out))
}

vecnorm = function(x) {
return(sqrt(sum(x^2)))
}

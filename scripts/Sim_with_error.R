sim.with.error=function(J,d=44,betasd=1,esd=0.11,n=400,rho=0.8,covmat){
  n=n
  
  covmat=lapply(seq(1:length(covmat)),function(x){covmat[[x]]/max(diag(covmat[[x]]))})
  
  
  
  library("mvtnorm")
  library("MASS")
  K=length(covmat)
  
  
  if(n!=0){
    z = sample(K,n,replace=TRUE)
    omega=abs(rnorm(n,mean=0,sd=betasd))##effect size variance can be big or small
    beta=t(sapply(seq(1:n),function(j){
      k=z[j]
      o=omega[j]
      mvrnorm(1,mu=rep(0,d),Sigma=o*covmat[[k]])
    }))
    beta=rbind(beta,matrix(rep(0,(J-n)*d),ncol=d))}
  if(n==0){
    beta=matrix(rep(0,(J-n)*d),ncol=d)
  }
  
  s.j.r=abs(rnorm(d,esd,0.001))##simulate with the same standard error for every J
  cormat=rep(1,d)%*%t(rep(1,d))*rho+(1-rho)*diag(1,d)#make the errors correlated
  v.mat=diag(s.j.r)%*%cormat%*%diag(s.j.r)##now v.j.r will be the same for every J
  cormat=cov2cor(v.mat)
  e=rmvnorm(J,mean=rep(0,d),sigma=v.mat)
  betahat = beta + e
  s.j=matrix(rep(as.matrix(s.j.r)),nrow(betahat),byrow=T,ncol=d)
  t.stat=betahat/abs(s.j)
  if(n!=0){
    return(list(beta=beta,betahat=betahat,component.mats=covmat,sebetahat=s.j,t.stat=t.stat,component.id=z,error=e,var.mat=v.mat,omega=omega,cormat=cormat))
  }
  if(n==0){
    return(list(beta=beta,betahat=betahat,sebetahat=s.j,t.stat=t.stat,error=e,var.mat=v.mat))
  }
}


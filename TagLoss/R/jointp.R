jointp <-
function(par)
{
  beta1=par[1]
  beta2=par[2]
  alpha=par[3]
  K=1+exp(beta1)+exp(beta2)+exp(beta1+beta2+alpha)
  jointps=matrix(c(1,exp(beta2),exp(beta1),exp(beta1+beta2+alpha))/K,nrow=2,byrow=T)
  rmarginal=c(exp(beta1)*(1+exp(beta2+alpha))/K)
  rmarginal=c(1-rmarginal,rmarginal)
  cmarginal=c(exp(beta2)*(1+exp(beta1+alpha))/K)
  cmarginal=c(1-cmarginal,cmarginal)
  conditional.1=matrix(c(1-plogis(beta1),1-plogis(beta1+alpha),plogis(beta1),plogis(beta1+alpha)),nrow=2,byrow=T)
  conditional.2=matrix(c(1-plogis(beta2),plogis(beta2),1-plogis(beta2+alpha),plogis(beta2+alpha)),nrow=2,byrow=T)  
  joint2= conditional.1*matrix(cmarginal,nrow=2,ncol=2,byrow=T)
  joint3= conditional.2*matrix(rmarginal,nrow=2,ncol=2,byrow=F)
  delta=1+(exp(alpha)-1)*plogis(beta1)*plogis(beta2)
 return(list(K=K,jointps=jointps,rmarginal=rmarginal,cmarginal=cmarginal,
       conditional.1=conditional.1,conditional.2=conditional.2,joint2=joint2,joint3=joint3,delta=delta))
}


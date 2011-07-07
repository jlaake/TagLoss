indep.loss <-
function(p,n01,n10,n11)
{
  n1=n01+n10
  n2=n11
  if(!is.matrix(p))
    return(n1*log(2*p*(1-p))+n2*log((1-p)^2)-(n1+n2)*log(1-p^2))
  else
    return(n01*log(p[,1]*(1-p[,2]))+n10*log(p[,2]*(1-p[,1]))+n11*log((1-p[,1])*(1-p[,2]))-(n1+n2)*log(1-p[,1]*p[,2]))
}


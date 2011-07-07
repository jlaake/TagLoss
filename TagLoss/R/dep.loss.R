dep.loss <-
function(par,n,nm)
{
  ps=as.vector(jointp(par)$jointps)
  lnl=sum(n*log(ps[1:3]))-sum(n)*log(sum(ps[1:3]))+sum(nm*log(ps))
  return(-lnl)
}


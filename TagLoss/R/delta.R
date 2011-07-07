delta <-
function(age,par)
{
      p1=p1(age,par)
      upper=1/p1
      lower=max(0,(2*p1-1)/p1^2)
      offset.v=log((1-lower)/(upper-1))
      delta=delta0(i,par,offset.v)*(upper-lower)+lower
      return(delta)
}


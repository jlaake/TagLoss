lnl.fi <-
function(xx,par)
{
    lnlike=0
    for (i in 0:15)
    {
      p1=p1(i,par)
      delta.v=1
      if(sum(xx[2:4,i+1])>0)
      lnlike=lnlike+(xx[2,i+1]+xx[3,i+1])*log(p1-p1^2*delta.v) +
          xx[4,i+1]*log(p1^2*delta.v) -
          (xx[2,i+1]+xx[3,i+1]+xx[4,i+1])*log(2*p1-p1^2*delta.v)
    }
    cat("\n par =",par," -lnl=",-lnlike)
    return(-lnlike)
 }


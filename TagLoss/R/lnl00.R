lnl00 <-
function(xx,par)
 {
    lnlike=0
    for (i in 0:16)
    {
      p1=p1(i,par)
      delta.v=delta(i,par)
      lnlike=lnlike+(xx[2,i+1]+xx[3,i+1])*log(p1-p1^2*delta.v) +
          xx[4,i+1]*log(p1^2*delta.v) + xx[1,i+1]*log(1-2*p1+p1^2*delta.v)
    }
    cat("\n par =",par," -lnl=",-lnlike)
    return(-lnlike)
 }


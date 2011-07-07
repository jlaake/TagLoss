chi.sq <-
function(x)
{
   p1=(x[3]+x[4])/sum(x)
   p2=(x[2]+x[4])/sum(x)
   expected=sum(x)*c((1-p1)*(1-p2),(1-p1)*p2,p1*(1-p2),p1*p2)
   return(sum((x-expected)^2/expected))
}


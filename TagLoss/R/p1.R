p1 <-
function(age,par)
{
  1/(1+exp(-par[1]))*(1 - exp( - (age/par[2])^( - par[4])))
}


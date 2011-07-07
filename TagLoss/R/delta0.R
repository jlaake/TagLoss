delta0 <-
function(age,par,offset.v)
{
    1/(1+exp(-par[3]*age^2-offset.v))
}


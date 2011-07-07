\name{jointp}
\alias{jointp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
jointp(par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
%%     ~~Describe \code{par} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

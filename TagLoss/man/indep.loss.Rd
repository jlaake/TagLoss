\name{indep.loss}
\alias{indep.loss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
indep.loss(p, n01, n10, n11)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{n01}{
%%     ~~Describe \code{n01} here~~
}
  \item{n10}{
%%     ~~Describe \code{n10} here~~
}
  \item{n11}{
%%     ~~Describe \code{n11} here~~
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
function(p,n01,n10,n11)
{
  n1=n01+n10
  n2=n11
  if(!is.matrix(p))
    return(n1*log(2*p*(1-p))+n2*log((1-p)^2)-(n1+n2)*log(1-p^2))
  else
    return(n01*log(p[,1]*(1-p[,2]))+n10*log(p[,2]*(1-p[,1]))+n11*log((1-p[,1])*(1-p[,2]))-(n1+n2)*log(1-p[,1]*p[,2]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

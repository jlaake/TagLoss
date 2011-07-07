\name{lnl.fi}
\alias{lnl.fi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lnl.fi(xx, par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{
%%     ~~Describe \code{xx} here~~
}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

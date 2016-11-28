\name{arp.mcmc}
\alias{arp.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Normal AR(p) MCMC
%%  ~~function to do ... ~~
}
\description{Performs a simple MCMC procedure for sampling from the posteriors of the parameters of a normal AR(p). 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
arp.mcmc(xdata, porder, n.iter = 1000, n.warmup = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xdata}{ The time series to be analyzed.
%%     ~~Describe \code{xdata} here~~
}
  \item{porder}{ Order of the AR. Must be less than 10.
%%     ~~Describe \code{porder} here~~
}
  \item{n.iter}{ Number of iterations used in the sample.
%%     ~~Describe \code{n.iter} here~~
}
  \item{n.warmup}{ Number of initial samples discarded. Must be at least 100 less than \code{n.iter}. 
%%     ~~Describe \code{n.warmup} here~~
}
}
\details{See Section 8.4 for details.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ A plot and a summary of the resulting samples for each parameter (including the variance) is returned.  In addition, the samples are returned invisibly as a vector of time series in parameter order and the variance last.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{D.S. Stoffer
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{}
\keyword{ ts }


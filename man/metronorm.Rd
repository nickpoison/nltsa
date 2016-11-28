\name{metronorm}
\alias{metronorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Metropolis Normal 
%%  ~~function to do ... ~~
}
\description{Metropolis algorithm for a multivariate normal.   
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
metronorm(CovProp, CovPi, nbiter, Xinit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ 
  \item{CovProp}{ Proposed normal covariance matrix.
%%     ~~Describe \code{CovProp} here~~
}
  \item{CovPi}{ Target normal covariance matrix.
%%     ~~Describe \code{CovPi} here~~
}
  \item{nbiter}{ Number of iterations.
%%     ~~Describe \code{nbiter} here~~
}
  \item{Xinit}{ Initial value of the algorithm.
%%     ~~Describe \code{Xinit} here~~
}
}
%%\details{
%%%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{X }{ Sample values, a vector of length \code{nbiter} - returned invisibly.}
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
\examples{
u = metronorm(3,1,1000,0)   # univariate standard normal
X = ts(u$X[-(1:500)])
plot(X)
acf(X)
plot(density(X))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }


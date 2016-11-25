\name{bi.coh}
\alias{bi.coh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Normalized Bispectrum
%%  ~~function to do ... ~~
}
\description{Produces a plot of the normalized bispectra of a series.  Larger probabilities indicate departure from linearity.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bi.coh(x.data, color = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.data}{ The time series.
%%     ~~Describe \code{x.data} here~~
}
  \item{color}{ If FALSE, the graphic is produced in gray scale.
%%     ~~Describe \code{color} here~~
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{\item{prob }{Matrix of probabilities returned invisibly}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Hinich, M. and Wolinsky, M. (2005). Normalizing bispectra. \emph{Journal of Statistical
Planning and Inference}, 130, 405--411.
%% ~put references to the literature/web site here ~
}
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
\examples{bi.coh(sunspots)}
\keyword{ ts }   % __ONLY ONE__ keyword per line

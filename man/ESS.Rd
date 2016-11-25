\name{ESS}
\alias{ESS}
\title{Effective Sample Size of importance weights}
\usage{
  ESS(weights)
}
\arguments{
  \item{weights}{Array of weights: they do not need to sum
  to 1}
}
\value{
  Effective Sample Size of importance weight
}
\description{
  Compute the efficient sample size of a set of possibly
  un-normalized weights.
}
\details{
  The weights will be automatically normalized by the
  function. If the sum of the weights is null or not
  finite, raise a warning and return \code{+Inf}.
}
\seealso{
  \code{\link{CV}}
}


\name{CV}
\alias{CV}
\title{Coefficient of Variation of importance weights}
\usage{
  CV(weights)
}
\arguments{
  \item{weights}{Array of weights: they do not need to sum
  to 1}
}
\value{
  Coefficient of variation of the weights
}
\description{
  Compute the coefficient of variation of a set of possibly
  un-normalized weights.
}
\details{
  The weights will be automatically normalized by the
  function. If the sum of the weights is null or not
  finite, raise a warning and return \code{+Inf}.
}
\seealso{
  \code{\link{ESS}}
}


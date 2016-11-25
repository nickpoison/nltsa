\name{normalized.exponential}
\alias{normalized.exponential}
\title{Normalized exponential of a vector}
\usage{
  normalized.exponential(x)
}
\arguments{
  \item{x}{Vector}
}
\value{
  exp(x)/sum(exp(x))
}
\description{
  Compute exp(x)/sum(exp(x)) in a numerically stable way
}
\details{
  Computing exp(x)/sum(exp(x)) can lead to numerical
  instability if some of the components of x are very large
  or very small. This procedure recenters by the median to
  ensure a bit more stability. It is not bullet-proof but
  is sufficient for practical cases.

  The computation of the median of a vector takes O(N log
  N), whereas taking the mean takes O(N). The mean can thus
  be used instead of the median if this function turns out
  to be a bottleneck, but the median is more robust to the
  few extreme components that typically cause numerical
  instability.
}


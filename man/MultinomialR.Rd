\name{MultinomialR}
\alias{MultinomialR}
\title{Multinomial resampling}
\usage{
  MultinomialR(p, n = length(p))
}
\arguments{
  \item{p}{Vector of probabilities, internally renormalized
  to sum to 1}

  \item{n}{Number of samples to take}
}
\value{
  Vector of length n containing i.i.d. samples taking
  values in 1:length(p) with probabilities proportional to
  p
}
\description{
  Multinomial resampling
}


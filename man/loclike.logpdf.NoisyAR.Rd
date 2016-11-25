\name{loclike.logpdf.NoisyAR}
\alias{loclike.logpdf.NoisyAR}
\title{Compute the local likelihood of a NoisyAR model}
\usage{
  loclike.logpdf.NoisyAR(self, x, y, ...)
}
\arguments{
  \item{self}{Instance of the NoisyAR class}

  \item{x}{Current state}

  \item{y}{Current observation}
}
\value{
  Logarithm of the local likelihood p(y|x)
}
\description{
  Compute the local likelihood of a NoisyAR model
}


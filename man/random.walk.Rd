\name{random.walk}
\alias{random.walk}
\title{Gaussian Random-walk proposal kernel}
\usage{
  random.walk(from, ...)
}
\arguments{
  \item{from}{State to start from}

  \item{...}{Other arguments passed on to
  \code{\link{rmvnorm}}}
}
\value{
  New state randomly sampled from a multivariate Gaussian
  centered on \code{from} with, by default, unit covariance
  matrix
}
\description{
  A Gaussian random-walk proposal kernel with identity
  covariance matrix
}
\details{
  This is an example of a proposal kernel for the
  Metropolis-Hastings \code{\link{MH}} and the
  Particle-Marginal Metropolis Hastings
  \code{\link{ParticleMarginalMH}} algorithms.

  All it does is take one random Gaussian step around the
  current state: \code{random.walk <- function(from, ...)
  rmvnorm(n=1, mean=from, ...)}

  The user is invited to write his own proposal kernel
  based on this simple example. At the very least, it is
  recommended to tune the covariance matrix used.
}


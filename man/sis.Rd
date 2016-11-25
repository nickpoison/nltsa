\name{sis}
\alias{sis}
\title{Sequential Importance Sampling with Prior kernel for 1-D NLSS}
\usage{
  sis(nlss, y, N)
}
\arguments{
  \item{nlss}{Non-linear state space model}

  \item{y}{Sequence of observations. Its length T is the
  number of timesteps.}

  \item{N}{Number of particles}

  \item{resampling}{Resampling scheme to use}
}
\value{
  A list with the following components:
  \item{particles}{Array (T, N) of the sampled particles}
  \item{logweights}{Array (T, N) of the \strong{logarithm}
  of the \strong{non-normalized} importance weights of the
  particles} \item{weights}{Array (T, N) of the
  \strong{normalized} importance weights of the particles}
  \item{t}{Indices 1 to T, included for ease of plotting}
}
\description{
  Runs sequential importance sampling \strong{without}
  resamplingon a given Non-Linear State Space models with
  prior kernel as proposal.
}
\details{
  This algorithm is not recommended and included for
  illustrative purposes only. This version is therefore a
  minimalistic and only supports NLLS with univariate
  states. Use \code{\link{sisr}} instead.

  The variant \code{\link{siskernel}} allows for arbitrary
  proposal kernel.
}
\seealso{
  \code{\link{siskernel}} \code{\link{sisr}}
}


\name{siskernel}
\alias{siskernel}
\title{Sequential Importance Sampling with arbitrary kernel for 1-D NLSS}
\usage{
  siskernel(nlss, y, N, proposal.rnd = prior.rnd,
    proposal.logpdf = prior.logpdf,
    initial.proposal.rnd = initial.rnd,
    initial.proposal.logpdf = initial.logpdf)
}
\arguments{
  \item{nlss}{Non-linear state space model}

  \item{y}{Sequence of observations. Its length T is the
  number of timesteps.}

  \item{N}{Number of particles}

  \item{proposal.rnd}{Function sampling from the proposal
  kernel to use}

  \item{proposal.logpdf}{Function computing the log-pdf of
  the proposal kernel}

  \item{initial.proposal.rnd}{Function sampling from the
  proposal kernel to use at initial timestep}

  \item{initial.proposal.logpdf}{Function computing the
  log-pdf of the proposal kernel at initial timestep}
}
\value{
  A list with the following components:
  \item{particles}{Array (T, N, D) of the sampled
  particles} \item{logweights}{Array (T, N) of the
  \strong{logarithm} of the \strong{non-normalized}
  importance weights of the particles} \item{weights}{Array
  (T, N) of the \strong{normalized} importance weights of
  the particles} \item{t}{Indices 1 to T, included for ease
  of plotting}
}
\description{
  Runs sequential importance sampling \strong{without}
  resampling on a given Non-Linear State Space models with
  user-specified kernel as proposal.
}
\details{
  This algorithm is a slightly more generic version of
  \code{\link{sis}}. It is not recommended and included for
  illustrative purposes only. This version is therefore a
  minimalistic and only supports NLLS with univariate
  states. Use \code{\link{sisr}} instead.
}
\seealso{
  \code{\link{sisr}}
}


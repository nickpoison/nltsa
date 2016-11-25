\name{sisr}
\alias{sisr}
\title{Particle filter via Sequential Importance Sampling with Resampling}
\usage{
  sisr(nlss, y, N, proposal.rnd = prior.rnd,
    proposal.logpdf = prior.logpdf,
    resampling = ResidualMultinomialR,
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

  \item{resampling}{Resampling scheme to use}
}
\value{
  A list with the following components:
  \item{particles}{Array (T, N, D) of the sampled
  particles} \item{logweights}{Array (T, N) of the
  \strong{logarithm} of the \strong{non-normalized}
  importance weights of the particles} \item{weights}{Array
  (T, N) of the \strong{normalized} importance weights of
  the particles} \item{ancestors}{Array (T, N) of the index
  of the ancestor of each particle in the previous
  generation} \item{t}{Indices 1 to T, included for ease of
  plotting} \item{success}{TRUE if the filtering succeeded}
}
\description{
  Runs a particle filter on a given Non-Linear State Space
  models
}


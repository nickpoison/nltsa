\name{ParticleMarginalMH}
\alias{ParticleMarginalMH}
\title{Run a Particle Marginal Metropolis-Hastings chain}
\usage{
  ParticleMarginalMH(initial.nlss, y, n.particles, n.mcmc,
    kernel.rnd = random.walk,
    kernel.logpdf = function(from, to) 0,
    prior.logpdf = function(nlss) 0,
    .progress = progress_text, ...)
}
\arguments{
  \item{initial.nlss}{NLSS whose parameters value will
  serve as starting point}

  \item{y}{Observations}

  \item{n.particles}{Number of particles}

  \item{n.mcmc}{Number of MCMC iterations}

  \item{kernel.rnd}{Sampling from proposal kernel on the
  parameter space, default to \code{\link{random.walk}}}

  \item{kernel.logpdf}{Evaluate the log-pdf of the proposal
  kernel}

  \item{prior.logpdf}{Evaluate the log-pdf of the prior
  distribution}

  \item{.progress}{Progress bar to use, from plyr package,
  default to progress_text}

  \item{...}{Extra parameters to pass to the SISR filter,
  see function \code{\link{sisr}}}
}
\value{
  A list with the same components as \code{\link{MH}}
}
\description{
  This function implements PMMH by calling a
  Metropolis-Hastings algorithm on the extended space.
}
\seealso{
  \code{\link{MH}}, \code{\link{sisr}},
  \code{\link{random.walk}}
}


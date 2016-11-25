\name{MH}
\alias{MH}
\title{Run a Metropolis-Hastings chain}
\usage{
  MH(initial.state, target.logpdf, n.mcmc,
    kernel.rnd = random.walk,
    kernel.logpdf = function(from, to) 0,
    .progress = progress_text)
}
\arguments{
  \item{initial.state}{Vector of starting point}

  \item{target.logpdf}{Evaluate the log-pdf of the target
  distribution}

  \item{n.mcmc}{Number of MCMC iterations}

  \item{kernel.rnd}{Sampling from proposal kernel on the
  parameter space, default to \code{\link{random.walk}}}

  \item{.progress}{Progress bar to use, from plyr package,
  default to \code{\link{progress_text}}}
}
\value{
  A list with the following components: \item{trace}{Matrix
  (\code{n.mcmc}, \code{length(initial.state)}) of accepted
  states} \item{logpdf}{Matrix (\code{n.mcmc}, 1) of the
  log-pdf of the accepted states} \item{proposed}{Matrix
  (\code{n.mcmc}, \code{length(initial.state)}) of the
  proposed states} \item{proposedpdf}{Matrix
  (\code{n.mcmc}, 1) of the log-pdf of the proposed states}
  \item{logmhratio}{Matrix (\code{n.mcmc}, 1) of the
  logarithm of Metropolish acceptance ratio}
}
\description{
  This function implements a Metropolis-Hastings MCMC
  algorithm.
}
\seealso{
  \code{\link{MH}}, \code{\link{progress_text}},
  \code{\link{random_walk}}
}


\name{CPFSGgrowth}
\alias{CPFSGgrowth}
\title{Run a Conditional Particle Stochastic Approximation SG chain}
\usage{
  CPFSGgrowth(initial.nlss, y, n.particles = function(i) {
        30 }, sa.rate = function(i) {
        min(1/(i + 1)^(2/3), 1) }, n.sg,
    .progress = progress_text, ...)
}
\arguments{
  \item{initial.nlss}{NLSS whose parameters value will
  serve as starting point}

  \item{y}{Observations}

  \item{n.particles}{Number of particles}

  \item{n.sg}{Number of SAEM iterations}

  \item{sa.rate}{Learning rate as a function of the
  iteration}

  \item{.progress}{Progress bar to use, from plyr package,
  default to progress_text}

  \item{...}{Extra parameters to pass to the SISR filter,
  see function \code{\link{sisr}}}
}
\value{
  A list with the same components as \code{\link{MH}}
}
\description{
  This function implements Particle Stochastic
  Approximation SG
}
\seealso{
  \code{\link{MH}}, \code{\link{sisr}},
  \code{\link{random.walk}}
}


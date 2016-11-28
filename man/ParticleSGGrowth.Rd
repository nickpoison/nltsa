\name{ParticleSGGrowth}
\alias{ParticleSGGrowth}
\title{Run a Particle Stochastic Gradient for the Growth model}
\usage{
  ParticleSGGrowth(initial.nlss, y,
    n.particles = function(i) {     30 },
    sa.rate = function(i) {     1/(i + 1)^(2/3) }, n.sg,
    .progress = progress_text, ...)
}
\arguments{
  \item{initial.nlss}{NLSS whose parameters value will
  serve as starting point}

  \item{y}{Observations}

  \item{n.particles}{Number of particles}

  \item{n.sg}{Number of SG iterations}

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
  This function implements Particle StochasticGradient for
  the Growth Model.
}
\seealso{
  \code{\link{MH}}, \code{\link{sisr}},
  \code{\link{random.walk}}
}


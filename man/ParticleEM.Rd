\name{ParticleEM}
\alias{ParticleEM}
\title{Run a Particle EM chain}
\usage{
  ParticleEM(initial.nlss, y, n.particles = function(i) {
        100 + i^2 }, n.em, .progress = progress_text, ...)
}
\arguments{
  \item{initial.nlss}{NLSS whose parameters value will
  serve as starting point}

  \item{y}{Observations}

  \item{n.particles}{Number of particles}

  \item{n.em}{Number of EM iterations}

  \item{.progress}{Progress bar to use, from plyr package,
  default to progress_text}

  \item{...}{Extra parameters to pass to the SISR filter,
  see function \code{\link{sisr}}}
}
\value{
  A list with the same components as \code{\link{MH}}
}
\description{
  This function implements Particle EM
}
\seealso{
  \code{\link{MH}}, \code{\link{sisr}},
  \code{\link{random.walk}}
}


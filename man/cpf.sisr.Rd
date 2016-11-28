\name{cpf.sisr}
\alias{cpf.sisr}
\title{Conditional Particle filter via Sequential Importance Sampling with Resampling}
\usage{
  cpf.sisr(x, nlss, y, N, proposal.rnd = prior.rnd,
    proposal.logpdf = prior.logpdf,
    resampling = MultinomialR,
    initial.proposal.rnd = initial.rnd,
    initial.proposal.logpdf = initial.logpdf)
}
\arguments{
  \item{x}{Trajectory on which to condition}

  \item{...}{same as classical sisr}
}
\value{
  A list with the same components as SISR, plus the
  components: \item{x}{Array (T, D) of the new selected
  trajectory}
}
\description{
  Runs a conditional particle filter on a given Non-Linear
  State Space models
}
\details{
  Warning: the resampling is Multinomial. Residual will
  bias the sampler, due to the way we overried the Nth
  particle by the conditioned trajectory.
}


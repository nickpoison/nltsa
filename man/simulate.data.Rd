\name{simulate.data}
\alias{simulate.data}
\title{Simulation of a dataset for a NLSS}
\usage{
  simulate.data(nlss, x1, T)
}
\arguments{
  \item{nlss}{An instance of a NLSS, e.g. ARCH or NoisyAR}

  \item{x1}{Initial hidden state}

  \item{T}{Number of time-steps including the initial one}
}
\value{
  \item{x}{Array (T, dimx) of hidden states} \item{y}{Array
  (T, dimy) of observations} \item{y}{Array (T, dimy) of
  observations}
}
\description{
  Generates a sequence of hidden states and of observations
  for a given NLSS. This NLSS can be any class that
  supports the methods loclike.rnd, prior.rnd, and has
  attributes spaces with fields dimx and dimy.
}


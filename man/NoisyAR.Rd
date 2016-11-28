\name{NoisyAR}
\alias{NoisyAR}
\title{Create a Noisy Auto-Regressive(1) NLSS model}
\usage{
  NoisyAR(phi, sigmaW, sigmaV)
}
\arguments{
  \item{phi}{State auto-regressive parameter}

  \item{sigmaW}{Dynamic noise standard deviation}

  \item{sigmaV}{Observation error standard deviation}
}
\value{
  S3 object of class NoisyAR
}
\description{
  Create a Noisy Auto-Regressive(1) NLSS model
}


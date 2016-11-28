\name{ARCH}
\alias{ARCH}
\title{Create an autoregressivemodel with ARCH(1) error}
\usage{
  ARCH(alpha0, alpha1, sigmaV)
}
\arguments{
  \item{alpha0}{Intercept of the variance as affine
  function of the square of the state}

  \item{alpha1}{Linear coefficient of the variance as
  affine function of the square of the state}

  \item{sigmaV}{Observation error standard deviation}
}
\value{
  S3 object of class ARCH
}
\description{
  Example NLSS model
}
\details{
  Xt = sigmaW(Xt-1) Wt

  Yt = Xt + sigmaV Vt,

  Wt and Vt are iid N(0,1) and sigmaW(x) = sqrt(alpha0 +
  alpha1*x^2) where alpha0 and alpha1 are positive.
}


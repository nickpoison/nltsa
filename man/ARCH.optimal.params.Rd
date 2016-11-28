\name{ARCH.optimal.params}
\alias{ARCH.optimal.params}
\title{Parameters for ARCH optimal kernel}
\usage{
  ARCH.optimal.params(arch, xpast, y)
}
\arguments{
  \item{arch}{Instance of the ARCH class}

  \item{xpast}{Value of the state at the former time-step}

  \item{y}{Observation at the current time-step}
}
\value{
  List with fields \item{mean}{Mean of the Gaussian kernel}
  \item{sd}{Standard deviation of the Gaussian kernel}
}
\description{
  Compute the mean and standard deviation of p(X_k|x_{k-1},
  y_k) for the ARCH model
}


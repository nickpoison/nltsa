\name{ARCH.optimal.rnd}
\alias{ARCH.optimal.rnd}
\title{Sample from ARCH optimal kernel}
\usage{
  ARCH.optimal.rnd(arch, xpast, y, ...)
}
\arguments{
  \item{arch}{Instance of the ARCH class}

  \item{xpast}{Value of the state at the former time-step}

  \item{y}{Observation at the current time-step}
}
\value{
  Sample from p(X_k|x_{k-1} = xpast, y_k = y)
}
\description{
  Sample from p(X_k|x_{k-1}, y_k) for the ARCH model
}


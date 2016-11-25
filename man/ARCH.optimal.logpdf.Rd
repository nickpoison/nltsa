\name{ARCH.optimal.logpdf}
\alias{ARCH.optimal.logpdf}
\title{Log-likelihood of  ARCH optimal kernel}
\usage{
  ARCH.optimal.logpdf(arch, xpast, xnew, y, ...)
}
\arguments{
  \item{arch}{Instance of the ARCH class}

  \item{xpast}{Value of the state at the former time-step}

  \item{xnew}{Value of the state at the current time-step}

  \item{y}{Observation at the current time-step}
}
\value{
  Logarithm of the optimal kernel density
  p(x_k=xnew|x_{k-1}=xpast, y_k=y)
}
\description{
  Evaluate p(X_|x_{k-1}, y_k) for the ARCH model
}


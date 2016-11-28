#' Coefficient of Variation of importance weights
#'
#' Compute the coefficient of variation of a set of possibly un-normalized weights.
#' 
#' The weights will be automatically normalized by the function.
#' If the sum of the weights is null or not finite, raise a warning and return \code{+Inf}.
#' 
#' @param weights Array of weights: they do not need to sum to 1
#'
#' @return Coefficient of variation of the weights
#' @seealso \code{\link{ESS}}
#' @export
CV <- function(weights) {
  s <- sum(weights)
  if(!is.finite(s) || s == 0) {
    warning('Sum of weights is null or NA')
    return(Inf)
  }
  N <- length(weights)
  sqrt(sum((N*weights/s - 1)^2)/N)
}

#' Effective Sample Size of importance weights
#'
#' Compute the efficient sample size of a set of possibly un-normalized weights.
#' 
#' The weights will be automatically normalized by the function.
#' If the sum of the weights is null or not finite, raise a warning and return \code{+Inf}.
#' 
#' @param weights Array of weights: they do not need to sum to 1
#'
#' @return Effective Sample Size of importance weight
#' @seealso \code{\link{CV}}
#' @export
ESS <- function(weights) {
  s <- sum(weights)  
  if(!is.finite(s) || s == 0) {
    warning('Sum of weights is null or NA')
    return(0)
  }
  1/sum((weights/sum(weights))^2)
}

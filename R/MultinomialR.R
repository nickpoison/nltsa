#' Multinomial resampling
#' 
#' @param p Vector of probabilities, internally renormalized to sum to 1
#' @param n Number of samples to take
#' @return Vector of length n containing i.i.d. samples taking values in 1:length(p) with probabilities proportional to p
#' @export
MultinomialR <- function (p, n=length(p)) {
  if(sum(p) == 0) stop("too few positive probabilities")
  rep(1:length(p), times=rmultinom(1, n, p))
}

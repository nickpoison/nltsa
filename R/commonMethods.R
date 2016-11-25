#' Log-pdf of the local likelihood of a Non linear state space model
#' 
#' @param nlss Instance of a nlss class
#' @param ... Passed to the specific implementation
#' @seealso \link{NoisyAR}
#' @export
loclike.logpdf <- function(nlss, ...) UseMethod("loclike.logpdf", nlss)


loclike.rnd <- function(nlss, ...) UseMethod("loclike.rnd", nlss)
prior.logpdf <- function(nlss, ...) UseMethod("prior.logpdf", nlss)
prior.f  <- function(nlss, ...) UseMethod("prior.f", nlss)
prior.rnd  <- function(nlss, ...) UseMethod("prior.rnd", nlss)
loclike.h  <- function(nlss, ...) UseMethod("loclike.h", nlss)
initial.rnd  <- function(nlss, ...) UseMethod("initial.rnd", nlss)
initial.logpdf  <- function(nlss, ...) UseMethod("initial.logpdf", nlss)
in.support  <- function(nlss, ...) UseMethod("in.support", nlss)

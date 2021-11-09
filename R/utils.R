#' Vector of frequencies of frequencies
#'
#' A function to calculate the abundance vector, or frequencies of frequencies of discrete or partly discrete
#' data vector `x`. The abundance vector is used as input in the functions `dPD()`, `MLEp()`, and `LMTp()`.
#' @param x Data vector `x`.
#' @keywords abundance
#' @details This function is equivalent to `table(table(x))`.
#' @return This function returns a named vector with the frequencies of the frequencies in the data vector x.
#' The function `base::table(x)` returns a contingency table with the frequencies in the input data vector `x` as
#' values. The `names(table(x))` are the unique values in data vector `x`. In `abundance(x)`,
#' the unique values in `table(x)` become the names of the values, while the values
#' themselves are the frequencies of the frequencies of data vector `x`.
#' @export
#' @examples
#' set.seed(111)
#' x<-rpois(10,10)
#' ## The frequency table of x:
#' print(table(x))
#' ## The frequency table of the frequency table of x:
#' abundance(x)


abundance<-function(x) {
  return(table(table(x)))
}




abundances<-function(freqs, freqs0=NULL) {
  if (!is.null(freqs0)) {
    L<-list(freqs, freqs0)
    freqs<-tapply(unlist(L), names(unlist(L)), sum)
  }
  return(table(freqs))
}

lognRF <- function(psi, n) {
  return(sum(log((rep(psi, n)+seq(0, n-1)))))
}

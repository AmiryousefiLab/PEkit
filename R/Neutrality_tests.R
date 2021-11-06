#' Test for the shape of the distribution
#'
#' This function performs a statistical test on the null hypothesis that a given
#' sample's underlying distribution is the Poisson-Dirichlet distribution. It
#' calculates a test statistic that is then used to gain a p-value from an
#' empirical distribution of the statistic from simulated samples from a
#' PD distribution.
#' @param x A discrete data vector.
#' @param rounds How many samples are simulated to obtain the empirical distribution.
#' @returns A p-value.
#' @export
#' @details The calculated test statistic is \deqn{W=\sum_{i=1}^n n_i^2 / n ,}
#' which is calculated from the sample. Here \eqn{n_i} are the frequencies of each unique value in the sample.
#' The MLE of \eqn{\psi} is then estimated from the sample with the function `MLEp()`, and an amount of samples
#' equal to the input parameter `rounds` are generated with that estimate of \eqn{psi}
#' and sample size \eqn{n}. The test statistic \eqn{W} is then calculated for each of thowe simulated samples.
#' The original \eqn{W} is then given a p-value based on what percentage of the simulated \eqn{W} it exceeds.
#' @references Watterson, G.A., (1978), The homozygosity test of neutrality. Genetics. 88(2):405-417.
#' @examples
#' ##Test whether a typical sample follows PD:
#' x<-rPD(1000,10)
#' tPD(x, 1000)
#'
#' ##Test whether a very atypical sample where frequencies of different values
#' ## are similar:
#'
#' x<-c(rep(1, 200), rep(2, 200), rep(3, 200), rep(4, 200), rep(5,200))
#' tPD(x,1000)

tPD<-function(x, rounds) {
  freq<-table(x)
  n<-length(x)
  W<-sum(freq**2)/n**2
  psi<-MLEp(table(table(x)))
  W.emp<- sapply(1:rounds, function(x) sum(table(rPD(n, psi))**2)/n**2 )

  return(sum(W>W.emp)/rounds)
}


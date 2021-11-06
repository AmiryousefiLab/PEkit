#' Maximum Likelihood Estimate of \eqn{\psi}
#'
#'
#' Numerically searches for the MLE of \eqn{\psi} given an abundance vector with a binary search algorithm.
#' @param abund An abundance vector.
#' @keywords maximum likelihood estimate \eqn{\psi}
#' @details Numerically searches for the MLE of \eqn{\psi} as the root of equation
#' \deqn{K=\sum_{i=1}^n\psi/(\psi+i-1),} where \eqn{K} is the observed number of
#' different species in the sample. The right side of the equation is monotonically
#' increasing when \eqn{\psi>0}, so a binary search is used to find the root.
#' An accepted \eqn{\psi} sets value of the right side
#' of the equation within R's smallest possible value of the actual value of \eqn{K}.
#' @return The MLE of \eqn{\psi}.
#' @export
#' @references W.J. Ewens, The sampling theory of selectively neutral alleles, Theoretical Population Biology, Volume 3, Issue 1,
#' 1972, Pages 87-112, ISSN 0040-5809, <https://doi.org/10.1016/0040-5809(72)90035-4>.
#' @examples
#' ##Find the MLE of psi of the vector (1,2,2).
#' ##The frequencies of the frequencies of the data vector are given as input:
#' MLEp(abundance(c(1,2,2)))
#'
#' ##Find the MLE of psi of a sample from the Poisson-Dirichlet distribution:
#' set.seed(1000)
#' x<-rPD(n=10000, psi=100)
#' MLEp(abundance(x))



MLEp<- function(abund) {
  n<-sum(as.integer(names(abund))*abund)
  k<-sum(abund)
  psi<-1
  asum<-0
  last<-psi/2

  while(abs(asum-k)>.Machine$double.xmin) {


    if (asum<k && last==psi/2) {
      last<-psi
      psi<-2*psi
    } else if (asum<k){
      psi1<-psi
      psi<- psi + abs(psi-last)/2
      last<-psi1
    } else if (asum>k && last==psi/2){
      last<-psi
      psi<-0.5*psi
    } else {
      psi1<-psi
      psi<- psi - abs(psi-last)/2
      last<-psi1
    }

    if(last==psi) {break}

    asum<-sum(psi/(psi+seq(0,n-1)))
  }
  return(psi)
}






#' Bootstrap confidence interval for the MLE of \eqn{\psi}
#'
#' A bootstrapped confidence interval for the Maximum Likelihood Estimate for
#' \eqn{\psi}.
#' @param x A categorical data vector.
#' @param level Level of confidence interval as number between 0 and 1.
#' @param rounds Number of bootstrap rounds. Default is 1000.
#' @param frac Percentage of data `x` used for each bootstrap round. 0.8 by default with accepted values between 0 and 1.
#' @export
#' @return The MLE of \eqn{\psi} as well as lower and upper bounds of the bootstrap
#' confidence interval.
#' @examples
#' ## Find a 95% -confidence interval for the MLE of psi given a sample from the
#' ## Poisson-Dirichlet distribution:
#' x<-rPD(n=10000, psi=100)
#' MLEp.bsci(x, 0.95)
#'


MLEp.bsci<-function(x, level=0.95, rounds=1000, frac=0.8) {
  if(frac<0.0 || frac>1.0) {
    print("frac must be a number between 0 and 1")
    break
  }

  n_bootstrap<- rounds
  n<-floor(frac*length(x))
  bootsrap_psis<-c()
  for (i in 1:n_bootstrap) {
    bootsrap_psis<-append(bootsrap_psis, MLEp(abundance(sample(x,n))))
  }
  bounds<-c((1-level)/2, 1-(1-level)/2)
  return(c("MLE"=MLEp(abundance(x)), stats::quantile(bootsrap_psis,bounds)))
}



#' Lagrange Multiplier Test for \eqn{\psi}
#'
#' Performs the Lagrange Multiplier test for the equality of the dispersion parameter \eqn{\psi} of a sample.
#' The null hypothesis of the test is \eqn{H_0: \psi = \psi_0}, where \eqn{\psi_0} is given as input here.
#' @param psi Target positive number \eqn{\psi_0} to be tested. Accepted values are "a" for absolute value 1,
#' "r" for relative value \eqn{n} (sample size), or any positive number.
#' @param abund An abundance vector of a sample.
#' @return The statistic \eqn{S} and a p-value of the two-sided test of the hypothesis.
#' @references Radhakrishna Rao, C, (1948), Large sample tests of statistical
#' hypotheses concerning several parameters with applications to problems of
#' estimation. Mathematical Proceedings of the Cambridge Philosophical Society,
#'  44(1), 50-57. <https://doi.org/10.1017/S0305004100023987>
#' @keywords score test
#' @details Calculates the Lagrange Multiplier test statistic \deqn{S\, = \,U(\psi_0)^2 / I(\psi_0),}
#' where \eqn{U} is the log-likelihood function of \eqn{\psi} and \eqn{I} is its Fisher information.
#' The statistic \eqn{S} follows \eqn{\chi^2}-distribution with 1 degree of freedom
#' when the null hypothesis \eqn{H_0:\psi=\psi_0} is true.
#' @export
#' @examples
#' ## Test the psi of a sample from the Poisson-Dirichlet distribution:
#' set.seed(10000)
#' x<-rPD(1000, 10)
#' ## Find the abundance of the data vector:
#' abund=abundance(x)
#' ## Test for the psi that was used, as well as a higher and a lower one:
#' LMTp(abund, 10)
#' LMTp(abund, 15)
#' LMTp(abund, 5)
#' LMTp(abund)       #test for psi=1
#' LMTp(abund, "r")  #test for psi=n

LMTp <- function(abund, psi="a") {
  n<-sum(as.integer(names(abund))*abund)
  k<-sum(abund)
  if (psi=="a") {
    psi<-1
  } else if (psi=="r") {
    psi<-n
  } else if (!is.numeric(psi) || psi<=0) {
    return("Psi must be a positive real number.")
  }
  asum<-sum(1/(psi + seq(0,n-1)))
  S<-(k/psi - asum)**2/(asum/psi - sum(1/(psi + seq(0,n-1))**2))
  return(c("p-value"=1-stats::pchisq(S,1), "S"=S))
}

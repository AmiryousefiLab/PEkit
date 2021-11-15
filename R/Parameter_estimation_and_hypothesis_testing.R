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
#' @param x A data vector.
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
    return(NULL)
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
#' sample.test(abund, 10)
#' sample.test(abund, 15)
#' sample.test(abund, 5)
#' sample.test(abund)       #test for psi=1
#' sample.test(abund, "r")  #test for psi=n

sample.test <- function(abund, psi="a") {
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


#' Two sample test for \eqn{\psi}
#'
#' Likelihood ratio test for the hypotheses \eqn{H_0: \: \psi_1=\psi_2} and
#' \eqn{H_1: \: \psi_1 \neq \psi_2}, where \eqn{\psi_1} and \eqn{\psi_2} are the
#' dispersal parameters of two input samples `s1` and `s2`.
#' @param s1,s2 The two data vectors to be tested.
#' @return Gives a vector with the Likelihood Ratio Test -statistic `Lambda`, as well as the
#' p-value of the test `p`.
#' @references Neyman, J., & Pearson, E. S. (1933). On the problem of the most
#' efficient tests of statistical hypotheses. Philosophical Transactions of the
#' Royal Society of London. Series A, Containing Papers of a Mathematical Or
#' Physical Character, 231(694-706), 289-337. <doi:10.1098/rsta.1933.0009>.
#' @details Calculates the Likelihood Ratio Test statistic
#' \deqn{-2log(L(\hat{\psi})/L(\hat{\psi}_1, \hat{\psi}_2)),}
#' where L is the likelihood function of observing the two input samples given
#' a single \eqn{\psi} in the numerator and two different parameters \eqn{\psi_1}
#' and \eqn{\psi_2} for each sample respectively in the denominator. According
#' to the theory of Likelihood Ratio Tests, this statistic converges in
#' distribution to a \eqn{\chi_d^2}-distribution under the null-hypothesis, where \eqn{d} is the
#' difference in the amount of parameters between the considered models, which
#' is 1 here. To calculate the statistic, the Maximum Likelihood Estimate for
#' \eqn{\psi_1,\: \psi_2} of \eqn{H_1} and the shared \eqn{\psi} of \eqn{H_0}
#' are calculated.
#' @export
#' @examples ##Create samples with different n and psi:
#' set.seed(111)
#' x<-rPD(500, 15)
#' y<-rPD(1000, 20)
#' z<-rPD(800, 30)
#' ##Run tests
#' two.sample.test(x,y)
#' two.sample.test(x,z)
#' two.sample.test(y,z)


two.sample.test<- function(s1 ,s2) {
  #estimate psi1 and psi2 for H_1 and psi for H_0

  #n1<-sum(as.integer(names(abund))*abund)
  n1<-length(s1)
  n2<-length(s2)
  #k<-sum(abund)
  k<-length(unique(s1)) + length(unique(s2))
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

    asum<-sum(psi/(psi+seq(0,n1-1))) + sum(psi/(psi+seq(0,n2-1)))
  }

  psi1<-MLEp(abundance(s1))
  psi2<-MLEp(abundance(s2))
  #psi1<-MLEp(s1) if we use abunds
  #psi2<-MLEp(s2)

  #test statistic
  L<- -2*( dlPD(abundance(s1), psi) + dlPD(abundance(s2), psi) -
             dlPD(abundance(s1), psi1) - dlPD(abundance(s2), psi2))
  p<-stats::pchisq(L, 1, lower.tail=F)

  return(c("Lambda"=L, "p-value"=p))
}




#' Test for \eqn{\psi} of multiple samples
#'
#' Likelihood ratio test for the hypotheses \eqn{H_0: \: \psi_1=\psi_2=...=\psi_d} and
#' \eqn{H_1: \: \psi_1 \neq \psi_2 \neq ... \neq \psi_d}, where \eqn{\psi_1,\psi_2,}...\eqn{,\psi_d} are the
#' dispersal parameters of the \eqn{d} samples in the columns of the input data array `x`.
#' @param x The data array to be tested. Each column of `x` is an independent sample.
#' @return Gives a vector with the Likelihood Ratio Test -statistic `Lambda`, as well as the
#' p-value of the test `p`.
#' @references Neyman, J., & Pearson, E. S. (1933). On the problem of the most
#' efficient tests of statistical hypotheses. Philosophical Transactions of the
#' Royal Society of London. Series A, Containing Papers of a Mathematical Or
#' Physical Character, 231(694-706), 289-337. <doi:10.1098/rsta.1933.0009>.
#' @details Calculates the Likelihood Ratio Test statistic
#' \deqn{-2log(L(\hat{\psi})/L(\hat{\psi}_1, \hat{\psi}_2, ..., \hat{\psi}_d)),}
#' where L is the likelihood function of observing the \eqn{d} input samples given
#' a single \eqn{\psi} in the numerator and \eqn{d} different parameters \eqn{\psi_1,\psi_2,}...\eqn{,\psi_d}
#' for each sample respectively in the denominator. According
#' to the theory of Likelihood Ratio Tests, this statistic converges in
#' distribution to a \eqn{\chi_{d-1}^2}-distribution when the null-hypothesis is true, where \eqn{d-1} is the
#' difference in the amount of parameters between the considered models. To
#' calculate the statistic, the Maximum Likelihood Estimate for
#' \eqn{\psi_1,\: \psi_2,\: ..., \: \psi_d} of \eqn{H_1} and the shared \eqn{\psi} of \eqn{H_0}
#' are calculated.
#' @export
#' @examples ##Create samples with different n and psi:
#' set.seed(111)
#' x<-rPD(1200, 15)
#' y<-c( rPD(1000, 20), rep(NA, 200) )
#' z<-c( rPD(800, 30), rep(NA, 400) )
#' samples<-cbind(cbind(x, y), z)
#' ##Run test
#' mult.sample.test(samples)

mult.sample.test<- function(x) {
  #estimate psi1 and psi2 for H_1 and psi for H_0

  #n1<-sum(as.integer(names(abund))*abund)
  d<-(dim(x)[2])
  N<-apply(x, 2, function(z) sum(!is.na(z)))
  #k<-sum(abund)
  k<-sum(apply(x, 2, function(z) length(unique(z[!is.na(z)]))))

  #MLE of shared psi
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

    asum<-sum(sapply(N, function(z) sum(psi/(psi+seq(0,z-1)))))
  }
  #shared psi is solved in variable psi

  #psi for each of d features
  psid<-apply(x, 2, function(z) MLEp(abundance(z[!is.na(z)])))
  #test statistic

  denominator<- c()
  for (i in 1:d) {
    denominator<-c(denominator, dlPD(abundance(x[,i]), psid[i]))
  }

  L<- -2* ( sum(apply(x, 2, function(z) dlPD(abundance(z[!is.na(z)]), psi))) - sum(denominator) )
  p<-stats::pchisq(L, d-1, lower.tail=F)

  return(c("Lambda"=L, "p-value"=p))
}










#' The Poisson-Dirichlet distribution
#'
#' Distribution function for the Poisson-Dirichlet distribution.
#' @usage dPD(abund, psi)
#' @param abund An abundance vector.
#' @param psi Dispersal parameter \eqn{\psi}. Accepted input values are positive real
#' numbers, "a" for absolute value \eqn{\psi}=1 by default, or "r" for relative
#' value \eqn{\psi=n}, where \eqn{n} is the size of the input sample.
#' @return The probability of the Poisson-Dirichlet distribution for the input
#' abundance vector, e.g. an exchangeable random partition,  and a dispersal parameter \eqn{\psi}.
#' @keywords Poisson-Dirichlet distribution
#' @export
#' @details Given an abundance vector `abunds`, calculates the probability
#' of a data vector `x` given by the Poisson-Dirichlet distribution. The higher the
#' dispersal parameter \eqn{\psi}, the higher the amount of distinct observed
#' species. In terms of the paintbox process, a high \eqn{\psi} increases the
#' size of the continuous part \eqn{p_0} of the process, while a low \eqn{\psi} will increase
#' the size of the discrete parts \eqn{p_{\neq 0}}.
#' @references W.J. Ewens, The sampling theory of selectively neutral alleles, Theoretical Population Biology, Volume 3, Issue 1,
#' 1972, Pages 87-112, ISSN 0040-5809, <\doi{10.1016/0040-5809(72)90035-4}>.
#' @examples
#' ## Get a random sample from the Poisson Dirichlet distribution, and
#' ## find the probability of such a sample with psi=5:
#' set.seed(111)
#' s <- rPD(n=100,psi=5)
#' a=abundance(s)
#' dPD(a, psi=5)
#'



dPD <- function(abund, psi="a") {
  n<-sum(as.integer(names(abund))*abund)
  if (psi=="a") {
    psi<-1
  } else if (psi=="r") {
    psi<-n
  } else if (!is.numeric(psi) || psi<=0) {
    return("Psi must be a positive number.")
  }
  return(exp(dlPD(abund, psi)))
}




dlPD <- function(abund, psi="a") {
  rho<-abund
  t<-as.integer(names(abund))
  n<-sum(t*rho)

  if (psi=="a") {
    psi<-1
  } else if (psi=="r") {
    psi<-n
  } else if (!is.numeric(psi) || psi<=0) {
    return("Psi must be a positive number.")
  }

  product<-sum(rho*log(psi/t) - lfactorial(rho))


  return(lfactorial(n) - lognRF(psi,n) + product)
}





#' Random sampling from the Poisson-Dirichlet Distribution
#'
#' rPD samples randomly from the PD distribution with a given \eqn{\psi} by simulating the Hoppe urn model.

#' @param n number of observations.
#' @param psi dispersal parameter.
#' @usage rPD(n, psi)
#' @return Returns a vector with a sample of size \eqn{n} from the Hoppe urn model with parameter \eqn{\psi}.
#' @keywords Poisson-Dirichlet distribution
#' @export
#' @references Hoppe, F.M. The sampling theory of neutral alleles and an urn model in population genetics.
#'  J. Math. Biology 25, 123–159 (1987). <\doi{10.1007/BF00276386}>.
#' @references W.J. Ewens, The sampling theory of selectively neutral alleles, Theoretical Population Biology, Volume 3, Issue 1,
#' 1972, Pages 87-112, ISSN 0040-5809, <\doi{10.1016/0040-5809(72)90035-4}>.
#' @details Samples random values with a given \eqn{\psi} from the Poisson-Dirichlet distribution by simulating the Hoppe urn model.
#' @examples
#' ## Get random sample from the PD distribution with different psi,
#' ## and estimate the psi of the samples:
#' s1<-rPD(1000, 10)
#' s2<- rPD(1000, 50)
#' print(c(MLEp(abundance(s1)), MLEp(abundance(s2))))
#'


rPD <- function(n, psi) {
  PDsample <- c(1)
  thesample<-c(1)
  for (i in 2:n) {
    newInt <- sample.int(length(PDsample)+1, 1, prob=c(PDsample/(i-1+psi), psi/(i-1+psi)))
    if (newInt==length(PDsample)+1) { PDsample<-append(PDsample,1) } else { PDsample[newInt]<-PDsample[newInt]+1}
    thesample<-append(thesample, newInt)
  }

  return(thesample)
}


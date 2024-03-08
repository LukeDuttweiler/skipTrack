#This document holds all of the functions that provide random draws from full conditional
#posteriors for the 'Duttweiler' algorithm.

#' Sample a value from the full conditional posterior of mu
#'
#' In our model the data are drawn from LogN(mu_i + log(c_{ij}), tau_i). The prior for mu_i
#' is given as N(mu, rho). This function draws from the conditional posterior of mu, given
#' that the prior on mu is N(priorMean, priorPre).
#'
#' @param muI Numeric vector, log of individuals mean values.
#' @param rho Numeric > 0, a sampled precision of the mu_i values
#' @param priorMean Numeric, prior mean of mu, default is log(30) to match common menstrual cycle lengths
#' @param priorPre Numeric > 0, prior precision of mu, default is 1 to have little influence
#'
#' @return Numeric
#' @export
#'
#' @examples
#' mui <- rep(log(31), 10)
#' postMu(mui, 10)
postMu <- function(muI, rho, priorMean = log(30), priorPre = 1){
  #n is the length of muI
  n <- length(muI)

  #Set Posterior Mean and precision
  postPre <- priorPre + n*rho
  postMean <- (priorPre*priorMean + rho*sum(muI))/postPre

  #Draw value for mu and return
  return(rnorm(1, mean = postMean, sd = sqrt(1/postPre)))
}

postBeta <- function(rho0, rho, Xi, muI){
  #Assuming Xi is (num Individuals)x(dimension of beta) matrix of covariates
  XTX <- t(Xi) %*% Xi
  Xmu <- t(Xi) %*% muI


}

#' Sample a value from the full conditional posterior of rho (UPDATED)
#'
#' In our model the data are drawn from LogN(mu_i + log(c_{ij}), tau_i). The prior for mu_i
#' is given as N(mu, rho). This function draws from the conditional posterior of rho, given
#' that the prior on rho is a uniform prior on the standard deviation.
#' Note that we parameterize with RATE, not SCALE.
#'
#' @param muI Numeric vector, log of individuals mean values.
#' @param mu Numeric, a sampled mean of the mu_i values
#'
#' @return Numeric
#' @export
#'
#' @examples
#' mui <- rep(log(31), 10)
#' postRho(mui, log(30))
postRho <- function(muI, xib){
  #n is the length of muI
  n <- length(muI)

  #Set posterior parameters
  postA <- (n-1)/2
  postB <- sum((muI - xib)^2)/2

  #Draw value for rho and return
  return(rgamma(1, shape = postA, rate = postB))
}

#' Sample a value from the full conditional posterior of mu_i (UPDATED)
#'
#' In our model the data are drawn from LogN(mu_i + log(c_{ij}), tau_i). The prior for mu_i
#' is given as N(x_i^T %*% beta, rho). This function draws from the conditional posterior of mu_i.
#'
#' Additionally, note that in order to vectorize the remainder of the MCMC algorithm
#' this function returns the sampled value repeated for length(yij)
#'
#' @param yij Numeric vector, cycle lengths for a single individual
#' @param cij Positive Integer vector, a sampled vector of length(yij) where the corresponding
#'  values in cij indicate a sampled number of TRUE cycles in each cycle length given by yij
#' @param taui Numeric > 0, A sampled precision for the yijs
#' @param xib Numeric, result of multiplying x_i^T %*% beta
#' @param rho Numeric > 0, sampled prior precision of mu_i
#'
#' @return Numeric vector, repeated sampled value of length(yij)
#' @export
postMui <- function(yij, cij, taui, xib, rho){
  #Ni is the length of yij
  Ni <- length(yij)

  #Set posterior mean and precision
  postPre <- Ni*taui + rho
  postMean <- (rho*xib + taui*sum(log(yij/cij)))/postPre

  #Draw from posterior and return
  dr <- rnorm(1, mean = postMean, sd = sqrt(1/postPre))
  return(rep(dr, Ni))
}

#' Sample a value from the full conditional posterior of tau_i
#'
#' In our model the data are drawn from LogN(mu_i + log(c_{ij}), tau_i). The prior for tau_i
#' is given as Gamma(priorA, priorB). This function draws from the conditional
#' posterior of tau_i. Note that we parameterize with RATE, not SCALE.
#'
#' Additionally, note that in order to vectorize the remainder of the MCMC algorithm
#' this function returns the sampled value repeated for length(yij)
#'
#' @param yij Numeric vector, cycle lengths for a single individual
#' @param cij Positive Integer vector, a sampled vector of length(yij) where the corresponding
#'  values in cij indicate a sampled number of TRUE cycles in each cycle length given by yij
#' @param mui Numeric, log of sampled mean of individuals yijs
#' @param priorA Numeric > 0, prior shape of tau_i
#' @param priorB Numeric > 0, prior rate of tau_i
#'
#' @return Numeric vector, repeated sampled value of length(yij)
#' @export
#'
#' @examples
#' ys <- rnorm(10, 30, 1)
#' cs <- rbinom(10, 2, .1) + 1
#' postTaui(ys, cs, log(30))
postTaui <- function(yij, cij, mui, priorA = .01, priorB = .01){
  #Ni is the length of yij
  Ni <- length(yij)

  #Set posterior parameters
  postA <- priorA + Ni/2
  postB <- priorB + sum((log(yij/cij) - mui)^2)/2

  #Draw from posterior and return
  dr <- rgamma(1, shape = postA, rate = postB)
  return(rep(dr, Ni))
}

#' Sample a vector of values from the full conditional posterior of the c_ij vector
#'
#' In our model the data are drawn from LogN(mu_i + log(c_{ij}), tau_i) The prior for
#'  c_ij is a categorical prior with category probabilities pi1, ..., pik, and c_ij can
#'  take values 1, ..., k where k is the length of pi. This function samples from the
#'  full conditional posterior of all c_ijs, given vectors of equal length yijs, muis, tauis
#'
#' @param yijs Numeric Vector, cycle lengths
#' @param pi Numeric vector, must sum to 1. Sampled probabilities for c_ijs
#' @param muis Numeric vector, log of sampled mean for all individuals yijs
#' @param tauis Numeric vector > 0, sampled precision for all individuals yijs
#'
#' @return Integer vector
#' @export
#'
#' @examples
#' tst <- simTrackData(1000, skipProb = c(.8, .15, .05))
#' postCij(tst$TrackedCycles, pi = c(.7, .2, .1), muis = tst$LogMean, tauis = tst$LogPrec)
postCij <- function(yijs, pi, muis, tauis){
  probs <- sapply(1:length(pi), function(j){
    #likelihood
    lik <- dlnorm(yijs, meanlog = muis + log(j), sdlog = sqrt(1/tauis))
    return(lik*pi[j])
  })
  probs <- pmax(probs, rep(0,length(probs)))

  return(glmnet::rmult(probs))
}

#' Sample a value from the full conditional posterior of pi
#'
#' In our model the data are drawn from LogN(mu_i + log(c_{ij}), tau_i) The prior for
#'  c_ij is a categorical prior with category probabilities pi1, ..., pik, and c_ij can
#'  take values 1, ..., k where k is the length of pi. This function samples from the
#'  posterior of pi = pi1, ..., pik, assuming that pi follows Dirichlet(priorAlphas)
#'
#' @param ci Integer vector, all of the sampled cij values for all individuals
#' @param priorAlphas Numeric vector, prior dirichlet parameters for pi
#'
#' @return Numeric vector
#' @export
#'
#' @examples
#' cs <- rbinom(1000, 2, .2) + 1
#' postPi(cs, priorAlphas = c(1,1,1))
postPi <- function(ci, priorAlphas){
  #Get posterior dirichlet parameters
  postAlphas <- sapply(1:length(priorAlphas), function(k){
    return(priorAlphas[k] + sum(ci == k))
  })

  #Sample from dirichlet and return
  return(as.numeric(LaplacesDemon::rdirichlet(1, alpha = postAlphas)))
}

#This document holds all of the functions that provide random draws from full conditional
#posteriors.

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
  #Check to make sure there are more than one mui values. If not, warn.
  if(length(muI) <= 1){
    warning('muI vector has only one observation. Did you mean to do this?')
  }

  #n is the length of muI
  n <- length(muI)

  #Set Posterior Mean and precision
  postPre <- priorPre + n*rho
  postMean <- (priorPre*priorMean + rho*sum(muI))/postPre

  #Draw value for mu and return
  return(rnorm(1, mean = postMean, sd = sqrt(1/postPre)))
}

#' Sample a value from the full conditional posterior of rho
#'
#' In our model the data are drawn from LogN(mu_i + log(c_{ij}), tau_i). The prior for mu_i
#' is given as N(mu, rho). This function draws from the conditional posterior of rho, given
#' that the prior on rho is Gamma(priorA, priorB). Note that we parameterize with RATE, not
#' SCALE.
#'
#' @param muI Numeric vector, log of individuals mean values.
#' @param mu Numeric, a sampled mean of the mu_i values
#' @param priorA Numeric > 0, prior shape of rho
#' @param priorB Numeric > 0, prior rate of rho
#'
#' @return Numeric
#' @export
#'
#' @examples
#' mui <- rep(log(31), 10)
#' postRho(mui, log(30))
postRho <- function(muI, mu, priorA = 1, priorB = 1){
  #Check to make sure there are more than one mui values. If not, warn.
  if(length(muI) <= 1){
    warning('muI vector has only one observation. Did you mean to do this?')
  }

  #n is the length of muI
  n <- length(muI)

  #Set posterior parameters
  postA <- priorA + n/2
  postB <- priorB + sum((muI - mu)^2)/2

  #Draw value for rho and return
  return(rgamma(1, shape = postA, rate = postB))
}

#' Sample a value from the full conditional posterior of mu_i
#'
#' In our model the data are drawn from LogN(mu_i + log(c_{ij}), tau_i). The prior for mu_i
#' is given as N(mu, rho). This function draws from the conditional posterior of mu_i.
#'
#' @param yij Numeric vector, cycle lengths for a single individual
#' @param cij Positive Integer vector, a sampled vector of length(yij) where the corresponding
#'  values in cij indicate a sampled number of TRUE cycles in each cycle length given by yij
#' @param taui Numeric > 0, A sampled precision for the yijs
#' @param mu Numeric, sampled prior mean of mu_i
#' @param rho Numeric > 0, sampled prior precision of mu_i
#'
#' @return Numeric
#' @export
#'
#' @examples
#' ys <- rnorm(10, 30, .1)
#' cs <- rep(1, 10)
#' postMui(ys, cs, 4, log(30), 1)
postMui <- function(yij, cij, taui, mu, rho){
  #If length of yij and cij do not match, throw error
  if(length(yij) != length(cij)){
    stop('Lengths of yij and cij must match.')
  }

  #If length of taui is more than 1, throw error
  if(length(taui) > 1){
    stop('taui cannot be a vector')
  }

  #Ni is the length of yij
  Ni <- length(yij)

  #Set posterior mean and precision
  postPre <- Ni*taui + rho
  postMean <- (rho*mu + taui*sum(log(yij/cij)))/postPre

  #Draw from posterior and return
  return(rnorm(1, mean = postMean, sd = sqrt(1/postPre)))
}

#' Sample a value from the full conditional posterior of tau_i
#'
#' In our model the data are drawn from LogN(mu_i + log(c_{ij}), tau_i). The prior for tau_i
#' is given as Gamma(priorA, priorB). This function draws from the conditional
#' posterior of tau_i. Note that we parameterize with RATE, not SCALE.
#'
#' @param yij Numeric vector, cycle lengths for a single individual
#' @param cij Positive Integer vector, a sampled vector of length(yij) where the corresponding
#'  values in cij indicate a sampled number of TRUE cycles in each cycle length given by yij
#' @param mui Numeric, log of sampled mean of individuals yijs
#' @param priorA Numeric > 0, prior shape of tau_i
#' @param priorB Numeric > 0, prior rate of tau_i
#'
#' @return Numeric
#' @export
#'
#' @examples
#' ys <- rnorm(10, 30, 1)
#' cs <- rbinom(10, 2, .1) + 1
#' postTaui(ys, cs, log(30))
postTaui <- function(yij, cij, mui, priorA = 1, priorB = .001){
  #If length of yij and cij do not match, throw error
  if(length(yij) != length(cij)){
    stop('Lengths of yij and cij must match.')
  }

  #If length of mui is more than 1, throw error
  if(length(mui) > 1){
    stop('mui cannot be a vector')
  }

  #Ni is the length of yij
  Ni <- length(yij)

  #Set posterior parameters
  postA <- priorA + Ni/2
  postB <- priorB + sum((log(yij/cij) - mui)^2)/2

  #Draw from posterior and return
  return(rgamma(1, shape = postA, rate = postB))
}

#' Sample a value from the full conditional posterior of c_ij
#'
#' In our model the data are drawn from LogN(mu_i + log(c_{ij}), tau_i) The prior for
#'  c_ij is a categorical prior with category probabilities pi1, ..., pik, and c_ij can
#'  take values 1, ..., k where k is the length of pi. This function samples from the
#'  full conditional posterior of c_ij
#'
#' @param yij1 Numeric, one cycle length for individual i
#' @param pi Numeric vector, must sum to 1. Sampled probabilities for c_ij
#' @param mui Numeric, log of sampled mean of individuals yijs
#' @param taui Numeric > 0, A sampled precision for the yijs
#'
#' @return Integer
#' @export
#'
#' @examples
#' postCij(yij1 = 60, pi = c(.7, .2, .1), mui = log(30), taui = 1)
postCij <- function(yij1, pi, mui, taui){
  #yij1 should be of length 1, otherwise stop
  if(length(yij1) != 1){
    stop('yij1 cannot be a vector')
  }

  #Get non-normalized probabilities
  probs <- sapply(1:length(pi), function(j){
    #likelihood
    lik <- dlnorm(yij1, meanlog = mui + log(j), sdlog = sqrt(1/taui))
    return(lik*pi[j])
  })
  probs <- pmax(probs, rep(0,length(probs)))

  if(all(probs == 0)){
    return(sample(1:length(pi), 1))
  }
  #Sample from probs and return
  return(sample(1:length(pi), 1, prob = probs))
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

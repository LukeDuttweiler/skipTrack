#This document holds all of the functions that provide random draws from full conditional
#posteriors for the 'Li' algorithm.

postLambdai <- function(yij, sij, priorK, priorG){
  #n is the number of sampled values
  n <- length(yij)

  #set posterior K and G
  postK <- priorK + sum(yij)
  postG <- priorG + n + sum(sij)

  #Return gamma draw
  return(rgamma(1, shape = postK, rate = postG))
}

postSij <- function(yijs, pii, lambdai, S){
  probs <- sapply(0:S, function(s){
    #Get truncated geometric probability
    p <- (pii^s)(1-pii)/(1-pii^(S+1))

    #Get likelihood
    lik <- dpois(yijs, lambda = lambdai*(s+1))

    return(p*lik)
  })
  probs <- pmax(probs, rep(0, length(probs)))

  return(glmnet::rmult(probs))
}

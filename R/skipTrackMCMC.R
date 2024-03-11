#' Perform MCMC sampling for identifying skips in cycle tracking.
#'
#' This function runs a Markov Chain Monte Carlo (MCMC) algorithm to update parameters in a hierarchical model (specified in...),
#' to identify skips in menstrual cycle tracking.
#'
#' @param cycleDat A data.frame with columns 'Individual' and 'TrackedCycles', representing the observed cycle lengths for each individual.
#' @param initialParams A list of initial parameter values for the MCMC algorithm.
#' Default values are provided for pi, muis, tauis, mu, rho, and cs.
#' @param reps The number of MCMC iterations (steps) to perform. Default is 1000.
#'
#' @return A list containing the MCMC draws for each parameter at each iteration. Each element
#' in the list is itself a list containing:
#' \describe{
#'   \item{ijDat}{A data.frame with updated parameters at the individual-observation level: Individual, ys, cs.}
#'   \item{iDat}{A data.frame with updated parameters at the individual level: Individual, mus, taus.}
#'   \item{mu}{Updated value of the global parameter mu.}
#'   \item{rho}{Updated value of the global parameter rho.}
#'   \item{pi}{Updated value of the global parameter pi.}
#'   \item{priorAlphas}{Vector of prior alpha values for updating pi.}
#' }
#'
#' @examples
#' # Example usage of skipTrackMCMC function
#' cycleDat <- simTrackData(10, skipProb = c(.7, .2, .1)) #Simulated Data
#' result <- skipTrackMCMC(cycleDat)
#'
#' #MORE REALISTIC VERSION
#' #COMMENTED TO SKIP DURING PACKAGE BUILD
#' #cycleDat <- simTrackData(1000, skipProb = c(.7, .2, .1)) #Simulated Data
#' #result <- skipTrackMCMC(cycleDat)
#'
#' @seealso \code{\link{gibbsStep}}
#'
#' @keywords mcmc hierarchical model skipped tracking cycles
#'
#' @export
#'
skipTrackMCMC <- function(Y,cluster,
                          X = matrix(1, nrow = length(unique(cluster))),
                          Z = matrix(1, nrow = length(unique(cluster))),
                          numSkips = 10,
                          reps = 1000,
                          initialParams = list(pi = rep(1/(numSkips+1), numSkips+1),
                                               muis = rep(log(30),
                                                          length(unique(cluster))),
                                               tauis = rep(5,
                                                          length(unique(cluster))),
                                               mu = log(30),
                                               rho = 1,
                                               cs = sample(1:3,
                                                           length(Y),
                                                           replace = TRUE),
                                               alphas = rep(1, numSkips +1),
                                               Beta = matrix(rep(0, ncol(as.matrix(X))),1),
                                               Gamma = matrix(rep(0, ncol(as.matrix(Z))),1),
                                               rhoBeta = 1,
                                               rhoGamma = 1,
                                               phi = .1)){
  #Checks for X and Z
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  if(nrow(X) != length(unique(cluster))){
    stop('X must be a num_individuals x num_covariates matrix')
  }
  if(nrow(Z) != length(unique(cluster))){
    stop('Z must be a num_individuals x num_covariates matrix')
  }


  #Organize data into initial list
  iDat <- data.frame('Individual' = unique(cluster),
                     'mus' = initialParams$muis,
                     'taus' = initialParams$tauis,
                     'thetas' = exp(Z %*% t(initialParams$Gamma)))
  ijDat <- data.frame('Individual' = cluster,
                      'ys' = Y,
                      'cs' = initialParams$cs,#)
                      'muis' = sapply(cluster, function(ind){iDat$mus[iDat$Individual == ind]}),
                      'tauis' = sapply(cluster, function(ind){iDat$taus[iDat$Individual == ind]}))
  fullDraws <- vector('list', reps + 1)
  fullDraws[[1]] <- list(ijDat = ijDat,
                         iDat = iDat,
                         rho = initialParams$rho,
                         pi = initialParams$pi,
                         Xi = X,
                         Zi = Z,
                         Beta = initialParams$Beta,
                         Gamma = initialParams$Gamma,
                         priorAlphas = initialParams$alphas,
                         indFirst = !duplicated(ijDat$Individual),
                         rhoBeta = initialParams$rhoBeta,
                         rhoGamma = initialParams$rhoGamma,
                         phi = initialParams$phi)
  #Progress bar
  pb <- utils::txtProgressBar(min = 0, max = reps, style = 3)

  #Do gibbs steps
  for(t in 1:reps){
    fullDraws[[t+1]] <- do.call('sampleStep', fullDraws[[t]])
    utils::setTxtProgressBar(pb, t)
  }
  return(fullDraws)
}

sampleStep <- function(ijDat, iDat, rho, pi,
                       Xi, Zi, Beta, Gamma,
                       priorAlphas, indFirst,
                       rhoBeta, rhoGamma, phi){
  #Start with high level (without Gamma and thetais as those are connected)
  newBeta <- postBeta(rhoBeta = rhoBeta, rho = rho, Xi = Xi, muI = iDat$mus)

  #Set xib based on newBeta
  xib <- Xi %*% t(newBeta)

  #Continue with high level information
  newRho <- postRho(muI = iDat$mus, xib = xib)
  newPi <- postPi(ci = ijDat$cs, priorAlphas = priorAlphas)

  #Now i level
  newMuis <- lapply(iDat$Individual, function(ind){
    postMui(yij = ijDat$ys[ijDat$Individual == ind],
            cij = ijDat$cs[ijDat$Individual == ind],
            taui = iDat$taus[iDat$Individual == ind],
            xib = xib, rho = newRho)
  })
  newMuis <- do.call('c', newMuis)
  newTauis <- lapply(iDat$Individual, function(ind){
    postTaui(yij = ijDat$ys[ijDat$Individual == ind],
             cij = ijDat$cs[ijDat$Individual == ind],
             mui = iDat$mus[iDat$Individual == ind],
             thetai = iDat$thetas[iDat$Individual == ind],
             phi = phi)
  })
  newTauis <- do.call('c', newTauis)

  #High level Gamma Things
  newGamList <- postGamma(taui = iDat$taus, Zi = Zi, currentGamma = Gamma, phi = phi,
                          rhoGamma = rhoGamma)
  newGamma <- newGamList$Gamma
  newThetas <- newGamList$thetai

  #Create new i level information
  iDatNew <- data.frame(Individual = iDat$Individual,
                        mus = newMuis[indFirst],
                        taus = newTauis[indFirst],
                        thetas = newThetas)

  ijDatNew <- ijDat
  ijDatNew$muis <- newMuis
  ijDatNew$tauis <- newTauis

  ijDatNew$cs <- postCij(ijDatNew$ys, pi = newPi,
                         muis = ijDatNew$muis, tauis = ijDatNew$tauis)

  return(list(ijDat = ijDatNew, iDat = iDatNew, rho = newRho,
              pi = newPi, Xi = Xi, Zi = Zi, Beta = newBeta,
              Gamma = newGamma, priorAlphas = priorAlphas, indFirst = indFirst,
              rhoBeta = rhoBeta, rhoGamma = rhoGamma, phi = phi))
}

#' Perform one sampling step for our skipTrackMCMC model
#'
#' This function updates parameters at three levels: global (\code{mu}, \code{rho}, \code{pi}),
#' individual (\code{mus}, \code{taus}), and individual-observation (\code{ys}, \code{cs}).
#'
#' @param ijDat A data.frame containing all data at the level of individual-observation (ij),
#' with columns: Individual, ys, cs.
#' @param iDat A data.frame containing all data at the level of individual (i),
#' with columns: Individual, mus, taus.
#' @param mu Current value of the global parameter mu.
#' @param rho Current value of the global parameter rho.
#' @param pi Current value of the global parameter pi.
#' @param priorAlphas Vector of prior alpha values for updating pi.
#' @param indFirst Logical vector indicating the first instance of each individual in the full data.
#'
#' @return A list containing the updated data and parameters:
#' \describe{
#'   \item{iDat}{A data.frame with updated parameters at the individual level: Individual, mus, taus.}
#'   \item{ijDat}{A data.frame with updated parameters at the individual-observation level: Individual, ys, cs.}
#'   \item{mu}{Updated value of the global parameter mu.}
#'   \item{rho}{Updated value of the global parameter rho.}
#'   \item{pi}{Updated value of the global parameter pi.}
#'   \item{priorAlphas}{Unchanged vector of prior alpha values for updating pi.}
#'   \item{indFirst}{Unchanged logical vector indFirst}
#' }
#'
#'
#' @seealso \code{\link{postMu}}, \code{\link{postRho}}, \code{\link{postPi}},
#' \code{\link{postMui}}, \code{\link{postTaui}}, \code{\link{postCij}}
#'
#' @keywords gibbs sampling hierarchical model update parameters
#'
gibbsStep <- function(ijDat, iDat, mu, rho, pi, priorAlphas, indFirst){
  #Start with high level
  newMu <- postMu(muI = iDat$mus, rho = rho)
  newRho <- postRho(muI = iDat$mus, mu = newMu)
  newPi <- postPi(ci = ijDat$cs, priorAlphas = priorAlphas)

  #Now i level
  newMuis <- lapply(iDat$Individual, function(ind){
    postMui(yij = ijDat$ys[ijDat$Individual == ind],
            cij = ijDat$cs[ijDat$Individual == ind],
            taui = iDat$taus[iDat$Individual == ind],
            mu = newMu, rho = newRho)
  })
  newMuis <- do.call('c', newMuis)
  newTauis <- lapply(iDat$Individual, function(ind){
    postTaui(yij = ijDat$ys[ijDat$Individual == ind],
             cij = ijDat$cs[ijDat$Individual == ind],
             mui = iDat$mus[iDat$Individual == ind])
  })
  newTauis <- do.call('c', newTauis)

  iDatNew <- data.frame(Individual = iDat$Individual,
                        mus = newMuis[indFirst],
                        taus = newTauis[indFirst])

  ijDatNew <- ijDat
  ijDatNew$muis <- newMuis
  ijDatNew$tauis <- newTauis

  ijDatNew$cs <- postCij(ijDatNew$ys, pi = newPi,
                         muis = ijDatNew$muis, tauis = ijDatNew$tauis)

  return(list(ijDat = ijDatNew, iDat = iDatNew, mu = newMu,
              rho = newRho, pi = newPi, priorAlphas = priorAlphas, indFirst = indFirst))
}


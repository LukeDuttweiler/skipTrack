#' Perform MCMC sampling to identify skips in cycle tracking.
#'
#' This function runs a Markov Chain Monte Carlo (MCMC) algorithm to update parameters in a hierarchical model
#' to identify skips in menstrual cycle tracking.
#'
#' @param Y A vector of observed cycle lengths.
#' @param cluster A vector indicating the individual cluster/group membership for each observation.
#' @param X A matrix of covariates for cycle length mean. Default is a matrix of 1's.
#' @param Z A matrix of covariates for cycle length precision. Default is a matrix of 1's.
#' @param numSkips The maximum number of skips to allow. Default is 10.
#' @param reps The number of MCMC iterations (steps) to perform. Default is 1000.
#' @param fixedSkips If TRUE cycle skip information (cijs) is not updated in sample steps and the inputs are instead assumed to be true.
#' @param initialParams A list of initial parameter values for the MCMC algorithm.
#' Default values are provided for pi, muis, tauis, mu, rho, cs, alphas, Beta, Gamma, rhoBeta, rhoGamma, and phi.
#'
#' @return A list containing the MCMC draws for each parameter at each iteration. Each element
#' in the list is itself a list containing:
#' \describe{
#'   \item{ijDat}{A data.frame with updated parameters at the individual-observation level: Individual, ys, cs.}
#'   \item{iDat}{A data.frame with updated parameters at the individual level: Individual, mus, taus.}
#'   \item{rho}{Updated value of the global parameter rho.}
#'   \item{pi}{Updated value of the global parameter pi.}
#'   \item{Xi}{Matrix of covariates for cycle length mean.}
#'   \item{Zi}{Matrix of covariates for cycle length precision.}
#'   \item{Beta}{Matrix of coefficients for cycle length mean.}
#'   \item{Gamma}{Matrix of coefficients for cycle length precision.}
#'   \item{priorAlphas}{Vector of prior alpha values for updating pi.}
#'   \item{indFirst}{A logical vector indicating the first occurrence of each individual.}
#'   \item{rhoBeta}{Updated value of the global parameter rhoBeta.}
#'   \item{rhoGamma}{Updated value of the global parameter rhoGamma.}
#'   \item{phi}{Value of the parameter phi.}
#'   \item{fixedSkips}{Logical. Indicates if skips were fixed.}
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
                          fixedSkips = FALSE,
                          initialParams = list(pi = rep(1/(numSkips+1), numSkips+1),
                                               muis = rep(log(30),
                                                          length(unique(cluster))),
                                               tauis = rep(5,
                                                          length(unique(cluster))),
                                               rho = 1,
                                               cs = sample(1:3, length(Y), replace = TRUE),
                                               alphas = rep(1, numSkips +1),
                                               Beta = matrix(rep(0, ncol(as.matrix(X))),1),
                                               Gamma = matrix(rep(0, ncol(as.matrix(Z))),1),
                                               rhoBeta = 1,
                                               rhoGamma = 1,
                                               phi = .001)){
  #Set initial params default list
  ip <- list(pi = rep(1/(numSkips+1), numSkips+1),
             muis = rep(log(30),
                        length(unique(cluster))),
             tauis = rep(5,
                         length(unique(cluster))),
             rho = 1,
             cs = sample(1:3, length(Y), replace = TRUE),
             alphas = rep(1, numSkips +1),
             Beta = matrix(rep(0, ncol(as.matrix(X))),1),
             Gamma = matrix(rep(0, ncol(as.matrix(Z))),1),
             rhoBeta = 1,
             rhoGamma = 1,
             phi = .001)

  #Replace anything that needs replacing
  for(i in 1:length(initialParams)){
    nm <- names(initialParams)[i]
    ip[nm] <- initialParams[nm]
  }

  #Replace whole thing
  initialParams <- ip

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
                         phi = initialParams$phi,
                         fixedSkips = fixedSkips)
  #Progress bar
  pb <- utils::txtProgressBar(min = 0, max = reps, style = 3)

  #Do gibbs steps
  for(t in 1:reps){
    fullDraws[[t+1]] <- do.call('sampleStep', fullDraws[[t]])
    utils::setTxtProgressBar(pb, t)
  }
  return(fullDraws)
}

#' Perform a single step of the MCMC sampling process.
#'
#' This function performs a single step of the Markov Chain Monte Carlo (MCMC) algorithm to update parameters
#' in a hierarchical model used for identifying skips in menstrual cycle tracking.
#'
#' @param ijDat A data.frame with individual-observation level parameters: Individual, ys, cs, muis, tauis.
#' @param iDat A data.frame with individual level parameters: Individual, mus, taus, thetas.
#' @param rho Updated value of the global parameter rho.
#' @param pi Updated value of the global parameter pi.
#' @param Xi A matrix of covariates for cycle length mean.
#' @param Zi A matrix of covariates for cycle length precision.
#' @param Beta Matrix of coefficients for cycle length mean.
#' @param Gamma Matrix of coefficients for cycle length precision.
#' @param priorAlphas Vector of prior alpha values for updating pi.
#' @param indFirst A logical vector indicating the first occurrence of each individual.
#' @param rhoBeta Updated value of the global parameter rhoBeta.
#' @param rhoGamma Updated value of the global parameter rhoGamma.
#' @param phi Value of the parameter phi.
#' @param fixedSkips Logical. If TRUE cycle skip information (cijs) is not updated in sample steps and the inputs are instead assumed to be true.
#'
#' @return A list containing updated parameters after performing a single MCMC step.
#' The list includes:
#' \describe{
#'   \item{ijDat}{A data.frame with updated parameters at the individual-observation level: Individual, ys, cs.}
#'   \item{iDat}{A data.frame with updated parameters at the individual level: Individual, mus, taus.}
#'   \item{rho}{Updated value of the global parameter rho.}
#'   \item{pi}{Updated value of the global parameter pi.}
#'   \item{Xi}{Matrix of covariates for cycle length mean.}
#'   \item{Zi}{Matrix of covariates for cycle length precision.}
#'   \item{Beta}{Matrix of coefficients for cycle length mean.}
#'   \item{Gamma}{Matrix of coefficients for cycle length precision.}
#'   \item{priorAlphas}{Vector of prior alpha values for updating pi.}
#'   \item{indFirst}{A logical vector indicating the first occurrence of each individual.}
#'   \item{rhoBeta}{Updated value of the global parameter rhoBeta.}
#'   \item{rhoGamma}{Updated value of the global parameter rhoGamma.}
#'   \item{phi}{Value of the parameter phi.}
#'   \item{fixedSkips}{Logical. Fixed skips input.}
#' }
#'
sampleStep <- function(ijDat, iDat, rho, pi,
                       Xi, Zi, Beta, Gamma,
                       priorAlphas, indFirst,
                       rhoBeta, rhoGamma, phi, fixedSkips){
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

  if(fixedSkips){
    ijDatNew$cs <- ijDat$cs
  }else{
    ijDatNew$cs <- postCij(ijDatNew$ys, pi = newPi,
                           muis = ijDatNew$muis, tauis = ijDatNew$tauis)
  }

  return(list(ijDat = ijDatNew, iDat = iDatNew, rho = newRho,
              pi = newPi, Xi = Xi, Zi = Zi, Beta = newBeta,
              Gamma = newGamma, priorAlphas = priorAlphas, indFirst = indFirst,
              rhoBeta = rhoBeta, rhoGamma = rhoGamma, phi = phi, fixedSkips = fixedSkips))
}

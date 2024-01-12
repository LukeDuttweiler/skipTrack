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
skipTrackMCMC <- function(cycleDat,
                          initialParams = list(pi = c(1/3, 1/3, 1/3),
                                               muis = rep(log(30),
                                                          length(unique(cycleDat$Individual))),
                                               tauis = rep(5,
                                                          length(unique(cycleDat$Individual))),
                                               mu = log(30), rho = 1,
                                               cs = sample(1:3,
                                                           nrow(cycleDat),
                                                           replace = TRUE)),
                          reps = 1000){
  #Set priorAlphas (currently) as 1s for each pi level
  priorAlphas <- rep(1, length(initialParams$pi))

  #Organize data into initial list
  ijDat <- data.frame('Individual' = cycleDat$Individual,
                      'ys' = cycleDat$TrackedCycles,
                      'cs' = initialParams$cs)
  iDat <- data.frame('Individual' = unique(cycleDat$Individual),
                     'mus' = initialParams$muis,
                     'taus' = initialParams$tauis)
  fullDraws <- vector('list', reps + 1)
  fullDraws[[1]] <- list(ijDat = ijDat, iDat = iDat,
                         mu = initialParams$mu, rho = initialParams$rho,
                         pi = initialParams$pi, priorAlphas = priorAlphas)
  #Progress bar
  pb <- utils::txtProgressBar(min = 0, max = reps, style = 3)

  #Do gibbs steps
  for(t in 1:reps){
    fullDraws[[t+1]] <- do.call('gibbsStep', fullDraws[[t]])
    utils::setTxtProgressBar(pb, t)
  }
  return(fullDraws)
}

#' Perform one Gibbs sampling step for our skipTrackMCMC model
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
#'
#' @return A list containing the updated data and parameters:
#' \describe{
#'   \item{iDat}{A data.frame with updated parameters at the individual level: Individual, mus, taus.}
#'   \item{ijDat}{A data.frame with updated parameters at the individual-observation level: Individual, ys, cs.}
#'   \item{mu}{Updated value of the global parameter mu.}
#'   \item{rho}{Updated value of the global parameter rho.}
#'   \item{pi}{Updated value of the global parameter pi.}
#'   \item{priorAlphas}{Unchanged vector of prior alpha values for updating pi.}
#' }
#'
#'
#' @seealso \code{\link{postMu}}, \code{\link{postRho}}, \code{\link{postPi}},
#' \code{\link{postMui}}, \code{\link{postTaui}}, \code{\link{postCij}}
#'
#' @keywords gibbs sampling hierarchical model update parameters
#'
gibbsStep <- function(ijDat, iDat, mu, rho, pi, priorAlphas){
  #Start with high level
  newMu <- postMu(muI = iDat$mus, rho = rho)
  newRho <- postRho(muI = iDat$mus, mu = newMu)
  newPi <- postPi(ci = ijDat$cs, priorAlphas = priorAlphas)

  #Now i level
  newMuis <- sapply(iDat$Individual, function(ind){
    postMui(yij = ijDat$ys[ijDat$Individual == ind],
            cij = ijDat$cs[ijDat$Individual == ind],
            taui = iDat$taus[iDat$Individual == ind],
            mu = newMu, rho = newRho)
  })
  newTauis <- sapply(iDat$Individual, function(ind){
    postTaui(yij = ijDat$ys[ijDat$Individual == ind],
             cij = ijDat$cs[ijDat$Individual == ind],
             mui = iDat$mus[iDat$Individual == ind])
  })

  iDatNew <- data.frame(Individual = iDat$Individual,
                        mus = newMuis,
                        taus = newTauis)

  #Now ij level
  #ijDatNew <- lapply(iDat$Individual, function(ind){
  #  newCijs <- sapply(ijDat$ys[ijDat$Individual == ind],
  #                    postCij,
  #                    pi = newPi,
  #                    mui = iDatNew$mus[iDatNew$Individual == ind],
  #                    taui = iDatNew$taus[iDatNew$Individual == ind])
  #  return(data.frame(Individual = ind,
  #                    ys = ijDat$ys[ijDat$Individual == ind],
  #                    cs = newCijs))
  #})
  #ijDatNew <- do.call('rbind', ijDatNew)
  ijDatNew <- ijDat
  for(k in 1:nrow(ijDatNew)){
    ijDatNew$cs[k] <- postCij(ijDat$ys[k], pi = newPi,
                              mui = iDatNew$mus[iDatNew$Individual == ijDat$Individual[k]],
                              taui = iDatNew$taus[iDatNew$Individual == ijDat$Individual[k]])
  }
  #ijDatNew$cs <- ijUpdate(ijDat, newPi, iDatNew)

  return(list(ijDat = ijDatNew, iDat = iDatNew, mu = newMu,
              rho = newRho, pi = newPi, priorAlphas = priorAlphas))
}

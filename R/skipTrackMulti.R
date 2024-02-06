#' Perform skipTrackMCMC on multiple chains using parallel or sequential computing
#'
#' This function runs skipTrackMCMC on multiple chains,
#' either in parallel or sequentially.
#'
#' @inheritParams skipTrackMCMC
#' @inheritDotParams skipTrackMCMC
#'
#' @param chains Number of chains to run in parallel or sequentially.
#' @param useParallel Logical, indicating whether to use parallel processing.
#'
#' @return A list containing the results of skipTrackMCMC for each chain.
#' @export
#'
#' @examples
#' #EXAMPLES SKIPPED TO AVOID 10 MINUTE PACKAGE BUILD
#' #tstData <- simTrackData(1000, skipProb = c(.8, .15, .05))
#'
#' # Run skipTrackMulti with parallel processing on 4 chains
#' #result_parallel <- skipTrackMulti(cycleDat = tstData, reps = 1000, chains = 4,
#' # useParallel = TRUE)
#'
#' # Run skipTrackMulti sequentially on 3 chains
#' #result_sequential <- skipTrackMulti(cycleDat = tstData, reps = 1000, chains = 3,
#' # useParallel = FALSE)
#'
#' # View the results
#' #View(result_parallel)
#' #View(result_sequential)
#'
#' @seealso \code{\link{skipTrackMCMC}}
#'
skipTrackMulti <- function(cycleDat, reps = 1000, chains, useParallel = TRUE, li = FALSE,
                           S = NULL,
                           liHyperparams = c(kappa = 180, gamma = 6, alpha = 2, beta = 20), ...){
  #if li, do hyperparameter inference first
  if(li){
    par <- tryCatch(liInference(cycleDat = cycleDat, S = S, startingParams = liHyperparams),
                    error = function(e){
                      print(e)
                      return(liHyperparams)
                    })
  }

  #Get number of cores
  numCores <- parallel::detectCores()

  #Split into different functionality based on useParallel
  if(useParallel){
    #Make sure number of cores is more than chains
    if(numCores < chains){
      stop('The number of chains must be <= the number of cores available if using parallel.')
    }

    #Set up cluster
    cl <- parallel::makeCluster(min(c(numCores, chains)))
    doParallel::registerDoParallel(cl)

    #Run skipTrackMCMC on each worker (or liMCMC if li == TRUE)
    res <- foreach::foreach(1:chains) %dopar% {
      if(li){
        liMCMC(cycleDat = cycleDat, reps = reps, hyperparams = par, S = S)
      }else{
        skipTrackMCMC(cycleDat = cycleDat, reps = reps)
      }
    }

    #Stop cluster
    parallel::stopCluster(cl)

    #Return
    return(res)
  }else{

    res <- foreach::foreach(1:chains) %do% {
      if(li){
        liMCMC(cycleDat = cycleDat, reps = reps, hyperparams = par, S = S)
      }else{
        skipTrackMCMC(cycleDat = cycleDat, reps = reps)
      }
    }

    #Return
    return(res)
  }
}

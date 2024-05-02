#' Runs the skipTrack Model, using multiple MCMC chains
#'
#' This function runs skipTrack.MCMC on multiple chains,
#' either in parallel or sequentially.
#' If li == TRUE can perform estimation using the model given in Li et al. (2022).
#'
#' @inheritParams skipTrack.MCMC
#' @inheritDotParams skipTrack.MCMC
#'
#' @param chains Number of chains to run in parallel or sequentially.
#' @param useParallel Logical, indicating whether to use parallel processing. Default is TRUE.
#' @param li Logical, if TRUE runs the estimation given the model provided in Li et al. (2022).
#' This model does not estimate covariates, so these will be ignored.
#' @param liHyperparams Named numeric vector, the initial hyperparameters named as in the Li algorithm. Defaults
#' are provided, will only be used if li == TRUE
#'
#' @return A list containing the results of skipTrack.MCMC for each chain.
#' @export
#'
#' @seealso \code{\link{skipTrack.MCMC}}
#'
skipTrack.fit <- function(Y,cluster,
                           X = matrix(1, nrow = length(cluster)),
                           Z = matrix(1, nrow = length(cluster)),
                           numSkips = 10,
                           reps = 1000, chains, useParallel = TRUE,
                           li = FALSE,
                           liHyperparams = c(kappa = 180, gamma = 6, alpha = 2, beta = 20), ...){
  #Checks to make sure dimensions are all good
  #FILL IN

  #Sort all inputs to make sure that individual's observations are all next to each other
  Y <- Y[order(cluster)]
  X <- X[order(cluster),]
  Z <- Z[order(cluster),]
  cluster <- cluster[order(cluster)]

  #if li, do hyperparameter inference first
  if(li){
    par <- tryCatch(liInference(Y = Y, cluster = cluster,
                                S = numSkips, startingParams = liHyperparams),
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

    #Run skipTrack.MCMC on each worker (or liMCMC if li == TRUE)
    res <- foreach::foreach(1:chains) %dopar% {
      if(li){
        liMCMC(Y = Y, cluster = cluster, reps = reps, hyperparams = par, S = numSkips)
      }else{
        skipTrack.MCMC(Y = Y, cluster = cluster, X = X, Z = Z, numSkips = numSkips, reps = reps,
                      ...)
      }
    }

    #Stop cluster
    parallel::stopCluster(cl)

    #Return
    return(res)
  }else{

    res <- foreach::foreach(1:chains) %do% {
      if(li){
        liMCMC(Y = Y, cluster = cluster, reps = reps, hyperparams = par, S = numSkips)
      }else{
        skipTrack.MCMC(Y = Y, cluster = cluster, X = X, Z = Z, numSkips = numSkips, reps = reps,
                      ...)
      }
    }

    #Return
    return(res)
  }
}

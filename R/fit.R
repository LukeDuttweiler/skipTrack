#' Fits the skipTrack Model using 1 or more MCMC chains
#'
#' This function fits the model using multiple instances of skipTrack.MCMC, either in parallel or sequentially.
#' If li == TRUE can perform estimation using the model given in Li et al. (2022).
#'
#' @inheritParams skipTrack.MCMC
#' @inheritDotParams skipTrack.MCMC
#'
#' @param chains Number of chains to run.
#' @param useParallel Logical, indicating whether to use parallel processing. Default is TRUE.
#' @param li Logical, if TRUE runs the estimation given the model provided in Li et al. (2022).
#' This model does not estimate covariates, so inputs for X and Z will be ignored.
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
  if(length(cluster) != length(Y)){
    stop('Length of cluster and Y should match.')
  }

  if(length(Y) != nrow(X) | length(Y) != nrow(Z)){
    stop('X and Z should be matrices with nrow == length(Y)')
  }

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
  }else{

    res <- foreach::foreach(1:chains) %do% {
      if(li){
        liMCMC(Y = Y, cluster = cluster, reps = reps, hyperparams = par, S = numSkips)
      }else{
        skipTrack.MCMC(Y = Y, cluster = cluster, X = X, Z = Z, numSkips = numSkips, reps = reps,
                      ...)
      }
    }
  }

  #Build return object
  retObj <- list('fit' = res,
                 'data' = list('Y' = Y,
                               'cluster' = cluster,
                               'X' = X,
                               'Z' = Z),
                 'useParallel' = useParallel,
                 'model' = ifelse(li, 'li', 'skipTrack'))
  if(li){
    retObj$liHyperparams <- liHyperparams
  }

  class(retObj) <- c('skipTrack.model', class(retObj))
  return(retObj)
}

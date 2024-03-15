#' Simulate user-tracked menstrual cycle data for multiple individuals using specified method.
#'
#' This function generates synthetic data for user tracked menstrual cycles given the specified method,
#' skip probabilities, and maximum cycles. It supports built-in methods ('dutt', 'li') and custom methods.
#'
#' @param n Number of individuals to simulate data for.
#' @param method Method for data simulation. Can be a character ('dutt', 'li') or a custom function.
#' @param skipProb Vector of probabilities for number of true cycles per tracked cycle. For
#' example, (.7, .2, .1) means that 70% of observed cycles will contain one true cycle, 20%
#' will contain 2 true cycles and 10% will contain 3 true cycles. Default is NULL. If
#' method == 'li', this should be ignored.
#' @param maxCycles Maximum number of cycles for generating skip cycles. Default is the length of skipProb.
#' If method == 'li' this must be specified, if method == 'dutt', leave as default.
#'
#' @return A data.frame with information dependent on the method.
#'
#' @examples
#' # Example usage of simTrackData function using the dutt method
#' resultDutt <- simTrackData(1000, method = 'dutt', skipProb = c(.7, .2, .1))
#'
#' #Example usage using the li method
#' resultLi <- simTrackData(1000, method = 'li', maxCycles = 3)
#'
#' @seealso \code{\link{duttSim}}, \code{\link{liSim}}
#'
#' @export
simTrackData <- function(n,
                         method = c('dutt', 'li'),
                         skipProb = NULL,
                         maxCycles = length(skipProb),
                         trueBetas = NULL,
                         trueGammas = trueBetas,
                         overlap = 0,
                         xCovF = NULL,
                         zCovF = xCovF){

  #Checks for built in methods, custom methods are on their own.
  if(is.character(method)){
    #If method is li, make sure length(skipProb) is 1
    if((method[1] == 'li') & !is.null(skipProb)){
      stop('with method li, skipProb should be NULL and maxCycles should be specified.')
    }

    #If method is dutt, make sure maxCycles is length(skipProb)
    if((method[1] == 'dutt') & maxCycles != length(skipProb)){
      stop(
        'for method dutt, maxCycles should be the length of skipProb'
      )
    }
  }

  #n gives the number of individuals, for each individual simulate tracked cycles
  #given method, skipProb, and maxCycles
  if(is.character(method)){
    if(method[1] == 'dutt'){
      simDat <- lapply(1:n, duttSim, skipProb = skipProb, maxCycles = maxCycles,
                       trueBetas = trueBetas,
                       trueGammas = trueGammas,
                       overlap = overlap,
                       xCovF = xCovF,
                       zCovF = zCovF)
      simDat <- do.call('rbind', simDat)
    }else if(method[1] == 'li'){
      simDat <- lapply(1:n, liSim, skipProb = skipProb, maxCycles = maxCycles)
      simDat <- do.call('rbind', simDat)
    }else{
      stop('specified method is unknown.')
    }
  }else if(is.function(method)){
    simDat <- lapply(1:n, method, skipProb = skipProb, maxCycles = maxCycles)
    simDat <- do.call('rbind', simDat)
  }else{
    stop('method must be one of the specified characters, or a function')
  }

  #Get X and Z
  X <- as.matrix(simDat[,grepl('X|Individual', names(simDat))])
  Z <- as.matrix(simDat[,grepl('Z|Individual', names(simDat))])
  X <- unique(X)[,-1]
  Z <- unique(Z)[,-1]

  #Return as list with specific components separated out
  return(list('Y' = simDat$TrackedCycles,
              'cluster' = simDat$Individual,
              'X' = X,
              'Z' = Z,
              'Beta' = trueBetas,
              'Gamma' = trueGammas,
              'NumTrue' = simDat$NumTrue,
              'Underlying' = simDat[, grepl('Mean|Prec', names(simDat))]))
}

#' Simulate user tracked menstrual cycle data for an individual using the dutt method.
#'
#' This function generates synthetic data for user tracked menstrual cycles for a
#' single individual using the dutt method.
#'
#' @param i Individual identifier. Character, numeric or integer.
#' @param skipProb Vector of probabilities for number of true cycles per tracked cycle. For
#' example, (.7, .2, .1) means that 70% of observed cycles will contain one true cycle, 20%
#' will contain 2 true cycles and 10% will contain 3 true cycles.
#' @param maxCycles Maximum number of true cycles per tracked cycle. Ignored.
#'
#' @return A data.frame with columns 'Individual', 'TrackedCycles', 'NumTrue',
#'  'LogMean', 'LogPrec'.
#'
#'
#' @seealso \code{\link{simTrackData}}
duttSim <- function(i, skipProb, maxCycles, trueBetas, trueGammas, overlap, xCovF, zCovF){
  #For each individual, generate the number of (tracked) cycles from poisson(7)
  #(restricted to > 0)
  numCycles <- max(rpois(1, 7), 1)

  #If TRUE xCov or zCov = 0, set mean/precision to parameters specifically,
  #otherwise, create the number of requested covariates and record effects
  if(is.null(trueBetas)){
    lm <- log(30)
    xi <- NULL
  }else{
    xi <- matrix(rnorm(length(trueBetas), .25), nrow = 1)
    lm <- log(30) + xi %*% trueBetas
  }
  if(is.null(trueGammas)){
    precm <- 1000
    zi <- NULL
  }else{
    #Which x to overlap?
    whichX <- ifelse(overlap == 0, 0, 1:overlap)
    zi <- matrix(c(xi[1,whichX],
                   rnorm(length(trueGammas)-overlap, .25)), nrow = 1)
    precm <- 1200 + exp(zi %*% trueGammas)
  }

  #For each individual sample a mean (on the log scale) and precision (on the log scale)
  phi0 <- .001 #Constant goes here for phi
  prec <- max(4,rgamma(1, precm*phi0, phi0)) #Don't let precision get absurdly low
  lmean <- rnorm(1, lm, .13)

  #Sample c (true cycles per tracked cycle) values for number of cycles
  cs <- sample(1:maxCycles, numCycles, replace = TRUE, prob = skipProb)

  #Sample tracked cycle lengths
  ys <- round(rlnorm(numCycles, meanlog = lmean + log(cs), sdlog = sqrt(1/prec)))

  #Create as many fake x and z covariates as requested
  xiF <- matrix(rnorm(length(xCovF), .25), nrow = 1)
  ziF <- matrix(rnorm(length(zCovF), .25), nrow = 1)

  xi <- as.data.frame(cbind(1, xi, xiF))
  zi <- as.data.frame(cbind(1, zi, ziF))
  names(xi) <- paste0('X', 0:(ncol(xi)-1))
  names(zi) <- paste0('Z', 0:(ncol(zi)-1))

  #Return as data.frame
  df <- data.frame('Individual' = i, 'TrackedCycles' = ys, 'NumTrue' = cs,
                   'LogMean' = lmean, 'LogPrec' = prec)
  df <- cbind(df, xi, zi)
  return(df)
}

#' Simulate user tracked menstrual cycle data for an individual using the li method.
#'
#' This function generates synthetic data for user tracked menstrual cycles for a
#' single individual using the li method.
#'
#' @param i Individual identifier. Character, numeric or integer.
#' @param skipProb Vector, ignored for this method.
#' @param maxCycles Integer, Maximum possible number of true cycles per tracked cycle.
#'
#' @return A data.frame with columns 'Individual', 'TrackedCycles',
#' 'NumTrue', 'Mean', 'SkipProb'.
#'
#' @seealso \code{\link{simTrackData}}
liSim <- function(i, skipProb, maxCycles){
  #For each individual, generate the number of (tracked) cycles from poisson(7)
  #(restricted to > 0)
  numCycles <- max(rpois(1, 7), 1)

  #For each individual sample a mean
  indMean <- rgamma(1, 180, 6)

  #For each individual sample a skip probability
  indSkip <- rbeta(1, 2, 20)

  #For each tracked cycle, sample the number of skipped cycles, cap at maxCycles
  numSkips <- rgeom(numCycles, 1-indSkip)
  numSkips <- pmin(numSkips, rep(maxCycles-1, numCycles))

  #For each tracked cycle, draw a length (dependent on numSkips)
  ys <- rpois(numCycles, indMean*(1+numSkips))

  #Return as data.frame
  df <- data.frame('Individual' = i, 'TrackedCycles' = ys, 'NumTrue' = numSkips + 1,
                   'Mean' = indMean, 'SkipProb' = indSkip)
  return(df)
}

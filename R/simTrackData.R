#' Simulate user-tracked menstrual cycle data for multiple individuals using specified method.
#'
#' This function generates synthetic data for user-tracked menstrual cycles given the specified method,
#' skip probabilities, and maximum cycles. It supports built-in methods ('duttweiler', 'li', 'mixture') and custom methods.
#'
#' @param n Number of individuals to simulate data for.
#' @param method Method for data simulation. Can be a character ('duttweiler', 'li', 'mixture') or a custom function.
#' @param skipProb Vector of probabilities for number of true cycles per tracked cycle. For
#' example, (.7, .2, .1) means that 70% of observed cycles will contain one true cycle, 20%
#' will contain 2 true cycles and 10% will contain 3 true cycles. Default is NULL. If
#' method == 'li', leave as default.
#' @param maxCycles Maximum number of cycles for generating skip cycles. Default is the length of skipProb.
#' If method == 'li', this must be specified; if method == 'duttweiler', leave as default.
#' @param trueBetas Optional. True values for the mean regression coefficients (not counting intercept which is automatic based on the method).
#' @param trueGammas Optional. True values for the precision regression coefficients (not counting intercept which is automatic based on the method).
#' @param overlap Optional. Number of shared (non-intercept) covariates between X and Z.
#'
#' @return A list with information dependent on the method.
#'
#' @examples
#' # Example usage of simTrackData function using the dutt method
#' resultDutt <- simTrackData(1000, method = 'duttweiler', skipProb = c(.7, .2, .1))
#'
#' # Example usage using the li method
#' resultLi <- simTrackData(1000, method = 'li', maxCycles = 3)
#'
#' @seealso \code{\link{duttSim}}, \code{\link{liSim}}, \code{\link{mixSim}}
#'
#' @export
simTrackData <- function(n,
                         method = c('duttweiler', 'li', 'mixture'),
                         skipProb = NULL,
                         maxCycles = length(skipProb),
                         trueBetas = NULL,
                         trueGammas = NULL,
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
    if((method[1] == 'duttweiler') & maxCycles != length(skipProb)){
      stop(
        'for method dutt, maxCycles should be the length of skipProb'
      )
    }
  }

  #n gives the number of individuals, for each individual simulate tracked cycles
  #given method, skipProb, and maxCycles
  if(is.character(method)){
    if(method[1] == 'duttweiler'){
      simDat <- lapply(1:n, duttSim, skipProb = skipProb, maxCycles = maxCycles,
                       trueBetas = trueBetas,
                       trueGammas = trueGammas,
                       overlap = overlap,
                       xCovF = xCovF,
                       zCovF = zCovF)
      simDat <- do.call('rbind', simDat)
    }else if(method[1] == 'li'){
      simDat <- lapply(1:n, liSim, skipProb = skipProb, maxCycles = maxCycles,
                       trueBetas = trueBetas, trueGammas = trueGammas)
      simDat <- do.call('rbind', simDat)
    }else if(method[1] == 'mixture'){
      simDat <- lapply(1:n, mixSim,
                       skipProb = skipProb,
                       maxCycles = maxCycles,
                       trueBetas = trueBetas,
                       trueGammas = trueGammas,
                       overlap = overlap)
      simDat <- do.call('rbind', simDat)
    }else{
      stop('specified method is unknown.')
    }
  }else if(is.function(method)){
    simDat <- lapply(1:n, method, skipProb = skipProb, maxCycles = maxCycles)
    simDat <- do.call('rbind', simDat)
  }else{
    stop('method must be one of the specified characters, or a function taking skipProb and maxCycles')
  }

  #Get X and Z
  X <- as.matrix(simDat[,grepl('X|Individual', names(simDat))])
  Z <- as.matrix(simDat[,grepl('Z|Individual', names(simDat))])
  X <- unique(X)[,-1]
  Z <- unique(Z)[,-1]

  #Get trueBetas and trueGammas based on settings
  if(is.null(trueBetas)){
    trueBetas <- simDat$Beta0[1]
  }else{
    trueBetas <- c(simDat$Beta0[1], trueBetas)
  }
  if(is.null(trueGammas)){
    trueGammas <- simDat$Gamma0[1]
  }else{
    trueGammas <- c(simDat$Gamma0[1], trueGammas)
  }

  #Return as list with specific components separated out
  return(list('Y' = simDat$TrackedCycles,
              'cluster' = simDat$Individual,
              'X' = X,
              'Z' = Z,
              'Beta' = trueBetas,
              'Gamma' = trueGammas,
              'NumTrue' = simDat$NumTrue,
              'Underlying' = simDat[, grepl('Mean|Prec', names(simDat)), drop = F]))
}

#' Simulate user tracked menstrual cycle data for an individual using the dutt method.
#'
#' This function generates synthetic data for user tracked menstrual cycles for a
#' single individual using the duttweiler method.
#'
#' @param i Individual identifier. Character, numeric or integer.
#' @param skipProb Vector of probabilities for number of true cycles per tracked cycle. For
#' example, (.7, .2, .1) means that 70% of observed cycles will contain one true cycle, 20%
#' will contain 2 true cycles and 10% will contain 3 true cycles.
#' @param maxCycles Maximum number of true cycles per tracked cycle. Ignored.
#' @param trueBetas Optional. True values for generated mean regression coefficients.
#' @param trueGammas Optional. True values for the generated precision regression coefficients.
#' @param overlap Optional. Number of shared covariates between X and Z.
#'
#' @return A data.frame with columns 'Individual', 'TrackedCycles', 'NumTrue',
#'  'LogMean', 'LogPrec', 'Beta0', 'Gamma0', 'X0',...,'Xn', 'Z0',...,'Zm'. Where n = length(trueBetas)
#'  and m = length(trueGammas).
#'
#' @seealso \code{\link{simTrackData}}
duttSim <- function(i, skipProb, maxCycles, trueBetas, trueGammas, overlap){
  #For each individual, generate the number of (tracked) cycles from poisson(7)
  #(restricted to > 0)
  numCycles <- max(rpois(1, 7), 1)

  #If trueBetas or trueGammas don't exist, set mean/precision to given average,
  #otherwise, create the number of requested covariates and record effects
  if(is.null(trueBetas)){
    lm <- log(30)
    xi <- NULL
  }else{
    xi <- matrix(rnorm(length(trueBetas), 0), nrow = 1)
    lm <- log(30) + xi %*% trueBetas
  }
  if(is.null(trueGammas)){
    precm <- exp(5.5)
    zi <- NULL
  }else{
    #Which x to overlap?
    if(overlap == 0){
      whichX <- 0
    }else{
      whichX <- 1:overlap
    }
    zi <- matrix(c(xi[1,whichX],
                   rnorm(length(trueGammas)-overlap, 0)), nrow = 1)
    precm <- exp(5.5 + zi %*% trueGammas)
  }

  #For each individual sample a mean (on the log scale) and precision (on the log scale)
  phi0 <- .01 #Constant goes here for phi
  prec <- max(4, rgamma(1, shape = precm*phi0, rate = phi0)) #Don't let precision get absurdly low
  lmean <- rnorm(1, lm, .13)

  #Sample c (true cycles per tracked cycle) values for number of cycles
  cs <- sample(1:maxCycles, numCycles, replace = TRUE, prob = skipProb)

  #Sample tracked cycle lengths
  ys <- round(rlnorm(numCycles, meanlog = lmean + log(cs), sdlog = sqrt(1/prec)))

  xi <- as.data.frame(cbind(1, xi))
  zi <- as.data.frame(cbind(1, zi))
  names(xi) <- paste0('X', 0:(ncol(xi)-1))
  names(zi) <- paste0('Z', 0:(ncol(zi)-1))

  #Return as data.frame
  df <- data.frame('Individual' = i, 'TrackedCycles' = ys, 'NumTrue' = cs,
                   'LogMean' = lmean, 'LogPrec' = prec, 'Beta0' = log(30), 'Gamma0' = 5.5)

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
#' @param trueBetas Optional. True values for generated mean regression coefficients.
#' @param trueGammas NULL, left for consistency. Will throw error if specified.
#'
#' @return A data.frame with columns 'Individual', 'TrackedCycles',
#' 'NumTrue', 'Mean', 'SkipProb', 'Z0', 'Beta0', 'Gamma0', 'X0', ..., 'Xn' where n = length(trueBetas).
#'
#' @seealso \code{\link{simTrackData}}
liSim <- function(i, skipProb, maxCycles, trueBetas, trueGammas = NULL){
  if(!is.null(trueGammas)){
    warning('Li data generation method does not take covs for precision, trueGamma input will be ignored')
  }

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

  #If there are trueBeta values, adjust individual mean for those
  if(is.null(trueBetas)){
    chnge <- 0
    xi <- NULL
  }else{
    xi <- matrix(rnorm(length(trueBetas), 0), nrow = 1)
    chnge <- xi %*% trueBetas
  }

  #Adjust based on change (on log scale)
  indMean <- as.numeric(exp(log(indMean) + chnge))

  #For each tracked cycle, draw a length (dependent on numSkips)
  ys <- rpois(numCycles, indMean*(1+numSkips))

  #Return as data.frame
  #NOTE: Beta0 is log(30) here as our model assumes a distribution on log(Y) not just Y,
  #Gamma0 is entered as log(60). NOT SURE WHY, but this seems to give close to consistent results. Need to
  #figure out details here.
  df <- data.frame('Individual' = i, 'TrackedCycles' = ys, 'NumTrue' = numSkips + 1,
                   'Mean' = indMean, 'SkipProb' = indSkip,
                   'Z0' = 1, Beta0 = log(30), Gamma0 = log(60))

  xi <- as.data.frame(cbind(1, xi))
  names(xi) <- paste0('X', 0:(ncol(xi)-1))


  df <- cbind(df, xi)

  return(df)
}

#' Simulate user tracked menstrual cycle data for an individual using the mixture method.
#'
#' This function generates synthetic data for user tracked menstrual cycles for a
#' single individual using the mixture method.
#'
#' @param i Individual identifier. Character, numeric, or integer.
#' @param skipProb Vector of probabilities for the number of true cycles per tracked cycle. For
#' example, (.7, .2, .1) means that 70% of observed cycles will contain one true cycle, 20%
#' will contain 2 true cycles, and 10% will contain 3 true cycles.
#' @param maxCycles Maximum number of true cycles per tracked cycle.
#' @param trueBetas Optional. True values for generated mean regression coefficients.
#' @param trueGammas Optional. True values for the generated precision regression coefficients.
#' @param overlap Optional. Number of shared covariates between X and Z.
#'
#' @return A data.frame with columns 'Individual', 'TrackedCycles', 'NumTrue',
#' 'Mean', 'Beta0', 'Gamma0', 'X0',...,'Xn', 'Z0',...,'Zm', where n = length(trueBetas)
#'  and m = length(trueGammas).
#'
#' @seealso \code{\link{simTrackData}}
mixSim <- function(i, skipProb, maxCycles, trueBetas, trueGammas, overlap){
  #For each individual, generate the number of (tracked) cycles from poisson(7)
  #(restricted to > 0)
  numCycles <- max(rpois(1, 7), 1)

  #Per cycle generate numTrue
  cs <- sample(1:maxCycles, numCycles, replace = TRUE, prob = skipProb)

  #If TRUE xCov or zCov = 0, set mean/precision to parameters specifically,
  #otherwise, create the number of requested covariates and record effects
  if(is.null(trueBetas)){
    m <- 30
    xi <- NULL
  }else{
    xi <- matrix(rnorm(length(trueBetas), 0), nrow = 1)
    m <- 30 + xi %*% trueBetas
  }

  #Get individual mean based on trueBetas
  indMean <- round(rgamma(1, m*5, 5))

  #Create probabilities affected by Gammas which sort individuals into 3 regularity categories
  if(is.null(trueGammas)){
    p <- .5
    zi <- NULL
  }else{
    #Which x to overlap?
    if(overlap == 0){
      whichX <- 0
    }else{
      whichX <- 1:overlap
    }
    zi <- matrix(c(xi[1,whichX],
                   rnorm(length(trueGammas)-overlap, 0)), nrow = 1)
    linComp <- zi %*% trueGammas
    p <- exp(linComp)/(1+exp(linComp))
  }

  #Draw category given p
  indCat <- rbinom(1, 2, p)

  #Different behavior depending on category
  if(indCat == 0){
    ys <- rpois(numCycles, lambda = indMean*cs)
  }else if(indCat == 1){
    ys <- sapply(cs, function(drC){
      indMean*drC + sum(sample((-2:2), drC, replace = TRUE, prob = c(.05, .25,.4, .25, .05)))
    })
  }else if(indCat == 2){
    ys <- indMean*cs
  }

  #Create X Z dfs
  xi <- as.data.frame(cbind(1, xi))
  zi <- as.data.frame(cbind(1, zi))
  names(xi) <- paste0('X', 0:(ncol(xi)-1))
  names(zi) <- paste0('Z', 0:(ncol(zi)-1))


  df <- data.frame('Individual' = i, 'TrackedCycles' = ys, 'NumTrue' = cs,
                   'Mean' = indMean, Beta0 = log(30), Gamma0 = 15)

  df <- cbind(df, xi, zi)

  return(df)
}

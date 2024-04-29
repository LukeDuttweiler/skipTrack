#' Get a table of Inference results from skipTrack.fit
#'
#' This function calculates inference results on Betas, Gammas, and cijs based on the provided MCMC results.
#' It returns summaries such as quantiles and confidence intervals for Betas, Gammas, and cijs, as well as Gelman-Rubin diagnostics.
#' Note that true values and converage are included in the output if trueVals is supplied, but otherwise not.
#'
#' @param stFit Object result of skipTrack.fit function.
#' @param trueVals Optional list containing true values for Betas, Gammas, and cijs.
#' @param burnIn Number of skipTrack iterations to discard as burn-in.
#'
#' @return A list containing the following elements:
#'   \item{Betas}{95% credible intervals and (if trueVals is supplied) true values for Betas and Coverage tag.}
#'   \item{Gammas}{95% credible intervals and (if trueVals is supplied) true values for Gammas and Coverage tag.}
#'   \item{GammaGR}{Gelman-Rubin diagnostics for Gammas.}
#'   \item{Cijs}{Wald confidence intervals and (if trueVals is supplied) true values for cijs and Coverage tags.}
#'   \item{CijGR}{Gelman-Rubin diagnostics for cijs.}
#'
#' @examples
#' # Example usage with simulated data (which includes access to ground truth):
#' #
#' # stResults(stFitObject, trueVals = simulatedData, burnIn = 750)
#' #
#' # Example usage with real data:
#' #
#' # stResults(stFitObject, burnIn = 750)
#'
#' @export
stResults <- function(stFit, trueVals = NULL, burnIn = 750){
  #Creates a dataframe with chain/draw specific betas and gammas
  betaDF <- lapply(1:length(stFit), function(chainI){
    #Extract and arrange
    chainIbetas <- sapply(stFit[[chainI]][burnIn:length(stFit[[chainI]])], getElement, 'Beta')
    if(is.matrix(chainIbetas)){
      chainIbetas <- t(chainIbetas)
    }else{
      chainIbetas <- matrix(chainIbetas, ncol = 1)
    }

    #Add iteration and chain info
    m <- cbind(matrix(1:nrow(chainIbetas), ncol = 1),
               chainIbetas, matrix(chainI, nrow = nrow(chainIbetas)))

    #Turn into df and name correctly
    df <- as.data.frame(m)
    names(df) <- c('t', paste0('Beta', 0:(ncol(chainIbetas)-1)), 'chain')
    return(df)
  })
  betaDF <- do.call('rbind', betaDF)

  gammaDF <- lapply(1:length(stFit), function(chainI){
    chainIgammas <- sapply(stFit[[chainI]][burnIn:length(stFit[[chainI]])], getElement, 'Gamma')

    if(is.matrix(chainIgammas)){
      chainIgammas <- t(chainIgammas)
    }else{
      chainIgammas <- matrix(chainIgammas, ncol = 1)
      }

    m <- cbind(matrix(1:nrow(chainIgammas), ncol = 1),
               chainIgammas, matrix(chainI, nrow = nrow(chainIgammas)))
    df <- as.data.frame(m)
    names(df) <- c('t', paste0('Gamma', 0:(ncol(chainIgammas)-1)), 'chain')
    return(df)
  })
  gammaDF <- do.call('rbind', gammaDF)

  betaQuants <- lapply(0:(ncol(betaDF)-3), function(i){
    ql <- quantile(betaDF[,paste0('Beta', i)], .025)
    mn <- mean(betaDF[,paste0('Beta', i)])
    qu <- quantile(betaDF[,paste0('Beta', i)], .975)
    return(data.frame('Lower' = ql, 'Mean' = mn, 'Upper' = qu, 'Beta' = i))
  })
  betaQuants <- do.call('rbind', betaQuants)

  gammaQuants <- lapply(0:(ncol(gammaDF)-3), function(i){
    ql <- quantile(gammaDF[,paste0('Gamma', i)], .025)
    mn <- mean(gammaDF[,paste0('Gamma', i)])
    qu <- quantile(gammaDF[,paste0('Gamma', i)], .975)
    return(data.frame('Lower' = ql, 'Mean' = mn, 'Upper' = qu, 'Gamma' = i))
  })
  gammaQuants <- do.call('rbind', gammaQuants)

  #Get cijs and confidence intervals
  cijs <- lapply(stFit, function(ch){
    chainCs <- sapply(ch, function(dr){
      return(dr$ijDat$cs)
    })
    return(chainCs[,burnIn:ncol(chainCs)])
  })
  cijs <- do.call('cbind', cijs)

  #cijQuants <- t(apply(cijs, 1, quantile, probs = c(.025, .5, .975)))
  cijMeans <- rowMeans(cijs)
  cijSDs <- apply(cijs, 1, sd)

  cijCIs <- data.frame('Lower' = cijMeans + qnorm(.025)*cijSDs,
                       'Mean' = cijMeans,
                       'Upper' = cijMeans + qnorm(.975)*cijSDs)


  #If supplied, get truth from trueVals
  if(!is.null(trueVals)){
    trueBetas <- trueVals$Beta
    trueGammas <- trueVals$Gamma
    trueCijs <- trueVals$NumTrue

    #Attach to quantiles
    betaQuants$TrueVals <- trueBetas[1:nrow(betaQuants)]
    gammaQuants$TrueVals <- trueGammas[1:nrow(gammaQuants)]
    cijCIs$TrueVals <- trueCijs

    #Calculate coverage
    betaQuants$Coverage <- betaQuants$Lower <= betaQuants$TrueVals & betaQuants$TrueVals <= betaQuants$Upper
    gammaQuants$Coverage <- gammaQuants$Lower <= gammaQuants$TrueVals & gammaQuants$TrueVals <= gammaQuants$Upper
    cijCIs$Coverage <- cijCIs$Lower <= cijCIs$TrueVals & cijCIs$TrueVals <= cijCIs$Upper
  }

  return(list('Betas' = betaQuants,
              'Gammas' = gammaQuants,
              'GammaGR' = stDiag(stFit = stFit, 'gammas', 'lanfear')@diagnostics$gelmanRubin$`Point est.`,
              'Cijs' = cijCIs,
              'CijGR' = stDiag(stFit = stFit, 'cijs', 'lanfear')@diagnostics$gelmanRubin$`Point est.`))
}

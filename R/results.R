#Function to get inference results on Betas and Gammas
stResults <- function(mcmcRes){
  #Creates a dataframe with chain/draw specific betas and gammas
  betaDF <- lapply(1:length(mcmcRes), function(chainI){
    #Extract and arrange
    chainIbetas <- sapply(mcmcRes[[chainI]], getElement, 'Beta')
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

  gammaDF <- lapply(1:length(mcmcRes), function(chainI){
    chainIgammas <- sapply(mcmcRes[[chainI]], getElement, 'Gamma')

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
    ql <- quantile(betaDF[,paste0('Beta', i)], .05)
    qu <- quantile(betaDF[,paste0('Beta', i)], .95)
    return(data.frame('Lower' = ql, 'Upper' = qu, 'Beta' = i))
  })
  betaQuants <- do.call('rbind', betaQuants)

  gammaQuants <- lapply(0:(ncol(gammaDF)-3), function(i){
    ql <- quantile(gammaDF[,paste0('Gamma', i)], .05)
    qu <- quantile(gammaDF[,paste0('Gamma', i)], .95)
    return(data.frame('Lower' = ql, 'Upper' = qu, 'Gamma' = i))
  })
  gammaQuants <- do.call('rbind', gammaQuants)

  View(gammaQuants)
}

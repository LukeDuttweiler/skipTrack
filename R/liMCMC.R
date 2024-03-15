#NEEDS DOCUMENTATION
liMCMC <- function(Y,
                   cluster,
                   S,
                   hyperparams = c(kappa = 180, gamma = 6, alpha = 2, beta = 20),
                   initialParams = list(pi = c(1/3, 1/3, 1/3),
                                               lambdais = rep(30,
                                                          length(unique(cycleDat$Individual))),
                                               piis = rep(.2,
                                                           length(unique(cycleDat$Individual))),
                                               ss = sample(0:S,
                                                           nrow(cycleDat),
                                                           replace = TRUE)),
                   reps = 1000){
  #Create cycleDat
  cycleDat <- data.frame('TrackedCycles' = Y, 'Individual' = cluster)

  #Organize data into initial list
  iDat <- data.frame('Individual' = unique(cycleDat$Individual),
                     'lambdas' = initialParams$lambdais,
                     'pis' = initialParams$piis)
  ijDat <- data.frame('Individual' = cycleDat$Individual,
                      'ys' = cycleDat$TrackedCycles,
                      'ss' = initialParams$ss,
                      'lambdais' = sapply(cycleDat$Individual, function(ind){iDat$lambdas[iDat$Individual == ind]}),
                      'piis' = sapply(cycleDat$Individual, function(ind){iDat$pis[iDat$Individual == ind]}))
  fullDraws <- vector('list', reps + 1)
  fullDraws[[1]] <- list(ijDat = ijDat, iDat = iDat,
                         kappa = as.numeric(hyperparams['kappa']),
                         gamma = as.numeric(hyperparams['gamma']),
                         alpha = as.numeric(hyperparams['alpha']),
                         beta = as.numeric(hyperparams['beta']),
                         S = S,
                         indFirst = !duplicated(ijDat$Individual))
  #Progress bar
  pb <- utils::txtProgressBar(min = 0, max = reps, style = 3)

  #Do gibbs steps
  for(t in 1:reps){
    fullDraws[[t+1]] <- do.call('gibbsStepLi', fullDraws[[t]])
    utils::setTxtProgressBar(pb, t)
  }
  return(fullDraws)
}

gibbsStepLi <- function(ijDat, iDat, kappa, gamma, alpha, beta, S, indFirst){
  #Now i level
  newLambdais <- lapply(iDat$Individual, function(ind){
    postLambdai(yij = ijDat$ys[ijDat$Individual == ind],
                sij = ijDat$ss[ijDat$Individual == ind],
                priorK = kappa, priorG = gamma)
  })
  newLambdais <- do.call('c', newLambdais)
  newPiis <- lapply(iDat$Individual, function(ind){
    postPii(sij = ijDat$ss[ijDat$Individual == ind],
            currentPii = iDat$pis[iDat$Individual == ind],
            priorA = alpha, priorB = beta, S = S)
  })
  newPiis <- do.call('c', newPiis)

  iDatNew <- data.frame(Individual = iDat$Individual,
                        lambdas = newLambdais[indFirst],
                        pis = newPiis[indFirst])

  ijDatNew <- ijDat
  ijDatNew$lambdais <- newLambdais
  ijDatNew$piis <- newPiis

  ijDatNew$ss <- postSij(ijDatNew$ys, pii = ijDatNew$piis,
                         lambdai = ijDatNew$lambdais, S = S)

  return(list(ijDat = ijDatNew, iDat = iDatNew, kappa = kappa,
              gamma = gamma, alpha = alpha, beta = beta, S = S, indFirst = indFirst))
}

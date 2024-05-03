#Contains utility functions, particularly methods for working with S3 object 'skipTrack.model'

#' Plot skipTrack.model objects
#'
#' @param x skipTrack.model object from the function skipTrack.fit
#' @param ... Needed for S3 consistency
#'
#' @return Invisible NULL. Prints plots from skipTrack.visualize
#' @export
#'
plot.skipTrack.model <- function(x, ...){
  plts <- skipTrack.visualize(x)
  plts$nrow <- 2
  do.call(gridExtra::grid.arrange, plts)
  return(invisible(NULL))
}

#' Print skipTrack.model to console
#'
#' @param x skipTrack.model object from the function skipTrack.fit
#' @param ... Needed for S3 consistency
#'
#' @return Invisible NULL. Prints info about skipTrack.model object
#' @export
#'
print.skipTrack.model <- function(x, ...){
  cat(paste0('-------------------------------------------\n',
             'SkipTrack Model Output\n',
             '-------------------------------------------\n',
             'Number of Observations:   ', length(x$data$Y), '\n',
             'Number of Individuals:    ', length(unique(x$data$cluster)), '\n',
             'Dimension of Beta:        ', ncol(x$data$X), '\n',
             'Dimension of Gamma:       ', ncol(x$data$Z), '\n',
             '-------------------------------------------\n',
             'Number of MCMC Chains:    ', length(x$fit), '\n',
             'Model:                    ', x$model, '\n',
             '-------------------------------------------\n'))

  return(invisible(NULL))
}

#' Report skipTrack.model structure to console
#'
#' @param x skipTrack.model object from the function skipTrack.fit
#'
#' @return Invisible NULL. Prints structure of skipTrack.model object
#' @export
#'
str.skipTrack.model <- function(x){
  cat(paste0('skipTrack.model S3 Object (also a list)\n\n',
             'NChains:    ', length(x$fit), '\n',
             'Model Type: ', x$model))
  return(invisible())
}

#' Report skipTrack.model results to console
#'
#' @param object skipTrack.model object from the function skipTrack.fit
#' @param ... arguments passed to skipTrack.results
#'
#' @return Invisible skipTrack.results. Prints info about skipTrack.model object
#'
#' @export
#'
summary.skipTrack.model <- function(object, ...){
  #Get results
  res <- skipTrack.results(object, ...)

  #Format tables
  betaTab <- res$Betas[c('Mean', 'Lower', 'Upper')]
  betaTab <- round(betaTab, 3)
  row.names(betaTab) <- c('(Intercept)', colnames(object$data$X)[-1])
  names(betaTab) <- c('Estimate', '      95% CI Lower', '95% CI Upper')

  gammaTab <- res$Gammas[c('Mean', 'Lower', 'Upper')]
  gammaTab <- round(gammaTab, 3)
  row.names(gammaTab) <- c('(Intercept)', colnames(object$data$Z)[-1])
  names(gammaTab) <- c('Estimate', '      95% CI Lower', '95% CI Upper')

  diagTab <- res$Diagnostics
  row.names(diagTab) <- diagTab$Parameter
  diagTab$Parameter <- NULL
  diagTab <- round(diagTab,2)
  names(diagTab) <- c('Effective Sample Size', '      Gelman-Rubin')

  #Print summary
  cat(paste0('----------------------------------------------------\n',
             'Summary of skipTrack.fit using ', object$model, ' model\n',
             '----------------------------------------------------\n',
             'Mean Coefficients: \n\n'))
  print(betaTab)
  cat(paste0('\n----------------------------------------------------\n',
             'Precision Coefficients: \n\n'))
  print(gammaTab)
  cat(paste0('\n----------------------------------------------------\n',
             'Diagnostics: \n\n'))
  print(diagTab)
  cat(paste0('\n----------------------------------------------------\n'))

  return(invisible(res))
}

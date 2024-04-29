#' Takes mcmcResults from skipTrack Multi and uses genMCMCDiag to get mcmc diagnostics.
#'
#' @param mcmcRes A list of MCMC results from the skipTrack Multi function.
#' @param param A character string specifying the parameter for which diagnostics are to be calculated.
#'   Must be one of: 'mu', 'rho', 'muis', 'tauis', or 'cijs'.
#' @param method An optional parameter specifying the method for calculating diagnostics. Default is NULL.
#' @param ... Additional parameters to be passed to the genDiagnostic function.
#' @inheritDotParams genMCMCDiag::genDiagnostic diagnostics distance verbose
#'
#' @return A mcmcDiag object of MCMC diagnostics for the specified parameter
#' @details If the parameter is 'rho' (the univariate parameters),
#'   the function extracts the specified parameter from the MCMC results and calculates
#'   diagnostics using the genDiagnostic function with the
#'   standard method. If the parameter is 'cijs', 'muis', or 'tauis', the
#'   function extracts the corresponding values and calculates diagnostics using the genDiagnostic
#'   function with the specified or default method ('ts') and hammingDist as the distance function.
#'
#' @examples
#' #Example usage:
#' #mcmcResults <- skipTrackMulti(...) #Need to fill this out
#' #stDiag(mcmcResults, 'mu')
#'
#' @seealso \code{\link{genDiagnostic}}, \code{\link{skipTrackMulti}}
#'
#' @export
stDiag <- function(mcmcRes, param, method = NULL, ...){
  if(param %in% c('mu', 'rho')){#If param is one of the univariate parameters, get diagnostics

    #Extract specified parameter
    mcmcExt <- lapply(mcmcRes, function(chain){
      #Get draws of specified parameter
      draws <- sapply(chain, getElement, name = param)

      #Return in expected format
      return(draws)
    })

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, method = 'standard', ...))

  }else if(param == 'cijs'){ #Continued alternative methods specific to skipTrack

    #Extract cijs
    mcmcExt <- lapply(mcmcRes, function(chain){
      #Get list of cij draws
      draws <- lapply(chain, function(d){
        return(d$ijDat$cs)
      })

      #Return in expected format
      return(draws)
    })

    #set method if not specified
    if(is.null(method)){
      method <- 'ts'
    }

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, method = method,
                                      distance = genMCMCDiag::hammingDist, ...))

  }else if(param == 'betas'){ #Continued alternative methods specific to skipTrack

    #Extract cijs
    mcmcExt <- lapply(mcmcRes, function(chain){
      #Get list of Beta draws
      draws <- lapply(chain, function(d){
        return(d$Beta)
      })

      #Return in expected format
      return(draws)
    })

    #set method if not specified
    if(is.null(method)){
      method <- 'ts'
    }

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, method = method,
                                      distance = genMCMCDiag::hammingDist, ...))

  }else if(param == 'gammas'){ #Continued alternative methods specific to skipTrack

    #Extract cijs
    mcmcExt <- lapply(mcmcRes, function(chain){
      #Get list of Beta draws
      draws <- lapply(chain, function(d){
        return(d$Gamma)
      })

      #Return in expected format
      return(draws)
    })

    #set method if not specified
    if(is.null(method)){
      method <- 'ts'
    }

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, method = method,
                                      distance = genMCMCDiag::hammingDist, ...))

  }else if(param == 'muis'){
    #Extract muis
    mcmcExt <- lapply(mcmcRes, function(chain){
      #Get list of mui draws
      draws <- lapply(chain, function(d){
        return(d$iDat$mus)
      })

      #Return in expected format
      return(draws)
    })

    #set method if not specified
    if(is.null(method)){
      method <- 'ts'
    }

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, method = method,
                                      distance = genMCMCDiag::hammingDist, ...))

  }else if(param == 'tauis'){
    #Extract tauis
    mcmcExt <- lapply(mcmcRes, function(chain){
      #Get list of taui draws
      draws <- lapply(chain, function(d){
        return(d$iDat$taus)
      })

      #Return in expected format
      return(draws)
    })

    #set method if not specified
    if(is.null(method)){
      method <- 'ts'
    }

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, method = method,
                                      distance = genMCMCDiag::hammingDist, ...))

  }else if(param == 'sijs'){ #Method for LI inference

    #Extract sijs
    mcmcExt <- lapply(mcmcRes, function(chain){
      #Get list of sij draws
      draws <- lapply(chain, function(d){
        return(d$ijDat$ss)
      })

      #Return in expected format
      return(draws)
    })

    #set method if not specified
    if(is.null(method)){
      method <- 'ts'
    }

    #Calculate diagnostics and return
    return(genMCMCDiag::genDiagnostic(mcmcExt, method = method,
                                      distance = genMCMCDiag::hammingDist, ...))

  }else{ #Throw error if param is not recognized
    stop("param must be character string 'mu', 'rho', 'muis', tauis', or 'cijs'")
  }
}

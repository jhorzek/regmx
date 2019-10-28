#' getFitMeasures
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#'
#' computes multiple fit measures
#' @author Jannik Orzek
#' @import OpenMx
#' @examples
#'
#' @export
#'
getFitMeasures <- function (regModel, regType, regOn, regIndicators, cvSample = NULL, zeroThresh = .001){

  # define return value:
  return_value <- data.frame("estimated_params" = NA,
                             "m2LL"= NA,
                             "AIC"= NA,
                             "BIC"= NA,
                             "CV m2LL" = NA,
                             "Boot m2LL" = NA
  )

  # get the number of estimated parameters:

  return_value$estimated_params <- getEstimatedParameters(regModel = regModel, regType = regType,
                                                          regOn = regOn, regIndicators = regIndicators,
                                                          zeroThresh = zeroThresh)$estimatedParameters

  ### compute Fit Indices:

  # get -2LogL:
  if(regModel$Submodel$data$type == "cov"){
    ExpCov <- mxGetExpected(regModel$Submodel, "covariance")
    ExpMean <- mxGetExpected(regModel$Submodel, "means")

    # get the data set:
    ObsCov <- regModel$Submodel$data$observed
    NObs <- regModel$Submodel$data$numObs
    NManif <- length(regModel$Submodel$manifestVars)

    # compute the m2LogL:

    return_value$m2LL <- getM2LogL_cov(ExpCov = ExpCov, ObsCov = ObsCov, NObs = NObs, NManif = NManif)

  }else if(regModel$Submodel$data$type == "raw"){
    # get the expected covariance/means
    ExpCov <- mxGetExpected(regModel$Submodel, "covariance")
    ExpMean <- mxGetExpected(regModel$Submodel, "means")

    # get the data set:
    ObsData <- regModel$Submodel$data$observed

    # compute the m2LogL:
    return_value$m2LL <- getM2LogL_FIML(ExpCov = ExpCov, ExpMean = ExpMean, ObsData = ObsData)
  } else{
    stop(paste("Error: Unknown data type: ", regModel$Submodel$data$type))
  }


  # AIC
  AIC <- return_value$m2LL + 2* return_value$estimated_params
  return_value$AIC <- AIC
  # BIC
  BIC <- return_value$m2LL + log(regModel$Submodel$data$numObs) * return_value$estimated_params
  return_value$BIC <- BIC

  #### if Cross -Validation is used ####
  ## CV m2LL
  if( !is.null(cvSample) ){ # if cvsample (as mxData) is provided:
    if(!class(cvSample) == "MxDataStatic"){
      stop("Provided cvSample data set is not an mxData file.")
    }

    if(regModel$Submodel$data$type == "cov"){
      if(!cvSample$type == "cov"){
        stop("Provided cvSample data set is not of the same type as the training set.")
      }
      ExpCov <- mxGetExpected(regModel$Submodel, "covariance")
      ExpMean <- mxGetExpected(regModel$Submodel, "means")

      # get the data set:
      ObsCov <- cvSample$observed
      NObs <- cvSample$numObs
      NManif <- length(regModel$Submodel$manifestVars)

      # compute the m2LogL:

      return_value$CV.m2LL <- getM2LogL_cov(ExpCov = ExpCov, ObsCov = ObsCov, NObs = NObs, NManif = NManif)

    }else if(regModel$Submodel$data$type == "raw"){
      if(!cvSample$type == "raw"){
        stop("Provided cvSample data set is not of the same type as the training set.")
      }

        #get the expected covariance/means

        ExpCov <- mxGetExpected(regModel$Submodel, "covariance")
        ExpMean <- mxGetExpected(regModel$Submodel, "means")

        # get the data set:
        ObsData <- cvSample$observed

        # compute the m2LogL:
        return_value$CV.m2LL <- getM2LogL_FIML(ExpCov = ExpCov, ExpMean = ExpMean, ObsData = ObsData)
    }
  }
  ret <- return_value

  return(ret)


}

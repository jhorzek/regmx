getCtFitMeasures <- function (regCtModel, regType, regOn, regIndicators, cvSample = NULL, zeroThresh = .001){

  # define return value:
  return_value <- data.frame("estimated_params" = NA,
                             "m2LL"= NA,
                             "AIC"= NA,
                             "BIC"= NA,
                             "CV m2LL" = NA,
                             "Boot m2LL" = NA
  )

  # get the number of estimated parameters:
  paramsAndModel <- getEstimatedParameters(regModel = regCtModel, regType = regType,
                                           regOn = regOn, regIndicators = regIndicators,
                                           zeroThresh = zeroThresh)
  return_value$estimated_params <- paramsAndModel$estimatedParameters

  ### compute Fit Indices:

  # get -2LogL:

  temp_ct_model <- mxRun(paramsAndModel$redefinedModel, useOptimizer = F, silent = T)
  return_value$m2LL <- temp_ct_model$fitfunction$result[[1]]

  # AIC
  AIC <- return_value$m2LL + 2* return_value$estimated_params
  return_value$AIC <- AIC
  # BIC
  BIC <- return_value$m2LL + log(regCtModel$Submodel$data$numObs) * return_value$estimated_params
  return_value$BIC <- BIC

  #### if Cross -Validation is used ####
  ## CV m2LL
  if( !is.null(cvSample) ){ # if cvsample (as mxData) is provided:
    if(!class(cvSample) == "MxDataStatic"){
      stop("Provided cvSample data set is not an mxData file.")
    }

    CVModel <- regCtModel$Submodel
    CVModel$data <- cvsample
    fit_CVModel <- mxRun(CVModel, useOptimizer = F, silent = T)
    return_value$CV.m2LL <- fit_CVModel$fitfunction$result[[1]]

  }

  ret <- return_value

  return(ret)


}

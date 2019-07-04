getEstimatedParameters <- function(regModel, regType, regOn, regIndicators, zeroThresh){
  matrices <- regModel$Submodel$matrices

  redefinedModel <- regModel$Submodel

  if(regType == "lasso"){
    ## set regularized parameters to zero if they are below the threshold:
    if(regModel$regValue$values > 0){
      for(matrix in regOn){
        mat <- matrices[[matrix]]
        if(any(abs(matrices[[matrix]]$values) < zeroThresh & matrices[[matrix]]$free & regIndicators[[matrix]] == 1)){
          # select parameters that are below the threshold
          selection <- abs(matrices[[matrix]]$values) < zeroThresh & matrices[[matrix]]$free & regIndicators[[matrix]] == 1
          # set to zero
          matrices[[matrix]]$values[selection] <- 0
          # set as fixed
          matrices[[matrix]]$free[selection] <- FALSE

          # replace matrix in redefinedModel
          redefinedModel <- mxModel(redefinedModel, matrices[[matrix]])
        }
      }
    }
  }

  # get estimated parameters:
  estimatedParameters <- 0
  for(matrix in matrices){
    if(any(!is.na(matrix$labels))){
      # elements without labels that are free:
      sum1 <- sum(is.na(matrix$labels)&& matrix$free)
      # unique elements with labels that are free
      sum2 <- length(unique(matrix$labels[matrix$free]))
      estimatedParameters <- estimatedParameters + sum1 + sum2
    }else{
      estimatedParameters <- estimatedParameters + sum(matrix$free)
    }

  }

  return(list("redefinedModel" = redefinedModel, "estimatedParameters" = estimatedParameters))
}

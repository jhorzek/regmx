#' getEstimatedParameters
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#'
#' computes number of estimated parameters
#' @author Jannik Orzek
#' @import OpenMx
#' @examples
#'
#' @export
#'

getEstimatedParameters <- function(regModel, alpha, gamma, regOn, regIndicators, zeroThresh){
  matrices <- regModel$Submodel$matrices

  redefinedModel <- regModel$Submodel

  if(alpha >0){
    ## set regularized parameters to zero if they are below the threshold:
    RegValuesNames <- paste("regValues",regOn, sep ="")
    for(matrix in regOn){
      if(regModel[[paste("regValues",matrix, sep ="")]]$values > 0){

        mat <- matrices[[matrix]]
        if(any(abs(matrices[[matrix]]$values) < zeroThresh & matrices[[matrix]]$free & regIndicators[[matrix]] == 1)){
          # select parameters that are below the threshold
          selection <- abs(matrices[[matrix]]$values) < zeroThresh & matrices[[matrix]]$free & regIndicators[[matrix]] == 1
          # set to zero
          matrices[[matrix]]$values[selection] <- 0
          # set as fixed
          matrices[[matrix]]$free[selection] <- FALSE

          # replace matrix in redefinedModel
          redefinedModel <- mxModel(redefinedModel, matrices[[matrix]])}
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

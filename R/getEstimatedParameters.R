getEstimatedParameters <- function(regModel, regType, regOn, regIndicators, zeroThresh){
  estimated <- 0
  for(matrix in regModel$Submodel$matrices){
    estimated <- estimated + sum(matrix$free)
  }
  ## subtract regularized parameters that fall below zeroThresh

  setToZero <- 0
  if(regType == "lasso"){
    for(matrix in 1:length(regOn)){
      setToZero <- setToZero + sum(abs(regModel$Submodel[[regOn[matrix]]]$values[regIndicators[[matrix]]==1 & regModel$Submodel[[regOn[matrix]]]$free]) < zeroThresh)
    }
  }
  newEstimated <- estimated - setToZero
  return(newEstimated)
}

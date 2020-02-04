#' createRegValueGrid
#'
#'
#' creates a grid of tuning parameters
#' @author Jannik Orzek
#'
#' @export
createRegValueGrid <- function(regValues, regOn){
  if(is.list(regValues)){
    # create a grid with all possible combinations of regValues:
    RegValuesGrid <- as.matrix(expand.grid(regValues))

  }else{
    RegValuesGrid <- matrix(NA, nrow = length(regValues), ncol = length(regOn))
    for(col in 1:length(regOn)){
      RegValuesGrid[,col] <- regValues
    }
    colnames(RegValuesGrid) <- regOn
  }
  return(RegValuesGrid)

}

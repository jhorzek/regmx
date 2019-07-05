#' Summary function for CvOptimRegModelObject
#'
#' @param object Object from optimRegModel
#' @param ... additional parameters.
#' @export
summary.CvOptimCtRegModelObject <- function(object,...){

  FinalModel <- object$`final Model`$Submodel
  bestPenalty <- object$`best penalty`
  if(object$call$autoCV){
    print(paste("The best penalty value based on", object$call$k,"fold cross-validation" ,"was:", bestPenalty))}
  print("Parameter values for the best model:")
  for(mat in FinalModel$matrices){
    # check if there are any free parameters in the matrix
    if(any(mat$free)){
      print(paste("Values for",mat$name,":"))
      print(mat$values)

    }
  }

  print("Fit measures for different penalty values:")
  print(object$`CV results`)

  retList <- list("FinalModel" = FinalModel, "bestPenalty" = bestPenalty, "fitMeasures" = object$`CV results`)
  class(retList) <- "summary.CvOptimRegModelObject"
}

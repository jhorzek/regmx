#' Summary function for OptimRegCtModelObject
#'
#' @param object Object from optimRegCtModel
#' @param ... additional parameters.
#' @export
summary.OptimRegCtModelObject <- function(object,...){

  FinalModel <- object$bestmodel$Submodel
  bestPenalty <- object$`best penalty`
  print(paste("The best penalty value based on", object$call$criterion, "was:", bestPenalty))
  print("Parameter values for the best model:")
  for(mat in FinalModel$matrices){
    # check if there are any free parameters in the matrix
    if(any(mat$free)){
      print(paste("Values for",mat$name,":"))
      print(mat$values)

    }
  }

  print("Fit measures for different penalty values:")
  print(object$`fit measures`)

  retList <- list("FinalModel" = FinalModel, "bestPenalty" = bestPenalty, "fitMeasures" = object$`fit measures`)
  class(retList) <- "summary.OptimRegCtModelObject"
}

#' computeFinalParameters
#'
#'
#' computes the finale parameters
#' @author Jannik Orzek
#'
#' @export
computeFinalParameters <- function(minimum_criterion, mxModelObject, alpha, gamma, regOn, regIndicators, results){

  if(nrow(minimum_criterion)>1){
    warning(paste("Multiple minima found. Will use the following penalty values:", paste(as.matrix(minimum_criterion)[nrow(minimum_criterion),],collapse = ","), sep = ":"))
    best_penalty = t(as.matrix(as.matrix(minimum_criterion)[nrow(minimum_criterion),]))
  }else{
    best_penalty = minimum_criterion
  }

  if(is.matrix(best_penalty)){
    regValue <- vector("list", length = ncol(best_penalty))
    names(regValue) <- colnames(best_penalty)
    regValue[colnames(best_penalty)] <- best_penalty
  }else{
    regValue <- best_penalty
  }

  reg_Model_criterion <- regModel(mxModelObject = mxModelObject, alpha = alpha, gamma = gamma,
                             regOn = regOn, regIndicators = regIndicators,
                             regValues = regValue)
  reg_Model_criterion <- mxOption(reg_Model_criterion, "Calculate Hessian", "No") # might cause errors; check
  reg_Model_criterion <- mxOption(reg_Model_criterion, "Standard Errors", "No") # might cause errors; check

  fit_reg_Model_criterion <- mxRun(reg_Model_criterion, silent = T)
  out <- list("best penalty" = minimum_criterion, "bestmodel" = fit_reg_Model_criterion, "fit measures" = results, "call" = call)
  return(out)
  }


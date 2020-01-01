computeFinalctParameters <- function(minimum_criterion, ctsemModelObject, alpha, gamma, regOn, regIndicators, link, dt, results){

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

  reg_CtModel_criterion <- regCtModel(ctsemModelObject = ctsemModelObject, alpha = alpha, gamma = gamma,
                                      regOn = regOn, regIndicators = regIndicators,
                                      regValues = regValue, link = link, dt = dt)
  reg_CtModel_criterion <- mxOption(reg_CtModel_criterion, "Calculate Hessian", "No") # might cause errors; check
  reg_CtModel_criterion <- mxOption(reg_CtModel_criterion, "Standard Errors", "No") # might cause errors; check

  fit_reg_CtModel_criterion <- mxRun(reg_CtModel_criterion, silent = T)
  out <- list("best penalty" = minimum_criterion, "bestmodel" = fit_reg_CtModel_criterion, "fit measures" = results, "call" = call)
  return(out)
}


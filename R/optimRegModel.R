#' optimRegModel
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#' optimRegModel creates a range of regularized models from an mxModel. It automatically tests different penalty values and returns the best model
#'
#' @param mxModelObject an already run mxModel
#' @param alpha alpha controls the type of penalty. For lasso regularization, set alpha = 1, for ridge alpha = 0. Values between 0 and 1 implement elastic net regularization
#' @param gamma gamma sets the power in the denominator of parameter specific weights when using adaptive lasso regularization. Make sure to set alpha to 1 when using a gamma other than 0.
#' @param regValue numeric value depicting the penalty size
#' @param regOn string vector with matrices that should be regularized. The matrices must have the same name as the ones provided in the mxModelObject (e.g., "A")
#' @param regIndicators list of matrices indicating which parameters to regularize in the matrices provided in regOn. The matrices in regIndicators must to have the same names as the matrices they correspond to (e.g., regIndicators = list("A" = diag(10))). 1 Indicates a parameter that will be regularized, 0 an unregularized parameter
#' @param regValue_start initial penalty value (recommended: 0)
#' @param regValue_end highest penalty value tested
#' @param regValue_stepsize increase in penValue between iterations
#' @param criterion Criterion for chosing the best final model. Possible are: AIC, BIC, CV.m2LL (only if manualCV is provided)
#' @param autoCV logical indicating if cross-validation should be performed automatically
#' @param k number of splits performed when autoCV = TRUE
#' @param Boot logical indicating if Bootstrap should be performed. Not yet implemented!
#' @param manualCV if cross-validation should be performed manually, provide a cross-validation sample (has to be of the same class as the data in the mxModelObject; e.g., mxData)
#' @param zeroThresh Threshold at which parameters are evaluated as being zero (necessary for AIC and BIC)
#' @param scaleCV indicate if the CV samples should be scaled automatically
#' @param cores Number of cores to use (>1 for parallel processing)
#'
#' @examples
#' # The following example is adapted from the regsem help to demonstrate the equivalence of both methods:
#'
#' library(lavaan)
#' library(OpenMx)
#' # put variables on same scale for regsem
#' HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
#'
#' # define variables:
#' latent = c("f1")
#' manifest = c("x1","x2","x3","x4","x5", "x6", "x7", "x8", "x9")
#'
#' # define paths:
#' loadings <- mxPath(from = latent, to = manifest, free = c(F,T,T,T,T,T,T,T,T), values = 1)
#' lcov <- mxPath(from = latent, arrows = 2, free = T, values = 1)
#' lmanif <- mxPath(from = manifest, arrows =2 , free =T, values = 1)
#'
#' # define model:
#' myModel <- mxModel(name = "myModel", latentVars = latent, manifestVars = manifest, type = "RAM",
#'                    mxData(observed = HS, type = "raw"), loadings, lcov, lmanif,
#'                    mxPath(from = "one", to = manifest, free = T)
#' )
#'
#' fit_myModel <- mxRun(myModel)
#'
#' # Show the names of the matrices in the model:
#' names(fit_myModel$matrices)
#'
#' # Penalize specific parameters from the A matrix (directional paths):
#' regOn <- c("A")
#'
#' selectedA <- matrix(0, ncol = ncol(fit_myModel$A$values), nrow = nrow(fit_myModel$A$values))
#' selectedA[c(2,3,7,8,9),10] <-1 # parameters that should be regularized have to be marked with 1
#' regIndicators <- list("A" = selectedA) # save in a list. Note the naming of the list element
#'
#' # Run the models:
#'
#' reg_model <- optimRegModel(mxModelObject = fit_myModel, alpha = 1, gamma = 0, regOn  = regOn, regIndicators = regIndicators)
#'
#' reg_model$`fit measures`
#'
#' reg_model$`best penalty`
#'
#' # Run the same model with 5-fold cross-validation
#'
#' CV_reg_model <- optimRegModel(mxModelObject = fit_myModel, alpha = 1, gamma = 0, regOn  = regOn, regIndicators = regIndicators,
#'                               autoCV = T, k = 5)
#' CV_reg_model$`CV results`
#'
#' @author Jannik Orzek
#' @import OpenMx
#' @export

optimRegModel <- function(mxModelObject, alpha = 1, gamma = 0, regOn, regIndicators,
                          regValue_start = 0, regValue_end = 1, regValue_stepsize = .01,
                          criterion = "BIC", autoCV = FALSE, k = 5, Boot = FALSE, manualCV = NULL, zeroThresh = .001, scaleCV = TRUE, cores = 1){

  if(cores == 1){

    ret <- SingleCoreOptimRegModel(mxModelObject = mxModelObject, alpha = alpha, gamma = gamma, regOn = regOn, regIndicators = regIndicators,
                                               regValue_start = regValue_start, regValue_end = regValue_end, regValue_stepsize = regValue_stepsize,
                                               criterion = criterion, autoCV = autoCV, k = k, Boot = Boot, manualCV = manualCV, zeroThresh = zeroThresh, scaleCV = scaleCV, cores = cores)
  }else{
    ret <- MultiCoreOptimRegModel(mxModelObject = mxModelObject, alpha = alpha, gamma = gamma, regOn = regOn, regIndicators = regIndicators,
                                   regValue_start = regValue_start, regValue_end = regValue_end, regValue_stepsize = regValue_stepsize,
                                   criterion = criterion, autoCV = autoCV, k = k, Boot = Boot, manualCV = manualCV, zeroThresh = zeroThresh, scaleCV = scaleCV, cores = cores)

  }

  return(ret)



}

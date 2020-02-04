#' optimRegCtModel
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#' optimRegCtModel creates a range of regularized models from a ctsem. It automatically tests different penalty values and returns the best model
#'
#' @param ctsemModelObject an already run ctsem object
#' @param alpha alpha controls the type of penalty. For lasso regularization, set alpha = 1, for ridge alpha = 0. Values between 0 and 1 implement elastic net regularization
#' @param gamma gamma sets the power in the denominator of parameter specific weights when using adaptive lasso regularization. Make sure to set alpha to 1 when using a gamma other than 0.
#' @param regValue numeric value depicting the penalty size
#' @param regOn string vector with matrices that should be regularized. The matrices must have the same name as the ones provided in the mxModelObject (e.g., "A")
#' @param regIndicators list of matrices indicating which parameters to regularize in the matrices provided in regOn. The matrices in regIndicators must to have the same names as the matrices they correspond to (e.g., regIndicators = list("A" = diag(10))). 1 Indicates a parameter that will be regularized, 0 an unregularized parameter
#' @param link list with functions that will be applied to the regularized matrices. For example, if the discrete time autoregressive and cross-lagged parameters should be regularized, use list("DRIFT" = "expm"). You can also define your own link function. DESCRIPTION NOT YET IMPLEMENTED
#' @param dt list with vectors depicting the time points for which the discrete time parameters should be regularized (e.g., list("DRIFT" = c(1)))
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
#' @param cores number of cores
#'
#'
#' @examples
#' library(ctsem)
#'
#' set.seed(175446)
#'
#' ## define the population model:
#'
#' # set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
#' ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)
#'
#' generatingModel<-ctModel(Tpoints=10,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                          MANIFESTVAR=diag(0,2),
#'                          LAMBDA=diag(1,2),
#'                          DRIFT=ct_drift,
#'                          DIFFUSION=matrix(c(.5,0,0,.5),2),
#'                          CINT=matrix(c(0,0),nrow=2),
#'                          T0MEANS=matrix(0,ncol=1,nrow=2),
#'                          T0VAR=diag(1,2))
#'
#' # simulate a training data set
#' traindata <- ctGenerate(generatingModel,n.subjects = 100)
#'
#' ## Build the analysis model. Note that drift eta1_eta2 is freely estimated
#' # although it is 0 in the population.
#'
#' myModel <- ctModel(Tpoints=10,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                    LAMBDA=diag(1,2),
#'                    MANIFESTVAR=diag(0,2),
#'                    CINT=matrix(c(0,0),nrow=2),
#'                    DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
#'                    T0MEANS=matrix(0,ncol=1,nrow=2),
#'                    T0VAR=diag(1,2))
#'
#' # fit the model using ctsem:
#' fit_myModel <- ctFit(traindata, myModel)
#' # discrete time parameters for deltaT = 1
#' expm(fit_myModel$mxobj$DRIFT$values)
#'
#' # select DRIFT values:
#' regOn = "DRIFT"
#' regIndicators = list("DRIFT" = matrix(c(0,1,0,0), byrow = T, ncol = 2))
#'
#' # regularize discrete time parameters for deltaT = 1
#' link = list("DRIFT" = "expm")
#' dt = list("DRIFT"= 1)
#'
#' # optimize regularized ctsem
#'
#' myRegCtModel <- optimRegCtModel(ctsemModelObject = fit_myModel, ralpha = 0, gamma = 1,
#'                                 regOn = regOn, regIndicators = regIndicators,
#'                                 link = link, dt = dt, autoCV = F, criterion = "BIC")
#'
#' # extract discrete time parameters for deltaT = 1
#' expm(myRegCtModel$bestmodel$Submodel$DRIFT$values)
#'
#' # optimize regularized ctsem with cross-validation
#'
#' CV_myRegCtModel <- optimRegCtModel(ctsemModelObject = fit_myModel, alpha = 0, gamma = 1,
#'                                    regOn = regOn, regIndicators = regIndicators,
#'                                    link = link, dt = dt, autoCV = T, k = 5, regValue_start = 0, regValue_end = .5, regValue_stepsize = .01)
#'
#' # extract discrete time parameters for deltaT = 1
#' expm(CV_myRegCtModel$`final Model`$Submodel$DRIFT$values)
#'
#' @author Jannik Orzek
#' @import OpenMx ctsem
#' @export

optimRegCtModel <- function(ctsemModelObject, alpha = 1, gamma = 0, regOn, regIndicators,
                            link = list("ident"), dt = NULL,
                            regValues,
                            criterion = "BIC", autoCV = FALSE, Boot = FALSE,
                            manualCV = NULL, k = 5, zeroThresh = .001, scaleCV = TRUE,
                            scaleFactors = scaleFactors, cores = 1){

  if(cores == 1){

    ret <- regmx::SingleCoreOptimRegCtModel(ctsemModelObject = ctsemModelObject, alpha = alpha, gamma = gamma, regOn=regOn, regIndicators = regIndicators,
                                     link = link , dt = dt,
                                     regValues = regValues,
                                     criterion = criterion, autoCV = autoCV, Boot = Boot, manualCV = manualCV, k = k, zeroThresh = zeroThresh, scaleCV = scaleCV, cores = cores)
  }else{
    ret <- regmx::MultiCoreOptimRegCtModel(ctsemModelObject = ctsemModelObject, alpha = alpha, gamma = gamma, regOn=regOn, regIndicators = regIndicators,
                                    link = link , dt = dt,
                                    regValues = regValues,
                                    criterion = criterion, autoCV = autoCV, Boot = Boot, manualCV = manualCV, k = k, zeroThresh = zeroThresh, scaleCV = scaleCV, cores = cores)

  }

  return(ret)



}

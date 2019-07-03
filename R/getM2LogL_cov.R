#' getM2LogL_cov
#'
#' Note: laremm is based on the R package \pkg{regsem}. Because of the early status of laremm, it is recommended to use regsem instead!
#'
#' This function computes the -2log(Likelihood) based on the formula reported in Bollen 1989, p. 133
#' @param ExpCov Expected covariance matrix
#' @param ObsCov Observed covariance matrix
#' @param NManif Number of manifest variables in the data set
#' @param NObs Sample size
#' @param useBiasedCov if TRUE, the biased covariance is used to compute the -2LogL
#' @author Jannik Orzek
#' @import OpenMx
#' @examples
#'
#' @export
#'
getM2LogL_cov <- function(ExpCov, ObsCov, NManif, NObs, useBiasedCov = T) {
  if(useBiasedCov == T){
    ObsCov = ObsCov*(NObs-1)/NObs

  }
  # Formula from Bollen 1989, p. 133
  LogL <- (-NObs/2)*((NManif) * log(2*pi) + log( det(ExpCov)) + tr(ObsCov%*%solve(ExpCov)))
  m2LogL <- -2*LogL
  return(m2LogL)
}

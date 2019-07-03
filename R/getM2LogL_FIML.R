#' getM2LogL_FIML
#'
#' Note: laremm is based on the R package \pkg{regsem}. Because of the early status of laremm, it is recommended to use regsem instead!
#'
#' This function computes the -2log(Likelihood) based on the formula reported in Hirose 2013 for full data
#' @param ExpCov Expected covariance matrix
#' @param ExpMean Expected means vector
#' @param ObsData observed raw data
#' @author Jannik Orzek
#' @import OpenMx
#' @examples
#'
#' @export
#'
getM2LogL_FIML <- function(ExpCov, ExpMean, ObsData) {
  # FIML based on https://raw.githubusercontent.com/OpenMx/OpenMx/master/docs/source/FIML_RowFit.rst

  m2LogL_comb <- 0
  NRow <- nrow(ObsData)
  for (i in 1:NRow){
    # k_rowi = 'number of non-missing observed variables in row *i* of the ObsData set' (https://raw.githubusercontent.com/OpenMx/OpenMx/master/docs/source/FIML_RowFit.rst)
    indice_isAN_rowi <- !is.na(ObsData[i,]) # vector indicating the missing values in row i
    k_rowi <- sum(!is.na(ObsData[i,])) #number of non-missing values in row i

    # ExpCov_rowi = '*filtered* model-implied manifest covariance matrix' (k_i x k_i matrix)
    ExpCov_rowi <- ExpCov[indice_isAN_rowi,indice_isAN_rowi]

    # X_rowi =  '*filtered* row *i* of the ObsData set' (1 x k_rowi Matrix)
    X_rowi <- as.numeric(ObsData[i,indice_isAN_rowi])

    # ExpMean_rowi = '*filtered* model-implied manifest means row vector' (1 x k_rowi Matrix)
    ExpMean_rowi <- ExpMean[indice_isAN_rowi]

    m2LogL_rowi <- k_rowi*log(2*pi) + log( det(ExpCov_rowi) ) + t(X_rowi - ExpMean_rowi) %*% solve(ExpCov_rowi) %*% (X_rowi - ExpMean_rowi) # Formel from Hiros 2013_FIML Formel

    m2LogL_comb <- m2LogL_comb + m2LogL_rowi
  }
  return(m2LogL_comb)
}

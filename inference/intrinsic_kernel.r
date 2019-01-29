
#' The 'Intrinsic_kernel' class object
#'
#' This a abstract class provide the estimated heat kernel matrices
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return an \code{\link{R6Class}} object which can be used for intrinsic GP regression
#' @format \code{\link{R6Class}} object.
#' @field cov_list vector containing the covariance matrices.  
#' @field cov_list_uf vector containing the covariance matrices for sparse GP. 
#' @export
#'
#' @author Mu Niu, \email{mu.niu@plymouth.ac.uk}

Intrinsic_kernel<-R6Class("Intrinsic_kernel",
  public = list(
    cov_list=NULL,
    cov_list_uf = NULL,

    initialize = function(cov_list = NULL, cov_list_uf = NULL) {
      self$cov_list <- cov_list
      self$cov_list_uf<-cov_list_uf
      self$greet()
    },
    
    greet = function() {
    }

  )
)





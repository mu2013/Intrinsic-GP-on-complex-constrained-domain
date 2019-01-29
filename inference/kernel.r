#' The 'Kernel' class object
#'
#' This a abstract class     provide the kernel function and the 1st order derivative of rbf kernel function. 
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return an  \code{\link{R6Class}} object which can be used for GP regression.
#' @format \code{\link{R6Class}} object.
#' @field k_par vector(of length n_hy) containing the hyper-parameter of kernel. n_hy is the length of kernel hyper parameters.
#' @section Methods:
#' \describe{
#'   \item{\code{kern(t1,t2)}}{This method is used to calculate the kernel function given two real inputs.}   
#'   \item{\code{dkd_kpar(t1,t2)}}{This method is used to calculate the gradient of kernel function against the kernel hyper parameters given two real inputs.}   
#'   \item{\code{dkdt(t1,t2)}}{This method is used to calculate the 1st order derivative of kernel function given two real inputs.} }
#' @export
#'
#' @author Mu Niu, \email{mu.niu@plymouth.ac.uk}

Kernel<-R6Class("Kernel",
  public = list(
    k_par=NULL,

    initialize = function(k_par = NULL) {
      self$k_par <- k_par
      self$greet()
    },
    
    greet = function() {
    },

    kern = function (t1,t2) {
    },

    dkd_kpar = function(t1,t2) {
    },

  dkdt = function (t1,t2) {  
    }

  )
)


#' The 'Heat_manif' class object
#'
#' This a R6 class. It inherits from 'Intrinsic_kernel' class. It provides the heat kernel cov list. 
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return an  \code{\link{R6Class}} object which can be used for GP regression.
#' @format \code{\link{R6Class}} object.
#' @export
#' @author Mu Niu, \email{mu.niu@plymouth.ac.uk}

Heat_manif <- R6Class("Heat_manif",
  inherit = Intrinsic_kernel,
  public = list(

    greet = function() {
      cat(paste0("length of cov is ", length(self$cov_list), ".\n"))
    },

    set_cov_list = function(val) {
      self$cov_list <- val
    }

  )
)

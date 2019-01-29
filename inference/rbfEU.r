#' The 'RBF' class object
#'
#' This a R6 class. It inherits from 'kernel' class. It provides the rbf kernel function and the 1st order derivative of rbf kernel function. 
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return an  \code{\link{R6Class}} object which can be used for GP regression.
#' @format \code{\link{R6Class}} object.
#' @export
#' @author Mu Niu, \email{mu.niu@plymouth.ac.uk}

RBF <- R6Class("RBF",
  inherit = Kernel,
  public = list(

    greet = function() {
      cat(paste0("RBF hyper par is", self$k_par, ".\n"))
    },

    set_k_par = function(val) {
      self$k_par <- val
    },

    kern = function (t1,t2) {
      x = sqrt(  sum( (t1-t2)^2 )  )
      self$k_par[1]*exp( -x^2/(2*self$k_par[2]) )
      #sig*exp( -r/(2*l) )
    },

    dkd_kpar = function(t1,t2) {
      x = sqrt(  sum( (t1-t2)^2 )  ) #x=t1-t2
      dkdl = self$k_par[1]*exp( -x^2/(2*self$k_par[2]) ) * x^2 / (2*self$k_par[2]^2)*self$k_par[2]
      dkdsig = self$k_par[1]*exp( -x^2/(2*self$k_par[2]) )
      c(dkdsig,dkdl)
    },

   dkdx = function (t1,t2) 
   {
     t1 = as.matrix(t1); t2 = as.matrix(t2)
     x = sqrt(  sum( (t1-t2)^2 )  ) #x=(t1-t2)
     nx = dim(t1)[1]
     grad=c()
     #if( nx>1 )
     #{
       for(i in 1:nx)
        {
        grad =cbind( grad, self$k_par[1]*exp( -x^2/(2*self$k_par[2]) )*(-(t1[i]-t2[i]))/self$k_par[2] )
        }
      #}else{
      # grad = self$k_par[1]*exp( -x^2/(2*self$k_par[2]) ) * (-x)/self$k_par[2]
      #}
     grad
    }

  )
)


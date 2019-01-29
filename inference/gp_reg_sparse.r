#' The 'sparse GP regression' class object
#'
#' This class provide the sparse GP regression model. 
#' @docType class
#' @importFrom R6 R6Class
#' @import pracma
#' @export
#' @keywords data
#' @return an  \code{\link{R6Class}} object which can be used for sparse GP regression.
#' @format \code{\link{R6Class}} object.
#' @field y matrix containing observation.
#' @field t vector containing time points for observation.
#' @field t vector containing inducing points.
#' @field ker kernel class object containing kernel.
#' @field nsig containing noise variance.
#' @section Methods:
#' \describe{
#'   \item{\code{upredict()}}{This method is used to make prediction on given time points}   
#'   \item{\code{optimiseM()}}{This method is used to optimise kernel hyper parameters.} }
#' @export
#' @examples
#' @author Mu Niu, \email{mu.niu@plymouth.ac.uk}
gp_reg_sparse <- R6Class("gp_reg_sparse",
  public = list(
    y = NULL,
    t = NULL,
    nsig=NULL,
    ker= NULL,
    jitter=NULL,
    mean_f = NULL,
    u= NULL,
    #initialize = function(y = NULL, t=NULL,b=NULL,lambda=NULL,ker=NULL) {
    initialize = function(y = NULL,u=NULL,t=NULL,nsig=NULL,ker=NULL,jitter=NULL,mean_f=NULL) 
    {
      self$y = y
      self$t = t
      self$u = u
      self$nsig = nsig
      self$ker = ker
      self$jitter = jitter
      if( missing(jitter) )
	    {  
          self$jitter = 1e-8
	    	} else{
	        self$jitter = jitter	
	    	}

      if( missing(mean_f) )
	    {  
          self$mean_f = 'off'
	    	} else{
	      self$mean_f = mean_f
	        }

      self$greet()
    },

    greet = function () {
      cat(paste0("ker is", self$ker$greet(), ".\n"))
    },

    k_m=function(X1,X2) 
    {
      nx1 = dim(X1)[2] ; nx2 = dim(X2)[2]
	  K<- matrix( rep(0, nx1*nx2 ), nrow= nx1 )
	    for (i in 1:nx1) 
	    {
	      for (j in 1:nx2) 
	      {
	        K[i,j] <- self$ker$kern(X1[,i], X2[,j])
	      }
	    }
	  return(K)
	},

   ulik=function(par,X1,y_n,u,jitter )  #par,X1,y_me,u,jitter
   {
   	y_n = matrix(y_n,nrow=1)
   	self$ker$k_par = exp(par[1:length(self$ker$k_par)])
   	self$nsig = exp(par[3])
	
	X2=X1
	n= max(dim(X1))
	nu=max(dim(u))

	Kfu = self$k_m(X1,u)
	Kuu = self$k_m(u,u)+jitter*diag(nu)
	kinv= chol2inv( chol(Kuu) )
	Ka= Kfu%*%kinv%*%t(Kfu) + self$nsig*diag(n)
    ## using woodbury matrix identity  (A + UCV)^{-1} = A^{-1} - A^{-1}U(C^{-1} + VA^-1 U)^-1 VA^-1
	## inverse of Ka (k_nm K_mm^-1 K_nm + Sigma)^-1  ## https://en.wikipedia.org/wiki/Woodbury_matrix_identity
	ins = 1/self$nsig*diag(n)
	kmed = Kuu +t(Kfu)%*%ins%*%Kfu
    
    ## use choleski for inversion
	iii = chol2inv( chol(kmed) )
	kainv = ins - ins%*%Kfu%*% iii %*% t(Kfu) %*% ins

	lnlik= -1/2*y_n%*%kainv%*%t(y_n) - 1/2*log(det(Ka))-n/2*log(2*pi)
	return(-lnlik)
   },

  optimisM=function( init,jitter) 
	{
	 if( missing(init) )
	    {  
          init = log( c( 1.5,1.5,1.5 ) )
	    	} else{
	        init = 	init	
	    	}
	 if( missing(jitter) )
	    {  
          jitter = 1e-5
	    	} else{
	      jitter = 	jitter	
	    	}

	 if(self$mean_f == 'off')
	 {
     scaled_y= self$y
     } else {
     scaled_y= self$y - mean(self$y )	
     }
     
     res1<-optim(c(init),self$ulik,,self$t,scaled_y,self$u,jitter,method="BFGS")
     self$ker$k_par = exp(res1$par)
	 return(res1$par )
	},

	upredict=function(X1,X2)
	{
	 if(self$mean_f == 'off')
	 {
     scaled_y= self$y
     } else {
     scaled_y= self$y - mean(self$y )	
     }

     k.ux<-self$k_m(self$u,X2)
     k.uf<-self$k_m(self$u,X1)
     k.uu<-self$k_m(self$u,self$u)
     kinv = solve(k.uu)
  
     ks = 1/self$nsig*k.uf%*%t(k.uf)+k.uu
     iks = chol2inv( chol(ks) )

     f.bar.star <- 1/self$nsig*t(k.ux)%*%iks%*%k.uf%*%t(scaled_y)
     cov.f.star <- t(k.ux)%*%iks%*%(k.ux)
     return(list("pre_m"=f.bar.star,"pre_v"=cov.f.star) )
	}

   )

)

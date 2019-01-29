#' The 'intrinsic sparse GP regression' class object
#'
#' This class provide the intrinsic sparse GP regression methods on manifold using heat kernel. 
#' @docType class
#' @importFrom R6 R6Class
#' @import pracma
#' @export
#' @keywords data
#' @return an  \code{\link{R6Class}} object which can be used for doing GP regression.
#' @format \code{\link{R6Class}} object.
#' @field y matrix containing observation.
#' @field t vector containing time points for observation.
#' @field u vector containing the inducing points.
#' @field Heat_manif contaning the heat kernel object of the manifold.
#' @section Methods:
#' \describe{
#'   \item{\code{predict()}}{This method is used to make prediction on given time points}  
#'   \item{\code{optimisMn()}}{This method is used to find the optimum cov matrices.} }
#' @export
#' @examples
#' @author Mu Niu, \email{mu.niu@plymouth.ac.uk}



gp_reg_manif_sparse <- R6Class("gp_reg_manif_sparse",
  public = list(
    y = NULL,
    u = NULL,
    x = NULL,
    nsig=NULL,
    Heat_manif= NULL,
    jitter=NULL,
    mean_f = NULL,
    initialize = function(y = NULL,u=NULL,x=NULL,nsig=NULL,Heat_manif=NULL,jitter=NULL,mean_f=NULL) 
    {
      self$y = y
      self$u = u
      self$x = x
      self$nsig = nsig
      self$Heat_manif = Heat_manif
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
      cat(paste0("length of cov is ", length(self$Heat_manif$cov_list), ".\n"))
    },

    lik_sparse = function(par,y_n,index,jitter)
    {
     if( missing(jitter) )
	    { jitter = 1e-5} else{
	    	jitter = jitter
	    }  

      y_n=matrix(y_n,nrow=1)
      sig= exp(par[1])
      nsig=exp(par[2])
      ksig = sig
      
      nu=max(dim(self$u))
      n= max(dim(self$x))

      mK<-self$Heat_manif$cov_list[[index]]
      eee<-eigen(mK, symmetric = TRUE)
	  lam=eee$values
	  if (lam[nu]<0)
	  {
		vvv=eee$vectors
		rlam=lam
		rlam[which(lam<0)]=lam[1]/1e+10
		rK = vvv%*%diag(rlam)%*%solve(vvv) 
	   } else {
		  rK=mK
	   }

      Kfu = ksig*t( self$Heat_manif$cov_list_uf[[index]] ) #rbf_k(X1,u,l,sig)#[,,1]

	  Kuu = ksig*rK +jitter*diag(nu)

	  kinv= chol2inv( chol(Kuu) ) #solve(Kuu)

	  Ka= Kfu%*%kinv%*%t(Kfu) + nsig*diag(n)
	  ## using woodbury matrix identity  (A + UCV)^{-1} = A^{-1} - A^{-1}U(C^{-1} + VA^-1 U)^-1 VA^-1
	  ## inverse of Ka (k_nm K_mm^-1 K_nm + Sigma)^-1  ## https://en.wikipedia.org/wiki/Woodbury_matrix_identity
	  ins = 1/nsig*diag(n)
	  kmed = Kuu +t(Kfu)%*%ins%*%Kfu

      ## use choleski for inversion
      iii = chol2inv( chol(kmed) )
      kainv = ins - ins%*%Kfu%*% iii %*% t(Kfu) %*% ins

      lnlik= -1/2*y_n%*%kainv%*%t(y_n) - 1/2*log(det(Ka))-n/2*log(2*pi)
      return(-lnlik)
    },

  optimisMn=function(len_range,ini_par)
  {
  lnlike = rep( -Inf, max(len_range) )
  opsig = c(0)
  opnsig = c(0)

  if(self$mean_f == 'off')
   {
     scaled_y= self$y  
     } else {
     scaled_y= self$y - mean(self$y ) 
     }
  
  for(i in len_range)
  {

  if(missing(ini_par)){
      iisig= 80
      nsig = 2
    }else{
      iisig = ini_par[1]
      nsig = ini_par[2]
    }

  par<- c( log(iisig), log(nsig) )
  resall<-optim(par,self$lik_sparse,,scaled_y,i,method="BFGS")
  lnlike[i] = -resall$value
  opsig[i] = exp(resall$par[1])
  opnsig[i] = exp(resall$par[2])
  }

  lindex=which(lnlike==max(lnlike)) 
  iisig = opsig[lindex] 
  nsig = opnsig[lindex]
  return(list('op_index'=lindex,'op_ksig'=iisig,'op_nsig'=nsig,'lnlike'=lnlike))
 },

 predict=function(op_par,cov_list_ux)
  {
  if(self$mean_f == 'off')
   {
     scaled_y= self$y  
     } else {
     scaled_y= self$y - mean(self$y ) 
     }

   nsig = op_par$op_nsig
   ksig = op_par$op_ksig

   k.uf<-ksig*self$Heat_manif$cov_list_uf[[op_par$op_index]]
   k.uu<-ksig*self$Heat_manif$cov_list[[op_par$op_index]]
   k.ux<-ksig*cov_list_ux

   kinv = solve(k.uu)

   ks = 1/nsig*k.uf%*%t(k.uf)+k.uu
   iks = chol2inv( chol(ks) )

   f.bar.star <- 1/nsig*(k.ux)%*%iks%*%k.uf%*%t(scaled_y)
   cov.f.star <- k.ux%*%iks%*%t(k.ux)

  return(list("pre_m"=f.bar.star,"pre_v"=cov.f.star))
  }

   )
 )
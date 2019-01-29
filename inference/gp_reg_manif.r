#' The 'intrinsic GP regression' class object
#'
#' This class provide the intrinsic GP regression methods on manifold using heat kernel. 
#' @docType class
#' @importFrom R6 R6Class
#' @import pracma
#' @export
#' @keywords data
#' @return an  \code{\link{R6Class}} object which can be used for doing intrinsic GP regression.
#' @format \code{\link{R6Class}} object.
#' @field y matrix containing observation.
#' @field t vector containing time points for observation.
#' @field Heat_manif contaning the heat kernel object of the manifold.
#' @section Methods:
#' \describe{
#'   \item{\code{predictM()}}{This method is used to make prediction on given time points} 	
#'   \item{\code{optimisM()}}{This method is used to find the optimum cov matrices.} }
#' @export
#' @examples
#' @author Mu Niu, \email{mu.niu@plymouth.ac.uk}

gp_reg_manif <- R6Class("gp_reg_manif",
  public = list(
    y = NULL,
    nsig=NULL,
    Heat_manif= NULL,
    jitter=NULL,
    mean_f = NULL,

    initialize = function(y = NULL,nsig=NULL,Heat_manif=NULL,jitter=NULL,mean_f=NULL) 
    {
      self$y = y
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

lik=function (par,y_n,index)
{
  y_n=matrix(y_n,nrow=1)
  sig= exp(par)
  ksig = sig#*2*pi*l #/( 1/sqrt(2*pi*l) ) 
  dn=length(y_n)

  mK<-self$Heat_manif$cov_list[[index]]

  eee<-eigen(mK, symmetric = TRUE)
  lam=eee$values

  if (lam[dn]<0)
  {
  vvv=eee$vectors
  rlam=lam
  rlam[which(lam<0)]=lam[1]/1e+10
  rK = vvv%*%diag(rlam)%*%solve(vvv) 
  } else {
    rK=mK
  }

  K=ksig*rK#+nsig*diag(dn) 
  kinv=solve(K)
  lnlik = -1/2*y_n%*%kinv%*%t(y_n) - 1/2*log(det(K))-dn/2*log(2*pi)

  return(-lnlik)
},

optimisM=function(len_range)
{
  lnlike = rep( -Inf, max(len_range) )
  opsig = c(0)

  if(self$mean_f == 'off')
   {
     scaled_y= self$y  
     } else {
     scaled_y= self$y - mean(self$y ) 
     }

  for(i in len_range)
  {
  iisig=1
  #nsig= 1e-10
  par<-c(log(iisig))
  resall<-optim(par,self$lik,,scaled_y,i,method="L-BFGS-B")
  lnlike[i] = -resall$value
  opsig[i] = exp(resall$par)
  }

  lindex=which(lnlike==max(lnlike)) 
  iisig = opsig[lindex] 
 return(list('op_index'=lindex,'op_ksig'=iisig,'lnlike'=lnlike))
},

likn=function (par,y_n,index)
{
  y_n=matrix(y_n,nrow=1)
  sig= exp(par[1])
  nsig=exp(par[2])
  ksig = sig#*2*pi*l #/( 1/sqrt(2*pi*l) ) 
  dn=length(y_n)

  mK<-self$Heat_manif$cov_list[[index]]

  eee<-eigen(mK, symmetric = TRUE)
  lam=eee$values

  if (lam[dn]<0)
  {
  vvv=eee$vectors
  rlam=lam
  rlam[which(lam<0)]=lam[1]/1e+10
  rK = vvv%*%diag(rlam)%*%solve(vvv) 
  } else {
    rK=mK
  }

  K=ksig*rK+nsig*diag(dn) 
  kinv=solve(K)
  lnlik = -1/2*y_n%*%kinv%*%t(y_n) - 1/2*log(det(K))-dn/2*log(2*pi)

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
      iisig= 1
      nsig = 5e-2^2
    }else{
      iisig = ini_par[1]
      nsig = ini_par[2]
    }

  #par<-c(log(iisig))
  #resall<-optim(par,self$lik,,scaled_y,i,method="L-BFGS-B")
  par<- c( log(iisig), log(nsig) )
  resall<-optim(par,self$likn,,scaled_y,i,method="L-BFGS-B")
  lnlike[i] = -resall$value
  opsig[i] = exp(resall$par[1])
  opnsig[i] = exp(resall$par[2])
  }

  lindex=which(lnlike==max(lnlike)) 
  iisig = opsig[lindex] 
  nsig = opnsig[lindex]
 return(list('op_index'=lindex,'op_ksig'=iisig,'op_nsig'=nsig,'lnlike'=lnlike))
},

predictM = function(op_par,k.xsx)
{
  k.xx = self$Heat_manif$cov_list[[op_par$op_index]]
  ksig = op_par$op_ksig
  if(self$mean_f == 'off')
   {
     scaled_y= self$y  
     } else {
     scaled_y= self$y - mean(self$y ) 
     }

#f.bar.star <- ksig*k.xsx%*%solve(ksig*k.xx + nsig*diag(1, ncol(k.xx)))%*%matrix(y_me, ncol=1)
  f.bar.star <- ksig*k.xsx%*%solve(ksig*k.xx)%*%matrix(scaled_y, ncol=1)

  return(list("pre_m"=f.bar.star) )
},

predictMn = function(op_par,k.xsx)
{
  k.xx = self$Heat_manif$cov_list[[op_par$op_index]]
  ksig = op_par$op_ksig
  nsig = op_par$op_nsig
  if(self$mean_f == 'off')
   {
     scaled_y= self$y  
     } else {
     scaled_y= self$y - mean(self$y ) 
     }

  f.bar.star <- ksig*k.xsx%*%solve(ksig*k.xx + nsig*diag(1, ncol(k.xx)))%*%matrix(scaled_y, ncol=1)
 #f.bar.star <- ksig*k.xsx%*%solve(ksig*k.xx)%*%matrix(scaled_y, ncol=1)

  return(list("pre_m"=f.bar.star) )
}

  )

)






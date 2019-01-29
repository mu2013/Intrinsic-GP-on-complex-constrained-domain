#' The 'GP regression' class object
#'
#' This class provide the GP regression model. 
#' @docType class
#' @importFrom R6 R6Class
#' @import pracma
#' @export
#' @keywords data
#' @return an  \code{\link{R6Class}} object which can be used for GP regression.
#' @format \code{\link{R6Class}} object.
#' @field y matrix containing observation.
#' @field t vector containing time points for observation.
#' @field ker kernel class object containing kernel.
#' @field nsig containing noise variance.
#' @section Methods:
#' \describe{
#'   \item{\code{predict()}}{This method is used to make prediction on given time points} 	
#'   \item{\code{optimiseM()}}{This method is used to optimise kernel hyper parameters.} }
#' @export
#' @examples
#' @author Mu Niu, \email{mu.niu@plymouth.ac.uk}

gp_reg <- R6Class("gp_reg",
  public = list(
  	
    y = NULL,
    t = NULL,
    nsig=NULL,
    ker= NULL,
    jitter=NULL,
    mean_f = NULL,

    initialize = function(y = NULL, t=NULL,nsig=NULL,ker=NULL,jitter=NULL,mean_f=NULL) 
    {
      self$y = y
      self$t = t
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
    
    lik=function (par,X1,y_n)  
    {
		self$ker$k_par = exp(par[1:length(self$ker$k_par)])
		X2=X1
		n = dim(y_n)[2]

		K1 = self$k_m(X1,X2) + self$nsig*diag(n)

        jitter = 1e-8
        varY = var( c(y_n)) #ij = 1

        repeat
         {
           K = K1+jitter*diag(n)

			invK = tryCatch({ chol2inv( chol(K,LINPACK =TRUE) ) 
	        }, warning = function(war)
	        { 
	          print(paste("log det: ",war))
	        },
              error = function(err) 
	        {# error handler picks up where error was generated
	          return( invK=diag(n) )
	        },finally = { } )
	        
	        jitter=jitter*10
	        #ij=ij+1
	        
	        if( !identical(invK,diag(n) ) )
		        {
		           break
		        } else {
		           if( jitter > varY )
		           {
		             break
		           }
		        }
          }

        self$jitter = jitter
		kinv = chol2inv( chol(K) )
		lnlik = -1/2*y_n%*%kinv%*%t(y_n) - 1/2*log(det(K))-n/2*log(2*pi)
		return(-lnlik)
	},

	likn=function (par,X1,y_n)  ## test
	{
		## l=exp(par[2])    # sig=exp(par[1])   # self$nsig=0.001
		self$nsig = exp(par[3])
		self$ker$k_par = exp(par[1:length(self$ker$k_par)])
		X2=X1
		n= dim(X1)[2]

        jitter = 1e-8
        varY = var( c(y_n) )

        repeat
         {
			K = self$k_m(X1,X2)+jitter*diag(n)

			invK = tryCatch({ chol2inv( chol(K, LINPACK =TRUE) ) 
	        }, warning = function(war)
	        { 
	          print(paste("log det: ",war))
	        },
	         error = function(err) 
	        {# error handler picks up where error was generated
	          return( invK=diag(n) )
	        },finally = { } )

	        jitter=jitter*10
	        self$jitter = jitter
	        #ij=ij+1
	        if( !identical(invK,diag(n) ) )
		        {
		           break
		        } else {
		        if( jitter > varY )
		           {
		             break
		           }
		        }
         }

		K = as.matrix( K + self$nsig*diag(n) ) 

		uK = chol(K)
		logdet = 2*sum( log( diag(uK) )  )
        
		#kinv = chol2inv( chol(K) )
	    #lnlik = -1/2*y_n%*%kinv%*%t(y_n) - 1/2*log( det(K) )-n/2*log(2*pi)

        pK = solve(K,t(y_n))
        lnlik = -1/2*y_n%*%pK - 1/2*logdet-n/2*log(2*pi)

		return(-lnlik)
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

	dk_m=function(X1,X2) 
	{
	  nx1 = dim(X1)[2] ; nx2 = dim(X2)[2]; xdim = dim(X1)[1]
      dK<- array(0,dim=c(nx1,nx2,xdim))
	  #dK<- matrix( rep(0, nx1*nx2 ), nrow= nx1 )
	    for (i in 1:nrow(dK)) 
	    {
	      for (j in 1:ncol(dK)) 
	      {
            dkm  <- self$ker$dkdx(X1[,i], X2[,j])
            for(d in 1:xdim)
            {
	          dK[i,j,d] <-dkm[d]
            }
	      }
	    }
	 return(dK)
	},

	dk_dTheta=function(X1,X2) 
	{
	  nx1 = dim(X1)[2] ; nx2 = dim(X2)[2]
	  dKdl<- matrix( rep(0, nx1*nx2 ), nrow= nx1 )
	  dKds<-dKdl
	    for (i in 1:nx1) 
	    {
	      for (j in 1:nx2) 
	      {
	      	dkdpar<- self$ker$dkd_kpar(X1[,i], X2[,j])
	        dKdl[i,j] <- dkdpar[1]
	        dKds[i,j] <- dkdpar[2]
	      }
	    }
	  return(list("dKdl"=dKdl,"dKds"=dKds))
	},

	predict=function(X2)
	{   
		X1 = self$t

		k.xsxs1 <- self$k_m(X2,X2)
        k.xsxs <- k.xsxs1+self$nsig*diag(1,ncol(k.xsxs1))

		k.xsx<- self$k_m(X2,X1)
		k.xx<- self$k_m(X1,X1)
		K = k.xx + (self$nsig+self$jitter)*diag(1, ncol(k.xx)) #solve(k.xx + self$nsig*diag(1, ncol(k.xx)))
		kinv = chol2inv( chol(K, LINPACK =TRUE) )
		
		if(self$mean_f == 'off')
		 {
	     scaled_y= self$y 
	     f.bar.star <- k.xsx%*%kinv%*%t(scaled_y)
	     } else {
	     scaled_y= self$y - mean(self$y )	
	     f.bar.star <- k.xsx%*%kinv%*%t(scaled_y) + mean(self$y )
	     }

		cov.f.star <- k.xsxs - k.xsx%*%kinv%*%t(k.xsx) #k.xxs

		dk.xsx <-  self$dk_m(X2,X1)   ##??? X2,X1
		xdim = dim(X2)[1]
		xnum = dim(X2)[2]
		f.dm = array(0, dim=c(xnum,1,xdim) );  f.dvar = array(0, dim=c(xnum,xnum,xdim) )

		for (dm in 1:xdim )
		{
		#f.dm <- dk.xsx%*%kinv%*%y_m
		  if(xnum == 1){
			dkmxsx<- as.matrix(dk.xsx[,,dm])
			 }else { 
			dkmxsx<- t( as.matrix(dk.xsx[,,dm])	)
			}
		 f.dm[,,dm] <- t(dkmxsx)%*%kinv%*%t(scaled_y)
		 f.dvar[,,dm] <- diag( -( t(dkmxsx)%*%kinv%*%t(k.xsx) + k.xsx%*%kinv%*%dkmxsx ) )
        }

	   return(list("pre_m"=f.bar.star,"pre_v"=cov.f.star,"pre_dm"=f.dm,"pre_dv"=f.dvar) )
	},

	minmax = function(val)
	 {
       if(val == 'min'){
         mpre = min(self$predict(self$t)$pre_m)
       	}else if(val=='max'){

         mpre = max(self$predict(self$t)$pre_m)
       	}
      mpre #self$t[mpre]
	 },

	num.grad.lik = function(par,X1,y)
	{
	  require(numDeriv)
      return( grad(self$lik,par,,,X1=X1,y_n=y) )
	},

    num.grad.likn = function(par,X1,y)
	{
	  require(numDeriv)
      return( grad(self$likn ,par,,,X1=X1,y_n=y) )
	},

    grad.hyp = function(par,X1,y_n)
    {
	     y_n=matrix(y_n,nrow=1)
		 self$ker$k_par = exp(par[1:length(self$ker$k_par)])
		 X2=X1
		 n = dim(y_n)[2]
		 K = self$k_m(X1,X2) + self$nsig*diag(n)
		 kinv = chol2inv( chol(K) )
	     dKdTheta = self$dk_dTheta(X1,X1) 
	     dLdl = 1/2*y_n%*%kinv%*%dKdTheta[[1]]%*%kinv%*%t(y_n) - 1/2*sum( diag(kinv%*%dKdTheta[[1]]) )
	     dLds = 1/2*y_n%*%kinv%*%dKdTheta[[2]]%*%kinv%*%t(y_n) - 1/2*sum( diag(kinv%*%dKdTheta[[2]]) )
	    # return negative loglikelihood
	    c(-dLdl,-dLds)
	},

    grad.hypn = function(par,X1,y_n)
    {
    	 self$nsig = exp(par[3])
	     y_n=matrix(y_n,nrow=1)
		 self$ker$k_par = exp(par[1:length(self$ker$k_par)])
		 X2=X1
		 n = dim(y_n)[2]

		 jitter = 1e-8
         varY = var( c(y_n))

        repeat
         {
			K = self$k_m(X1,X2)+jitter*diag(n)

			invK = tryCatch({ chol2inv( chol(K, LINPACK =TRUE) ) 
	        }, warning = function(war)
	        { 
	          print(paste("log det: ",war))
	        },
	         error = function(err) 
	        {# error handler picks up where error was generated
	          return( invK=diag(n) )
	        },finally = { } )

	        jitter=jitter*10
	        #ij=ij+1
	        
	        if( !identical(invK,diag(n) ) )
		        {
		           break
		        } else {
		        if( jitter > varY )
		           {
		             break
		           }
		        }
         }

		 K = as.matrix( K + self$nsig*diag(n) ) 
		 #K = self$k_m(X1,X2) + self$nsig*diag(n)

		 kinv = chol2inv( chol(K) )
	     dKdTheta = self$dk_dTheta(X1,X1)
	     dKdn =  self$nsig*diag(n)

	     dLdl = 1/2*y_n%*%kinv%*%dKdTheta[[1]]%*%kinv%*%t(y_n) - 1/2*sum( diag(kinv%*%dKdTheta[[1]]) )
	     dLds = 1/2*y_n%*%kinv%*%dKdTheta[[2]]%*%kinv%*%t(y_n) - 1/2*sum( diag(kinv%*%dKdTheta[[2]]) )
	     dLdn = 1/2*y_n%*%kinv%*%dKdn%*%kinv%*%t(y_n) - 1/2*sum( diag( kinv%*%dKdn ) )
	    # return negative loglikelihood
	    c(-dLdl,-dLds,-dLdn)
	},

	optimisM = function( init) 
	{
	 if( missing(init) )
	    {  
          init = log( c( 1.5,1.5 ) )
	    	} else{
	        init = 	init	
	    	}

	 if(self$mean_f == 'off')
	 {
     scaled_y= self$y
     } else {
     scaled_y= self$y - mean(self$y )	
     }
     
     res1<-optim(log(c(init)),self$lik,self$grad.hyp,self$t,scaled_y,method="BFGS")
     self$ker$k_par = exp(res1$par)
	 return(self$ker$k_par )
	},

	optimisMn = function( init,bounded ) 
	{
	 tbd = array(c(1), c(length(self$ker$k_par)) )
	 if( missing(bounded) )
	    {  
         bounded= c(-Inf,Inf)
         lbound = bounded[1]*tbd
         ubound = bounded[2]*tbd
	    } else {
	    lbound = log(bounded[1]*tbd)
	    ubound = log(bounded[2]*tbd)
	    }
	 if( missing(init) )
	    {  
          init = c( 1.5,1.5,1e-2 )
	    	} else{
	        init = 	init	
	    	}
	 
	 if(self$mean_f == 'off')
	 {
     scaled_y= self$y  
     } else {
     scaled_y= self$y - mean(self$y )	
     }
     #res1<-optim(log(c(init)),self$likn,,self$t,scaled_y,method="BFGS")
     res1<-optim(log(c(init)),self$likn,self$grad.hypn,self$t,scaled_y,method="BFGS")
     self$ker$k_par = exp(res1$par)
	 return(self$ker$k_par )
	},

    optimisMn_l = function( init ) 
	{
	 tbd = array(c(1), c(length(self$ker$k_par)) )

	 if( missing(init) )
	    {  
          init = c( 1.5,1.5,1e-2 )
	    	} else{
	        init = 	init	
	    	}
	 
	 if(self$mean_f == 'off')
	 {
     scaled_y= self$y  
     } else {
     scaled_y= self$y - mean(self$y )	
     }
     #res1<-optim(log(c(init)),self$likn,,self$t,scaled_y,method="BFGS")
     res1<-optim(log(c(init)),self$likn,self$grad.hypn,self$t,scaled_y,method="L-BFGS-B")
     self$ker$k_par = exp(res1$par)
	 return(self$ker$k_par )
	}


  )

)



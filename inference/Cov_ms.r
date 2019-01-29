
## the functions in this file are used to combine rows of heat kernel estimation and produce covariance matrices.

Cov_ms<- function(data_length,directory,T_len)
{
rowlist<-list()  ## store many cov matrices for each lengthscale 
    load(file=paste(directory ,1,'Krow.RData',sep=""))
    ppp= sKrow  ## number of lengthscale * length of data points (length of row of cov matrix)

  for (ilen in 1:T_len) ### startlen:endlen  loop through possible lengthscale  800 in this case 
    {
     rowlist<-c(rowlist,list(sKrow[ilen,]) )
     }

  for( ipoint in 2:data_length )
  {
     load(file=paste(directory ,ipoint,'Krow.RData',sep=""))
     ppp= sKrow  ## number of lengthscale * length of data points (length of row of cov matrix)

   ###  go through all lengthscale and combined row vectors.
      for(ilenr in 1:T_len )
      {
       rowlist[[ilenr]] = rbind(rowlist[[ilenr]],ppp[ilenr,])
      }
  }

## use lower off diagonal only
rowlistl = rowlist
      for(ilenr in 1:T_len )
      {
      rowlistl[[ilenr]] = rowlist[[ilenr]][1:data_length,1:data_length] 
      rowlistl[[ilenr]][upper.tri(rowlistl[[ilenr]])]= t(rowlistl[[ilenr]])[upper.tri(rowlistl[[ilenr]])]
      }
return(list('cov_list'=rowlistl))
}


Cov_fixl<-function(directory,data_length)
{
  gmatx=c()
 for( ipoint in 1:data_length )
  {
    load(file=paste(directory,ipoint,'gKrow.RData',sep=""))
     ppp= gKr
    gmatx=rbind(gmatx,ppp)
  }
return(list('cov_fix'=gmatx))
}


Cov_ms_uf<- function(u_length,directory,T_len)
{
  fulist<-list() 
    load(file=paste(directory,1,'Krow.RData',sep=""))
    ppp= sKrow 

  for (ilen in 1:T_len) ### startlen:endlen  loop through possible lengthscale 
    {
     fulist<-c(fulist,list(sKrow[ilen,]) )
    }

  for( ipoint in 2:u_length )
  {
     load(file=paste(directory,ipoint,'Krow.RData',sep=""))
     ppp= sKrow  
   ###  go through all lengthscale and combined row vectors.
      for(ilenr in 1:T_len )
      {
       fulist[[ilenr]] = rbind(fulist[[ilenr]],ppp[ilenr,])
      }
  }

return( list('cov_list_uf'= fulist ) )
}


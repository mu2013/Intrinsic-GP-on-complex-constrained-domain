

## horseshoe example
library(Matrix)
require(mgcv)

U_gen_data<-function(){
## plot the function, and its boundary...
fsb <- fs.boundary()
bb=1.5
## m grid in x direction      n grid in y direction
m<-300;n<-150
xm <- seq(-1,4,length=m)
yn<-seq(-1,1,length=n)
xx <- rep(xm,n);yy<-rep(yn,rep(m,n))
tru <- matrix(fs.test(xx,yy,b=bb),m,n)  ## truth

dax1<-c(-0.2 , 0.6 , 1.4 , 2.2 , 3.0)   
dax2<-c(-0.7,-0.2,0.3,0.7)

dmx1<-length(dax1)
dnx2<-length(dax2)

daxx1<-rep(dax1,dnx2)
daxx2<-rep(dax2,rep(dmx1,dnx2))
#dmat <- matrix(fs.test(daxx1,daxx2),dmx1,dnx2)
dmat <- matrix(fs.test(daxx1,daxx2,b=bb),ncol=2)

start_list = rbind( daxx1, daxx2)

mg<-30;ng<-15
xmg <- seq(-1,4,length=mg)
yng<-seq(-1,1,length=ng)
xxg <- rep(xmg,ng);
yyg<-rep(yng,rep(mg,ng))
listrue<-fs.test(xxg,yyg,b=bb)
trug <- matrix(listrue,mg,ng) 

grid_list = rbind(xxg,yyg)

y_n=c(dmat[,1],dmat[,2])

fsb1 <- list(fs.boundary())
names(fsb1[[1]]) <- c("xxg","yyg")
ind <- inSide(fsb1,x=xxg,y=yyg) 

return(list('start_list'=start_list,'grid_list'=grid_list,'y_n'=y_n,'fsb'=fsb,'dmat'=dmat,'truth'=tru,'xs'=xm,'ys'=yn,'xg'=xmg,'yg'=yng,'outsider'=ind))
}


## cross product of two vector which is described as four points
crosspt<-function(A,B,C,D)
{
return( (B[1]-A[1])*(D[2]-C[2])-(D[1]-C[1])*(B[2]-A[2]) )
}

## detect wheter the propsed point cross a line (gap)
intersect<-function(A,B,C,P)
{
 inout1 = crosspt(P,A,P,B) < 0
 inout2 = crosspt(C,A,C,B) < 0
 inout3 = crosspt(A,C,A,P) < 0
 inout4 = crosspt(B,C,B,P) < 0

 inout0=( !(inout1 == inout2)  &&  !(inout3 == inout4) )
 return (inout0);
}


## BM on U shape
U_BM<-function(sta,win,Tmax,Tlen,N,step)
{
Tgrid = seq(step,Tmax,length.out=Tlen);
xdis = rep(sta[1],Tlen); 
x2dis = rep(sta[2],Tlen); 
ydis = matrix(c(0),nrow=N,ncol=Tlen);
y2dis = matrix(c(0),nrow=N,ncol=Tlen);
crossbds = cbind(c(-0.05,3.7),c(0,0))

for (j in 1:N)
{
  for(i in 1:(Tlen-1))
  {
   repeat 
   {
    xdis[i+1]=rnorm(1,xdis[i],sd=sqrt(step) )
    x2dis[i+1]=rnorm(1,x2dis[i],sd=sqrt(step) )
    prpoint=c(xdis[i+1],x2dis[i+1])
    
    if( !intersect(  prpoint,c(xdis[i],x2dis[i]),crossbds[1,],crossbds[2,] ) )
     {
      if ( !any(is.na( fs.test(xdis[i+1],x2dis[i+1],r0=.04,r=.52,l=3.15) )  )    ) break
     }
   }
  }

ydis[j,]=xdis
y2dis[j,]=x2dis
}
return(list("ydis"=ydis,"y2dis"=y2dis))
}


## construct row of cov matrices for different lengthscale 
mcCovrow<-function(index,endlist,N,win,xall,yall)
{
  dn = dim(endlist)[2]
  rK = c(0)
    for(i in 1:dn)
  {
    end = endlist[,i]
    tranp = length(  which(  sqrt( (xall[,index]-end[1])^2+(yall[,index]-end[2])^2 ) < win )   )
    rK[i] = tranp/N/(pi*win^2)
  }
 return(rK)
}



Cov_row_list<-function(start_list,N,T_len,win,BM_dir,Density_dir)
{
  dn<-dim(start_list)[2]

  for (ipoint in 1:dn)    ## loop through data points and set each as starting point of BM
  {
  ##############################  load data  for one single point
  simbatch = 1# 50
  bloop=c(1:simbatch)   ### loop over batches of simulated data points simbatch is defined by simulation

  rrr=matrix(,nrow=1,ncol=T_len)
   for(ibatch in bloop)
   {
   load(file=paste(BM_dir,ipoint,'b',ibatch,'xlist.RData',sep=""))
   ppp= ylist
   rrr= rbind(rrr,ppp[[1]])
   }
   xall=rrr[-1,]

  rrr=matrix(,nrow=1,ncol=T_len)  
   for(ibatch in bloop)
   {
   load(file=paste(BM_dir,ipoint,'b',ibatch,'ylist.RData',sep=""))
   ppp= y2list
   rrr= rbind(rrr,ppp[[1]])
   }
   yall=rrr[-1,]

  end_list = start_list
  Krow=matrix(,nrow=1,ncol=dn)

  for(i in 1:T_len)
  {
   Kr= mcCovrow(i,end_list,N,win,xall,yall)
   Krow= rbind(Krow,Kr)
  }

  sKrow = Krow[-1,]
  save(sKrow, file = paste(Density_dir,ipoint,'Krow.RData',sep="") )
  }

}



Cov_row_list_fix<-function(indexlen,start_list,end_list,N,Tlen,win,BM_dir,Density_dir)
{
  dn <- dim(start_list)[2]
  gn <- dim(end_list)[2]

  for (ipoint in 1:dn)    ## loop through grid points and set each as starting point of BM
  {
  ##############################  load data  for one single point 
  simbatch = 1 # 50
  bloop=c(1:simbatch)   ### loop over batches of simulated data points simbatch is defined by simulation

  rrr=matrix(,nrow=1,ncol=T_len)
   for(ibatch in bloop)
   {
   load(file=paste(BM_dir,ipoint,'b',ibatch,'xlist.RData',sep=""))
   ppp= ylist
   rrr= rbind( rrr,ppp[[1]] )
   }
   xall=rrr[-1,]

  rrr=matrix(,nrow=1,ncol=T_len)  
   for(ibatch in bloop)
   {
   load(file=paste(BM_dir,ipoint,'b',ibatch,'ylist.RData',sep=""))
   ppp= y2list
   rrr= rbind(rrr,ppp[[1]])
   }
   yall=rrr[-1,]

   gKr= mcCovrow(indexlen,end_list,N,win,xall,yall)  ## endlist   gridlist
   
  save(gKr, file = paste(Grid_Density_dir,ipoint,'gKrow.RData',sep="") )

  }

}







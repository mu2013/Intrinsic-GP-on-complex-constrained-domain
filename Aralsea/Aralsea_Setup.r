
set.seed(12345)

require(gamair);require(mgcv)
data(aral); data(aral.bnd)

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

## detect whether point is within a polygon
AnglePoly <- function(bounds,prpoint)
{
	a = bounds[,1]
	b = bounds[,2]
	xt = prpoint[1]
	yt = prpoint[2]
	pn=max(dim(bounds))
	theta = 0
	for(i in 1: (pn-1) )
	  {
	  	U = (a[i]-xt)*(a[i+1]-xt) + (b[i]-yt)*(b[i+1]-yt)
	  	V = (a[i]-xt)*(b[i+1]-yt) - (b[i]-yt)*(a[i+1]-xt)
	  if(V>0||V==0)
	  	{
	    	dtheta = acos( U/sqrt(U^2+V^2) )
	  	}else{
	    	dtheta = (-1)*acos( U/sqrt(U^2+V^2) )
	  	}
		theta = theta+dtheta
		#cat(dtheta,"\n")
	  }
    return(theta)
}

Aralsea_gen_data<-function()
{
## make boundary
lon = aral.bnd[[1]]
lat = aral.bnd[[2]]
lon[8:36] = aral.bnd[[1]][8:36]+0.08
lon[47:52] = aral.bnd[[1]][47:52]-0.08
lon[58:60] = aral.bnd[[1]][58:60]+0.08
lon[64:65] = aral.bnd[[1]][64:65]+0.05
lon[72:93] = aral.bnd[[1]][72:93]-0.08
lon[100:107] = aral.bnd[[1]][100:107]-0.08
lat[1:7] = aral.bnd[[2]][1:7]+0.08
lat[37:42] = aral.bnd[[2]][37:42]-0.08
lat[53:57] = aral.bnd[[2]][53:57]-0.08
lat[94:99] = aral.bnd[[2]][94:99]+0.08
lat[67:71] = aral.bnd[[2]][67:71]-0.08

rmp = -c(4, 8,9,11,13,15,17,20,23,25,28,29,32,33,35,38,40,43,45,46,48,49,51,53,57,61,64,65,
	68,71,73,76,78,79,81,83,84,87,89,91,93,95,97,98,99,101,103:106 )
tbounds = cbind(lon[rmp],lat[rmp])

## rescale the boundary and data points
shrinke = 10

tbounds[,1] = (tbounds[,1]-mean(aral[,1]) ) * shrinke
tbounds[,2] = (tbounds[,2]-mean(aral[,2]) ) * shrinke
bounds = rbind(tbounds,tbounds[1,]) ## add first boundary point as the last point in bounds vector

olon = aral.bnd[[1]]
olat = aral.bnd[[2]]
otbounds = cbind(olon,olat)
otbounds[,1] = (otbounds[,1]-mean(aral[,1]) ) * shrinke
otbounds[,2] = (otbounds[,2]-mean(aral[,2]) ) * shrinke
obounds = rbind(otbounds,otbounds[1,])

crossbds = cbind(c(bounds[31,1],bounds[38,1]),c(bounds[31,2],min(bounds[,2]) )  )
## picked inducing points
pick = c(36,41,47,52, 110,114,120,125,129, 218,227,232,236, 358,365,370,374, 420,423,429,434, 481)

## normalise data points
taa = aral
taa[,1] = (taa[,1]-mean(aral[,1]) ) * shrinke
taa[,2] = (taa[,2]-mean(aral[,2]) ) * shrinke

x= sort(unique(taa[,1]))
y= sort(unique(taa[ ,2]))
sepr = diff( sort(unique(taa[,1])) )[1]    

## na data  113 195 361
    knt <- list(lon=c(58.55,59.09,59.36,59.64,59.91,60.18,58.27,58.55,59.09,
    59.36,59.64,59.91,60.18,60.45,58.27,58.55,58.82,59.09,59.36,59.64,59.91,
    60.18,60.45,58.27,58.55,59.36,59.64,59.91,60.18,58.55,59.36,59.64,59.91,
    60.18,58.55,58.82,59.36,59.64,59.91,60.18,60.45,58.82,59.09,59.64,59.91,
    60.18,59.64),
    lat=c(44.27,44.27,44.27,44.27,44.27,44.27,44.55,44.55,44.55,44.55,44.55,
    44.55,44.55,44.55,44.82,44.82,44.82,44.82,44.82,44.82,44.82,44.82,44.82,
    45.09,45.09,45.09,45.09,45.09,45.09,45.36,45.36,45.36,45.36,45.36,45.64,
    45.64,45.64,45.64,45.64,45.64,45.64,45.91,45.91,45.91,45.91,45.91,46.18))

ind_lon = (knt[[1]]-mean(aral[,1]) ) * shrinke
ind_lat = (knt[[2]]-mean(aral[,2]) ) * shrinke
rmind=c(7,17,15,24,35)

## inducing points
ind=rbind(ind_lon,ind_lat)[,-rmind]

## data points
use = which( !is.na(taa[,3]) )
X1 = rbind(taa[use,1],taa[use,2])
y_n = as.numeric(taa[use,3])
y_me = scale(y_n, center=TRUE, scale=FALSE)
mean_y= mean(y_n)


###### plot grid and data points
xs2 =  sort(c( x[1]-sepr ,x-sepr/2,x ,x[30]+sepr/2,x[30]+sepr  ))    #round( seq(x[1]-sepr,x[length(x)]+sepr,by=sepr/2) , digits=5)
ys2 =  sort(c( y[1]-sepr,y-sepr/2,y ,y[27]+sepr/2,y[27]+sepr ))   #round( seq(y[1]-sepr,y[length(y)]+sepr,by=sepr/2) , digits=5)
nxs2 = length(xs2)
nys2 = length(ys2)

ndat = length(taa[,3])
damx2 = matrix(c(NA),nrow=nxs2,ncol=nys2 )

for(i in 1:nxs2)
{  
  xinde = which( taa[,1] == xs2[i] )
  taa[xinde,4] = i
}
for(j in 1:nys2)
{  
  yinde = which( taa[,2] == ys2[j] )
  taa[yinde,5] = j
}
for(i in 1:ndat)
{
  damx2[ taa[i,4], taa[i,5] ]= taa[i,3]
}

#############################
xxg <- rep(xs2,nys2);yyg<-rep(ys2,rep(nxs2,nys2)) ## grid location
grids = rbind(xxg,yyg,NA)
for(i in 1:(nxs2*nys2) )
{
 if( abs( AnglePoly(obounds,grids[1:2,i]) ) > 1e-8 ) 
 {
   grids[3,i]=1
 }
}
inbound=(which(grids[3,]==1))

return( list('inducing'=ind,'y_me'=y_me,'x_n'=X1,'grids'=grids,'lon'=xs2,'lat'=ys2,'data_col'=damx2,'in_bound'=inbound,'crossbds'=crossbds,'bounds'=bounds ) )

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
#### make cov matrix for data points or inducing points kmm
dn<-dim(start_list)[2]
simbatch=1

for (ipoint in 1:dn)    ## loop through data points and set each as starting point of BM
{
########################## load data  for one single point  # ipoint=8   ## ith data point  
simbatch = 1
bloop=c(1:simbatch)   ### loop over batches of simulated data points simbatch is defined by simulation

rrr=matrix(,nrow=1,ncol=T_len)  ## combine all BM sample path into one array  50*1000*200  batch*N*Tlen
 for(ibatch in bloop)
 {
 load(file=paste(BM_dir,ipoint,'b',ibatch,'xlist.RData',sep=""))
 ppp= xlist
 rrr= rbind(rrr,ppp[[1]])
 }
 xall=rrr[-1,]

rrr=matrix(,nrow=1,ncol=T_len)  
 for(ibatch in bloop)
 {
 load(file=paste(BM_dir,ipoint,'b',ibatch,'ylist.RData',sep=""))
 ppp= ylist
 rrr= rbind(rrr,ppp[[1]])
 }
 yall=rrr[-1,]

  endlist = start_list
  Krow=matrix(,nrow=1,ncol=dn)
  for(i in 1:T_len)
  {
   Kr= mcCovrow(i,endlist,N,win,xall,yall)
   Krow= rbind(Krow,Kr)
  }

  sKrow = Krow[-1,]
  save(sKrow, file = paste(Density_dir,ipoint,'Krow.RData',sep="") )
 }
  
}


Cov_row_list_uf<-function(start_list,f_list,N,T_len,win,BM_dir,Density_dir)
{
  ### kmn
dn=max(dim(start_list))
X1 = f_list

 for (ipoint in 1:dn)    ## loop through data points and set each as starting point of BM
 {
  ##############################  load data  for one single point  # ipoint=8   ## ith data point  
  simbatch = 1 
  bloop=c(1:simbatch)   ### loop over batches of simulated data points simbatch is defined by simulation

  rrr=matrix(,nrow=1,ncol=T_len)  ## combine all BM sample path into one array  50*1000*200  batch*N*Tlen
   for(ibatch in bloop)
   {
   load(file=paste(BM_dir,ipoint,'b',ibatch,'xlist.RData',sep=""))
   ppp= xlist
   rrr= rbind(rrr,ppp[[1]])
   }
   xall=rrr[-1,]

  rrr=matrix(,nrow=1,ncol=T_len)  
   for(ibatch in bloop)
   {
   load(file=paste(BM_dir,ipoint,'b',ibatch,'ylist.RData',sep=""))
   ppp= ylist
   rrr= rbind(rrr,ppp[[1]])
   }
   yall=rrr[-1,]

  win=1
  endlist = X1#stalist
  Krow=matrix(,nrow=1,ncol=max(dim(endlist)) )
    for(i in 1:T_len)
    {
     Kr= mcCovrow(i,endlist,N,win,xall,yall)
     Krow= rbind(Krow,Kr)
    }

    sKrow = Krow[-1,]
    save(sKrow, file = paste(Density_dir,ipoint,'Krow.RData',sep="") )
 }


}



Cov_row_list_fix<-function(indexlen,grid_list,start_list,N,T_len,win,BM_dir,Density_dir)
{

  dn=max(dim(start_list))
  endlist = grid_list  ##use grid list rather than start list

  for (ipoint in 1:dn)    ## loop through inducing points and set each as starting point of BM
  {
  ##############################  load data  for one single point 
  simbatch = 1
  bloop=c(1:simbatch)   ### loop over batches of simulated data points simbatch is defined by simulation

  rrr=matrix(,nrow=1,ncol=T_len)
   for(ibatch in bloop)
   {
   load(file=paste(BM_dir,ipoint,'b',ibatch,'xlist.RData',sep=""))
   ppp= xlist
   rrr= rbind( rrr,ppp[[1]] )
   }
   xall=rrr[-1,]

  rrr=matrix(,nrow=1,ncol=T_len)  
   for(ibatch in bloop)
   {
   load(file=paste(BM_dir,ipoint,'b',ibatch,'ylist.RData',sep=""))
   ppp= ylist
   rrr= rbind(rrr,ppp[[1]])
   }
   yall=rrr[-1,]

   gKr= mcCovrow(indexlen,endlist,N,win,xall,yall)  ## endlist   gridlist
   save(gKr, file = paste(Density_dir,ipoint,'gKrow.RData',sep="") )
  }

}


aralBM<-function(sta,Tmax,Tlen,N,step,crossbds,bounds)
{
  xdis = rep(sta[1],Tlen); 
  ydis = rep(sta[2],Tlen); 

  xlis = array(c(0),c(N,Tlen) );
  ylis = array(c(0),c(N,Tlen) );

  for (j in 1:N)
  {
  for(i in 1:(Tlen-1))
   {
     repeat 
     {
     xdis[i+1]=rnorm(1,xdis[i],sd=sqrt(step) )
     ydis[i+1]=rnorm(1,ydis[i],sd=sqrt(step) )
     
     prpoint=c(xdis[i+1],ydis[i+1])

     if( !intersect(  prpoint,c(xdis[i],ydis[i]),crossbds[1,],crossbds[2,] ) ){
     if ( abs( AnglePoly(bounds,prpoint)) > 1e-10 ) break ## if angle !=0 inside
      }

     }
   }
  xlis[j,]=xdis
  ylis[j,]=ydis  
  }
return(list("xlis"=xlis,"ylis"=ylis))
}



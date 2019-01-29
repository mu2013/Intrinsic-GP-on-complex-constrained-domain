
library(lattice)
library(rgl)
library(pracma)
library(Matrix)
library(R6)

#source('kernel.r')
#source('rbfEU.r')
#source('gp_reg.r')

Torus_transform_coord<-function(R,r,phi,theta)
{
 if(length(phi)==length(theta))
 {
   xp = cos(phi)*(R+r*cos(theta))
   yp = sin(phi)*(R+r*cos(theta))
   zp = r*sin(theta)
   bmp = rbind(xp,yp,zp)
   vbmp = rbind(  as.vector(xp) , as.vector(yp) , as.vector(zp)  )
  }else{
   bmp = NULL
   vbmp =NULL
  }

 xm<-outer(phi,theta,function(phi, theta)cos(phi)*(R+r*cos(theta)))
 ym<-outer(phi,theta,function(phi, theta)sin(phi)*(R+r*cos(theta)))
 zm<-outer(phi,theta,function(phi, theta) r*sin(theta))

 return(list('coord'=bmp,'vector'=vbmp,'x'=xm,'y'=ym,'z'=zm))
}

Torus_true_fun<-function(phi,theta)
{
   wtheta = 0.1
   fun_tor = as.vector( phi + wtheta*sin(theta) ) 
   fun_tor
}


Torus_color_fun<-function(phi,theta,R,r)
{
   igrid = meshgrid(phi,theta)
   tphi = igrid$X
   utheta = igrid$Y

   ## grid points as vector
   eu_grid = Torus_transform_coord(R,r,tphi,utheta)$vector

   ## true function on grid points
   fun_tor = Torus_true_fun(tphi,utheta)

   ## color coding 
   nmp = length(theta)             
   cut_itev=cut( fun_tor,breaks=nmp*(nmp-1) ) 
   crp.rg  <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) # jet.colors
   cols = crp.rg( nmp*(nmp-1))
   mcols= t(matrix( c(cols[cut_itev]),ncol=nmp-1 ))

 return(list('true_fun'=fun_tor,'grids'=eu_grid,'true_cols'=mcols))
}

################################ generate data points ##########################
Torus_gen_data<-function()
{
   phi_te = seq( pi/6, 6.5*pi/4, length.out =6 ) 
   theta_te = c(0,pi/2,3*pi/2,pi)
   x1_grid = meshgrid(phi_te,theta_te)
   x1_phi = as.vector( x1_grid$X )
   x1_theta = as.vector( x1_grid$Y )
   use_list = rbind( x1_phi, x1_theta)[,-c(4,8,12,16,20)]  
   use_list[1,1:3] = use_list[1,1:3]+0.15
   
   start_list = use_list
   Radius = 6
   radius = 4
   test_p = Torus_transform_coord(Radius,radius,start_list[1,],start_list[2,])

   noise = 0.05
   #y_n = matrix( c( c( start_list[1,]+wtheta*sin(start_list[2,]) + rnorm( max(dim(start_list)) ,0,noise) ) )
   y_n = matrix( c( Torus_true_fun(start_list[1,],start_list[2,]) + rnorm( max(dim(start_list)) ,0,noise) ) )

   nmp = 25 
   phi <- seq(0, 5.5/3.2*pi, length.out= nmp)[-1]; 
   theta <- seq(0, 2*pi, length.out= nmp); 
   eu_coord = Torus_transform_coord(Radius,radius,phi,theta)
########################### generate true function on Torus ################################
## define color coding based on true function on Torus
   True_col = Torus_color_fun(phi,theta,Radius,radius)

   return( list('points'=test_p$coord,'y_n'=y_n,'start_list'=start_list,'grid_list'=True_col$grids , 'color'=True_col$true_cols, 'eu_coord'=eu_coord ,'True_fun'=True_col$true_fun) )
}


##########  BM on Bitten Torus
Torus_BM<-function(sta,Tmax,N,step)
{
R <- 6; 
r <- 4; 

Tgrid = seq(step,Tmax,by=step); 
Tlen = length(Tgrid);
thetaDis = rep(sta[2],Tlen); 
phiDis = rep(sta[1],Tlen); 

thetaM = matrix(c(0),nrow=N,ncol=Tlen);
phiM = matrix(c(0),nrow=N,ncol=Tlen);

for (j in 1:N)
 {
   for(i in 1:(Tlen-1))
   {
	  repeat 
	  {
	  thetaDis[i+1] = thetaDis[i]- 1/2*r^(-1)*(R+r*cos( thetaDis[i] ) )^(-1)* sin( thetaDis[i] )*step + r^(-1)*rnorm(1,0,sd=sqrt(step) )
	  phiDis[i+1] = rnorm( 1, phiDis[i], sd= (R+r*cos(thetaDis[i]))^(-1)*sqrt(step) )
      if(phiDis[i+1]<2*pi & phiDis[i+1]>pi/40) break
      }
   }
   thetaM[j,] = thetaDis
   phiM[j,] = phiDis
 }

 return( list( "thetaM" = thetaM, "phiM" = phiM) )
}


## construct row of cov matrices for different lengthscale 
mcCovrow<-function(index,endlist,N,win,phiAll,thetaAll)
{
  R <- 6; 
  r <- 4; 

  dn = dim(endlist)[2]
  rK = c(0)
    for(i in 1:dn)
  {
    end = endlist[,i]

    x_bm = cos(phiAll[,index])*( R+r*cos(thetaAll[,index]) ) 
    y_bm = sin(phiAll[,index])*( R+r*cos(thetaAll[,index]) )
    z_bm = r*sin(thetaAll[,index])

    tranp = length(  which(  sqrt( (x_bm-end[1])^2+(y_bm-end[2])^2 +(z_bm-end[3])^2 ) < win )   )
    rK[i] = tranp/N/(pi*win^2)
  }
 return(rK)
}



Cov_row_list<-function(point_list,N,T_len,win,BM_dir,Density_dir)
{
  dn<-dim(point_list)[2]

  for (ipoint in 1:dn)    ## loop through data points and set each as starting point of BM
  {
  ##############################  load data  for one single point
  simbatch = 1# 50
  bloop=c(1:simbatch)   ### loop over batches of simulated data points simbatch is defined by simulation

  rrr=matrix(,nrow=1,ncol=T_len)
   for(ibatch in bloop)
   {
   load(file=paste(BM_dir,ipoint,'b',ibatch,'thetalist.RData',sep=""))
   ppp= thetalist
   rrr= rbind(rrr,ppp[[1]])
   }
   theta_all=rrr[-1,]

  rrr=matrix(,nrow=1,ncol=T_len)  
   for(ibatch in bloop)
   {
   load(file=paste(BM_dir,ipoint,'b',ibatch,'philist.RData',sep=""))
   ppp= philist
   rrr= rbind(rrr,ppp[[1]])
   }
   phi_all=rrr[-1,]

  end_list = point_list
  Krow=matrix(,nrow=1,ncol=dn)

  for(i in 1:T_len)
  {
   Kr= mcCovrow(i,end_list,N,win,phi_all,theta_all)
   Krow= rbind(Krow,Kr)
  }

  sKrow = Krow[-1,]
  save(sKrow, file = paste(Density_dir,ipoint,'Krow.RData',sep="") )
  }

}



Cov_row_list_fix<-function(indexlen,point_list,end_list,N,Tlen,win,BM_dir,Density_dir)
{
  dn <- dim(point_list)[2]
  gn <- dim(end_list)[2]

  for (ipoint in 1:dn)    ## loop through grid points and set each as starting point of BM
  {
  ##############################  load data  for one single point 
  simbatch = 1 # 50
  bloop=c(1:simbatch)   ### loop over batches of simulated data points simbatch is defined by simulation

  rrr=matrix(,nrow=1,ncol=T_len)
   for(ibatch in bloop)
   {
   load(file=paste(BM_dir,ipoint,'b',ibatch,'thetalist.RData',sep=""))
   ppp= thetalist
   rrr= rbind( rrr,ppp[[1]] )
   }
   theta_all=rrr[-1,]

  rrr=matrix(,nrow=1,ncol=T_len)  
   for(ibatch in bloop)
   {
   load(file=paste(BM_dir,ipoint,'b',ibatch,'philist.RData',sep=""))
   ppp= philist
   rrr= rbind(rrr,ppp[[1]])
   }
   phi_all=rrr[-1,]

   gKr= mcCovrow(indexlen,end_list,N,win,phi_all,theta_all)  ## endlist   gridlist
   
  save(gKr, file = paste(Grid_Density_dir,ipoint,'gKrow.RData',sep="") )

  }

}

## set the color coding for the predictin on grid
pre_color=function( pre_y, Torus_data)
{
## color coding 
nmp = dim(Torus_data$color)[2]
pre_y_m = t(matrix( c(pre_y),ncol=nmp-1 ))

oooo<-as.vector(pre_y) ; 
oooo[ which(oooo==min(oooo))]=min(Torus_data$True_fun); 
oooo[ which(oooo==max(oooo))]=max(Torus_data$True_fun); 

crp.rg  <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) # jet.colors
cols = crp.rg( nmp*(nmp-1)) #cols<-rainbow(nmp*(nmp-1)) #heat.colors(nmp*(nmp-1)) #topo.colors(nmp*(nmp-1)) 

cut_pre_y=cut( oooo,breaks=nmp*(nmp-1) ) 
pre_cols= t(matrix( c(cols[cut_pre_y]),ncol=nmp-1 ))
pre_cols
}

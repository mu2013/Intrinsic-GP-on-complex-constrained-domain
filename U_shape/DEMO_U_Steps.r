
set.seed(12345)
########################### initialisatino and plot true function #######################################
source('U_setup.r')
source('../inference/Cov_ms.r')
library(R6)
source('../inference/rbfEU.r')
source('../inference/kernel.r')
source('../inference/gp_reg.r')
source('../inference/intrinsic_kernel.r')
source('../inference/Heat_manif.r')
source('../inference/gp_reg_manif.r')

U_data = U_gen_data()
y_n = U_data$y_n
data_length = length(y_n)

true_fun = U_data$truth
fsb = U_data$fsb
start_list = U_data$start_list
grid_list = U_data$grid_list


############################################################################################################
### plot U shape and data points with true function
dev.new()
image(U_data$xs,U_data$ys,true_fun,col=heat.colors(100),xlab="x",ylab="y",main='True Function')
lines(fsb$x,fsb$y,lwd=3)
contour(U_data$xs,U_data$ys,true_fun,levels=seq(-10,10,by=.5),add=TRUE,labcex=1)  ## b>1
points(start_list[1,],start_list[2,],lwd=2,pch=4)


############################### classic GP inference and plot GP prediction #######################################
## use classic normal RBF kernel 
ker = RBF$new(c(1,1))   

## GP regression with RBF kernel
reg_m = gp_reg$new( t(y_n),start_list,nsig= 1,ker,jitter=0,mean_f='on')  
reg_m$optimisMn_l(c(1,1,1))

## GP predictin on grid points
res = reg_m$predict( grid_list )
pre_y = res$pre_m

## plot prediction as a heat map
pre_y[!U_data$outsider]= NA
plpg = matrix(pre_y,length(U_data$xg),length(U_data$yg))

dev.new()
image(U_data$xg,U_data$yg,plpg,col=heat.colors(100),xlab="x",ylab="y",main='normal GP')
lines(fsb$x,fsb$y,lwd=3)
contour(U_data$xg,U_data$yg,plpg,levels=seq(-6,6,by=.5),add=TRUE,labcex=1,lwd=1) ## b<=1
points(start_list[1,],start_list[2,],lwd=2,pch=4)


############################ Intrinsic GP on U shape domain  ########################################

######################### Create Cov matrices from BM paths #########################################
Tmax = 0.4;   ## maximum  lengthscale 
T_len= 80     ## step = Tmax/T_len 
N = 100000    ## monte carlo trajectories of BM
directory = paste('./Cov_row_list/d',sep="")
cov_list = Cov_ms(data_length,directory,T_len)$cov_list

############### find the optimum cov matrix and corresponding parameters ########################
## use intrinsic heat kernel on manifold
in_ker = Heat_manif$new(cov_list) 
## intrinsic GP regression 
in_reg_m = gp_reg_manif$new( t(y_n),nsig= 1,in_ker,jitter=0,mean_f='on')  
len_range=c(60:80) ## we can ignore small lengtscales for numerical stability
op_par = in_reg_m$optimisM(len_range)


########################################################################################
## The following block of codes does not need to be run. Because the results have been generated.
## If the user run the BM simulation and generate the sample paths, then the user can run
## the follwoing block to get the cross covariance of grid points and data points given
## the optimised lengthscale or diffusion time.
#########################################################################################
# ####### construct cov matrices of grid points given a chosen lengthscale 
# indexlen = op_par$op_index ## this is defined by running GP regression 
# start_list = U_data$start_list
# end_list = U_data$grid_list
# BM_dir = paste('./BM_Paths/d',sep="")
# Grid_Density_dir = paste('./Cov_row_list_fix/d',sep="")
# win = 0.1
# Cov_row_list_fix(indexlen,start_list,end_list,N,Tlen,win,BM_dir,Density_dir)
##################################################################


############### combine the rows of the cross covariance to make the cross covariance matrix
directory = paste('./Cov_row_list_fix/d',sep="")
cov_fix = Cov_fixl(directory,data_length)$cov_fix


### intrinsic GP prediction on grid points on U shape domain.
gridpr = in_reg_m$predictM(op_par,t(cov_fix) )$pre_m

gplpg = matrix(gridpr,length(U_data$xg),length(U_data$yg))
gplpg[!U_data$outsider]= NA

dev.new()
image(U_data$xg,U_data$yg,gplpg,col=heat.colors(100),xlab="x",ylab="y",main='Intrisic GP on manifold')
lines(fsb$x,fsb$y,lwd=3)
contour(U_data$xg,U_data$yg,gplpg,levels=seq(-6,6,by=.5),add=TRUE,labcex=1) ## b<=1
points(start_list[1,],start_list[2,],lwd=2,pch=4)


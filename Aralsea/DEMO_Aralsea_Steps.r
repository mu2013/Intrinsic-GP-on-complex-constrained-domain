
set.seed(12345)
########################### initialisation and plot true function #######################################

source('../inference/Cov_ms.r')
library(R6)
source('../inference/rbfEU.r')
source('../inference/kernel.r')
source('../inference/gp_reg.r')
source('../inference/gp_reg_sparse.r')
source('../inference/intrinsic_kernel.r')
source('../inference/Heat_manif.r')
source('../inference/gp_reg_manif_sparse.r')
source('../inference/gp_reg_manif.r')

source('Aralsea_Setup.r')

Aral_data = Aralsea_gen_data()

############################################################################################################
### plot Aral sea data points 
image(Aral_data$lon,Aral_data$lat,Aral_data$data_col,xlab='longitude',ylab='lattitude')


############################### classic Sparse GP inference and plot GP prediction #######################################
## use classic normal RBF kernel 
ker = RBF$new(c(1,1))   
## set the inducing point
u = Aral_data$inducing
## sparse GP regression with RBF kernel
reg_m = gp_reg_sparse$new( t(matrix(Aral_data$y_me) ),u,Aral_data$x_n,nsig= 1,ker,jitter=0,mean_f='off')  
reg_m$optimisM( log(c(8,5,2)) )

## sparse GP prediction
grids = Aral_data$grids[1:2,Aral_data$in_bound]
ppu<-reg_m$upredict(Aral_data$x_n,grids)

## plot prediction
grid_bound = Aral_data$grids
grid_bound[3,Aral_data$in_bound]=ppu$pre_m
pdamxg = matrix(grid_bound[3,], length(Aral_data$lon), length(Aral_data$lat))

image(Aral_data$lon,Aral_data$lat,pdamxg,xlab='x',ylab='y',main='normalGP')
contour(Aral_data$lon,Aral_data$lat,pdamxg,levels=seq(-6,6,by=1),add=TRUE)



############################### Intrinsic sparse GP inference and plot GP prediction #######################################
######################### Create Cov matrices from BM paths #########################################
 Tmax = 8   ## maximum  lengthscale or diffusion time
 N = 100000  ## number of BM paths
 T_len =80   ## length of BM steps
 win = 1
 data_length = dim(Aral_data$inducing)[2]

## covariance between inducing points k_uu
 directory = paste('./Cov_row_list/d',sep="")
 cov_list = Cov_ms(data_length,directory,T_len)$cov_list

## covariance between inducing points and data points k_uf
 directory = paste('./Cov_row_list_uf/d',sep="")
 cov_list_uf = Cov_ms_uf(data_length,directory,T_len)$cov_list_uf
 
############### find the optimum cov matrix and corresponding parameters ########################
## use intrinsic heat kernel on manifold
in_ker = Heat_manif$new(cov_list,cov_list_uf) 
## intrinsic GP regression 
in_reg_m = gp_reg_manif_sparse$new( t(Aral_data$y_me),Aral_data$inducing,Aral_data$x_n,nsig= 1,in_ker,jitter=0,mean_f='off')  
len_range=c(60:80) ## we can ignore small lengtscales for numerical stability
op_par = in_reg_m$optimisMn(len_range)


########################################################################################
## The following block of codes does not need to be run. Because the results have been generated.
## If the user run the BM simulation and generate the sample paths, then the user can run
## the follwoing block to get the cross covariance of grid points and data points given
## the optimised lengthscale or diffusion time.
######################################################################################### 
# ####### construct cov matrices of grid points given a chosen lengthscale  
#Aral_data = Aralsea_gen_data()
#Tmax = 8;   ## maximum  lengthscale 
#N = 100000    ## monte carlo trajectories of BM
#T_len=80      ## Tlen=(Tmax-Tmin)/0.1  
#win=1
#indexlen= op_par$op_index 
#grids = Aral_data$grids[1:2,Aral_data$in_bound]
#Density_dir= paste('./Cov_row_list_fix/d',sep="")
#Cov_row_list_fix(indexlen,grids,start_list,N,T_len,win,BM_dir,Density_dir)


############### combine the rows of the cross covariance to make the cross covariance matrix
directory = paste('./Cov_row_list_fix/d',sep="")
cov_list_ux = Cov_fixl(directory,data_length)$cov_fix

### intrinsic GP prediction on grid points
pre = in_reg_m$predict(op_par,t(cov_list_ux) )

grid_bound = Aral_data$grids
grid_bound[3,Aral_data$in_bound]=pre$pre_m
pdamxg = matrix(grid_bound[3,], length(Aral_data$lon), length(Aral_data$lat))

image(Aral_data$lon,Aral_data$lat,pdamxg,xlab='x',ylab='y',main='normalGP')
contour(Aral_data$lon,Aral_data$lat,pdamxg,levels=seq(-6,6,by=1),add=TRUE)




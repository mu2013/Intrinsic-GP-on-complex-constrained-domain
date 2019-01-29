
set.seed(123456)
########################### initialisatino and plot true function #######################################
source('Torus_setup.r')
source('../inference/Cov_ms.r')
library(R6)
source('../inference/rbfEU.r')
source('../inference/kernel.r')
source('../inference/gp_reg.r')
source('../inference/intrinsic_kernel.r')
source('../inference/Heat_manif.r')
source('../inference/gp_reg_manif.r')


SEED = 123456
set.seed(SEED)

Torus_data = Torus_gen_data()


############################################################################################################
### plot the bitten Torus and data points with true function
bg3d(color="#887777")
surface3d(Torus_data$eu_coord$x, Torus_data$eu_coord$y, Torus_data$eu_coord$z,col= Torus_data$color, specular="black")
points3d(Torus_data$points[1,],Torus_data$points[2,],Torus_data$points[3,],size=4)


############################### classic GP inference and plot GP prediction #######################################
## use classic normal RBF kernel 
ker = RBF$new(c(1,1))

## GP regression with RBF kernel
reg_m = gp_reg$new( t(Torus_data$y_n),Torus_data$points,nsig= 1,ker,jitter=0,mean_f='on')
reg_m$optimisMn(c(10,10,1e-1))

## GP prediction on Bitten Torus
res = reg_m$predict( Torus_data$grid_list )

## color coding
pre_y = t(res$pre_m)
pre_cols = pre_color( pre_y, Torus_data)

bg3d(color="#887777")
surface3d(Torus_data$eu_coord$x[-1,], Torus_data$eu_coord$y[-1,], Torus_data$eu_coord$z[-1,],col= pre_cols[-1,] , specular="black") ## ,back="lines"



############################ Intrinsic GP on U shape domain  ########################################

######################### Create Cov matrices from BM paths #########################################
 Tmax = 70   ## maximum  lengthscale or diffusion time
 N = 100000  ## number of BM paths
 T_len =350  ## length of BM steps
 win = 2
 directory = paste('./Cov_row_list/d',sep="")
 data_length = length(Torus_data$y_n)
 cov_list = Cov_ms(data_length,directory,T_len)$cov_list

############### find the optimum cov matrix and corresponding parameters ########################
## use intrinsic heat kernel on manifold
in_ker = Heat_manif$new(cov_list) 
## intrinsic GP regression 
in_reg_m = gp_reg_manif$new( t(Torus_data$y_n),nsig= 1,in_ker,jitter=0,mean_f='on')  
len_range=c(200:350) ## we can ignore small lengtscales for numerical stability
op_par = in_reg_m$optimisMn(len_range)

# ####### construct cov matrices of grid points given a chosen lengthscale 
# Torus_data = Torus_gen_data()
# point_list = Torus_data$points
# end_list = Torus_data$grid_list
# Tmax = 70   ## maximum  lengthscale or diffusion time
# N = 100000  ## number of BM paths
# T_len =350  ## length of BM steps
# win = 2
# indexlen = op_par[[1]] ## this is defined by running intrinsic GP regression in file DEMO_Torus_steps.r
# BM_dir = paste('./BM_Paths/d',sep="")
# Grid_Density_dir = paste('./Cov_row_list_fix/d',sep="")
# Cov_row_list_fix(indexlen,point_list,end_list,N,Tlen,win,BM_dir,Density_dir)
########################################################################################

############### combine the rows of the cross covariance to make the cross covariance matrix
directory = paste('./Cov_row_list_fix/d',sep="")
cov_fix = Cov_fixl(directory,data_length)$cov_fix

### intrinsic GP prediction on Bitten Torus
gridpr = in_reg_m$predictMn(op_par,t(cov_fix) )$pre_m

## color coding of in-GP prediction
pre_y_in = t(gridpr) + mean(Torus_data$y_n)
pre_cols_in = pre_color( pre_y_in, Torus_data)

bg3d(color="#887777")
surface3d(Torus_data$eu_coord$x[-1,], Torus_data$eu_coord$y[-1,], Torus_data$eu_coord$z[-1,],col= pre_cols_in[-1,] , specular="black") ## ,back="lines"


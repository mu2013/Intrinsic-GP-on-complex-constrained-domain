########################################################################################
## The following codes do not need to be run unless the user want to run the BM simulation 
## and generate the sample paths. These codes are used to produce lists of rows of covariance
## matrices given different lengthscale or diffusion time.
#########################################################################################

set.seed(12345)
source('U_setup.r')

U_data = U_gen_data()

start_list = U_data$start_list
Tmax = 0.4;  ## 0.4 maximum  lengthscale 
N = 100000   ## 50000; ## monte carlo trajectories of BM
T_len = 80   ## Tlen=(Tmax-Tmin)/step 

win=0.1


BM_dir = paste('./BM_Paths/d',sep="")
Density_dir= paste('./COV_row_list/d',sep="")

Cov_row_list(start_list,N,T_len,win,BM_dir,Density_dir)




# ####### construct cov matrices of grid points given a chosen lengthscale 

# source('U_setup.r')
# U_data = U_gen_data()

# start_list = U_data$start_list
# end_list = U_data$grid_list
# Tmax = 0.4  ## maximum  lengthscale or diffusion time
# N = 100000  ## number of BM paths
# T_len = 80  ## length of BM steps

# win = 0.1
# indexlen = 78 ## this is defined by running intrinsic GP regression in file DEMO_U_steps.r

# BM_dir = paste('./BM_Paths/d',sep="")
# Grid_Density_dir = paste('./Cov_row_list_fix/d',sep="")

# Cov_row_list_fix(indexlen,start_list,end_list,N,Tlen,win,BM_dir,Density_dir)


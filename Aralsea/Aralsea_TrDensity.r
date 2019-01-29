########################################################################################
## The following codes do not need to be run unless the user want to run the BM simulation 
## and generate the sample paths. These codes are used to produce lists of rows of covariance
## matrices given different lengthscale or diffusion time.
#########################################################################################

set.seed(12345)
source('Aralsea_setup.r')

Aral_data = Aralsea_gen_data()

Tmax = 8;    ## maximum  lengthscale 
N = 100000   ## monte carlo trajectories of BM
T_len=80     ## Tlen=(Tmax-Tmin)/0.1  
win=1


BM_dir = paste('./BM_Paths/d',sep="")
start_list = Aral_data$inducing

## row of covariance between inducing points  k_uu
Density_dir= paste('./Cov_row_list/d',sep="")
Cov_row_list(start_list,N,T_len,win,BM_dir,Density_dir)

## orw of covariance between inducing points and data points k_uf
Density_dir= paste('./Cov_row_list_uf/d',sep="")
f_list = Aral_data$x_n
Cov_row_list_uf(start_list,f_list,N,T_len,win,BM_dir,Density_dir)


# ####### construct cov matrices of grid points given a chosen lengthscale 
#indexlen=67
#grids = Aral_data$grids[1:2,Aral_data$in_bound]
#Density_dir= paste('./Cov_row_list_fix/d',sep="")
#Cov_row_list_fix(indexlen,grids,start_list,N,T_len,win,BM_dir,Density_dir)


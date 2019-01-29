########################################################################################
## The following codes do not need to be run unless the user want to run the BM simulation 
## and generate the sample paths. These codes are used to produce lists of rows of covariance
## matrices given different lengthscale or diffusion time.
#########################################################################################

source('Torus_Setup.r')

Torus_data = Torus_gen_data()
point_list = Torus_data$points

Tmax = 70    ## maximum  lengthscale 
N = 100000   ## monte carlo trajectories of BM
T_len = 350  ## Tlen=(Tmax-Tmin)/step 
win = 2

BM_dir = paste('./BM_Paths/d',sep="")
Density_dir= paste('./Cov_row_list/d',sep="")

Cov_row_list(point_list,N,T_len,win,BM_dir,Density_dir)


# ####### construct cov matrices of grid points given a chosen lengthscale 

# Torus_data = Torus_gen_data()
# point_list = Torus_data$points
# end_list = Torus_data$grid_list

# Tmax = 70   ## maximum  lengthscale or diffusion time
# N = 100000  ## number of BM paths
# T_len =350  ## length of BM steps

# win = 2
# indexlen = 344 ## this is defined by running intrinsic GP regression in file DEMO_Torus_steps.r

# BM_dir = paste('./BM_Paths/d',sep="")
# Grid_Density_dir = paste('./Cov_row_list_fix/d',sep="")

# Cov_row_list_fix(indexlen,point_list,end_list,N,Tlen,win,BM_dir,Density_dir)



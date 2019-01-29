
source('Torus_Setup.r')

Torus_data = Torus_gen_data()

### we need to start BM from different starting points
start_list = Torus_data$start_list
SEED = 12345

ptm <- proc.time()
######## simulating BM code  ##############################################
for( idat in 1:dim(start_list)[2] )
{
SEED = SEED +100*idat
set.seed(SEED)
sta= start_list[,idat] ## starting point for phi and theta

N = 100000     ## number of sample paths
step = 0.2     ##  delta t
Tmax = 70      ##  max diffusion time

spath<-Torus_BM(sta,Tmax,N,step)

thetalist<-list()
philist<-list()

## combine monte carlo list for different starting points.
thetalist[[1]]<-spath$thetaM
philist[[1]]<-spath$phiM

save(thetalist,file=paste('./BM_Paths/','d',idat,'b',1,'thetalist.RData',sep=""))
save(philist,file=paste('./BM_Paths/','d',idat,'b',1,'philist.RData',sep=""))
}
proc.time() - ptm

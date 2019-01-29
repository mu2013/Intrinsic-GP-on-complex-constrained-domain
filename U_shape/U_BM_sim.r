

source('U_setup.r')

### we need to start BM from 20 different starting points
U_data = U_gen_data()

start_list = U_data$start_list

Tmax = 0.4;      ## maximum  lengthscale 
N = 100000       ## monte carlo trajectories of BM
Tlen = 80        ## Tlen=(Tmax)/step  every two lengthscale sample have difference Tmax/Tlen
step = Tmax/Tlen


ptm <- proc.time()

for( iid in 1:dim(start_list)[2])
{
idat = iid
SEED = 12345
batch = 1 
set.seed(SEED)

ylist<-list()
y2list<-list()
## set stating point and simulate N BM trajectoriers with choise of lenghtscale from 0.01 to Tmax (number of steps Tmax)
sta = start_list[,idat]

spath<-U_BM(sta,win,Tmax,Tlen,N,step)

ydis<-spath$ydis
y2dis<-spath$y2dis

## combine monte carlo list for different starting points.
ylist[[length(ylist)+1]] <- ydis
y2list[[length(y2list)+1]] <- y2dis

save(ylist,file=paste('./BM_Paths/','d',idat,'b',batch,'xlist.RData',sep=""))
save(y2list,file=paste('./BM_Paths/','d',idat,'b',batch,'ylist.RData',sep=""))
}

proc.time() - ptm

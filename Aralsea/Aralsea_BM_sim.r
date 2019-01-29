


source('Aralsea_Setup.r')
### we need to start BM from different inducing points
Aral_data = Aralsea_gen_data()

ptm <- proc.time()

SEED=  54321

ptm <- proc.time()

for(idat in 1:42)
{
 N = 100000        ## number of seeds * N * number of starting points
 Tlen = 80         ## biggest lengthscale is Tlen*step
 sta= Aral_data$inducing[,idat]   ## starting location (inducing points location)
 step = 0.1

SEED =SEED+idat*100
set.seed(SEED)
## simulate N BM trajectoriers with choise of lenghtscale from 0.01 to Tmax (number of steps Tmax)
bml = aralBM(sta,Tmax,Tlen,N,step,Aral_data$crossbds,Aral_data$bounds)

xdis = bml$xlis
ydis = bml$ylis

xlist<-list()
ylist<-list()
xlist[[length(xlist)+1]] <- xdis
ylist[[length(ylist)+1]] <- ydis

save(xlist,file=paste('./BM_Paths/','d',idat,'b',1,'xlist.RData',sep=""))
save(ylist,file=paste('./BM_Paths/','d',idat,'b',1,'ylist.RData',sep=""))

}

proc.time() - ptm





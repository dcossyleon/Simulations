#Cycle through V's

library(reshape2)
library(MCMCglmm)
library(MCMCvis)
library(MASS)

n=100
mzcor <- 0.4
dzcor <- 0.2
mz <- as.data.frame(mvrnorm(n, mu=c(0, 0), Sigma=matrix(c(1,mzcor,mzcor,1), nrow=2), empirical=T)) #n is # of families
dz <- as.data.frame(mvrnorm(n, mu=c(0, 0), Sigma=matrix(c(1,dzcor,dzcor,1), nrow=2), empirical=T))

mz$fam <- seq(1:dim(mz)[1])
mz <- melt(mz, id='fam')
mz$zyg=1

dz$fam <- c((dim(dz)[1]+1):(dim(dz)[1]*2))
dz <- melt(dz, id='fam')
dz$zyg=2

tw <- rbind(mz,dz)

tw <- tw[order(tw$fam),]
family <- as.data.frame(tw$fam)
tw$id <- c(1:dim(tw)[1])

a <- apply(family, c(1), function(x) ifelse(tw$fam==x & tw$zyg==1, 1, ifelse(tw$fam==x & tw$zyg==2, 0.5, 0)))
diag(a) <- 1
rownames(a) <-  tw$id
colnames(a) <-  tw$id

N <- dim(tw)[1]
a <- a + diag(0.00001, N)
i <- rep(1:N,rep(N,N))
j <- rep(1:N,N)
s <- Matrix::spMatrix(N,N,i,j,as.vector(a))
Ginv <- Matrix::solve(s)
class(Ginv) <- 'dgCMatrix'
rownames(Ginv) <- Ginv@Dimnames[[1]] <- with(tw,id)



Vs <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

start <- Sys.time()
m1 <- mclapply(Vs, FUN=function(x) MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, nitt=100000, burnin=10000, thin=100,
                                      prior=list(R=list(V=1, nu=0.002), G=list(G1=list(V=x, nu=0.2)))), mc.cores=4)
end <- Sys.time()
end-start

m1.HPD<- as.data.frame(do.call(rbind, mclapply(1:10, function(x) HPDinterval((m1[[x]]$VCV[,'id'])/(m1[[x]]$VCV[,'id']+m1[[x]]$VCV[,"units"]))[c(1,2)], mc.cores=1)), row.names=Vs)
  m1.HPD$CI <- m1.HPD$V2-m1.HPD$V1
  #m1.HPD
  plot(m1.HPD$CI)


#raw scale
m1.HPD.r<- as.data.frame(do.call(rbind, mclapply(1:10, function(x) HPDinterval(m1[[x]]$VCV)[c(1,3)], mc.cores=1)), row.names=Vs)
  m1.HPD.r$CI <- m1.HPD.r$V2-m1.HPD.r$V1
  #m1.HPD.r
  plot(m1.HPD.r$CI, sub="m1.HPD.r$CI, raw")

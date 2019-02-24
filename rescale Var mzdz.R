#rescale
mzcor <- 0.8
dzcor <- 0.4
mz <- as.data.frame(mvrnorm(100, mu=c(0, 0), Sigma=matrix(c(2,mzcor,mzcor,2), nrow=2), empirical=T)) #n is # of families
dz <- as.data.frame(mvrnorm(100, mu=c(0, 0), Sigma=matrix(c(2,dzcor,dzcor,2), nrow=2), empirical=T))

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

test <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=0.4, nu=2)))

fitAE <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=test, burnin=10000, nitt=100000, thin=100)
plot(fitAE$VCV)
data.frame('a2'=mean(fitAE$VCV[,'id']/(fitAE$VCV[,'id']+fitAE$VCV[, 'units'])), 'e2'=mean(fitAE$VCV[,'units']/(fitAE$VCV[,'id']+fitAE$VCV[, 'units'])))

library(MASS)
library(mvdalab)
library(psych)

loads1=c(.4,.4,.2,.2,.2,0,0,0,0,0) #pick how you want your variables to load onto your factor
loads2=c(0,0,0,0,0,.2,.2,.2,.2,.2) #how they load onto factor 2
corfactors <-0.0 #how correlated should both factors be
resid <- 0.01 #R doesn't like it when it's perfectly 0, but how well do you want it to fit?
N=900 #sample size

#loads <- as.matrix(cbind(loads1, loads2))
#factorvar <- matrix(c(corfactors,corfactors),ncol(loads),ncol(loads),byrow=F)
#diag(factorvar) <- 1
#sigma <- loads%*%factorvar%*%t(loads)

sigma<-round(cor(beh[,c(7,8,11,12,13,15,17,18,19,20)], use="pairwise.complete.obs", method="spearman"),2)



diag(sigma) <- 1
data=mvrnorm(N,mu=rep(0,10),Sigma=sigma, empirical=T)

#data=mvrnorm.svd(n = N, mu=rep(0,6), Sigma = sigm <a, empirical = T, Dist = "poisson", poisson.mean = 1)


factanal(data, factors=4, rotation="promax") #maximum likelihood extraction EFA

fa(data, nfactors=3, rotate="oblimin",fm="pa")  #principal axis, for non-normal
 



#Notes
#When loads are al .7, corfactors=0, and resid=0:
  #BIC=-18.42
  #RMSEA index=0 (should be less than .1, ideally less than 0.05)
  #Fit based off diag=1
#When corfactors are .5, stuff loads heavily onto factor 1, but weakly or split onto 2
  #mostly similar between ML and PA extraction methods
  #BIC=-22.82
  #RMSEA =0
  #fit based off diag, 1
  #Tucker Lewis 1.036

#If 3 factors, with .2,.7.,7 loadings, 
  #misspecify 2 factors: drops the .2, says is sufficient
  #misspecify 4 factors: takes the .2 loading and splits into factor3 (.47 loading) and factor 4(.38 loading)

#Cumulative Variance always seems low. Even when properly specified...why?



#Scree Plot
library(nFactors)

ev <- eigen(cor(data, use="pairwise.complete.obs"))
ap <- parallel(subject=nrow(data),var=ncol(data),rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)

plotnScree(nS)

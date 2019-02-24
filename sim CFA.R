library(MASS)
library(mvdalab)
library(psych)
library(lavaan)

loads1=c(.7,.7,.7,0,0,0) #pick how you want your variables to load onto your factor
loads2=c(0,0,0,.5,.6,.7) #how they load onto factor 2
corfactors <- 0 #how correlated should both factors be
resid <- 0.5 #R doesn't like it when it's perfectly 0, but how well do you want it to fit?
N=300 #sample size

loads <- as.matrix(cbind(loads1, loads2))
factorvar <- matrix(c(1,corfactors,corfactors,1),ncol(loads),ncol(loads),byrow=T)
sigma <- loads%*%factorvar%*%t(loads)
diag(sigma) <- 1
data=mvrnorm(N,mu=rep(0,6),Sigma=sigma,empirical=T)

colnames(data) <- c("var1", "var2", "var3", "var4", "var5", "var6")

model <- "
lv1=~var1+var2+var3
lv2=~var4+var5+var6
"
fit<- cfa(model=model, data)
summary(fit, standardized=T)

fitMeasures(fit)
#factor: factor loadings lv vs std
#The std covar is the correlations between the two factors
#Then I'm getting the residual variance and the LV variance
#fit measures-- which I have to eval.

#cfi above 95?
#rmsea.pvalue -- want something high as well?
#aic vs bic - lower is better.....

cor(predict(fit)[,1], data[,1])


data=mvrnorm.svd(n = N, mu=rep(0,6), Sigma = sigma, empirical = T, Dist = "poisson")


factanal(data, factors=2, rotation="promax")







#Scree Plot
library(nFactors)

ev <- eigen(cor(data, use="pairwise.complete.obs"))
ap <- parallel(subject=nrow(data),var=ncol(data),rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)

plotnScree(nS)

#Jan 29, 2018 ZANB Model Simulation
#Two-part model

library(aster)
library(VGAM)

#1) Generate data
n=1000
X1 <- runif(n, min=0, max=1)
size=3 #DISPERSION PARAMETER
k=0 #truncation limit
B0 <- 4
B1 <- 2
C0 <- 1.6
C1 <- 5 #is hurdle estimating properly?

Y=ifelse(rbinom(n=n, size = 1, prob=(exp(C0+C1*X1)/(1+(exp(C0+C1*X1)))))==0, 0, rktnb(n=n, size=size, k=k,  mu=exp(B0+B1*X1), xpred=1)) #lots of zeros, plus the zero-truncated neg.binomial 
d <- data.frame(x=X1,y=Y)
hist(d$y)

#data manip
d.no0 <- d
d.no0$y <- ifelse(d.no0$y>0, d.no0$y, NA)
d.bin <- d
d.bin$y <- ifelse(d.bin$y>0, 1, 0)



#glm for ZANB? 
summary(hurdle(y~x, dist="negbin", link="logit", data=d))

#----
#GLM for Binomial
summary(glm(y~x, family = binomial, data=d.bin))
glm.bin<- glm(y~x, family = binomial, data=d.bin)
logLik(glm.bin)

#GLM for truncated
summary(vglm(y~x, family=posnegbinomial, data=d.no0)) #*Intercept2 is the log of the dispersion






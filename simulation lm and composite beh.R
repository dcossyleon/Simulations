#March 28, 2018

perc=.01 #set how much variance explained I want
         # note: 1=100%, .1=10% and .01=1%, etc
         # NOTE: effect size is defined for p=q=.5

geffect=(perc*1)/(1-perc) #not sure what this one is really doing
a1=sqrt(geffect/(.5)) #a=sqrt(beta/2pq), derived from the equation you'd use to find the mean of the population, when working with hardy weinberg equation. 
b1=c(-a1,0,a1) #allows us to translate R2 into a beta for a given Variance


d1 <- mvrnorm(n=250, mu=b1[1], Sigma=1, empirical=TRUE)
d2 <- mvrnorm(n=500, mu=b1[2], Sigma=1, empirical=TRUE)
d3 <- mvrnorm(n=250, mu=b1[3], Sigma=1, empirical=TRUE)

d <- rbind(d1,d2,d3)

df <- data.frame(pheno=d, geno=c(rep(-1,250), rep(0, 500), rep(1,250)))
                 

summary(lm(pheno~geno, data=df))



#Second trial
perc=.2
geffect=(perc*1)/(1-perc)
a1=sqrt(geffect/(.5))
b1=c(-a1,0,a1)

d1 <- mvrnorm(n=250, mu=b1[1], Sigma=1, empirical=TRUE)
d2 <- mvrnorm(n=500, mu=b1[2], Sigma=1, empirical=TRUE)
d3 <- mvrnorm(n=250, mu=b1[3], Sigma=1, empirical=TRUE)

d <- rbind(d1,d2,d3)

df <- data.frame(pheno=d, geno=c(rep(-1,250), rep(0, 500), rep(1,250)))


summary(lm(pheno~geno, data=df))



#Now with phenos that are correlated

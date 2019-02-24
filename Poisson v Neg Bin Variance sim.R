#1/23/2018 #Poisson vs Neg Binomial Variance Sim

#Binomial "studies"
bi.glm <- function(x){
  n=1000
  size=1
  X1 <- runif(n, min=0, max=1)
  B0 <- 1
  B1 <- 2
  Y <- rbinom(n=n, size=size, prob=(exp(B0+B1*X1)/(1+(exp(B0+B1*X1))))) 
  int<- summary(glm(Y~X1, family = binomial))$coefficient[1]
  b1<- summary(glm(Y~X1, family = binomial))$coefficient[2]
  d <- data.frame(int,b1)
  return(d)
}
b.sim <- do.call(rbind, lapply(1:1000, function(x) bi.glm()))
b.var <- data.frame(V.int=var(b.sim$int), V.b1=var(b.sim$b1))
b.var

#Poisson "studies"
p.glm <- function(x){
  n=1000
  X1 <- runif(n, min=0, max=1)
  B0 <- 3
  B1 <- 2
  Y <- rpois(n=n, lambda=exp(B0+B1*X1))
  int<-(glm(Y~X1, family="poisson")$coefficient[1])
  b1<- (glm(Y~X1, family="poisson")$coefficient[2]) 
  d <- data.frame(int,b1)
}
p.sim <- do.call(rbind, lapply(1:1000, function(x) p.glm()))
p.var <- data.frame(V.int=var(p.sim$int), V.b1=var(p.sim$b1))
p.var






#____Binomial Y equations
#Binomial "studies"

#Y: prob=(exp(B0+B1*X1)/(1+(exp(B0+B1*X1))))
bi.glm <- function(x){
  n=1000
  size=1
  X1 <- runif(n, min=0, max=1)
  B0 <- 1
  B1 <- 2
  Y <- rbinom(n=n, size=size, prob=(exp(B0+B1*X1)/(1+(exp(B0+B1*X1))))) 
  int<- summary(glm(Y~X1, family = binomial))$coefficient[1]
  b1<- summary(glm(Y~X1, family = binomial))$coefficient[2]
  d <- data.frame(int,b1)
  return(d)
}
b.sim <- do.call(rbind, lapply(1:1000, function(x) bi.glm()))
b.mu <- data.frame(mu.int=mean(b.sim$int), mu.b1=mean(b.sim$b1))
b.mu


Y2.bi.glm <- function(x){
  n=1000
  size=1
  X1 <- runif(n, min=0, max=1)
  B0 <- 1
  B1 <- 2
  Y2 <-rbinom(n=n, size=size, prob=(1/(1+exp(-(B0+B1*X1)))))
  int<- summary(glm(Y2~X1, family = binomial))$coefficient[1]
  b1<- summary(glm(Y2~X1, family = binomial))$coefficient[2]
  d <- data.frame(int,b1)
  return(d)
}
Y2b.sim <- do.call(rbind, lapply(1:1000, function(x) Y2.bi.glm()))
Y2b.mu <- data.frame(Y2mu.int=mean(Y2b.sim$int), Y2mu.b1=mean(Y2b.sim$b1))
Y2b.mu


Y3.bi.glm <- function(x){ #NOT GOOD
  n=1000
  size=1
  X1 <- runif(n, min=0, max=1)
  B0 <- 1
  B1 <- 2
  Z <- B0+B1*X1
  Y3 <-rbinom(n=n, size=size, prob=(1/(1+exp(-(Z)))))
  int<- summary(glm(Y3~X1, family = binomial))$coefficient[1]
  b1<- summary(glm(Y3~X1, family = binomial))$coefficient[2]
  d <- data.frame(int,b1)
  return(d)
}
Y3b.sim <- do.call(rbind, lapply(1:1000, function(x) Y3.bi.glm()))
Y3b.mu <- data.frame(Y3mu.int=mean(Y3b.sim$int), Y3mu.b1=mean(Y3b.sim$b1))
Y3b.mu


Y1.bi.glm <- function(x){ 
  n=1000
  size=1
  X1 <- runif(n, min=0, max=1)
  B0 <- 1
  B1 <- 2
  Z <- B0+B1*X1
  PR <- 1/(1+(exp(-Z)))
  Y1 <- rbinom(n, 1, PR)
  int<- summary(glm(Y1~X1, family = binomial))$coefficient[1]
  b1<- summary(glm(Y1~X1, family = binomial))$coefficient[2]
  d <- data.frame(int,b1)
  return(d)
}
Y1b.sim <- do.call(rbind, lapply(1:1000, function(x) Y1.bi.glm()))
Y1b.mu <- data.frame(Y1mu.int=mean(Y1b.sim$int), Y1mu.b1=mean(Y1b.sim$b1))
Y1b.mu



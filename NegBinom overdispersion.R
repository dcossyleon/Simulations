#Feb 16, 2018


#proof of principle, Var=mu+(mu^2)/k
test<-rnbinom(n=10000, size=3, mu=20)
v <- var(test)
m <- mean(test)
k <- m^2/(v-m) #solves for k using formula




h <- t(do.call(cbind(lapply(1:10000, function (x) do.call(rbind(lapply(seq(1,1000,1), function(x) tryCatch(NB(), warning = function(w) NA))))))))
        #why t?



#1) two neg binom distributions, mu1=10, mu2=20
n=5
size=5
mu1=10
mu2=20

NB <- function(x){
  g1 <-rnbinom(n=n, size=size, mu=mu1)
  g2 <- rnbinom(n=n, size=size, mu=mu2)
  df<- data.frame(Y=c(g1,g2), X=c(rep("g1",n), rep("g2", n)))
  pvalue<- summary(glm.nb(Y~X, data=df))$coefficient[8]
  return(pvalue)
}

results <- do.call(rbind,lapply(1:10000, function(x) tryCatch(NB(), warning = function(w) NA)))
t <- table(results<0.05)
t[[2]]/(t[[1]]+t[[2]]) #power (%) ~21.4%



#2) mu1=100, mu2=200
n=5
size=5
mu1=100
mu2=200

NB2 <- function(x){
  g1 <-rnbinom(n=n, size=size, mu=mu1)
  g2 <- rnbinom(n=n, size=size, mu=mu2)
  df<- data.frame(Y=c(g1,g2), X=c(rep("g1",n), rep("g2", n)))
  pvalue<- summary(glm.nb(Y~X, data=df))$coefficient[8]
  return(pvalue)
}

results2 <- do.call(rbind,lapply(1:10000, function(x) tryCatch(NB2(), warning = function(w) NA)))
t2 <- table(results2<0.05)
t2[[2]]/(t2[[1]]+t2[[2]]) #power (%) ~21.9%




#3) mu1=1000, mu2=2000
n=5
size=5
mu1=1000
mu2=2000

NB3 <- function(x){
  g1 <-rnbinom(n=n, size=size, mu=mu1)
  g2 <- rnbinom(n=n, size=size, mu=mu2)
  df<- data.frame(Y=c(g1,g2), X=c(rep("g1",n), rep("g2", n)))
  pvalue<- summary(glm.nb(Y~X, data=df))$coefficient[8]
  return(pvalue)
}

results3 <- do.call(rbind,lapply(1:10000, function(x) tryCatch(NB2(), warning = function(w) NA)))
t3 <- table(results3<0.05)
t3[[2]]/(t3[[1]]+t3[[2]]) #power (%) ~ 21.7%


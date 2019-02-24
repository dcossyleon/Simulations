#May 31, 2017
#Extreme Groups Simulation
  library(MASS)
  library(parallel)
  library(ggplot2)
  
  par(mfrow=c(2,3))
  
  c=0.05
  sigma<-matrix(c(1,c,c,1), nrow=2) #covariance matrix
  size <- c(24,60,120,240,360,480,600,800,1000) #Sample sizes for each "study"
  qsize <- 2*size #bottom+top quartiles still equal same N as no-split
  tsize  <- (3/2)*size
  dsize <- 5*size
  sample <- function(g){
    data=as.data.frame(mvrnorm(g, mu=c(0,10), Sigma=sigma, empirical=F))
    colnames(data) <- c("X", "Y") 
    p <- summary(lm(Y~X, data=data))$coefficient[2,4]
    return(p)
  }

#Run sample() 1000x for each size
  results<- t(do.call(cbind, mclapply(1:1000, function(x) do.call(rbind, lapply(size, function(x) sample(x))), mc.cores = 4, mc.cleanup = TRUE))) #df of 1000 pvalues, columns=diff samples sizes

#Calculate Power
    prop <- function(x){
      pr <- length(x[x<0.05])/1000 #1000--> number of rows
      return(pr)
    }
    power <- apply(results, 2, prop)
    power



#Quartile splits
  Qsample <- function(g){
    data=as.data.frame(mvrnorm(g, mu=c(0,10), Sigma=sigma, empirical=F))
    colnames(data) <- c("X", "Y") 
    quant <- data[data$X <= quantile(data$X, 0.25)| data$X>= quantile(data$X, 0.75), ]
    p <- summary(lm(Y~X, data=quant))$coefficient[2,4]
    return(p)
  }

  Qresults<- t(do.call(cbind, mclapply(1:1000, function(x) do.call(rbind, lapply(qsize, function(x) Qsample(x))), mc.cores = 4, mc.cleanup = TRUE)))
  
  Qpwr <- apply(Qresults, 2, prop)
  Qpwr
  
  
#Tertile splits
  Tsample <- function(g){
    data=as.data.frame(mvrnorm(g, mu=c(0,10), Sigma=sigma, empirical=F))
    colnames(data) <- c("X", "Y") 
    quant <- data[data$X <= quantile(data$X, 1/3)| data$X>= quantile(data$X, 2/3), ]
    p <- summary(lm(Y~X, data=quant))$coefficient[2,4]
    return(p)
  }
  
  Tresults<- t(do.call(cbind, mclapply(1:1000, function(x) do.call(rbind, lapply(tsize, function(x) Tsample(x))), mc.cores = 4, mc.cleanup = TRUE)))
  
  Tpwr <- apply(Tresults, 2, prop)
  Tpwr
  
  
#Decile splits
  Dsample <- function(g){
    data=as.data.frame(mvrnorm(g, mu=c(0,10), Sigma=sigma, empirical=F))
    colnames(data) <- c("X", "Y") 
    quant <- data[data$X <= quantile(data$X, .1)| data$X>= quantile(data$X, .9), ]
    p <- summary(lm(Y~X, data=quant))$coefficient[2,4]
    return(p)
  }
  
  Dresults<- t(do.call(cbind, mclapply(1:1000, function(x) do.call(rbind, lapply(dsize, function(x) Dsample(x))), mc.cores = 4, mc.cleanup = TRUE)))
  
  Dpwr <- apply(Dresults, 2, prop)
  Dpwr
  
#Plots
  plot(power~size, ylim=c(0.0,1.0), ylab="Power", xlab="Sample Size", xlim=c(0,1000), pch=16) #no split
  par(new=TRUE)
  plot(Qpwr~size, ylim=c(0.0,1.0), xlim=c(0,1000), ylab="", xlab="", col="blue", pch=24) #quartile split
  par(new=TRUE)
  plot(Tpwr~size, ylim=c(0.0,1.0), xlim=c(0,1000), ylab="", xlab="", col="red") #tertile split
  par(new=TRUE)
  plot(Dpwr~size, ylim=c(0.0,1.0), xlim=c(0,1000), ylab="", xlab="", col="green") #decile split
  par(new=TRUE)
  abline(h = 0.80)
  par(new=FALSE)
  
  
  
#1) to reach 80% power. what are the different samples sizes you would need in both cases. the "traditional" way

#2) with this sample size, what's the difference in power?
  #- 15, 15
  #- more relevant to the Karen Parker paper
  

  
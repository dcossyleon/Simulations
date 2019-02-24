#Feb 19, 2018: see if you can automate the power thing to automatically stop when it's powered at 80%

#1) two neg binom distributions, mu1=10, mu2=20
size=5
mu1=10
mu2=20

NB <- function(n){
  g1 <-rnbinom(n=n, size=size, mu=mu1)
  g2 <- rnbinom(n=n, size=size, mu=mu2)
  df<- data.frame(Y=c(g1,g2), X=c(rep("g1",n), rep("g2", n)))
  pvalue<- summary(glm.nb(Y~X, data=df))$coefficient[8]
  return(pvalue)
}





g <- t(do.call(cbind, lapply(1:20, function (x) do.call(rbind, lapply(seq(8,10,1), function(x) tryCatch(NB(x), warning=function(w) NA))))))
        #why t?


results <- do.call(rbind,lapply(1:10000, function(x) tryCatch(NB(), warning = function(w) NA)))
t <- table(results<0.05)
power <- t[[2]]/(t[[1]]+t[[2]]) #power (%) ~21.4%


library(dplyr)

big.function <- function(x){
  
}










#Sequence until t <- table(g<0.05) and 
  ifelse(power>=0.795, return(n), "x")
  
  
  power.reached <- function(x) {
    repeat {
      # do something
      i <- sample(nrow(df), 1)
      x <- df[sample(nrow(df), 1), ]
      # exit if the condition is met
      if (x$SCORE > 0) break
    }
    return(x)
  }
  
  
  df <- data.frame(NAME=c(rep('Frank',10),rep('Mary',10)), SCORE=rnorm(20))
  df
  
  random.sample <- function(x) {
    x <- df[sample(nrow(df), 1), ]
    if (x$SCORE > 0) return(x)
    #if (x$SCORE <= 0) run the function again
  }
  
  random.sample(df)
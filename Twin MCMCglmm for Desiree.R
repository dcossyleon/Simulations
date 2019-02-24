library(reshape2)
library(MCMCglmm)
library(MCMCvis)
library(MASS)

mzcor <- 0.30
dzcor <- 0.15
mz <- as.data.frame(mvrnorm(100, mu=c(0, 0), Sigma=matrix(c(1,mzcor,mzcor,1), nrow=2), empirical=T)) #n is # of families
dz <- as.data.frame(mvrnorm(100, mu=c(0, 0), Sigma=matrix(c(1,dzcor,dzcor,1), nrow=2), empirical=T))

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

priorAE <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
  priorAE.v1 <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=0.4, nu=0.002)))
  priorAE.v04 <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=.2)))
  
n1 <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.0002)))
n2 <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.2)))
n3 <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=2.0)))
n4 <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=20)))
  n4.v1<- list(R=list(V=1, nu=0.002), G=list(G1=list(V=0.4, nu=20)))
n5 <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=2000)))
  n5.v1 <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=0.4, nu=2000)))
n.TEST <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=0.4, nu=2000)))
  


  fitAE <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=priorAE, burnin=10, nitt=150, thin=1)
  #plot(fitAE$VCV)
  #data.frame('a2'=mean(fitAE$VCV[,'id']/(fitAE$VCV[,'id']+fitAE$VCV[, 'units'])), 'e2'=mean(fitAE$VCV[,'units']/(fitAE$VCV[,'id']+fitAE$VCV[, 'units'])))
  h2.fitAE<- (fitAE$VCV[,'id'])/(fitAE$VCV[,'id']+fitAE$VCV[,"units"])
  HPDinterval(h2.fitAE) #
  plot(h2.fitAE)


#Low H2: G nu=0.002, V=1
  #low1 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=priorAE, burnin=10000, nitt=100000, thin=100)
  (df.low1<- data.frame('a2'=mean(low1$VCV[,'id']/(low1$VCV[,'id']+low1$VCV[, 'units'])), 'e2'=mean(low1$VCV[,'units']/(low1$VCV[,'id']+low1$VCV[, 'units']))))
  h2.low1<- (low1$VCV[,'id'])/(low1$VCV[,'id']+low1$VCV[,"units"])
  HPDinterval(h2.low1)#0.0002637536 0.1255228
  HPDinterval(low1$VCV)
  plot(h2.low1, xlim=c(0,1), ylim=c(0,50), trace = FALSE, sub="h2.low1")

#Low H2: G nu=2, V=0.04
  low2 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n.TEST, burnin=10000, nitt=100000, thin=100)
  plot(low2$VCV, sub="low2")
  (df.low2<- data.frame('a2'=mean(low2$VCV[,'id']/(low2$VCV[,'id']+low2$VCV[, 'units'])), 'e2'=mean(low2$VCV[,'units']/(low2$VCV[,'id']+low2$VCV[, 'units']))))
  h2.low2<- (low2$VCV[,'id'])/(low2$VCV[,'id']+low2$VCV[,"units"])
  HPDinterval(h2.low2) #
  plot(h2.low2, xlim=c(0,1), ylim=c(0,50), trace = FALSE, sub="h2.low2")  
  

#COMPARE
#R G nu=0.002, V=1 for both
#Changing N: 100
  #fit1 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=priorAE, burnin=10000, nitt=100000, thin=100)
    plot(fit1$VCV, sub="fit1")
    (df1<- data.frame('a2'=mean(fit1$VCV[,'id']/(fit1$VCV[,'id']+fit1$VCV[, 'units'])), 
                      'e2'=mean(fit1$VCV[,'units']/(fit1$VCV[,'id']+fit1$VCV[, 'units']))))
    h2.fit1<- (fit1$VCV[,'id'])/(fit1$VCV[,'id']+fit1$VCV[,"units"])
    e2.fit1<- (fit1$VCV[,'units'])/(fit1$VCV[,'id']+fit1$VCV[,"units"])
    HPDinterval(h2.fit1) #23-52
    plot(h2.fit1, xlim=c(0,1), ylim=c(0,10),trace = FALSE, sub="h2.fit1")
    plot(e2.fit1, xlim=c(0,1), ylim=c(0,10),trace = FALSE, sub="e2.fit1")

          #Changing G nu: 0.0002, leaving R
          #fit1.a <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n1, burnin=10000, nitt=100000, thin=100)
          plot(fit1.a$VCV, sub="fit1.a")
          (df1.a<- data.frame('a2'=mean(fit1.a$VCV[,'id']/(fit1.a$VCV[,'id']+fit1.a$VCV[, 'units'])), 
                            'e2'=mean(fit1.a$VCV[,'units']/(fit1.a$VCV[,'id']+fit1.a$VCV[, 'units']))))
          h2.fit1.a<- (fit1.a$VCV[,'id'])/(fit1.a$VCV[,'id']+fit1.a$VCV[,"units"])
          e2.fit1.a<- (fit1.a$VCV[,'units'])/(fit1.a$VCV[,'id']+fit1.a$VCV[,"units"])
          HPDinterval(h2.fit1.a) #24.5-54.1
          par(mfrow=c(2,1))
          plot(h2.fit1.a, xlim=c(0,1), ylim=c(0,10),trace = FALSE, sub="h2.fit1.a")
          plot(e2.fit1.a, xlim=c(0,1),ylim=c(0,10), trace = FALSE, sub="e2.fit1.a")
          
                #Changing G nu: 0.2, leaving R
                #fit1.a1 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n2, burnin=10000, nitt=100000, thin=100)
                plot(fit1.a1$VCV, sub="fit1.a1")
                (df1.a1<- data.frame('a2'=mean(fit1.a1$VCV[,'id']/(fit1.a1$VCV[,'id']+fit1.a1$VCV[, 'units'])), 
                                     'e2'=mean(fit1.a1$VCV[,'units']/(fit1.a1$VCV[,'id']+fit1.a1$VCV[, 'units']))))
                h2.fit1.a1<- (fit1.a1$VCV[,'id'])/(fit1.a1$VCV[,'id']+fit1.a1$VCV[,"units"])
                e2.fit1.a1<- (fit1.a1$VCV[,'units'])/(fit1.a1$VCV[,'id']+fit1.a1$VCV[,"units"])
                HPDinterval(h2.fit1.a1) #22-53
                plot(h2.fit1.a1, xlim=c(0,1),ylim=c(0,10), trace = FALSE, sub="h2.fit1.a1")
                plot(e2.fit1.a1, xlim=c(0,1),ylim=c(0,10), trace = FALSE, sub="e2.fit1.a1")
                abline(v=0.4)
                
                      #Changing G nu: 2.0, leaving R
                      #fit1.a2 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n3, burnin=10000, nitt=100000, thin=100)
                      plot(fit1.a2$VCV, sub="fit1.a2")
                      (df1.a2<- data.frame('a2'=mean(fit1.a2$VCV[,'id']/(fit1.a2$VCV[,'id']+fit1.a2$VCV[, 'units'])), 
                                           'e2'=mean(fit1.a2$VCV[,'units']/(fit1.a2$VCV[,'id']+fit1.a2$VCV[, 'units']))))
                      h2.fit1.a2<- (fit1.a2$VCV[,'id'])/(fit1.a2$VCV[,'id']+fit1.a2$VCV[,"units"])
                      e2.fit1.a2<- (fit1.a2$VCV[,'units'])/(fit1.a2$VCV[,'id']+fit1.a2$VCV[,"units"])
                      HPDinterval(h2.fit1.a2) #28-54
                      plot(h2.fit1.a2, xlim=c(0,1),ylim=c(0,10), trace = FALSE, sub="h2.fit1.a2")
                      plot(e2.fit1.a2, xlim=c(0,1), ylim=c(0,10),trace = FALSE, sub="e2.fit1.a2")
                      
                            #Changing G nu: 20, leaving R
                           # fit1.a3 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n4, burnin=10000, nitt=100000, thin=100)
                            plot(fit1.a3$VCV, sub="fit1.a3")
                            (df1.a3<- data.frame('a2'=mean(fit1.a3$VCV[,'id']/(fit1.a3$VCV[,'id']+fit1.a3$VCV[, 'units'])), 
                                                 'e2'=mean(fit1.a3$VCV[,'units']/(fit1.a3$VCV[,'id']+fit1.a3$VCV[, 'units']))))
                            h2.fit1.a3<- (fit1.a3$VCV[,'id'])/(fit1.a3$VCV[,'id']+fit1.a3$VCV[,"units"])
                            e2.fit1.a3<- (fit1.a3$VCV[,'units'])/(fit1.a3$VCV[,'id']+fit1.a3$VCV[,"units"])
                            HPDinterval(h2.fit1.a3) #40-60
                            plot(h2.fit1.a3, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit1.a3")
                            plot(e2.fit1.a3, xlim=c(0,1), ylim=c(0,10),trace = FALSE, sub="e2.fit1.a3")
                            
                                  #Changing G nu: 2000, leaving R
                                  #fit1.a4 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n5, burnin=10000, nitt=100000, thin=100)
                                  plot(fit1.a4$VCV, sub="fit1.a4")
                                  (df1.a4<- data.frame('a2'=mean(fit1.a4$VCV[,'id']/(fit1.a4$VCV[,'id']+fit1.a4$VCV[, 'units'])), 
                                                       'e2'=mean(fit1.a4$VCV[,'units']/(fit1.a4$VCV[,'id']+fit1.a4$VCV[, 'units']))))
                                  h2.fit1.a4<- (fit1.a4$VCV[,'id'])/(fit1.a4$VCV[,'id']+fit1.a4$VCV[,"units"])
                                  e2.fit1.a4<- (fit1.a4$VCV[,'units'])/(fit1.a4$VCV[,'id']+fit1.a4$VCV[,"units"])
                                  HPDinterval(h2.fit1.a4) #62-72
                                  plot(h2.fit1.a4, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit1.a4")
                                  plot(e2.fit1.a4, xlim=c(0,1), ylim=c(0,10),trace = FALSE, sub="e2.fit1.a4")
                                  
                                  #Changing G nu: 2000, making G V= correct h2
                                  #fit1.a4.v1 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n5.v1, burnin=10000, nitt=100000, thin=100)
                                  plot(fit1.a4.v1$VCV, sub="fit1.a4.v1")
                                  (df1.a4.v1<- data.frame('a2'=mean(fit1.a4.v1$VCV[,'id']/(fit1.a4.v1$VCV[,'id']+fit1.a4.v1$VCV[, 'units'])), 
                                                       'e2'=mean(fit1.a4.v1$VCV[,'units']/(fit1.a4.v1$VCV[,'id']+fit1.a4.v1$VCV[, 'units']))))
                                  h2.fit1.a4.v1<- (fit1.a4.v1$VCV[,'id'])/(fit1.a4.v1$VCV[,'id']+fit1.a4.v1$VCV[,"units"])
                                  e2.fit1.a4.v1<- (fit1.a4.v1$VCV[,'units'])/(fit1.a4.v1$VCV[,'id']+fit1.a4.v1$VCV[,"units"])
                                  HPDinterval(h2.fit1.a4.v1) #34-45
                                  plot(h2.fit1.a4.v1, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit1.a4.v1")
                                  plot(e2.fit1.a4.v1, xlim=c(0,1), ylim=c(0,10),trace = FALSE, sub="e2.fit1.a4.v1")
                                  
                                      #Changing anything
                                      fit1.TEST <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n.TEST, burnin=10000, nitt=100000, thin=100)
                                      plot(fit1.a4.v1$VCV, sub="fit1.a4.v1")
                                      (df1.a4.v1<- data.frame('a2'=mean(fit1.TEST$VCV[,'id']/(fit1.TEST$VCV[,'id']+fit1.TEST$VCV[, 'units'])), 
                                                              'e2'=mean(fit1.TEST$VCV[,'units']/(fit1.TEST$VCV[,'id']+fit1.TEST$VCV[, 'units']))))
                                      h2.fit1.TEST<- (fit1.TEST$VCV[,'id'])/(fit1.TEST$VCV[,'id']+fit1.TEST$VCV[,"units"])
                                      HPDinterval(h2.fit1.TEST)
                                      plot(h2.fit1.TEST, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit1.TEST")

#N: 50  
  #fit2 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=priorAE, burnin=10000, nitt=300000, thin=100) #nu=0.002
    (df2<- data.frame('a2'=mean(fit2$VCV[,'id']/(fit2$VCV[,'id']+fit2$VCV[, 'units'])), 'e2'=mean(fit2$VCV[,'units']/(fit2$VCV[,'id']+fit2$VCV[, 'units']))))
    h2.fit2<- (fit2$VCV[,'id'])/(fit2$VCV[,'id']+fit2$VCV[,"units"])
    HPDinterval(h2.fit2) #12.4-59.6
    plot(h2.fit2, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit2")

      #Changing G nu: 0.0002, leaving R
      #fit2.a <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n1, burnin=10000, nitt=300000, thin=100)
      (df2.a<- data.frame('a2'=mean(fit2.a$VCV[,'id']/(fit2.a$VCV[,'id']+fit2.a$VCV[, 'units'])), 'e2'=mean(fit2.a$VCV[,'units']/(fit2.a$VCV[,'id']+fit2.a$VCV[,'units']))))
      h2.fit2.a<- (fit2.a$VCV[,'id'])/(fit2.a$VCV[,'id']+fit2.a$VCV[,"units"])
      HPDinterval(h2.fit2.a) #11.98-59.1
      plot(h2.fit2.a, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit2.a")

          #Changing G nu: 0.2, leaving R
          #fit2.a1 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n2, burnin=10000, nitt=300000, thin=100)
          (df2.a1<- data.frame('a2'=mean(fit2.a1$VCV[,'id']/(fit2.a1$VCV[,'id']+fit2.a1$VCV[,'units'])),'e2'=mean(fit2.a1$VCV[,'units']/(fit2.a1$VCV[,'id']+fit2.a1$VCV[,'units']))))
          h2.fit2.a1<- (fit2.a1$VCV[,'id'])/(fit2.a1$VCV[,'id']+fit2.a1$VCV[,"units"])
          HPDinterval(h2.fit2.a1) #15-57
          plot(h2.fit2.a1, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit2.a1")
          
              #Changing G nu: 2000, leaving R
              #fit2.a5 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n5, burnin=10000, nitt=300000, thin=100)
              (df2.a1<- data.frame('a2'=mean(fit2.a5$VCV[,'id']/(fit2.a5$VCV[,'id']+fit2.a5$VCV[,'units'])),'e2'=mean(fit2.a5$VCV[,'units']/(fit2.a5$VCV[,'id']+fit2.a5$VCV[,'units']))))
              h2.fit2.a5<- (fit2.a5$VCV[,'id'])/(fit2.a5$VCV[,'id']+fit2.a5$VCV[,"units"])
              HPDinterval(h2.fit2.a5) #61-74
              plot(h2.fit2.a5, xlim=c(0,1), ylim=c(0,20), trace = FALSE, sub="h2.fit2.a5")
    
#N: 25
 # fit25 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=priorAE, burnin=10000, nitt=100000, thin=100)
  (df25<- data.frame('a2'=mean(fit25$VCV[,'id']/(fit25$VCV[,'id']+fit25$VCV[, 'units'])), 'e2'=mean(fit25$VCV[,'units']/(fit25$VCV[,'id']+fit25$VCV[, 'units']))))
  h2.fit25<- (fit25$VCV[,'id'])/(fit25$VCV[,'id']+fit25$VCV[,"units"])
  HPDinterval(h2.fit25) #0-58
  plot(h2.fit25, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit25")
      
      #nu=0.002, V=0.4
  #fit25a <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=priorAE.v1, burnin=10000, nitt=100000, thin=100)
      (df25.v1<- data.frame('a2'=mean(fit25a$VCV[,'id']/(fit25a$VCV[,'id']+fit25a$VCV[, 'units'])), 'e2'=mean(fit25a$VCV[,'units']/(fit25a$VCV[,'id']+fit25a$VCV[, 'units']))))
      h2.fit25.v1<- (fit25a$VCV[,'id'])/(fit25a$VCV[,'id']+fit25a$VCV[,"units"])
      HPDinterval(h2.fit25.v1) #0-58
      plot(h2.fit25.v1, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit25.v1")
      
          #nu=0.2, V=0.5
          fit25.test <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n.TEST, burnin=10000, nitt=100000, thin=100)
          plot(fit25.test$VCV, sub="fit25.test")
          (df25.test<- data.frame('a2'=mean(fit25.test$VCV[,'id']/(fit25.test$VCV[,'id']+fit25.test$VCV[, 'units'])), 'e2'=mean(fit25.test$VCV[,'units']/(fit25.test$VCV[,'id']+fit25.test$VCV[, 'units']))))
          h2.fit25.test<- (fit25.test$VCV[,'id'])/(fit25.test$VCV[,'id']+fit25.test$VCV[,"units"])
          #HPDinterval(h2.fit25.test) #0.02-60, v=0.5
          #HPDinterval(h2.fit25.test) #0.056-0.61, v=0.4
          HPDinterval(h2.fit25.test) #0.051-0.62, v=0.3
          plot(h2.fit25.test, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit25.test")
  
              
#N: 10  
  #fit3 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=priorAE, burnin=10000, nitt=100000, thin=100)
  plot(fit3$VCV, sub="fit3")
  (df3<- data.frame('a2'=mean(fit3$VCV[,'id']/(fit3$VCV[,'id']+fit3$VCV[, 'units'])),'e2'=mean(fit3$VCV[,'units']/(fit3$VCV[,'id']+fit3$VCV[, 'units']))))
  h2.fit3<- (fit3$VCV[,'id'])/(fit3$VCV[,'id']+fit3$VCV[,"units"])
  HPDinterval(h2.fit3) #.0005-0.66
  plot(h2.fit3, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit3")
  
    #fit3.v1 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=priorAE.v1, burnin=10000, nitt=100000, thin=100)
    plot(fit3.v1$VCV, sub="fit3.v1")
    (df3.v1<- data.frame('a2'=mean(fit3.v1$VCV[,'id']/(fit3.v1$VCV[,'id']+fit3.v1$VCV[, 'units'])),'e2'=mean(fit3.v1$VCV[,'units']/(fit3.v1$VCV[,'id']+fit3.v1$VCV[, 'units']))))
    h2.fit3.v1<- (fit3.v1$VCV[,'id'])/(fit3.v1$VCV[,'id']+fit3.v1$VCV[,"units"])
    HPDinterval(h2.fit3.v1) #.0005-0.66
    plot(h2.fit3.v1, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit3.v1")
  
      #Changing G nu: 0.0002, leaving R
      #fit3.a <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n1, burnin=10000, nitt=100000, thin=100)
      plot(fit3.a$VCV, sub="fit3.a")
      (df3.a<- data.frame('a2'=mean(fit3.a$VCV[,'id']/(fit3.a$VCV[,'id']+fit3.a$VCV[, 'units'])),'e2'=mean(fit3.a$VCV[,'units']/(fit3.a$VCV[,'id']+fit3.a$VCV[, 'units']))))
      h2.fit3.a<- (fit3.a$VCV[,'id'])/(fit3.a$VCV[,'id']+fit3.a$VCV[,"units"])
      HPDinterval(h2.fit3.a) #.0005-0.66
      plot(h2.fit3.a, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit3.a")
      
      
        #Changing G nu: 20, leaving R
        #fit3.a3 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n4, burnin=10000, nitt=100000, thin=100)
        plot(fit3.a3$VCV, sub="fit3.a3")
        (df3.a3<- data.frame('a2'=mean(fit3.a3$VCV[,'id']/(fit3.a3$VCV[,'id']+fit3.a3$VCV[, 'units'])),'e2'=mean(fit3.a3$VCV[,'units']/(fit3.a3$VCV[,'id']+fit3.a3$VCV[, 'units']))))
        h2.fit3.a3<- (fit3.a3$VCV[,'id'])/(fit3.a3$VCV[,'id']+fit3.a3$VCV[,"units"])
        HPDinterval(h2.fit3.a3) #.0005-0.66
        plot(h2.fit3.a3, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit3.a3")
        
            #Changing G nu: 20 with V=0.4, leaving R
            #fit3.a3.v1 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n4.v1, burnin=10000, nitt=100000, thin=100)
            plot(fit3.a3.v1$VCV, sub="fit3.a3.v1")
            (df3.a3.v1<- data.frame('a2'=mean(fit3.a3.v1$VCV[,'id']/(fit3.a3.v1$VCV[,'id']+fit3.a3.v1$VCV[, 'units'])),
                                'e2'=mean(fit3.a3.v1$VCV[,'units']/(fit3.a3.v1$VCV[,'id']+fit3.a3.v1$VCV[, 'units']))))
            h2.fit3.a3.v1<- (fit3.a3.v1$VCV[,'id'])/(fit3.a3.v1$VCV[,'id']+fit3.a3.v1$VCV[,"units"])
            HPDinterval(h2.fit3.a3.v1) #20-64
            plot(h2.fit3.a3.v1, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit3.a3.v1")
    
#N: 5 
    #fit4 <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=priorAE, burnin=10000, nitt=100000, thin=100)
    plot(fit4$VCV, sub="fit4")
    (df4<- data.frame('a2'=mean(fit4$VCV[,'id']/(fit4$VCV[,'id']+fit4$VCV[, 'units'])), 'e2'=mean(fit4$VCV[,'units']/(fit4$VCV[,'id']+fit4$VCV[, 'units']))))
    h2.fit4<- (fit4$VCV[,'id'])/(fit4$VCV[,'id']+fit4$VCV[,"units"])
    plot(h2.fit4, xlim=c(0,1), ylim=c(0,10), trace = FALSE, sub="h2.fit4")
    
            #Changing G nu: 0.0002, leaving R
            #fit4.a <- MCMCglmm(value~1, random=~id, ginverse=list(id=Ginv), data=tw, prior=n1, burnin=10000, nitt=100000, thin=100)
            plot(fit4.a$VCV, sub="fit4.a")
            (df4.a<- data.frame('a2'=mean(fit4.a$VCV[,'id']/(fit4.a$VCV[,'id']+fit4.a$VCV[, 'units'])),'e2'=mean(fit4.a$VCV[,'units']/(fit4.a$VCV[,'id']+fit4.a$VCV[, 'units'])))) 
                  
    
    
    
    

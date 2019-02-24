#Aug 23, 2018
#Ordinal MCMCglmm

library(MCMCglmm)
library(pedigreemm)
library(parallel)


#Priors
prior <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
C1prior <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002), G2 = list(V = 1, nu = 0.002)))


#for ordinal traits
ph$ord.tn <- round((ph$tn/ph$Time)*60,0) #counts/hour
ph$ord.sc <- round((ph$sc/ph$Time)*60,0) #counts/hour
beh.ord <- c("ord.sc", "ord.tn")

ph.2x$ord.tn <- round((ph.2x$tn/ph.2x$Time)*60,0) #counts/hour
ph.2x$ord.sc <- round((ph.2x$sc/ph.2x$Time)*60,0) #counts/hour

#MCMCglmm N=214
    start.time <- Sys.time()
    mcmc.ord <- mclapply(beh.ord, function(x) MCMCglmm(as.formula(paste0(x,"~1")), random=~animal, data=ph, family="ordinal", prior=prior, pedigree=pedigree, nitt=10000000, burnin=1000000, thin=10000), mc.cores=4) 
    end.time <- Sys.time()
    end.time-start.time

  mclapply(1:2, FUN=function(x) autocorr.diag(mcmc.ord[[x]]$VCV), mc.cores=4) #autocorr is horrible. 

ord.sc <- mcmc.ord[[1]]
ord.tn <- mcmc.ord[[2]]

#H2 
    h2.tn.ord <- (ord.tn$VCV[,'animal'])/(ord.tn$VCV[,'animal']+ord.tn$VCV[,"units"])
    plot(h2.tn.ord)
    hist(h2.tn.ord, sub="nu= regular, n=214") #getting back priors
    h2.sc.ord <- (ord.sc$VCV[,'animal'])/(ord.sc$VCV[,'animal']+ord.sc$VCV[,"units"])
    plot(h2.sc.ord)
    hist(h2.sc.ord, sub="nu= regular, n=214")
    

#increasing N by 2x
  animal2 <- paste(pedigree$animal,2)
  sire2 <- ifelse(pedigree$sire != "<NA>", paste(pedigree$sire,2), "<NA>")
  dam2 <- ifelse(pedigree$dam != "<NA>", paste(pedigree$dam,2), "<NA>")
    pedigree2 <- editPed(animal2, sire = sire2, dam = dam2)[,1:3]
      colnames(pedigree2) <- c("animal", "sire", "dam")
      
  pedigree.2x <- rbind(pedigree, pedigree2)
    
  ph2 <- cbind("animal"=paste(ph$animal,2), ph[,c(1,3:33)]) #ord.tn=32, ord.sc=33
  ph.2x <- rbind(ph, ph2)
  

#MCMCglmm N= 428  
  start.time <- Sys.time()
  mcmc.ord.2x <- mclapply(beh.ord, function(x) MCMCglmm(as.formula(paste0(x,"~1")), random=~animal, data=ph.2x, family="ordinal", prior=prior, pedigree=pedigree.2x, nitt=10000000, burnin=700000, thin=9000), mc.cores=4) 
  end.time <- Sys.time()
  end.time-start.time
  
  mclapply(1:2, FUN=function(x) autocorr.diag(mcmc.ord.2x[[x]]$VCV), mc.cores=4) #autocorr

  ord.sc.2x <- mcmc.ord.2x[[1]]
  ord.tn.2x <- mcmc.ord.2x[[2]]
  
    #H2: SC
      h2.sc.ord.2x <- (ord.sc.2x$VCV[,'animal'])/(ord.sc.2x$VCV[,'animal']+ord.sc.2x$VCV[,"units"])
      plot(h2.sc.ord.2x, main="h2 scratch ordinal, regular nu", sub="ordinal, n=428, R&G nu=0.002") #increasing burnin to 1mill may help
      hist(h2.sc.ord.2x)
      effectiveSize(h2.sc.ord.2x) #20
        mean(h2.sc.ord.2x) #21
        HPDinterval(h2.sc.ord.2x) #0-53
    #H2: TN
      h2.tn.ord.2x <- (ord.tn.2x$VCV[,'animal'])/(ord.tn.2x$VCV[,'animal']+ord.tn.2x$VCV[,"units"])
      plot(h2.tn.ord.2x)
      hist(h2.tn.ord.2x)
      effectiveSize(h2.tn.ord.2x) #608
        mean(h2.tn.ord.2x) #30
        HPDinterval(h2.tn.ord.2x) #0, 0.62
      
    
      #same, but with stronger priors (for both R and G)
        start.time <- Sys.time()
        m2 <- mclapply(beh.ord, function(x) MCMCglmm(as.formula(paste0(x,"~1")), random=~animal, data=ph.2x, family="ordinal", pedigree=pedigree.2x, nitt=9000000, burnin=10000, thin=5000,
                                                     prior=list(R=list(V=1, nu=0.02), G=list(G1=list(V=1, nu=0.02)))), mc.cores=3) 
        end.time <- Sys.time()
        end.time-start.time
        
        mclapply(1:2, FUN=function(x) autocorr.diag(m2[[x]]$VCV), mc.cores=4) #autocorr
        
        m2.sc <- m2[[1]]
        m2.tn <- m2[[2]]
        
        #H2
        h2.m2.sc <- (m2.sc$VCV[,'animal'])/(m2.sc$VCV[,'animal']+m2.sc$VCV[,"units"])
        plot(h2.m2.sc, main="h2 scratch ordinal, stronger nu", sub="ordinal, n=428, R&G nu=0.02")
        hist(h2.m2.sc, sub="n=428, RG nu=0.02", 30)
        effectiveSize(h2.m2.sc) #24
        mean(h2.m2.sc) #23
        HPDinterval(h2.m2.sc) #0-67
        
        h2.m2.tn <- (m2.tn$VCV[,'animal'])/(m2.tn$VCV[,'animal']+m2.tn$VCV[,"units"])
        plot(h2.m2.tn)
        hist(h2.m2.tn, sub="n=428, RG nu=0.02")
        effectiveSize(h2.m2.tn) #85 ---> very diff from w/ weak prior (previous run)
        mean(h2.m2.tn) #29.5, similar to weak prior
        HPDinterval(h2.m2.tn) #0, 0.62
    

#increasing N by 4x
        animal4 <- paste(pedigree$animal,4)
        sire4 <- ifelse(pedigree$sire != "<NA>", paste(pedigree$sire,4), "<NA>")
        dam4 <- ifelse(pedigree$dam != "<NA>", paste(pedigree$dam,4), "<NA>")
        pedigree4 <- editPed(animal4, sire = sire4, dam = dam4)[,1:3]
        colnames(pedigree4) <- c("animal", "sire", "dam")
        
        pedigree.4x <- rbind(pedigree.2x, pedigree2)
        
        ph4 <- cbind("animal"=paste(ph$animal,4), ph[,c(1,3:33)]) #ord.tn=32, ord.sc=33
        ph.4x <- rbind(ph, ph4)
        
        
        start.time <- Sys.time()
        mcmc.ord.4x <- mclapply(beh.ord, function(x) MCMCglmm(as.formula(paste0(x,"~1")), random=~animal, data=ph.4x, family="ordinal", prior=prior, pedigree=pedigree.4x, nitt=10000000, burnin=700000, thin=9000), mc.cores=4) 
        end.time <- Sys.time()
        end.time-start.time
        
        
        start.time <- Sys.time()
        mcmc.ord.4x.1 <- mclapply(beh.ord, function(x) MCMCglmm(as.formula(paste0(x,"~1")), random=~animal, data=ph.4x, family="ordinal", pedigree=pedigree.4x, nitt=10000000, burnin=700000, thin=9000, prior=list(R=list(V=1, nu=0.02), G=list(G1=list(V=1, nu=0.02)))), mc.cores=4) 
        end.time <- Sys.time()
        end.time-start.time

        
        start.time <- Sys.time()
        mcmc.ord.4x.1 <- mclapply(beh.ord, function(x) MCMCglmm(as.formula(paste0(x,"~1")), random=~animal, data=ph.4x, family="ordinal", pedigree=pedigree.4x, nitt=10000000, burnin=700000, thin=9000, prior=list(R=list(V=1, nu=0.02), G=list(G1=list(V=1, nu=0.02)))), mc.cores=4) 
        end.time <- Sys.time()
        end.time-start.time

library(MASS)
library(tibble)
#FULL SIBS
   
#figure out the 3 member fam
#simulations: with larger families, setting sigma=1, empirical T: end goal: hist var that is a single bar. 

#reduced variance in the family as ICC

#take var of the "does it make sense" (plot against kinship  values), ins
  #how is correlation actually calculated



# Matrix > 3
fam5 <- matrix(c(1.0, 0.5, 0.5, 0.5, 0.5,
                 0.5, 1.0, 0.5, 0.5, 0.5,
                 0.5, 0.5, 1.0, 0.5, 0.5,
                 0.5, 0.5, 0.5, 1.0, 0.5,
                 0.5, 0.5, 0.5, 0.5, 1.0), nrow=5)

        t(chol(fam5))

dz5 <- matrix(c(1.00001, 1.0, 1.0, 1.0, 1.0,
                 1.0, 1.00001, 1.0, 1.0, 1.0,
                 1.0, 1.0, 1.00001, 1.0, 1.0,
                 1.0, 1.0, 1.0, 1.00001, 1.0,
                 1.0, 1.0, 1.0, 1.0, 1.00001), nrow=5) 
        t(chol(dz5))



fam1000 <- matrix(rep(0.5,1e+06), nrow=1000) 
diag(fam1000) <- 1
      f <- as.tibble (t(chol(fam1000)))
      fam1000.ch <- t(chol(fam1000))
      hist(diag(fam1000.ch), df)
      


#FULL
   sib.cor=0.5
       vars <- function(x) {
         sib <- matrix(c(1, sib.cor,
                         sib.cor, 1), nrow=2)
         Y <- mvrnorm(2, mu=0, Sigma=1, empirical=F)
         Y.ch <- t(chol(sib))%*%Y
         var.Y <- var(Y)
         var.Y.ch <- var(Y.ch)
         df <- data.frame(var.Y, var.Y.ch)
         return(df)
       }
      vars.full <- do.call(rbind, lapply(1:10000, function(x) vars(x)))
      hist(vars.full$var.Y, xlim=c(0,15),  ylim=(c(0,1000)))
      hist(vars.full$var.Y.ch, xlim=c(0,15),  ylim=(c(0,1000)))
      mean(vars.full$var.Y)
      mean(vars.full$var.Y.ch)
      
      
#HALF SIBS
      sib.cor=0.125
      vars <- function(x) {
        sib <- matrix(c(1, sib.cor,
                        sib.cor, 1), nrow=2)
        Y <- mvrnorm(2, mu=0, Sigma=1, empirical=F)
        Y.ch <- t(chol(sib))%*%Y
        var.Y <- var(Y)
        var.Y.ch <- var(Y.ch)
        df <- data.frame(var.Y, var.Y.ch)
        return(df)
      }
      #half
      vars.half <- do.call(rbind, lapply(1:10000, function(x) vars(x)))
      hist(vars.half$var.Y,  xlim=c(0,15), ylim=c(0,1000))
      hist(vars.half$var.Y.ch, xlim=c(0,15), ylim=c(0,1000))
      mean(vars.half$var.Y)
      mean(vars.half$var.Y.ch)

#Family of 4
  fam <- matrix(c(1.0, 0.0, 0.25, 0.25, 
                  0.0, 1.0, 0.25, 0.25,
                  0.25, 0.25, 1.0, 0.25,
                  0.25, 0.25, 0.25, 1.0), nrow=4)
  t(chol(fam))
  Z <- matrix(c(2,3,3,3), nrow=4)
  t(chol(fam))%*%Z


#___Family of 3________
#EMPIRICAL = T
vars3 <- function(x) {
  fam3 <- matrix(c(1.0, 0.5, 0.5, 
                   0.5, 1.0, 0.5,
                   0.5, 0.5, 1.0), nrow=3)
  Y3 <- mvrnorm(3, mu=0, Sigma=1, empirical=T)
  Y3.ch <- t(chol(fam3))%*%Y3 #Why does this look this way?
  var.Y3 <- var(Y3)
  var.Y3.ch <- var(Y3.ch)
  df3 <- data.frame(var.Y3, var.Y3.ch)
  return(df3)
}

vars3.full <- do.call(rbind, lapply(1:10000, function(x) vars3(x)))
head(vars3.full)
  hist(vars3.full$var.Y3)
  hist(vars3.full$var.Y3.ch)
  mean(vars3.full$var.Y3)
  mean(vars3.full$var.Y3.ch)

#EMPIRICAL = FALSE
  vars3 <- function(x) {
    fam3 <- matrix(c(1.0, 0.5, 0.5, 
                     0.5, 1.0, 0.5,
                     0.5, 0.5, 1.0), nrow=3)
    Y3 <- mvrnorm(3, mu=0, Sigma=1, empirical=FALSE)
    Y3.ch <- t(chol(fam3))%*%Y3 #Why does this look this way?
    var.Y3 <- var(Y3)
    var.Y3.ch <- var(Y3.ch)
    df3 <- data.frame(var.Y3, var.Y3.ch)
    return(df3)
  }
  
  vars3.full <- do.call(rbind, lapply(1:10000, function(x) vars3(x)))
  head(vars3.full)
  hist(vars3.full$var.Y3, sub="EMP FALSE")
  hist(vars3.full$var.Y3.ch, sub="EMP FALSE")
  mean(vars3.full$var.Y3)
  mean(vars3.full$var.Y3.ch)



# t(chol(sib))%*%(chol(sib))
# 
# id <- matrix(c(1,0,0,1), nrow=2)
# B <- matrix(rnorm(2), ncol=1) #pheno
# 
# 
# t(chol(sib)) %*% B
# B %*% t(chol(sib)) #doesn't work

kinship=kk2
kinship_chol <- t(chol(kinship))
N=214
P=1
#genBgShared
    B <- matrix(mvrnorm(N * P, mu=0, Sigma=1, empirical=T), ncol = P) 
    A <- matrix(rep(1, P * P), ncol = P) # only changes the scale, but nothing else....
    #A[, 1] <- rnorm(P)
    genBgShared <- kinship_chol %*% (B %*% t(A))
    genBgShared2 <- data.frame(animal=rownames(genBgShared), "Trait_1"=genBgShared)
     hist(genBgShared, sub="genBgShared, empirical=TRUE")

# #Bg Shared noise
#     Bn <- matrix(rnorm(N * P, mean = 0, sd = 0.0001), ncol = P)
#     An <- matrix(rep(1, P * P), ncol = P) #How does changing A rep value affect final estimates? Only factor/scale of values. Distr is the same
#     #An[, 1] <- rnorm(P, mean = mean, sd = sd)
#     noiseBgShared <- Bn %*% t(An)
#     hist(noiseBgShared)
# 
# # Rescale variance: Takes average column variance and rescales to the user-specified proportion of variance
#     genBg_shared_scaled <- rescaleVariance(genBg$shared, propvar=shared * genVar)
#     genBg_independent_scaled <- rescaleVariance(genBg$independent, propvar=independent * genVar) #Null when genBg all shared
#     noiseBg_shared_scaled <- rescaleVariance(noiseBg$shared, propvar= shared * noiseVar)
#     noiseBg_independent_scaled <- rescaleVariance(noiseBg$independent, propvar= independent *noiseVar) #Null
# 
# Y <- scale(genBg_shared_scaled$component + noiseBg_shared_scaled$component)
# # + genBg_independent_scaled$component + noiseBg_independent_scaled$component)
# Y2 <- data.frame(animal=rownames(Y), Y)



#H2
m1 <- MCMCglmm(Trait_1~1, random=~animal, pedigree=pedigree, data=genBgShared2, burnin=10000, nitt=1000000, thin=100,
               prior=list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002))))
plot(m1$VCV)
(m1.val<- data.frame('a2'=mean(m1$VCV[,'animal']/(m1$VCV[,'animal']+m1$VCV[, 'units'])), 
                     'e2'=mean(m1$VCV[,'units']/(m1$VCV[,'animal']+m1$VCV[, 'units']))))




#rescaleVariance

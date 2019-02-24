#Phenotype Simulator

#requires installing "snpStats" from bioconductor. To install:
#source("https://bioconductor.org/biocLite.R")
#biocLite("snpStats")


library(pbmcapply)
library(PhenotypeSimulator)
library(pedigreemm)

ph<- read.csv(file="Documents/Behavioral Genetics Study/R Projects/SNPs/final geno pheno_2018.07.18.csv",header=TRUE, na.strings = c("", "NA")) #pheno file with 'animal' IDs
p1<-read.csv(file="Documents/Behavioral Genetics Study/R Projects/Heritability/ped_dd.csv") #pedigree
p1[p1==""] <- NA 
p1$Dam[is.na(p1$Sire)] <- NA 
pedit <- editPed(p1$Sire, p1$Dam, p1$ego)[,1:3]
kk2 <- kinship(pedit[,1], pedit[,2], pedit[,3])*2
l <- levels(ph$animal)
kk2 <- kk2[l,l]



genVar <- 0.40
noiseVar <- 1- genVar
shared <- 1 
independent <- 1 - shared

genBg <- geneticBgEffects(N=214, kinship=kk2, P=1)
noiseBg <- noiseBgEffects(N = 214, P = 1)

# Rescale variance: Takes average column variance and rescales to the user-specified proportion of variance
genBg_shared_scaled <- rescaleVariance(genBg$shared, propvar=shared * genVar)
genBg_independent_scaled <- rescaleVariance(genBg$independent, propvar=independent * genVar) #Null when genBg all shared
noiseBg_shared_scaled <- rescaleVariance(noiseBg$shared, propvar= shared * noiseVar)
noiseBg_independent_scaled <- rescaleVariance(noiseBg$independent, propvar= independent *noiseVar) #Null


# Total variance proportion shave to add up to 1
total <- independent * noiseVar + independent * genVar + shared * noiseVar + shared * genVar
total == 1

# Confirm that this reflects the h2 I set
var(genBg_shared_scaled$component)/(var(genBg_shared_scaled$component)+var(noiseBg_shared_scaled$component)) 
# Y1 <- genBg_shared_scaled$component+noiseBg_shared_scaled$component
# Y1a <- data.frame(animal=rownames(Y1), Y1)
# 
#           #Do kin show expected correlations?
#           inds= which(kk2==0.5, arr.ind=TRUE)
#           cor(Y1[inds[,1]],Y1[inds[,2]]) #47
# 
#           inds2= which(kk2==0.25, arr.ind=TRUE)
#           cor(Y1[inds2[,1]],Y1[inds2[,2]]) #39
# 
#           inds3= which(kk2==0.125, arr.ind=TRUE)
#           cor(Y1[inds3[,1]],Y1[inds3[,2]]) #19

# combine components into final phenotype
Y <- scale(genBg_shared_scaled$component + noiseBg_shared_scaled$component)
# + genBg_independent_scaled$component + noiseBg_independent_scaled$component)
Y2 <- data.frame(animal=rownames(Y), Y)

#example of full-sibs
Y["RDf8",]
Y["RWe9",]


#Do kin show expected correlations?
inds= which(kk2==0.5, arr.ind=TRUE)
cor(Y[inds[,1]],Y[inds[,2]]) #

inds2= which(kk2==0.25, arr.ind=TRUE)
cor(Y[inds2[,1]],Y[inds2[,2]]) #

inds3= which(kk2==0.125, arr.ind=TRUE)
cor(Y[inds3[,1]],Y[inds3[,2]]) #


# Model Pheno with MCMCglmm
m1 <- MCMCglmm(Trait_1~1, random=~animal, pedigree=pedigree, data=Y2, burnin=10000, nitt=100000, thin=100,
               prior=list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002))))
plot(m1$VCV)
(m1.val<- data.frame('a2'=mean(m1$VCV[,'animal']/(m1$VCV[,'animal']+m1$VCV[, 'units'])), 
                     'e2'=mean(m1$VCV[,'units']/(m1$VCV[,'animal']+m1$VCV[, 'units']))))
h2.m1 <- (m1$VCV[,'animal'])/(m1$VCV[,'animal']+m1$VCV[,"units"])
plot(h2.m1) 
posterior.mode(h2.m1)


#400 sims

ps.run <- function(){
  
  genBg <- geneticBgEffects(N=214, kinship=kk2, P=1) #this is what incluldes the sampling, i think?
  noiseBg <- noiseBgEffects(N = 214, P = 1)
  genBg_shared_scaled <- rescaleVariance(genBg$shared, propvar=shared * genVar)
  noiseBg_shared_scaled <- rescaleVariance(noiseBg$shared, propvar= shared * noiseVar)
  
  Y <- scale(genBg_shared_scaled$component + noiseBg_shared_scaled$component) #or maybe this?
  Y2 <- data.frame(animal=rownames(Y), Y)
  
  m1 <- MCMCglmm(Trait_1~1, random=~animal, pedigree=pedigree, data=Y2, burnin=10000, nitt=100000, thin=100,
                 prior=list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002))))
  
  m1.val<- data.frame('a2'=mean(m1$VCV[,'animal']/(m1$VCV[,'animal']+m1$VCV[, 'units'])), 
                      'e2'=mean(m1$VCV[,'units']/(m1$VCV[,'animal']+m1$VCV[, 'units'])))
  return(m1.val$a2)
}


do.call(rbind, pbmclapply(1:2, function(x) {
  ps.run()
  
}, ignore.interactive = TRUE, mc.cores=2))

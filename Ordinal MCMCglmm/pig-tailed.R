# Change ... to your desired working directory
setwd("~/Documents/Behavioral Genetics Study/R Projects/simulations/Ordinal MCMCglmm/")


# Load required packages
library(MCMCglmm)

# Load the pedigree file and data
mcmc.ped<-read.delim("ordered-ped.txt",stringsAsFactors=F) #542
mcmc.dat<-read.delim("datafile.txt",stringsAsFactors=F)

# Set priors on the variance components
prior1<-list(G=list(G1=list(V=1,nu=0.2)),R=list(V=1,nu=0.2))

# Run the threshold model (this will take a while)
output<-MCMCglmm(trt~age+adult.mass, random=~animal,pedigree=mcmc.ped,data=mcmc.dat,prior=prior1,
                 family="ordinal",nitt=1000000,thin=1000,burnin=100000,verbose=F,pl=TRUE)

# Heritability estimate
h2<-mean(output$VCV[,"animal"]/(output$VCV[,"animal"]+output$VCV[,"units"]))

# 95% HPD Interval
CI.95<-HPDinterval(output$VCV[,"animal"]/(output$VCV[,"animal"]+output$VCV[,"units"]))

# Trace plots
plot(output$VCV)

# Histogram of h2 posterior distribution
hist(output$VCV[,"animal"]/(output$VCV[,"animal"]+output$VCV[,"units"]))


# As described in the manuscript, this histogram will indicate that there
# are not enough data here to estimate heritability, since it will
# reflect the prior distribution. Also, the trace plots of the variance
# components will likely show them to ``blow up,'' but this in and of itself
# is not necessarily a problem, as also discussed in the text.



###############################################
# Simulation on WaNPRC Pedigree
#
# Ordinal variable simulated according to the threshold model,
# showing h^2 estimability with a phenotyped sample size of 250 as described in the text.
#
# While there still may be some degree of stickiness at h^2=1 for this sample
# size, the MCMC generations do appear to contain reasonable information regarding
# h^2. Specifically, the posterior means of the MCMC chain do tend towards the true simulated
# value of h^2. 

library(MASS)		# for mvrnorm function
library(kinship2)	# for kinship function


# Create covariance matrix skeleton based on kinship coefficients calculated from pedigree
phi<-kinship(mcmc.ped$EGO,mcmc.ped$SIRE,mcmc.ped$DAM)
cov.mat<-6*2*phi+diag(4,dim(mcmc.ped)[1],dim(mcmc.ped)[1])

# Create dataframe with all 542 monkeys
mcmc.dat.sim<-cbind(read.csv("full-ped.txt",header=F,stringsAsFactors=F)[,1:2],rep(NA,542),rep(NA,542))
colnames(mcmc.dat.sim)<-c("fam","animal","trt","x1")


## make age same as real data, or rnorm(18,4)
for(i in 1:dim(mcmc.dat.sim)[1]){
   if(is.element(mcmc.dat.sim$animal[i],mcmc.dat$animal)){
      mcmc.dat.sim$x1[i]<-mcmc.dat[which(mcmc.dat$animal==mcmc.dat.sim$animal[i]),]$age
   }
   else{
      mcmc.dat.sim$x1[i]<-rnorm(1,18,4)
   }
}


# keep the most related to 173 up to 250 (greedy search)
orig.173<-mcmc.dat$animal
keep<-orig.173

# first, kick out the ones of the 173 who aren't related to anyone (why are they even in here??)
kins.orig<-rep(NA,length(orig.173))
for(i in 1:length(orig.173)){
   kins.orig[i]<-sum(phi[orig.173[i],orig.173])
}
keep<-orig.173[-which(kins.orig==0.5)]

# Now add in the most related ones (in a greedy manner) until we reach 250 monkeys
kins<-rep(NA,dim(phi)[1])
j<-length(keep)
while(length(keep)<250){
   for(i in 1:dim(phi)[1]){
      kins[i]<-sum(phi[i,keep])
   }
   kins[which(is.element(row.names(phi),keep))]<-0	# so that it doesn't consider the ones already in
   k<-which(kins==max(kins))				# find the max... could be more than one
   keep<-c(keep,row.names(phi)[k])			# keep it/them
}

# Now generate trait according to threshold model
mu<-(-10)+1.5*mcmc.dat.sim$x1 # mean vector according to parameter estimates from data
temp<-as.numeric(t(mvrnorm(1,mu,cov.mat))) 
mcmc.dat.sim[,3]<-ifelse(temp<7,1,ifelse(temp<8,2,ifelse(temp<9,3,ifelse(temp<10,4,ifelse(temp<12,5,ifelse(temp<14,6,
   ifelse(temp<17,7,ifelse(temp<20,8,ifelse(temp<28,9,10)))))))))

for(i in 1:dim(mcmc.dat.sim)[1]){
   if(!is.element(mcmc.dat.sim$animal[i],keep)){
      mcmc.dat.sim$trt[i]<-NA
   }
}


# Run MCMCglmm
output.250<-MCMCglmm(trt~x1, random=~animal,pedigree=mcmc.ped,data=mcmc.dat.sim,prior=prior1,family="ordinal",
                nitt=210000000,thin=100000,burnin=10000000,verbose=F)

# Plot histogram of estimated posterior distribution of h^2
hist(output.250$VCV[,"animal"]/(output.250$VCV[,"animal"]+output.250$VCV[,"units"]),
   xlab="heritability",main="estimated posterior distribution")



#Van der sluis Loop

rm(list=ls(all=TRUE))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# data simulation script used in
#
# "Phenotypic complexity, measurement bias, and poor phenotypic resolution 
#  contribute to the missing heritability problem in genetic association studies"
#  van der Sluis, Verhage, Posthuma, & Dolan, PLOS ONE, 2010
#
# Please refer to this paper if you use the scripts
#
#
#  SIMULATE FACTOR MODEL IN TWO POPULATIONS
#  e.g. men vs women or lab1 versus lab2
#  the latent factor is regressed on one QTL, 
#  which explains some variance in the latent factor
#  the parameters in the factor model can differ across the two sub pops
#  effect size QTL depends on the exact parameters estimates in the two samples
#  but is defined for group1 !! 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

library(norm) # norm has to be present!
library(MASS) # library MASS should be present!

# file name: first number refers to lambda,
#            second to allele freq p
#            third to perc variance explained
#            so set3_1_05=lambda=.3,p=.1,%varexpl=1/2
namefile=as.character("set39_5_1")            # nr of items

rwrite=T       # write results to external files
check=T        # check the results


N=1200   # anticipated total sample N (may change slightly due to rounding below)
pc=.5    # percentage of total sample group1 
N1=round(pc*N)   # group1
N2=round((1-pc)*N)   # group2
ngroup=2
ng=3    # nr subgroups [genotype groups]
nvar=6
ngene=1
nfac=1

# -------------------------- DEFINE SAMPLE ------------------------------------------

# assign effect size [see also below]
perc=.01 # percentage explained variance in group1 
# note: 1=100%, .1=10% and .01=1%, etc
# NOTE: effect size is defined for p=q=.5

# DEFINE GENE FREQS
# allele freqs entire sample (group1 + group2)
p=.50   # allele freq sample
q=1-p
allfreq=c(p^2,2*p*q,q^2)

# allele freqs group1
p1=p
q1=q
allfreq1=allfreq

# allele freqs group2
p2=p
q2=q
allfreq2=allfreq

ns1=round(N1*allfreq1)
ns2=round(N2*allfreq2)
Ntot1=sum(ns1)
Ntot2=sum(ns2)
Ntot=rbind(ns1,ns2)
Ntot

group=c(rep(1,Ntot1),rep(2,Ntot2))

gentot=c(rep(-1,ns1[1]),
         rep(0,ns1[2]),
         rep(1,ns1[3]),
         rep(-1,ns2[1]),
         rep(0,ns2[2]),
         rep(1,ns2[3]))

# variance of the gene
x1=gentot[1:Ntot1]
x2=gentot[(Ntot1+1):(Ntot1+Ntot2)]
varg1=var(x1)*((Ntot1-1)/Ntot1)  #adjust variance because Mx divides by N, and not N-1
varg2=var(x2)*((Ntot2-1)/Ntot2)

# -------------------------- DEFINE GENETIC EFFECTS ------------------------------------------
# GENETIC EFFECT ON LATENT FACTOR
# for convenience, variance of the factor excl. gene effect =1
# how large should the gene effect x then be to explain y % of the variance?
# to calculate percentage explained by x
# y = x / var(tot) = x / (var(eta)+x)
# so: x = (y*var(eta)) / (1-y)
# where y is coded as .01 if 1% and 1 if 100%
# e.g., if QTL explains 3% (i.e., y=.03) and var(eta)=3:
# x = .03*3 / 1-.03 = .09/.97=.093
#    indeed: .093 / 3+.093 = .03

geffect=(perc*1)/(1-perc)

# expected phenotypic mean in the population is  
# mpop= a(p-q) +2pqd, which reduces to a(p-q) when dominance is absent
# expected variance in the population:
# spop= 2pq(a-(q-p)d)^2 +2pqd^2, which reduces to 2pqa^2 when dominance is absent
# so beta=2pqa^2 => a^2=(beta/2pq) => a=sqrt(beta/2pq)
# this is used below 

# note that the genetic effect is defined for p=q=.5
# which is why we always divide by .5, and not by 2pq

a1=sqrt(geffect/(.5))
a2=sqrt(geffect/(.5))

b1=c(-a1,0,a1)
b2=c(-a2,0,a2)


# -------------------------- DEFINE PHENOTYPIC DATA ------------------------------------------

# define factor loadings
L1=matrix(c(.3,.3,.3,.9,.9,.9),nvar,nfac,byrow=F)
L2=matrix(c(.3,.3,.3,.9,.9,.9),nvar,nfac,byrow=F)

# define factor variances
Y1=matrix(c(1),nfac,nfac)            # factor variance=1 (unconditional)
Y2=matrix(c(1),nfac,nfac)

# define residuals
# calc residuals such that L^2 + E^2 = 1
E1=diag(1,nvar)
for (i in 1:nvar) {
  for (j in 1:nfac) {
    E1[i,i]=.6
  } }

E2=diag(1,nvar)
for (i in 1:nvar) {
  for (j in 1:nfac) {
    E2[i,i]=.6
  } }

# calculate sigma
sigma1=L1%*%Y1%*%t(L1)+(E1)             
sigma2=L2%*%Y2%*%t(L2)+(E2)

# observed means / intercepts (observed scores)

omean1=matrix(0,ng,nvar)
omean2=matrix(0,ng,nvar)

# latent means 
lmean1=matrix(0,ng,nfac)
lmean2=matrix(0,ng,nfac)

# -------------------------- CREATE DATA ------------------------------------------

# group1 (Ntot1)
ii=0
dat1=matrix(0,N1,nvar)

for (k in 1:ng) {
  corn=(Ntot[1,k]-1)/Ntot[1,k] # adjust for unbiased estimator: exact fit in Mx
                                #First row (means group1), Kth column (means geno A, B, or C)
  for (i in 1:nvar) {                     # Mx divides by N, rather then N-1
    for (j in 1:nfac) {
      omean1[k,i]=omean1[k,i]+L1[i,j]*(b1[k]+lmean1[k])
    }}
  tmp1=mvrnorm(n=Ntot[1,k],mu=omean1[k,],Sigma=(sigma1/corn),empirical=T) 
  for (ip in 1:Ntot[1,k]) {
    ii=ii+1
    dat1[ii,]=tmp1[ip,]
  }
}

# group1 (Ntot2)
ii=0
dat2=matrix(0,N2,nvar)

for (k in 1:ng) {
  corn=(Ntot[2,k]-1)/Ntot[2,k]            # adjust for unbiased estimator: exact fit in Mx
  for (i in 1:nvar) {
    for (j in 1:nfac) {
      omean2[k,i]=omean2[k,i]+L2[i,j]*(b2[k]+lmean2[k])
    }}
  tmp2=mvrnorm(Ntot[2,k],omean2[k,],(sigma2/corn),empirical=T) 
  for (ip in 1:Ntot[2,k]) {
    ii=ii+1
    dat2[ii,]=tmp2[ip,]
  }
}



# -------------------- WRITE -------------------------------------

tmpx=rbind(dat1,dat2)
var1=tmpx[,1]
var2=tmpx[,2]
var3=tmpx[,3]
var4=tmpx[,4]
var5=tmpx[,5]
var6=tmpx[,6]

som=c()
for (i in 1:sum(Ntot)) {
  som[i]=var1[i]+var2[i]+var3[i]+var4[i]+var5[i]+var6[i]
}

dat=cbind(group,gentot,(rbind(dat1,dat2)),som,gentot)

if (rwrite==T) {
  write(t(dat),file=namefile,nc=nvar+4) }

# -------------------------- CHECK -----------------------------


group=dat[,1]
gen=dat[,2]
var1=dat[,3]
var2=dat[,4]
var3=dat[,5]
var4=dat[,6]
var5=dat[,7]
var6=dat[,8]


if (check==T) {
  # CHECK N1
  check1=matrix(0,ng, nvar)
  check1[1,1]=mean(var1[group==1&gen==-1],digits=6)
  check1[1,2]=mean(var2[group==1&gen==-1],digits=6)
  check1[1,3]=mean(var3[group==1&gen==-1],digits=6)
  check1[1,4]=mean(var4[group==1&gen==-1],digits=6)
  check1[1,5]=mean(var5[group==1&gen==-1],digits=6)
  check1[1,6]=mean(var6[group==1&gen==-1],digits=6)
  check1[2,1]=mean(var1[group==1&gen==0],digits=6)
  check1[2,2]=mean(var2[group==1&gen==0],digits=6)
  check1[2,3]=mean(var3[group==1&gen==0],digits=6)
  check1[2,4]=mean(var4[group==1&gen==0],digits=6)
  check1[2,5]=mean(var5[group==1&gen==0],digits=6)
  check1[2,6]=mean(var6[group==1&gen==0],digits=6)
  check1[3,1]=mean(var1[group==1&gen==1],digits=6)
  check1[3,2]=mean(var2[group==1&gen==1],digits=6)
  check1[3,3]=mean(var3[group==1&gen==1],digits=6)
  check1[3,4]=mean(var4[group==1&gen==1],digits=6)
  check1[3,5]=mean(var5[group==1&gen==1],digits=6)
  check1[3,6]=mean(var6[group==1&gen==1],digits=6)
  
  min1=check1-omean1 ;  round(min1)
}

if (check==T) {
  # CHECK N2
  check2=matrix(0,ng, nvar)
  check2[1,1]=mean(var1[group==2&gen==-1],digits=6)
  check2[1,2]=mean(var2[group==2&gen==-1],digits=6)
  check2[1,3]=mean(var3[group==2&gen==-1],digits=6)
  check2[1,4]=mean(var4[group==2&gen==-1],digits=6)
  check2[1,5]=mean(var5[group==2&gen==-1],digits=6)
  check2[1,6]=mean(var6[group==2&gen==-1],digits=6)
  check2[2,1]=mean(var1[group==2&gen==0],digits=6)
  check2[2,2]=mean(var2[group==2&gen==0],digits=6)
  check2[2,3]=mean(var3[group==2&gen==0],digits=6)
  check2[2,4]=mean(var4[group==2&gen==0],digits=6)
  check2[2,5]=mean(var5[group==2&gen==0],digits=6)
  check2[2,6]=mean(var6[group==2&gen==0],digits=6)
  check2[3,1]=mean(var1[group==2&gen==1],digits=6)
  check2[3,2]=mean(var2[group==2&gen==1],digits=6)
  check2[3,3]=mean(var3[group==2&gen==1],digits=6)
  check2[3,4]=mean(var4[group==2&gen==1],digits=6)
  check2[3,5]=mean(var5[group==2&gen==1],digits=6)
  check2[3,6]=mean(var6[group==2&gen==1],digits=6)
  
  min2=check2-omean2 ; round(min2)
}


# ---------------------------- ANALYIS --------------------------------
# see analysis_continuous.r
# and group1.mx

b1
b2
p1
p2

Ntot




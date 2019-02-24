#CpG: How does R2 change if I used individual beh vs Noisy composite?
#March 21, 2018

library(pscl)
library(lmtest)
library(BaylorEdPsych)
setwd("~/Documents/Behavioral Genetics Study/R Projects/simulations/")

#1) Read file with Methylation and Behavior
Data3<- read.csv(file ="~/Documents/Behavioral Genetics Study/Shota Methylation Study/OXTR epi file_2018.02.08.csv")
str(Data3)

#1a)
Data2 <- Data3[!is.na(Data3$CpG1),]  #Removes Missing Values
Data<- Data2[!is.na(Data2$sc),]  #Removes Missing Values
str(Data)


#1) Prox
Data$resid_px<-residuals(lm(px_adultfem~compound, data=Data, na.action=na.exclude))
x1<-summary(lm(px_adultfem~CpG1+compound, data=Data))
x2<-summary(lm(resid_px~CpG1, data=Data))$r.squared

Data$resid_px1<-residuals(lm(px_all_notownix~compound, data=Data, na.action=na.exclude))
x3<-summary(lm(px_all_notownix~CpG1+compound, data=Data))
x4<-summary(lm(resid_px1~CpG1, data=Data))$r.squared

x2
x4


#3) cpg1 and Anxiety
SS<-summary(glm.nb((sc)~CpG1+compound+offset(log(Time)), data=Data))
#Pseudo R2
pR2(glm.nb((sc)~CpG1+compound+offset(log(Time)), data=Data))
pR2(glm.nb((sum_anx)~CpG1+compound+offset(log(Time)), data=Data))




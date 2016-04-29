#### PEER_CM.RData

##(a) Full Model

# Get the data

options(stringsAsFactors = FALSE)
rm(list=ls())
setwd("C:/Users/jillh/Dropbox/DHGLab/commonmind/")

load("FinalProcData_CM.RData")

datMeta$Dx <- factor(datMeta$Dx, levels=c("Control","BP","SCZ"))
datMeta$Ethnicity <- factor(datMeta$Ethnicity, levels=c("Caucasian","African-American","Hispanic","Asian"))

## Get the covariate data
mDx <- model.matrix(~datMeta$Dx)
mDx <- mDx[,-1]
age <- as.numeric(as.factor(datMeta[,"Age_of_Death"]))
sex <- as.numeric(as.factor(datMeta[,"Gender"]))-1
RIN <- as.numeric(as.character(datMeta[,"RIN"]))
race <- model.matrix(~datMeta$Ethnicity)
race <- race[,-1]
PMI <- as.numeric(as.factor(datMeta$PMI_hrs))
seqStatPC1 <- as.numeric(datMeta$`SeqPC1 (98.7%) - Depth`)
seqStatPC2 <- as.numeric(datMeta$`SeqPC2 (<2%) - GC/Length`)
regvars <- as.data.frame(cbind(mDx,age,sex,RIN,race,seqStatPC1,seqStatPC2))
colnames(regvars)[1:2]=c("Bipolar","Schizo")

mod = regvars

## Start PEER
library(peer)

model = PEER()

PEER_setPhenoMean(model,as.matrix(t(datExpr)))
dim(PEER_getPhenoMean(model))
PEER_setCovariates(model,as.matrix(mod))
dim(PEER_getCovariates(model))

PEER_setNk(model,100)
PEER_getNk(model)

PEER_update(model)

factors = PEER_getX(model)
dim(factors)
weights = PEER_getW(model)
dim(weights)
precision = PEER_getAlpha(model)
dim(precision)
residuals = PEER_getResiduals(model)
dim(residuals)
#plot(precision)

save(factors,mod,file="PEER_factors_CM100.RData")




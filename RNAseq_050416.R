## Surrogate Variable Analysis and Collinearity Presentation 
## RNAseq_050416_1-4.R
## Use the Common Mind dataset to first illustrate what to do with collinear data.
## Then, apply PEER and HCP to the Common Mind data and compare the obtained beta values to 
## the microarray beta values.
## Finally, use the Brain GVEX data (a non-confounded dataset) to illustrate what PEER and HCP
## do again by compring to the same microarray beta values
## In a nutshell:
## (1) Common Mind Collinearity Adjustment, obtain beta
## (2) Common Mind betas v. Microarray betas
## (3) Common Mind PEER betas 
## (4) Common Mind HCP betas 
## (5) Common Mind SVA betas

options(stringsAsFactors = FALSE)
setwd("C:/Users/jillh/Dropbox/GitHub/RNAseq/")
rm(list=ls())

#### (0) Check that Brain GVEX is not confounded, in case I want to use it later

datMeta = read.csv("C:/Users/jillh/Dropbox/DHGLab/CDRNASeqExpr/BrainGVEX/datMeta.csv")[1:153,]
rownames(datMeta)=datMeta[,1]
datMeta = datMeta[,-4]
datMeta$Diagnosis = factor(datMeta$Diagnosis,levels=c("Control","BP","SCZ"))

colnames(datMeta)

pdf(file="CovariatesANOVA_BrainGVEX.pdf")

par(mfrow=c(2,2),mar=c(3,2.5,2,1),cex=0.75)
plot(datMeta$Diagnosis, ylab="Number", main="Subjects")

for(i in colnames(datMeta)[-c(1,2,11)]){
  
  if( i == "BrainBank" || i == "Sex"  || i == "YearAutopsy" || i == "Ethnicity" || i == "CauseDeath" 
      || i == "Hemisphere" || i == "LibraryBatch"){
    print(paste(i,"Character Graph",sep=" "))
    A = anova(lm(as.numeric(as.factor(datMeta[,i])) ~ datMeta$Diagnosis)); p = A$"Pr(>F)"[1]   
    plot(as.factor(datMeta[,i]) ~ datMeta$Diagnosis, main=paste(i," p=", signif(p,2)), ylab="", xlab="")
  }
  else{
    print(paste(i,"Number Graph",sep=" "))
    A = anova(lm(as.numeric(datMeta[,i]) ~ datMeta$Diagnosis)); p = A$"Pr(>F)"[1]   
    plot(as.numeric(as.character(datMeta[,i])) ~ datMeta$Diagnosis, main=paste(i," p=", signif(p,2)), ylab="", xlab="")
  }
}

dev.off()

## Need to remove some males from Brain GVEX

idx = which(datMeta$Sex=="M" & datMeta$Diagnosis=="Control") 
datMetaCM= datMeta[idx,]

## Remove the top 10 with the highest brain weight

datMetaCM = datMetaCM[order(datMetaCM$BrainWeight),]
to_remove = rownames(datMetaCM)[35:44]

idx = which(rownames(datMeta)%in%to_remove)

datMeta_sub = datMeta[-idx,]

## Also remove 5 bipolar females with lowest RIN

idx = which(datMeta$Sex=="F" & datMeta$Diagnosis=="BP")
datMetaBF = datMeta[idx,]
datMetaBF = datMetaBF[order(datMetaBF$RIN),]
to_remove = rownames(datMetaBF)[1:5]

idx = which(rownames(datMeta_sub)%in%to_remove)
datMeta_sub = datMeta_sub[-idx,]

## Redo the covariate analysis

pdf(file="CovariatesANOVA_BrainGVEX_Subset.pdf")

par(mfrow=c(2,2),mar=c(3,2.5,2,1),cex=0.75)
plot(datMeta_sub$Diagnosis, ylab="Number", main="Subjects")

for(i in colnames(datMeta)[-c(1,2,11)]){
  
  if( i == "BrainBank" || i == "Sex"  || i == "YearAutopsy" || i == "Ethnicity" || i == "CauseDeath" 
      || i == "Hemisphere" || i == "LibraryBatch"){
    print(paste(i,"Character Graph",sep=" "))
    A = anova(lm(as.numeric(as.factor(datMeta_sub[,i])) ~ datMeta_sub$Diagnosis)); p = A$"Pr(>F)"[1]   
    plot(as.factor(datMeta_sub[,i]) ~ datMeta_sub$Diagnosis, main=paste(i," p=", signif(p,2)), ylab="", xlab="")
  }
  else{
    print(paste(i,"Number Graph",sep=" "))
    A = anova(lm(as.numeric(datMeta_sub[,i]) ~ datMeta_sub$Diagnosis)); p = A$"Pr(>F)"[1]   
    plot(as.numeric(as.character(datMeta_sub[,i])) ~ datMeta_sub$Diagnosis, main=paste(i," p=", signif(p,2)), ylab="", xlab="")
  }
}

dev.off()

write.csv(datMeta_sub,file="C:/Users/jillh/Dropbox/DHGLab/CDRNASeqExpr/BrainGVEX/datMeta_subset_nocoll_BrainGVEX.csv")

## You want all of your relevant datMeta categories to have p > 0.01. 
## Especially age, ethnicity, sex, RIN, and batch.
## Some things (like cause of death) are going to be too variable to ever balance
## Those types of covariates can be excluded from your model

## Brain GVEX model will include: 
## Diagnosis,PMI, pH, Sex, Ethnicity, AgeDeath, RIN, SeqStats
## LibraryBatch/BrainBank will be removed by ComBat prior to linear model

#### (1) Common Mind Collinearity Adjustment, obtain beta

## See commonMind.R for pre-processing

rm(list=ls())
load("C:/Users/jillh/Dropbox/DHGLab/commonmind/FinalProcData_CM.RData")

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
seqStatPC1 <- as.numeric(datMeta$`SeqPC1 (98.7%) - Depth`)
seqStatPC2 <- as.numeric(datMeta$`SeqPC2 (<2%) - GC/Length`)
regvars <- as.data.frame(cbind(mDx,age,sex,RIN,race,seqStatPC1,seqStatPC2))
colnames(regvars)[1:2]=c("Bipolar","Schizo")

## Make object for ttables
ttable = list("Bipolar"=matrix(NA,nrow=nrow(datExpr),ncol=5),
              "Schizo"=matrix(NA,nrow=nrow(datExpr),ncol=5))
rownames(ttable$Bipolar) <- rownames(ttable$Schizo) <- rownames(datExpr)
colnames(ttable$Bipolar) <- colnames(ttable$Schizo) <- c("Estimate","Std.Error","t.value","Pr(>|t|)","p.adj")

## Run the regression
datExpr.reg <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(datExpr))
rownames(datExpr.reg) <- rownames(datExpr)
colnames(datExpr.reg) <- colnames(datExpr)

lmmod <- apply(as.matrix(datExpr),1,function(y) lm(y~age+sex+RIN+race+seqStatPC1+seqStatPC2,data=regvars))

for (i in 1:nrow(datExpr)) {
  if (i%%1000 == 0) {print(i)}
  
  coef <- coef(lmmod[[i]])
  datExpr.reg[i,] <- coef[1] + residuals(lmmod[[i]]) 
  
  lmmod2 <- lm(datExpr.reg[i,]~Bipolar+Schizo,data=regvars)
  lmmod2$df.residual <- lmmod2$df.residual - 8   ## 8 the number of coefficients previously removed
  summary = summary(lmmod2)
  
  ttable$Bipolar[i,1:4] <- coef(summary)[2,]
  ttable$Schizo[i,1:4] <- coef(summary)[3,]
  
  ## The full data minus all of the covariates but condition 
  ## Also equivalent to <- datExpr - coef*var (unwanted variates) model
  ## So we have a datExpr object with only the contributions of Dx
  
  ## Then adjust the degrees of freedom for the residual before we do our t tests
  ## Ignore the warning
}


## Get adjusted p values
for(i in c(1:2)){
  ttable[[i]][,5] <- p.adjust(ttable[[i]][,4],method="fdr")
}

write.csv(ttable$Bipolar,file="DiffExpr_Bipolar_CM.csv")
write.csv(ttable$Schizo,file="DiffExpr_Schizo_CM.csv")

quantile(datExpr[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))
quantile(datExpr.reg[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))

save(datMeta,datSeq,datExpr.reg,file="datExpr.reg_CM.RData")

all_betas_lm = cbind(ttable$Bipolar[,1],ttable$Schizo[,1]) 
colnames(all_betas_lm) = c("BP","SCZ")
rownames(all_betas_lm) = rownames(datExpr.reg)

library(corrplot)

par(mfrow=c(1,1)); par(mar=c(5,4,2,2)); par(oma=c(1,1,1,1))
pdf(file="CM_BetaCorr.pdf")
corrplot.mixed(cor(all_betas_lm, method="spearman"), lower="ellipse", upper= "number")
title(main="Correlation CM - RLog Covariates Removed",line=1)
dev.off()
## 0.36

#### (2) Compare pure betas to microarray betas
## Now compare to the 'gold standard' microarray data

home = getwd()
setwd("C:/Users/jillh/Dropbox/DHGLab/CDRNASeqExpr/Final/NonLieberData/")

mSCZ = read.csv(file="scz_meta-analysis_Chen_Maycox_Lanz_Iwa_Narayan_062615.csv")
mBP = read.csv(file="bad_meta-analysis_Chen_Lanz_Iwa_062615.csv")

rSCZ = all_betas_lm[,"SCZ"]
rBP = all_betas_lm[,"BP"]

setwd(home)

## SCZ

intGenes = intersect(mSCZ[,1],rownames(all_betas_lm))
idx = match(intGenes,mSCZ[,1])
mSCZ = mSCZ[idx,]
idx = match(intGenes,rownames(all_betas_lm))
rSCZ = rSCZ[idx]

## BP

intGenes = intersect(mBP[,1],rownames(all_betas_lm))
idx = match(intGenes,mBP[,1])
mBP = mBP[idx,]
idx = match(intGenes,rownames(all_betas_lm))
rBP = rBP[idx]

##pdf(file="RNASeqVMicroarray_ComparisonGraphs_CM_regOut.pdf")
pdf(file="RNASeqVMicroarray_ComparisonGraphs_CM_HCPb_regOut.pdf")

par(mfrow=c(1,1)); par(mar=c(5,5,3,5)); par(oma=c(1,1,1,1))
plot(mSCZ[,2] ~ rSCZ,main="RNASeq v. Microarray SCZ",xlab="RNASeq",ylab="Microarray")
abline(fit <- lm(mSCZ[,2] ~ rSCZ),col="red")
legend("bottomright",bty="n",legend=paste("R =",format(sqrt(summary(fit)$adj.r.squared),digits=4)))

plot(mBP[,2] ~ rBP,main="RNASeq v. Microarray BD",xlab="RNASeq",ylab="Microarray")
abline(fit <- lm(mBP[,2] ~ rBP),col="red")
legend("bottomright",bty="n",legend=paste("R =",format(sqrt(summary(fit)$adj.r.squared),digits=4)))

dev.off()

save(all_betas_lm,file="CM_Betas_pure.RData")

#### (3) Common Mind PEER betas v. Microarray betas

# See PEER_CM.R

load("PEERreg_CM.RData")
datExpr = datExpr.res
## Regress out the PEER factors and then...
## (a) Try full model, accounting for 8 removed covariates
## (b) Try regressing out everything but Dx, accounting for 16 removed covariates

datMeta$Dx <- factor(datMeta$Dx, levels=c("Control","BP","SCZ"))
datMeta$Ethnicity <- factor(datMeta$Ethnicity, levels=c("Caucasian","African-American","Hispanic","Asian"))
mDx <- model.matrix(~datMeta$Dx)
mDx <- mDx[,-1]
race <- model.matrix(~datMeta$Ethnicity)
race <- race[,-1]

## Now linear regression (a) full

## Make object for ttables
ttable = list("Bipolar"=matrix(NA,nrow=nrow(datExpr),ncol=5),
              "Schizo"=matrix(NA,nrow=nrow(datExpr),ncol=5))
rownames(ttable$Bipolar) <- rownames(ttable$Schizo) <- rownames(datExpr)
colnames(ttable$Bipolar) <- colnames(ttable$Schizo) <- c("Estimate","Std.Error","t.value","Pr(>|t|)","p.adj")

lmmod <- apply(as.matrix(datExpr),1,function(y) lm(y~mDx+age+sex+RIN+race+seqStatPC1+seqStatPC2,data=mod))

for (i in 1:nrow(datExpr)) {
  if (i%%1000 == 0) {print(i)}
  
  lmmod[[i]]$df.residual <- lmmod[[i]]$df.residual - 8   ## 8 the number of PEER coefficients previously removed
  summary = summary(lmmod[[i]])
  
  ttable$Bipolar[i,1:4] <- coef(summary)[2,]
  ttable$Schizo[i,1:4] <- coef(summary)[3,]
}


## Get adjusted p values
for(i in c(1:2)){
  ttable[[i]][,5] <- p.adjust(ttable[[i]][,4],method="fdr")
}

all_betas_lm = cbind(ttable$Bipolar[,1],ttable$Schizo[,1]) 
colnames(all_betas_lm) = c("BP","SCZ")
rownames(all_betas_lm) = rownames(datExpr)

library(corrplot)

par(mfrow=c(1,1)); par(mar=c(5,4,2,2)); par(oma=c(1,1,1,1))
corrplot.mixed(cor(all_betas_lm, method="spearman"), lower="ellipse", upper= "number")
title(main="Correlation CM - PEERa factors removed",line=1)
##0.28
save(all_betas_lm,file="CM_Betas_PEERa.RData")

## Now linear regression (b) regress out

## Make object for ttables
ttable = list("Bipolar"=matrix(NA,nrow=nrow(datExpr),ncol=5),
              "Schizo"=matrix(NA,nrow=nrow(datExpr),ncol=5))
rownames(ttable$Bipolar) <- rownames(ttable$Schizo) <- rownames(datExpr)
colnames(ttable$Bipolar) <-  colnames(ttable$Schizo) <- c("Estimate","Std.Error","t.value","Pr(>|t|)","p.adj")

## Run the regression
datExpr.reg <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(datExpr))
rownames(datExpr.reg) <- rownames(datExpr)
colnames(datExpr.reg) <- colnames(datExpr)

lmmod <- apply(as.matrix(datExpr),1,function(y) lm(y~age+sex+RIN+race+seqStatPC1+seqStatPC2,data=mod))

for (i in 1:nrow(datExpr)) {
  if (i%%1000 == 0) {print(i)}
  
  coef <- coef(lmmod[[i]])
  datExpr.reg[i,] <- coef[1] + residuals(lmmod[[i]]) 
  
  lmmod2 <- lm(datExpr.reg[i,]~Bipolar+Schizo,data=mod)
  lmmod2$df.residual <- lmmod2$df.residual - 16   ## 16 the number of coefficients previously removed
  summary = summary(lmmod2)
  
  ttable$Bipolar[i,1:4] <- coef(summary)[2,]
  ttable$Schizo[i,1:4] <- coef(summary)[3,]
  
  ## The full data minus all of the covariates but condition 
  ## Also equivalent to <- datExpr - coef*var (unwanted variates) model
  ## So we have a datExpr object with only the contributions of Dx
  
  ## Then adjust the degrees of freedom for the residual before we do our t tests
  ## Ignore the warning
}


## Get adjusted p values
for(i in c(1:2)){
  ttable[[i]][,5] <- p.adjust(ttable[[i]][,4],method="fdr")
}

all_betas_lm = cbind(ttable$Bipolar[,1],ttable$Schizo[,1]) 
colnames(all_betas_lm) = c("BP","SCZ")
rownames(all_betas_lm) = rownames(datExpr)

library(corrplot)

par(mfrow=c(1,1)); par(mar=c(5,4,2,2)); par(oma=c(1,1,1,1))
corrplot.mixed(cor(all_betas_lm, method="spearman"), lower="ellipse", upper= "number")
title(main="Correlation CM - PEERb factors removed",line=1)
##0.33
save(all_betas_lm,file="CM_Betas_PEERb.RData")

# use LM_ConfoundRemove.R, modify, for (a) and as is for (b), then compare to each other and microarray

#### (4) Common Mind HCP betas 

## Regress out the HCP factors and then...
## (a) Try full model, accounting for 20 removed covariates
## (b) Try regressing out everything but Dx, accounting for 28 removed covariates

## get data

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
seqStatPC1 <- as.numeric(datMeta$`SeqPC1 (98.7%) - Depth`)
seqStatPC2 <- as.numeric(datMeta$`SeqPC2 (<2%) - GC/Length`)
regvars <- as.data.frame(cbind(mDx,age,sex,RIN,race,seqStatPC1,seqStatPC2))
colnames(regvars)[1:2]=c("Bipolar","Schizo")

mod = regvars

## Get HCP functions

mod.std = standardize(mod)

datExpr.std = standardize(t(datExpr))

hcp = hidden_covariate_linear(mod.std,datExpr.std,k=20)

datExpr.reg = datExpr.std - hcp$Z%*%hcp$B

datExpr = t(datExpr.reg) ## hidden factors now regressed out

save(datExpr,file="datExpr_HCP_CM.RData")


## Now linear regression (a) full
load("datExpr_HCP_CM.RData")

## Make object for ttables
ttable = list("Bipolar"=matrix(NA,nrow=nrow(datExpr),ncol=5),
              "Schizo"=matrix(NA,nrow=nrow(datExpr),ncol=5))
rownames(ttable$Bipolar) <- rownames(ttable$Schizo) <- rownames(datExpr)
colnames(ttable$Bipolar) <- colnames(ttable$Schizo) <- c("Estimate","Std.Error","t.value","Pr(>|t|)","p.adj")

lmmod <- apply(as.matrix(datExpr),1,function(y) lm(y~mDx+age+sex+RIN+race+seqStatPC1+seqStatPC2,data=regvars))

for (i in 1:nrow(datExpr)) {
  if (i%%1000 == 0) {print(i)}
  
  lmmod[[i]]$df.residual <- lmmod[[i]]$df.residual - 20   ## 20 the number of HCP coefficients previously removed
  summary = summary(lmmod[[i]])
  
  ttable$Bipolar[i,1:4] <- coef(summary)[2,]
  ttable$Schizo[i,1:4] <- coef(summary)[3,]
}


## Get adjusted p values
for(i in c(1:2)){
  ttable[[i]][,5] <- p.adjust(ttable[[i]][,4],method="fdr")
}

all_betas_lm = cbind(ttable$Bipolar[,1],ttable$Schizo[,1]) 
colnames(all_betas_lm) = c("BP","SCZ")
rownames(all_betas_lm) = rownames(datExpr)

library(corrplot)

par(mfrow=c(1,1)); par(mar=c(5,4,2,2)); par(oma=c(1,1,1,1))
pdf(file="CM_BetaCorrHCPa.pdf")
corrplot.mixed(cor(all_betas_lm, method="spearman"), lower="ellipse", upper= "number")
title(main="Correlation CM - HCP factors removed",line=1)
dev.off()
##0.23
save(all_betas_lm,file="CM_Betas_HCPa.RData")

## Now linear regression (b) regress out

## Make object for ttables
ttable = list("Bipolar"=matrix(NA,nrow=nrow(datExpr),ncol=5),
              "Schizo"=matrix(NA,nrow=nrow(datExpr),ncol=5))
rownames(ttable$Bipolar) <- rownames(ttable$Schizo) <- rownames(datExpr)
colnames(ttable$Bipolar) <-  colnames(ttable$Schizo) <- c("Estimate","Std.Error","t.value","Pr(>|t|)","p.adj")

## Run the regression
datExpr.reg <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(datExpr))
rownames(datExpr.reg) <- rownames(datExpr)
colnames(datExpr.reg) <- colnames(datExpr)

lmmod <- apply(as.matrix(datExpr),1,function(y) lm(y~age+sex+RIN+race+seqStatPC1+seqStatPC2,data=regvars))

for (i in 1:nrow(datExpr)) {
  if (i%%1000 == 0) {print(i)}
  
  coef <- coef(lmmod[[i]])
  datExpr.reg[i,] <- coef[1] + residuals(lmmod[[i]]) 
  
  lmmod2 <- lm(datExpr.reg[i,]~Bipolar+Schizo,data=regvars)
  lmmod2$df.residual <- lmmod2$df.residual - 28   ## 28 the number of coefficients previously removed
  summary = summary(lmmod2)
  
  ttable$Bipolar[i,1:4] <- coef(summary)[2,]
  ttable$Schizo[i,1:4] <- coef(summary)[3,]
  
  ## The full data minus all of the covariates but condition 
  ## Also equivalent to <- datExpr - coef*var (unwanted variates) model
  ## So we have a datExpr object with only the contributions of Dx
  
  ## Then adjust the degrees of freedom for the residual before we do our t tests
  ## Ignore the warning
}


## Get adjusted p values
for(i in c(1:2)){
  ttable[[i]][,5] <- p.adjust(ttable[[i]][,4],method="fdr")
}

all_betas_lm = cbind(ttable$Bipolar[,1],ttable$Schizo[,1]) 
colnames(all_betas_lm) = c("BP","SCZ")
rownames(all_betas_lm) = rownames(datExpr)

library(corrplot)

par(mfrow=c(1,1)); par(mar=c(5,4,2,2)); par(oma=c(1,1,1,1))
pdf(file="CM_BetaCorrHCPb.pdf")
corrplot.mixed(cor(all_betas_lm, method="spearman"), lower="ellipse", upper= "number")
title(main="Correlation CM - HCP factors removed b",line=1)
dev.off()
##0.28
save(all_betas_lm,file="CM_Betas_HCPb.RData")

#### (5) Common Mind SVA betas

## get data

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
seqStatPC1 <- as.numeric(datMeta$`SeqPC1 (98.7%) - Depth`)
seqStatPC2 <- as.numeric(datMeta$`SeqPC2 (<2%) - GC/Length`)
regvars <- as.data.frame(cbind(mDx,age,sex,RIN,race,seqStatPC1,seqStatPC2))
colnames(regvars)[1:2]=c("Bipolar","Schizo")

mod = cbind(rep(1,567),regvars)
colnames(mod)[1]="Intercept"
mod0 = mod[,-c(2,3)]

## sva removal

library(sva)

n.sv = num.sv(as.matrix(datExpr),as.matrix(mod))
n.sv #34

svobj = sva(as.matrix(datExpr),as.matrix(mod),as.matrix(mod0),n.sv=n.sv)

factors = svobj$sv

###(a) with full model

## Make object for ttables
ttable = list("Bipolar"=matrix(NA,nrow=nrow(datExpr),ncol=5),
              "Schizo"=matrix(NA,nrow=nrow(datExpr),ncol=5))
rownames(ttable$Bipolar) <- rownames(ttable$Schizo) <- rownames(datExpr)
colnames(ttable$Bipolar) <-  colnames(ttable$Schizo) <- c("Estimate","Std.Error","t.value","Pr(>|t|)","p.adj")

## Run the regression
datExpr.reg <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(datExpr))
rownames(datExpr.reg) <- rownames(datExpr)
colnames(datExpr.reg) <- colnames(datExpr)

lmmod <- apply(as.matrix(datExpr),1,function(y) lm(y~Bipolar+Schizo+age+sex+RIN+race+seqStatPC1+seqStatPC2+factors,data=regvars))

for (i in 1:nrow(datExpr)) {
  if (i%%1000 == 0) {print(i)}
  
  coef <- coef(lmmod[[i]])
  datExpr.reg[i,] <- coef[1] + coef[2]*regvars[,"Bipolar"] + coef[3]*regvars[,"Schizo"] + 
    coef[4]*regvars[,"age"] + coef[5]*regvars[,4] + coef[6]*regvars[,5] + coef[7]*regvars[,6] + 
    coef[8]*regvars[,7] + coef[9]*regvars[,8] +coef[10]*regvars[,9] + coef[11]*regvars[,10] +residuals(lmmod[[i]]) 
  
  lmmod2 <- lm(datExpr.reg[i,]~Bipolar+Schizo+age+sex+RIN+race+seqStatPC1+seqStatPC2o,data=regvars)
  lmmod2$df.residual <- lmmod2$df.residual - n.sv   ## n.sv the number of coefficients previously removed
  summary = summary(lmmod2)
  
  ttable$Bipolar[i,1:4] <- coef(summary)[2,]
  ttable$Schizo[i,1:4] <- coef(summary)[3,]
  
  ## The full data minus all of the covariates but condition 
  ## Also equivalent to <- datExpr - coef*var (unwanted variates) model
  ## So we have a datExpr object with only the contributions of Dx
  
  ## Then adjust the degrees of freedom for the residual before we do our t tests
  ## Ignore the warning
}


## Get adjusted p values
for(i in c(1:2)){
  ttable[[i]][,5] <- p.adjust(ttable[[i]][,4],method="fdr")
}

all_betas_lm = cbind(ttable$Bipolar[,1],ttable$Schizo[,1]) 
colnames(all_betas_lm) = c("BP","SCZ")
rownames(all_betas_lm) = rownames(datExpr)

library(corrplot)

par(mfrow=c(1,1)); par(mar=c(5,4,2,2)); par(oma=c(1,1,1,1))
pdf(file="CM_BetaCorrSVAa.pdf")
corrplot.mixed(cor(all_betas_lm, method="spearman"), lower="ellipse", upper= "number")
title(main="Correlation CM - SVA factors removed a",line=1)
dev.off()

save(all_betas_lm,file="CM_Betas_SVAa.RData")





###(b) with regressed model



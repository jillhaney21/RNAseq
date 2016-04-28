### PEER testing on example data

setwd("/geschwindlabshares/CrossDisorder_transcriptome_comparison/PEER/")

library(peer)
library(qtl)

#setwd("./examples/data/")

#cross <- read.cross(format="csvs",genfile="genotype.csv",phefile = "expression.csv",genotypes=c(0,1))

#geno = read.csv(file="genotype.csv")

expr = read.csv(file="expression.csv")
dim(expr)

model = PEER()

PEER_setPhenoMean(model,as.matrix(expr))
dim(PEER_getPhenoMean(model))

PEER_setNk(model,10)
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
plot(precision)


### Apply to our CDRNA expression data

rm(list=ls())
library(WGCNA)
options(stringsAsFactors=FALSE)

## Load Sequencing Data

setwd("/geschwindlabshares/CrossDisorder_transcriptome_comparison/PEER/")

library(peer)
library(qtl)

load("FinalProcData.rda")

datMeta$Dx = factor(datMeta$Dx, levels= c("Control", "Schizo", "Bipolar", "MDD"))
datMeta$SeqPC1= datMeta$"SeqPC1 (98%)\nDepth"
datMeta$SeqPC2= datMeta$"SeqPC2 (2%)\nGC/Length"

design = "~Dx+age+Sex+RIN+Race+SeqPC1+SeqPC2"
mod = model.matrix(as.formula(design),data = datMeta)
colnames(mod) = make.names(colnames(mod))
mod = mod[,-1]

model = PEER()

PEER_setPhenoMean(model,as.matrix(t(datExpr)))
dim(PEER_getPhenoMean(model))
PEER_setCovariates(model,as.matrix(mod))
dim(PEER_getCovariates(model))

PEER_setNk(model,100)
PEER_getNk(model)

PEER_setAdd_mean(model,TRUE)  ## an additional covariate will be added to account for mean expression

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

modMat = cbind(rep(1,109),factors)
facs = c(1:31)
facNames = paste("F",facs,sep="")
names = c("Intercept",colnames(mod),facNames)
colnames(modMat) = names

library(limma)

#design = "~Dx+age+Sex+RIN+Race+SeqPC1+SeqPC2"
#mod = model.matrix(as.formula(design),data = datMeta)
#colnames(mod) = make.names(colnames(mod))
lm = lmFit(datExpr,modMat)
lm = eBayes(lm)

datExpr.reg = datExpr
for (i in 1:nrow(datExpr)){
  datExpr.reg[i,] = datExpr[i,] - lm$coefficients[i,5]*mod[,5] - lm$coefficients[i,6]*mod[,6] - lm$coefficients[i,7]*mod[,7] - lm$coefficients[i,10]*mod[,10]
}

quantile(datExpr[1,])
quantile(datExpr.reg[1,])
## datExpr now only has the contributions of Dx, race, and seqStatPC1 (bc of confoundment by these two last covariates)
save(datExpr.reg,file="datExpr_regressed.rda")

scz = topTable(lm,coef = 2,number = Inf,sort.by = "none")
bd = topTable(lm,coef = 3,number = Inf,sort.by = "none")
mdd = topTable(lm,coef = 4,number = Inf,sort.by = "none")

all_betas_lm = cbind(scz[,1],bd[,1],mdd[,1])
colnames(all_betas_lm) = c("SCZ","BP","MDD")
rownames(all_betas_lm) = rownames(datExpr)

library(corrplot)

plot.new(); par(mfrow=c(1,1)); par(mar=c(5,4,2,2)); par(oma=c(1,1,1,1))
corrplot.mixed(cor(all_betas_lm, method="spearman"), lower="ellipse", upper= "number")
title(main="Corr LM - 30 PEER factors + sum + noSSPCs ",line=-2)

### Try..
### 10 factors
### 10 factors + sum factor
### 10 factors + sum factor + without seq stats
### 30 factors
### 30 factors + sum factor
### 30 factors + sum factor + without seq stats

save(model_10,model_10sum,model_10sum_noSPCs,file="PEERmods_10F.RData")

load("PEERmods_30F.RData")
save(model_30,model_30sum,model_30sum_noSPCs,file="PEERmods_30F.RData")


## Basic Linear Model for RNAseq Data
## How to do your analysis when condition is confounded with any other covariates

options(stringsAsFactors = FALSE)
setwd("C:/Users/jillh/Dropbox/DHGLab/CDRNASeqExpr/Final/")
load("FinalProcData_cqnNoQN.rda")

datMeta$Dx <- factor(datMeta$Dx, levels=c("Control","Bipolar","MDD","Schizo"))

## Get the covariate data
mDx <- model.matrix(~datMeta$Dx)
mDx <- mDx[,-1]
age <- as.numeric(as.factor(datMeta[,"age"]))
sex <- as.numeric(as.factor(datMeta[,"Sex"]))-1
RIN <- as.numeric(as.character(datMeta[,"RIN"]))
race <- as.numeric(datMeta[,"Race"])-1
seqStatPC1 <- as.numeric(datMeta[,"SeqPC1 (98%)\nDepth"])
seqStatPC2 <- as.numeric(datMeta[,"SeqPC2 (2%)\nGC/Length"])
regvars <- as.data.frame(cbind(mDx,age,sex,RIN,race,seqStatPC1,seqStatPC2))
colnames(regvars)[1:3]=c("Bipolar","MDD","Schizo")

## Run the regression
datExpr.reg <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(datExpr))
rownames(datExpr.reg) <- rownames(datExpr)
colnames(datExpr.reg) <- colnames(datExpr)

lmmod <- apply(as.matrix(datExpr),1,function(y) lm(y~age+sex+RIN+race+seqStatPC1+seqStatPC2,data=regvars))

for (i in 1:nrow(datExpr)) {
  if (i%%1000 == 0) {print(i)}
  
  coef <- coef(lmmod[[i]])
  datExpr.reg[i,] <- coef[1] + residuals(lmmod[[i]]) 
  
  ## The full data minus all of the covariates but condition 
  ## Also equivalent to <- datExpr - coef*var (unwanted variates) model
  ## So we have a datExpr object with only the contributions of Dx
}

ttest <- function(reg, coefnum,resid){
  co <- coef(summary(reg))
  tstat <- (co[coefnum,1]-0)/co[coefnum,2]
  pval <- (2 * pt(abs(tstat), (reg$df.residual-resid), lower.tail = FALSE))
  se <- co[coefnum,2]
  coef <- co[coefnum,1]
  stats <- list("coef"=coef,"error"=se,"t"=tstat,"p"=pval)
  return(stats)
}

## We need to have a function that uses the correct degrees of freedom for the residuals

## Bipolar only as our example

BPstats = matrix(NA,nrow=nrow(datExpr.reg),ncol=5)
colnames(BPstats) = c("Estimate","Std.Error","t value","Pr(>|t|)","p.adj")
rownames(BPstats) = rownames(datExpr.reg)

for (i in 1:dim(BPstats)[1]){
  if (i%%1000 == 0) {print(i)}
  stats=ttest(lm(datExpr.reg[i,]~Bipolar+MDD+Schizo,data=regvars),2,6) 
  ## 2 = Biplar, 6 = number of coef in lmmod
  BPstats[i,1] = stats$coef
  BPstats[i,2] = stats$error
  BPstats[i,3] = stats$t
  BPstats[i,4] = stats$p
}

BPstats[,5] = p.adjust(BPstats[,4],method = "fdr")

write.csv(BPstats,file="DiffExpr_Bipolar.RData")

quantile(datExpr[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))
quantile(datExpr.reg[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))

save(datMeta,datExpr.reg,coefmat.reg,file="datExpr.reg_exapmle_RData")

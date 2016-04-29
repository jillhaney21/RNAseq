## LM_ConfoundRemove_Short.R

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

## Make object for ttables
ttable = list("Bipolar"=matrix(NA,nrow=nrow(datExpr),ncol=5),
              "MDD"=matrix(NA,nrow=nrow(datExpr),ncol=5),
              "Schizo"=matrix(NA,nrow=nrow(datExpr),ncol=5))
rownames(ttable$Bipolar) <- rownames(ttable$MDD) <- rownames(ttable$Schizo) <- rownames(datExpr)
colnames(ttable$Bipolar) <- colnames(ttable$MDD) <- colnames(ttable$Schizo) <- c("Estimate","Std.Error","t.value","Pr(>|t|)","p.adj")

## Run the regression
datExpr.reg <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(datExpr))
rownames(datExpr.reg) <- rownames(datExpr)
colnames(datExpr.reg) <- colnames(datExpr)

lmmod <- apply(as.matrix(datExpr),1,function(y) lm(y~age+sex+RIN+race+seqStatPC1+seqStatPC2,data=regvars))

for (i in 1:nrow(datExpr)) {
  if (i%%1000 == 0) {print(i)}
  
  coef <- coef(lmmod[[i]])
  datExpr.reg[i,] <- coef[1] + residuals(lmmod[[i]]) 

  lmmod2 <- lm(datExpr.reg[i,]~Bipolar+MDD+Schizo,data=regvars)
  lmmod2$df.residual <- lmmod2$df.residual - 6   ## 6 the number of coefficients previously removed
  summary = summary(lmmod2)
  
  ttable$Bipolar[i,1:4] <- coef(summary)[2,]
  ttable$MDD[i,1:4] <- coef(summary)[3,]
  ttable$Schizo[i,1:4] <- coef(summary)[4,]
  
  ## The full data minus all of the covariates but condition 
  ## Also equivalent to <- datExpr - coef*var (unwanted variates) model
  ## So we have a datExpr object with only the contributions of Dx
  
  ## Then adjust the degrees of freedom for the residual before we do our t tests
  ## Ignore the warning
}


## Get adjusted p values
for(i in c(1:3)){
  ttable[[i]][,5] <- p.adjust(ttable[[i]][,4],method="fdr")
}


write.csv(ttable$Bipolar,file="DiffExpr_Bipolar.csv")
write.csv(ttable$MDD,file="DiffExpr_MDD.csv")
write.csv(ttable$Schizo,file="DiffExpr_Schizo.csv")

quantile(datExpr[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))
quantile(datExpr.reg[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))

save(datMeta,datExpr.reg,file="datExpr.reg_example_RData")
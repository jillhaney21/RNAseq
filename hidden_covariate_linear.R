## HCP for R
## Transcribed by Michael Gandal

setwd("~/GeschwindLab/R_scripts/HCP/")
Y = read.csv("Y.csv",head=F)
F = read.csv("F.csv",head=F)

standardize<- function(X)
{
  X = as.matrix(X)
  n = dim(X)[1]
  p = dim(X)[2]
  
  X = scale(X, center = TRUE, scale=FALSE)
  X = scale(X,center=FALSE, scale=sqrt(apply(X^2,2,sum)))

  # m = apply(X,2,mean)
  # st = sqrt(apply(X^2,2,sum));
  # st_mat = matrix(st, nrow = length(st), ncol = dim(X)[2], byrow=FALSE)
  # X2 = X / st_mat
  
}




 hidden_covariate_linear <- function(F,Y,k = 10 ,lambda = 5,lambda2 = 1, lambda3 =1 , iter = 100) {
   #
   #
   # function [Z,B,U,o,error,error1,error2,dz,db,du] = hidden_covariate_linear(F,Y,k,lambda,lambda2,lambda3,iter);
   # input:
   #      F: a matrix nxd of known covariates, where n is the number of
   #      subjects and d is the number of known covariates. *must be standardize (columns have 0 mean and constant SS).
   #      Y: a matrix of nxg of expression data (must be standardized (columns
   #      scaled to have constant SS and mean 0). ** use standardize function to standardize F and Y.
   #      k: number of inferred hidden components (k is an integer)
   #      lambda, lambda2, lambda3 are model parameters
   #      (optional) iter: number of iterations (default = 100);
   #
   #      note: k>0, lambda>0, lambda2>0, lambda3>0 must be set by the user based on the data at hand. one can set these values
   #      using cross-validation, by evaluating the "performance" of the  resulting residual data on a desired task.
   #      typically, if lambda>5, then hidden factors match the known covariates closely. 
   #
   # objective:
   #
   # this function solves the following problem:
   # argmin_{Z,B,U}   ||Y-Z*B||_2 + \lambda*||Z-F*U||_2 + \lambda2*||B||_2 + \lambda_3||U||_2
   #
   # output:
   #      Z: matrix of hidden components, dimensionality: nxk
   #      B: matrix of effects of hidden components, dimensionality: kxg
   #      o: value of objective function on consecutive iterations.
   #
   # to use the residual data: Residual = Y - Z*B
   library(MASS)
   library(pracma)
   
   tol = 1e-6;
   
   U = matrix(0, nrow=dim(F)[2],k)
   Z = matrix(0, nrow=dim(F)[1],k)
   B = matrix(runif(dim(Z)[2]*dim(Y)[2]), nrow=dim(Z)[2], ncol=dim(Y)[2])
   
   F = as.matrix(F)
   
   n1 = dim(F)[1];
   d1 = dim(F)[2]
 
   n2 = dim(Y)[1]
   d2 = dim(Y)[2]
   
   if(n1!=n2)    stop("number of rows in F and Y must agree")
   
   
   if (k<1 | lambda<1e-6 | lambda2<1e-6 | lambda3<1e-6 ) {
     stop("lambda, lambda2, lambda3 must be positive and/or k must be an integer");
   }
   
   o = vector(length=iter)
   
   for (ii in 1:iter) {
     o[ii] = sum((Y - Z%*%B)^2) + sum((Z -  F%*%U)^2)*lambda + (sum(B^2))*lambda2 + lambda3*(sum(U^2));
     Z = (Y %*% t(B) + lambda * F %*%U) %*% ginv(B %*% t(B) + lambda * diag(dim(B)[1]))
     B = mldivide(t(Z) %*% Z + lambda2 * diag(dim(Z)[2]), (t(Z) %*% Y))
     U = mldivide(t(F) %*% F * lambda + lambda3 * diag(dim(U)[1]), lambda * t(F) %*% Z)
     
     if(ii > 1 &&  (abs(o[ii]-o[ii-1])/o[ii]) < tol)  break
   }
         
   error =  sum((Y - Z%*%B)^2) / sum(Y^2)  + sum((Z - F%*%U)^2)/sum((F%*%U)^2)
   error1 = sum((Y - Z%*%B)^2) / sum(Y^2);
   error2 = sum((Z - F%*%U)^2) / sum((F%*%U)^2);
   
   dz = Z%*%(B%*%t(B) + lambda*diag(dim(B)[1]))-(Y%*%t(B) + lambda*F%*%U);
   db = (t(Z)%*%Z + lambda2*diag(dim(Z)[2]))%*%B - t(Z)%*%Y;
   du = (t(F)%*%F*lambda + lambda3*diag(dim(U)[1]))%*%U-lambda*t(F)%*%Z;
                    
           
   dataout = list(Z = Z, B = B, U = U)
   return(dataout)  
   
 }


################################################################## Try with the Lieber dataset (pure log2 normalize then standardize)

setwd("~/DHGLab/CDRNASeqExpr/Final/")

load("cuff_ht_filtered_unmatched.rda")
datSeq = read.csv("PicardToolsQC.csv")

## Summarize median expression

medhtg <- apply(htg,1,median)
medhte <- apply(hte,1,median)
medcuff <- apply(cuff,1,median)

### We will ultimately use HTSC union gene as our expression data (datExpr), since this gives 'whole gene' expression
### We use cufflinks and HTSC exon as filters, Cufflinks as a 'lower bound' for expression and a 'replicate' for capturing
## gene expression, and HTSC exon to ensure that the genes we find have significant exonic expression

### Filter for gene length, remove outliers, and take intersection of all datasets for datExpr to feed into hcp

library(biomaRt)
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position",
             "end_position","strand","band","gene_biotype","percentage_gc_content")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl",
                host="sep2013.archive.ensembl.org") ## Gencode v19, otherwise use ##mart <- useMart(biomart="ensembl",dataset="hsapiens <- gene <- ensembl")
#a = listAttributes(mart)

datExpr = htg

geneAnno1 <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values=rownames(datExpr),mart=mart)

geneAnno <- geneAnno1[match(rownames(datExpr),geneAnno1[,1]),]
geneAnno <- geneAnno[!is.na(geneAnno[,1]),]
datExpr <- datExpr[!is.na(geneAnno[,1]),]
datDepth <- apply(datExpr,2,sum)

geneAnno$length = geneAnno$end_position - geneAnno$start_position 
gc18unionAnno <- geneAnno[geneAnno$length>200,]
rownames(gc18unionAnno) = gc18unionAnno[,1]

keepvec <- intersect(substr(rownames(gc18unionAnno),1,15),rownames(htg))
geneAnno <- gc18unionAnno[match(keepvec,substr(rownames(gc18unionAnno),1,15)),]
rownames(geneAnno) <- substr(rownames(geneAnno),1,15)

htg_raw = htg
htg <- htg[match(rownames(geneAnno),rownames(htg)),]

# get intersection of all datasets

htscint <- intersect(rownames(hte),intersect(rownames(htg),rownames(cuff)))

### remove outliers

datExpr = log2(htg + 1)

## For all samples
library(WGCNA)
normadj <- (0.5+0.5*bicor(datExpr))^2 ## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
par(mfrow=c(2,1))
C <- netsummary$ClusterCoef
z.C <- (C-mean(C))/sqrt(var(C))
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))

datLabel <- paste(datMeta[,"Dx"],datMeta[,"age"],datMeta[,"Sex"])
plot(1:length(z.ku),z.ku)
text(1:length(z.ku),z.ku,label=datLabel,pos=2,cex=0.6)
abline(h= -2, col="red")
flag <- z.ku < -2

plot(z.ku,z.C)
text(z.ku,z.C,label=datLabel,pos=2,cex=0.6)

#datMeta <- datMeta[!flag,]
datSeq <- datSeq[!flag,]
datExpr.htg <- htg[,!flag]
datExpr.cuff <- cuff[,!flag]

idx = match(colnames(datExpr.htg),rownames(datMeta))
datMeta = datMeta[idx,]

# Assemble final raw dataset to feed into hcp

datExpr = log2(datExpr.htg[htscint,]+1)

sex = as.numeric(as.factor(datMeta[,6]))
race = as.numeric(as.factor(datMeta[,7]))
datMeta$Dx = factor(datMeta$Dx, levels= c("Control", "Schizo", "Bipolar", "MDD"))
Dx = datMeta$Dx

mod = cbind(datMeta[,4:5],sex,race,Dx,datSeq[,c(2,4,14,16,18)])
design = "~Dx+age+sex+RIN+race+TOTAL_READS+PF_HQ_ALIGNED_READS+MEDIAN_5PRIME_TO_3PRIME_BIAS+AT_DROPOUT+READ_PAIR_OPTICAL_DUPLICATES"
model = model.matrix(as.formula(design),data = mod)
colnames(model) = make.names(colnames(model))
model = model[,-1]

################################################################################### Run hcp

setwd("~/DHGLab/CDRNASeqExpr/Final/")
load("FinalProcData_vst.rda")

datMeta$Dx = factor(datMeta$Dx, levels= c("Control", "Schizo", "Bipolar", "MDD"))
datMeta$SeqPC1= datMeta$"SeqPC1 (98%)\nDepth"
datMeta$SeqPC2= datMeta$"SeqPC2 (2%)\nGC/Length"

design = "~Dx+age+Sex+RIN+Race+SeqPC1+SeqPC2"
mod = model.matrix(as.formula(design),data = datMeta)
colnames(mod) = make.names(colnames(mod))
mod = mod[,-1]

mod.std = standardize(mod)

datExpr.std = standardize(t(datExpr))

hcp = hidden_covariate_linear(mod.std,datExpr.std,k=20)

datExpr.res = datExpr.std - hcp$Z%*%hcp$B

datExpr = t(datExpr.res)   ## hidden factors now regressed out

# Now try linear regression on hcp output

#datMeta$Dx = factor(datMeta$Dx, levels= c("Control", "Schizo", "Bipolar", "MDD"))

design = "~Dx+age+Sex+RIN+Race+SeqPC1+SeqPC2"
model = model.matrix(as.formula(design),data = datMeta)
colnames(model) = make.names(colnames(model))
library(limma)
lm = lmFit(datExpr,model)
lm = eBayes(lm)

#datExpr.reg = datExpr
#for (i in 1:nrow(datExpr)){
#  datExpr.reg[i,] = datExpr[i,] - lm$coefficients[i,5]*mod[,5] - lm$coefficients[i,6]*mod[,6] - lm$coefficients[i,7]*mod[,7] - lm$coefficients[i,10]*mod[,10]
#}

#quantile(datExpr[1,])
#quantile(datExpr.reg[1,])
## datExpr now only has the contributions of Dx, race, and seqStatPC1 (bc of confoundment by these two last covariates)
#save(datExpr.reg,file="datExpr_regressed_cqnQuantNorm.rda")

scz = topTable(lm,coef = 2,number = Inf,sort.by = "none")
bd = topTable(lm,coef = 3,number = Inf,sort.by = "none")
mdd = topTable(lm,coef = 4,number = Inf,sort.by = "none")

all_betas_lm = cbind(scz[,1],bd[,1],mdd[,1])
colnames(all_betas_lm) = c("SCZ","BP","MDD")
rownames(all_betas_lm) = rownames(datExpr)

library(corrplot)

pdf(file="HCP_k20_CorrPlot.pdf")
par(mfrow=c(1,1)); par(mar=c(5,4,2,2)); par(oma=c(1,1,1,1))
corrplot.mixed(cor(all_betas_lm, method="spearman"), lower="ellipse", upper= "number")
title(main="Correlation LM - log2 HCP All CoV's k=20",line=-2)
dev.off()

## Add in ETOH and ASD to RNA Seq LM Correlation

setwd("~/DHGLab/CDRNASeqExpr/Final/NonLieberData/")

ASD_stats = read.csv(file="ASD_RNAseq_DEX_N106_with_permutation.csv")
ASD_stats = ASD_stats[-20999,]
ETOH_stats = read.csv(file="Mayfield_ETOH_limma_fullmodel.csv")

asd_beta = as.data.frame(ASD_stats[,"beta.condition"])
rownames(asd_beta) = ASD_stats[,2]

etoh_beta = as.data.frame(ETOH_stats[,"logFC"])
rownames(etoh_beta) = ETOH_stats[,1]

intGenes = intersect(rownames(asd_beta),intersect(rownames(etoh_beta),rownames(all_betas_lm)))
idx = match(intGenes,rownames(asd_beta))
asd_beta = as.data.frame(asd_beta[idx,],row.names = rownames(asd_beta)[idx])
idx = match(intGenes,rownames(etoh_beta))
etoh_beta = as.data.frame(etoh_beta[idx,],row.names = rownames(etoh_beta)[idx])
idx = match(intGenes,rownames(all_betas_lm))
all_betas_lm = as.data.frame(all_betas_lm[idx,],row.names = rownames(all_betas_lm)[idx])

all_betas_lm = cbind(all_betas_lm,asd_beta,etoh_beta)
colnames(all_betas_lm)[c(4,5)] = c("ASD","ETOH")

beta_graph = cbind(all_betas_lm[,"ASD"],all_betas_lm[,"SCZ"],all_betas_lm[,"BP"],all_betas_lm[,"MDD"],all_betas_lm[,"ETOH"])
colnames(beta_graph) = c("ASD","SCZ","BD","MDD","ETOH")

setwd("../")
pdf(file="HCP_k20_CorrPlot_AllDx.pdf")
par(mfrow=c(1,1)); par(mar=c(5,5,3,5)); par(oma=c(1,1,1,1))
corrplot.mixed(cor(beta_graph, method="spearman"), lower="ellipse", upper= "number",
               main="RNASeq Validation Data Corr log2 HCP k=20",sub="ETOH = Microarray",line=-1.25,cex.main=1.5,cex.sub=0.75)
dev.off()


#### Now look for correlations between microarray and RNASeq datasets

setwd("~/DHGLab/CDRNASeqExpr/Final/NonLieberData/")

mSCZ = read.csv(file="scz_meta-analysis_Chen_Maycox_Lanz_Iwa_Narayan_062615.csv")
mBP = read.csv(file="bad_meta-analysis_Chen_Lanz_Iwa_062615.csv")
mMDD = read.csv(file="mdd_meta-analysis_Sibille_Lanz_0626215.csv")
mASD = read.csv(file="asd_meta-analysis_Voin_Chow_Garbett_062615.csv")
mETOH = read.csv(file="Mayfield_ETOH_limma_fullmodel.csv")

rSCZ = all_betas_lm[,"SCZ"]
rBP = all_betas_lm[,"BP"]
rMDD = all_betas_lm[,"MDD"]
rASD = all_betas_lm[,"ASD"]

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

## MDD

intGenes = intersect(mMDD[,1],rownames(all_betas_lm))
idx = match(intGenes,mMDD[,1])
mMDD = mMDD[idx,]
idx = match(intGenes,rownames(all_betas_lm))
rMDD = rMDD[idx]

## ASD

intGenes = intersect(mASD[,1],rownames(all_betas_lm))
idx = match(intGenes,mASD[,1])
mASD = mASD[idx,]
idx = match(intGenes,rownames(all_betas_lm))
rASD = rASD[idx]

## Now graph

setwd("../")

pdf(file="RNASeqVMicroarray_ComparisonGraphs_hcpFROMcqnNoQN.pdf")

par(mfrow=c(1,1)); par(mar=c(5,5,3,5)); par(oma=c(1,1,1,1))
plot(mSCZ[,2] ~ rSCZ,main="RNASeq v. Microarray SCZ",xlab="RNASeq",ylab="Microarray")
abline(fit <- lm(mSCZ[,2] ~ rSCZ),col="red")
legend("bottomright",bty="n",legend=paste("R =",format(sqrt(summary(fit)$adj.r.squared),digits=4)))

plot(mBP[,2] ~ rBP,main="RNASeq v. Microarray BD",xlab="RNASeq",ylab="Microarray")
abline(fit <- lm(mBP[,2] ~ rBP),col="red")
legend("bottomright",bty="n",legend=paste("R =",format(sqrt(summary(fit)$adj.r.squared),digits=4)))

plot(mMDD[,2] ~ rMDD,main="RNASeq v. Microarray MDD",xlab="RNASeq",ylab="Microarray")
abline(fit <- lm(mMDD[,2] ~ rMDD),col="red")
legend("bottomright",bty="n",legend=paste("R =",format(sqrt(summary(fit)$adj.r.squared),digits=4)))

plot(mASD[,2] ~ rASD,main="RNASeq v. Microarray ASD",xlab="RNASeq",ylab="Microarray")
abline(fit <- lm(mASD[,2] ~ rASD),col="red")
legend("bottomright",bty="n",legend=paste("R =",format(sqrt(summary(fit)$adj.r.squared),digits=4)))

dev.off()




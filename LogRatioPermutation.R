library(TCGAbiolinks)
library(dplyr)
library(survival)
library(matrixStats)
library(fdrci)
##### Download survival and data and make it for further analysis ######
clinical_paad <- TCGAquery_clinic("paad","clinical_patient")
clinical_paad$new_death <- c()
for (i in 1:length(as.numeric(as.character(clinical_paad$days_to_death)))){
  clinical_paad$new_death[i] <- ifelse(is.na(as.numeric(as.character(clinical_paad$days_to_death))[i]),
                                       as.numeric(as.character(clinical_paad$days_to_last_followup))[i],as.numeric(as.character(clinical_paad$days_to_death))[i])
}
clinical_paad$death_event <- ifelse(clinical_paad$vital_status == "Alive", 0,1)
rownames(clinical_paad) <- clinical_paad$bcr_patient_barcode
########### This part of code is still not tested #########
# cox1 <- coxph(Surv(time,cens)~tgrade,data=GBSG2)
# summary(cox1)
# coxph(formula = Surv(time, cens) ~ tgrade, data = GBSG2)
####Both log-rank p-value are same ################
# t <- c(1:15)
# m <- median(t)
# sd <- sd(t)
# ifelse(t >= m+sd, "2", ifelse(t < m-sd, "1","0"))
# ifelse(t >= m+sd, "High", ifelse(t < m-sd, "Low","Modr"))
###########################################################
## X ==Survival age, A== Censor, Y== MPS data for each gene, B=MPS-class level
#attach(clinical_paad)
myanalysis = function(X,Y){
  attach(X)
  ntests = ncol(Y)
  rslts = as.data.frame(matrix(NA,nrow=ntests,ncol=2))
  names(rslts) = c("ID","pvalue")
  rslts[,"ID"] = 1:ntests
  for(i in 1:ntests){
    MPS.med <- colMedians(Y) ### Median for particular gene
    MPS.median <- MPS.med[i]
    MPS.class <- ifelse(Y[,i] >= MPS.median, "High","Low")
    #fit = survdiff(Surv(as.numeric(clinical_paad$new_death/30.4167),clinical_paad$death_event)~tmp3)
    fit = survdiff(Surv(as.numeric(new_death/30.4167), death_event)~ MPS.class)
    p.value <- round(1 - pchisq(fit$chisq, length(fit$n) - 1),3)
    rslts[i,"pvalue"] = p.value
  }
  return(rslts)
} # End myanalysis
## Generate observed results
MPS1 <- t(MPS)
rownames(MPS1) <- gsub("-01","", rownames(MPS1))
rownames(clinical_paad) <- clinical_paad$bcr_patient_barcode
clinical_paad.1 <- clinical_paad[rownames(MPS1),]
obs = myanalysis(clinical_paad.1,MPS1)
## Generate permuted results
nperm = 10
X <- MPS1[,1:25]
ncol_ <- ncol(MPS) ### Number column in X1
perml = vector('list',nperm)
for(p_ in 1:nperm){
  #X1 = MPS1[order(runif(ncol_)),
  X1 = MPS1[order(runif(ncol_)),]
  ##I think in my case both "sample" and "order" command work well
  perml[[p_]] = myanalysis(clinical_paad.1,X1)
}
## FDR results table
fdrTbl(obs$pvalue,perml,"pvalue",ncol_,1, 1.5)
#log10(0.01)=-2,log10(0.05)=-1.30103-- We can change lower and upper cutoff

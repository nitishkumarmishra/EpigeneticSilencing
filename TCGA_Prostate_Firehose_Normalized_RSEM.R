## R code for Firehose Gene Normalized RSEM plot for MALATA1 gene
## We use GDC harmozized clinical (from TCGABiolinks) as a latest clinical file
setwd("F:/OneDrive - University of Nebraska Medical Center/Manish Tripathi UTHSC/Prostate Cancer/")
load("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/GDC-Harmonized/PRAD/PRAD.RData")
#rm(list= ls()[!(ls() %in% c('Prostate.Clinical'))])
GDC_RNASeq.GTF <- data.table::fread("../TCGA_GDC_Harmonided_RNASeq-GFT.txt", header = TRUE)
#GDC_RNASeq.GTF[which(GDC_RNASeq.GTF$Gene_Symbol=="SLC2A1"),] ## It's give idea ENSG00000117394.18 is SLC2A1/GLUT1

exp.gene <- PRAD.HTSeq.Counts
exp.miR <- PRAD.miRNA.Counts
colnames(exp.gene) <- substr(colnames(exp.gene), 1,16)
colnames(exp.miR) <- substr(colnames(exp.miR), 1,16)
GLUT1 <- log2(exp.gene[which(rownames(exp.gene)=="ENSG00000117394.18"),])
mir132 <- log2(exp.miR[which(rownames(exp.miR)=="hsa-mir-132"),])
Expression <- merge(t(GLUT1), t(mir132), by="row.names")
rownames(Expression) <- Expression$Row.names; Expression$Row.names <- NULL
Expression.t <- Expression[substr(rownames(Expression), 14,15)=="01",]
Expression.t$submitter_id <- substr(rownames(Expression.t), 1,12)
# ## Use RSEM gene normalized data, as they have data for MALAT1 gene

#INT(N days * 0.032854884083862)  = N months
library(dplyr)
PRAD.Clinical.Select <- PRAD.Clinical %>% select(submitter_id, vital_status, days_to_death,days_to_last_follow_up)
rownames(PRAD.Clinical.Select) <- PRAD.Clinical.Select$submitter_id
PRAD.Clinical.Select$Days <- as.integer(ifelse(PRAD.Clinical.Select$vital_status == "alive", PRAD.Clinical.Select$days_to_last_follow_up, PRAD.Clinical.Select$days_to_death))
PRAD.Clinical.Select$Time <- PRAD.Clinical.Select$Days*0.032854884083862
PRAD.Clinical.Select$Status <- ifelse(PRAD.Clinical.Select$vital_status=="alive", 0, 1)

data.PRAD <- merge(Expression.t, PRAD.Clinical.Select, by="submitter_id")
data.PRAD$GLUT1 <- ifelse( data.PRAD$ENSG00000117394.18 >= mean(data.PRAD$ENSG00000117394.18), "High", "Low")
data.PRAD$miR132 <- ifelse( data.PRAD$`hsa-mir-132` >= mean(data.PRAD$`hsa-mir-132`), "High", "Low")

library(survival)
library(survminer)
#MALAT1.clinic$status <- ifelse(MALAT1.clinic$VITAL_STATUS=="Alive", 0, 1)
#Survival_data <- MALAT1.clinic[-which(MALAT1.clinic$OS_MONTHS=="[Not Available]"),]
#Survival_data$OS_MONTHS <- as.numeric(as.character(Survival_data$OS_MONTHS))
#Survival_data$Status <- ifelse(Survival_data$MALAT1 >= (mean(Survival_data$MALAT1)+sd(Survival_data$MALAT1)), "High", ifelse(Survival_data$MALAT1 <= (mean(Survival_data$MALAT1)-sd(Survival_data$MALAT1)), "Low", "Medium"))
data.PRAD$GLUT1 <- ifelse(data.PRAD$ENSG00000117394.18 >= (mean(data.PRAD$ENSG00000117394.18)+sd(data.PRAD$ENSG00000117394.18)), "High", ifelse(data.PRAD$ENSG00000117394.18 <= (mean(data.PRAD$ENSG00000117394.18)-sd(data.PRAD$ENSG00000117394.18)), "Low", "Medium"))
data.PRAD$miR132 <- ifelse(data.PRAD$`hsa-mir-132` >= (mean(data.PRAD$`hsa-mir-132`)+sd(data.PRAD$`hsa-mir-132`)), "High", ifelse(data.PRAD$`hsa-mir-132` <= (mean(data.PRAD$`hsa-mir-132`)-sd(data.PRAD$`hsa-mir-132`)), "Low", "Medium"))

sfit <- survfit(Surv(Time, Status)~GLUT1, data=data.PRAD)
ggsurvplot(sfit, pval = TRUE, risk.table = "abs_pct", risk.table.height=0.15, risk.table.fontsize=4, surv.median.line = "hv", pval.size = 5, pval.coord = c(120,1),  main="Kaplan-Meier Overall Survival Curves", ylab = "Probability of survival", palette=c("blue", "red", "green"))
#dev.print(pdf, 'MALAT1 Survival Plot.pdf', width = 10, height = 10)

sfit <- survfit(Surv(Time, Status)~miR132, data=data.PRAD)
ggsurvplot(sfit, pval = TRUE, risk.table = "abs_pct", risk.table.height=0.15, risk.table.fontsize=4, surv.median.line = "hv", pval.size = 5, pval.coord = c(120,1),  main="Kaplan-Meier Overall Survival Curves", ylab = "Probability of survival", palette=c("blue", "red", "green"))


data.PRAD.1 <- data.PRAD[which(data.PRAD$GLUT1 != "Medium"),]
sfit <- survfit(Surv(Time, Status)~GLUT1, data=data.PRAD.1)
ggsurvplot(sfit, pval = TRUE, risk.table = "abs_pct", risk.table.height=0.15, risk.table.fontsize=4, surv.median.line = "hv", pval.size = 5, pval.coord = c(100,1),  main="Kaplan-Meier Overall Survival Curves", ylab = "Probability of survival", palette=c("blue", "red", "green"))
#dev.print(pdf, 'MALAT1 Survival Plot.pdf', width = 10, height = 10)


data.PRAD.1 <- data.PRAD[which(data.PRAD$miR132 != "Medium"),]
sfit <- survfit(Surv(Time, Status)~miR132, data=data.PRAD.1)
ggsurvplot(sfit, pval = TRUE, risk.table = "abs_pct", risk.table.height=0.15, risk.table.fontsize=4, surv.median.line = "hv", pval.size = 5, pval.coord = c(100,1),  main="Kaplan-Meier Overall Survival Curves", ylab = "Probability of survival", palette=c("blue", "red", "green"))
#dev.print(pdf, 'MALAT1 Survival Plot.pdf', width = 10, height = 10)

save.image("GLUT1_miR132.RData")

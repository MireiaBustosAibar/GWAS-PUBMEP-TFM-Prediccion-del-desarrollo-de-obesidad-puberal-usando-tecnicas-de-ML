######################################################################3
#
#                 LOAD MAIN DIRECTORY USING BASH AND DATA IN R
#
######################################################################3

#a=/media/mireia/HardDiskHDD/Documentos/Universidad/Master_Bioinformatica_UOC/Curso_21_22/SEMESTRE_2/TFM/PERFORMING_ANALYSIS/2_DATABASES_ML
#cd $a
#R
setwd("/media/mireia/HardDiskHDD/Documentos/Universidad/Master_Bioinformatica_UOC/Curso_21_22/SEMESTRE_2/TFM/PERFORMING_ANALYSIS/2_DATABASES_ML")
#setwd("/media/mireia/HardDiskHDD/Documentos/Universidad/Master_Bioinformatica_UOC/Curso_21_22/SEMESTRE_2/TFM/PERFORMING_ANALYSIS/2_DATABASES_ML/2022_05_26_Datasets_Pubmep_TFM_mireia")
library("readxl")

lepData <- as.data.frame(read_excel("/media/mireia/HardDiskHDD/Documentos/Universidad/Master_Bioinformatica_UOC/Curso_21_22/SEMESTRE_2/TFM/PERFORMING_ANALYSIS/2_DATABASES_ML/2022_05_26_Datasets_Pubmep_TFM_mireia/2022_05_26_DB_pubmep_bioquimica_T1_imputada.xlsx"));
rownames(lepData) <- lepData$Code_new_T2; lepData <- lepData[,-26]

iniData <- as.data.frame(read_excel("PUBMEP_FENOTIPOS_RAW.xlsx"));
rownames(iniData) <- iniData$Code_new_T2

ms <- as.data.frame(read_excel("PUBMEP_FENOTIPOS_RAW_MS.xlsx"));
rownames(ms) <- ms[,1]; ms <- ms [,-1]

acelImp <- as.data.frame(read_excel("PUBMEP_Acelerometria_IMPUTADA.xlsx"));
rownames(acelImp) <- acelImp[,1]; acelImp <- acelImp[,-1]
acelnoImp <-as.data.frame(read_excel("PUBMEP_Acelerometria_NO_IMPUTADA.xlsx"));
rownames(acelnoImp) <- acelnoImp[,1]; acelnoImp <- acelnoImp[,-1]

bq <- lepData
#bq <- as.data.frame(read_excel("PUBMEP_BQ_IMPUTADA.xlsx"));
#rownames(bq) <- bq[,1]; bq <- bq[,-1]

prs <- as.data.frame(read_excel("ALL_PRSs.xlsx"))
rownames(prs) <- prs[,1]; prs <- prs [,-1]

table(rownames(prs)==rownames(ms))
table(rownames(prs)==rownames(acelImp))
table(rownames(prs)==rownames(acelnoImp))
table(rownames(prs)==rownames(bq))

#ms <- apply(ms,2,as.factor)
head(ms,3)
table(ms)
table(is.na(ms))
#ms_clean <- as.data.frame(ms[-which(is.na(ms)),])
ms_clean <- na.omit(ms)
ms_clean[,c(1,3)] <- apply(ms_clean[,c(1,3)],2,as.factor)
ms_clean <- as.data.frame(ms_clean)
apply(ms_clean,2,table)
#cat(capture.output(print(ms_clean), file="Groups_Frequence_table_PUBMEP.txt"))

# METABOLIC SYNDROME
ob_sm <- ms_clean[which(ms_clean$Obesity_WC == "Obese" & (ms_clean$Altered_components_Metabolic_Health >= 2)),]
#ob_sm <- ms_clean[which(ms_clean$Obesity_WC == "Obese" & ms_clean$Metabolic_Health_Status == "Unhealthy"),]

ob_sm$Metabolic_Syndrome_2COMP <- rep(1,ncol(ob_sm))
apply(ob_sm,2,table)

no_sm <- ms_clean[-which(ms_clean$Obesity_WC == "Obese" & (ms_clean$Altered_components_Metabolic_Health >= 2)),]
#no_sm <- ms_clean[-which(ms_clean$Obesity_WC == "Obese" & ms_clean$Metabolic_Health_Status == "Unhealthy"),]
no_sm$Metabolic_Syndrome_2COMP <- rep(0,ncol(no_sm))
apply(no_sm,2,table)

ms_new <- rbind(ob_sm,no_sm)
ms_new <- ms_new[ order(row.names(ms_new)), ]
table(rownames(ms_new) == rownames(ms_clean))
apply(ms_new,2,table)

Metabolic_Syndrome <- as.data.frame(ms_new[,4])
rownames(Metabolic_Syndrome) <- rownames(ms_new)
names(Metabolic_Syndrome) <- "Metabolic_Syndrome_2COMP"
#cat(capture.output(print(ms_ob), file="Only_obesity_Frequence_tables_PUBMEP.txt"))
#cat(capture.output(print(ms_nw), file="Only_normal_weight_Frequence_tables_PUBMEP.txt"))



################################################################################
#                     DEFINIG IMPUTED DATABASES
################################################################################

ml1 <- merge(Metabolic_Syndrome,prs,by=0)
rownames(ml1) <- ml1$Row.names; ml1 <- ml1[,-1]
names(ml1)

iniml3 <- merge(Metabolic_Syndrome,bq,by=0)
rownames(iniml3) <-iniml3$Row.names; iniml3 <- iniml3[,-1]
ml3 <- merge(iniml3[],acelImp,by=0)
rownames(ml3) <- ml3$Row.names; ml3 <- ml3[,-1]
names(ml3)

ml2 <- ml3[,-c(31:37)]
rownames(ml2) <- rownames(ml3)
names(ml2)

ml4 <- merge(ml1,ml2[,-1],by=0)
rownames(ml4) <- ml4$Row.names; ml4 <- ml4[,-1]
names(ml4)

ml5 <- merge(ml1,ml3[,-1],by=0)
rownames(ml5) <- ml5$Row.names; ml5 <- ml5[,-1]
names(ml5)

ml6 <- iniml3
## SAVE

#write.csv(ml1,"1_ML_Only_GRS.csv")
#write.csv(ml2,"2_ML_BQpre_AcePre_Imputada.csv")
#write.csv(ml3,"3_ML_BQpre_Acepre_AcePub_Imputada.csv")
#write.csv(ml4,"4_ML_GRSs_BQpre_Acepre_Imputada.csv")
#write.csv(ml5,"5_ML_GRSs_BQpre_Acepre_AcePub_Imputada.csv")
#write.csv(ml6,"6_ML_GRSs_BQpre.csv")


################################################################################
# CORR PLOT
################################################################################
#BiocManager::install(c('faraway','corrplot'))
library(faraway)
library(corrplot)
library(ggcorrplot)
msData <- iniData[which(iniData$Code_new_T2%in%rownames(ms_new)),]

tocorr <- merge(ml5,msData[,c(3:13,15,17:19)],by=0)
rownames(tocorr) <- tocorr$Row.names; tocorr <- tocorr[,-1]
tocorr <- tocorr[-which("Z421"==rownames(tocorr)),]


corrAll <- cor(tocorr)
p.mat <- cor_pmat(tocorr)

pdf("SIG_Corelation Plot of ML included variables and PUBER variables.pdf")
corrplot(corrAll, tl.col = "red", tl.srt = 45, bg = "White",
         title = "\n\n Correlation Plot Of GRSs + BQpre + AcePre + AcePub",
         type = "lower",
         p.mat=p.mat,
         tl.cex=0.5,
       pch.cex=0.3)
dev.off()

################################################################################
#                     HOMA ALTERED DATASETS
################################################################################

msData <- iniData
table(msData$Code_new_T2==rownames(prs))

HOMA_R <- data.frame(HOMA_R=msData$HOMA)
rownames(HOMA_R) <- rownames(msData);
HOMA <- na.omit(HOMA)
dim(HOMA)

ml1 <- merge(HOMA,prs[,c(1,2,5,9)],by=0)
rownames(ml1) <- ml1$Row.names; ml1 <- ml1[,-1]
names(ml1)
dim(ml1)

iniml3 <- merge(HOMA,bq,by=0)
rownames(iniml3) <-iniml3$Row.names; iniml3 <- iniml3[,-1]
ml3 <- merge(iniml3[],acelImp,by=0)
rownames(ml3) <- ml3$Row.names; ml3 <- ml3[,-1]
names(ml3)
dim(ml3)

ml2 <- ml3[,-c(31:37)]
rownames(ml2) <- rownames(ml3)
names(ml2)
dim(ml2)

ml4 <- merge(ml1,ml2[,-1],by=0)
rownames(ml4) <- ml4$Row.names; ml4 <- ml4[,-1]
names(ml4)
dim(ml4)

ml5 <- merge(ml1,ml3[,-1],by=0)
rownames(ml5) <- ml5$Row.names; ml5 <- ml5[,-1]
names(ml5)
dim(ml5)

ml6 <- iniml3
names(ml6)
## SAVE

write.csv(ml1,"3.1.1_HOMA_ML_Only_GRS.csv")
write.csv(ml2,"3.1.2_HOMA_ML_BQpre_AcePre_Imputada.csv")
write.csv(ml3,"3.1.3_HOMA_ML_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml4,"3.1.4_HOMA_ML_GRSs_BQpre_Acepre_Imputada.csv")
write.csv(ml5,"3.1.5_HOMA_ML_GRSs_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml6,"3.1.6_HOMA_ML_GRSs_BQpre.csv")

################################################################################
#                     TRYGLICERIDES ALTERED DATASETS
################################################################################

msData <- iniData
table(msData$Code_new_T2==rownames(prs))

TRYGLICERIDES <- data.frame(TRYGLICERIDES_R=msData$TG_R)
rownames(TRYGLICERIDES) <- rownames(msData);
TRYGLICERIDES <- na.omit(TRYGLICERIDES)
dim(TRYGLICERIDES)

prsTG <- data.frame(PGS002197_Triglycerides=prs[,8])
rownames(prsTG) <- rownames(prs)

ml1 <- merge(TRYGLICERIDES,prsTG,by=0)
rownames(ml1) <- ml1$Row.names; ml1 <- ml1[,-1]
names(ml1)
dim(ml1)

iniml3 <- merge(TRYGLICERIDES,bq,by=0)
rownames(iniml3) <-iniml3$Row.names; iniml3 <- iniml3[,-1]
ml3 <- merge(iniml3[],acelImp,by=0)
rownames(ml3) <- ml3$Row.names; ml3 <- ml3[,-1]
names(ml3)
dim(ml3)

ml2 <- ml3[,-c(31:37)]
rownames(ml2) <- rownames(ml3)
names(ml2)
dim(ml2)

ml4 <- merge(ml1,ml2[,-1],by=0)
rownames(ml4) <- ml4$Row.names; ml4 <- ml4[,-1]
names(ml4)
dim(ml4)

ml5 <- merge(ml1,ml3[,-1],by=0)
rownames(ml5) <- ml5$Row.names; ml5 <- ml5[,-1]
names(ml5)
dim(ml5)

ml6 <- iniml3
names(ml6)
dim(ml6)
## SAVE

write.csv(ml1,"3.2.1_TRYGLICERIDES_ML_Only_GRS.csv")
write.csv(ml2,"3.2.2_TRYGLICERIDES_ML_BQpre_AcePre_Imputada.csv")
write.csv(ml3,"3.2.3_TRYGLICERIDES_ML_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml4,"3.2.4_TRYGLICERIDES_ML_GRSs_BQpre_Acepre_Imputada.csv")
write.csv(ml5,"3.2.5_TRYGLICERIDES_ML_GRSs_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml6,"3.2.6_TRYGLICERIDES_ML_GRSs_BQpre_Acepre_AcePub_Imputada.csv")

################################################################################
#                     HDL ALTERED DATASETS
################################################################################

msData <- iniData
table(msData$Code_new_T2==rownames(prs))

HDL <- data.frame(HDL_R=msData$HDL_R)
rownames(HDL) <- rownames(msData);
HDL <- na.omit(HDL)
dim(HDL)

prsHDL <- data.frame(PGS000686_HDL=prs[,4])
rownames(prsHDL) <- rownames(prs)

ml1 <- merge(HDL,prsHDL,by=0)
rownames(ml1) <- ml1$Row.names; ml1 <- ml1[,-1]
names(ml1)
dim(ml1)

iniml3 <- merge(HDL,bq,by=0)
rownames(iniml3) <-iniml3$Row.names; iniml3 <- iniml3[,-1]
ml3 <- merge(iniml3[],acelImp,by=0)
rownames(ml3) <- ml3$Row.names; ml3 <- ml3[,-1]
names(ml3)
dim(ml3)

ml2 <- ml3[,-c(31:37)]
rownames(ml2) <- rownames(ml3)
names(ml2)
dim(ml2)

ml4 <- merge(ml1,ml2[,-1],by=0)
rownames(ml4) <- ml4$Row.names; ml4 <- ml4[,-1]
names(ml4)
dim(ml4)

ml5 <- merge(ml1,ml3[,-1],by=0)
rownames(ml5) <- ml5$Row.names; ml5 <- ml5[,-1]
names(ml5)
dim(ml5)

ml6 <- iniml3
names(ml6)

## SAVE

write.csv(ml1,"3.3.1_HDL_ML_Only_GRS.csv")
write.csv(ml2,"3.3.2_HDL_ML_BQpre_AcePre_Imputada.csv")
write.csv(ml3,"3.3.3_HDL_ML_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml4,"3.3.4_HDL_ML_GRSs_BQpre_Acepre_Imputada.csv")
write.csv(ml5,"3.3.5_HDL_ML_GRSs_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml6,"3.3.6_HDL_ML_GRSs_BQpre.csv")

################################################################################
#                     OBESITY ALTERED DATASETS
################################################################################

msData <- iniData

OBESITY <- data.frame(OBESITY_R=msData$Obesity_WC)
rownames(OBESITY) <- rownames(msData);
OBESITY <- na.omit(OBESITY)
dim(OBESITY)

prsOBESITY <- data.frame(PGS002033_Obesity=prs[,7])
rownames(prsOBESITY) <- rownames(prs)
head(prsOBESITY) == head(prs$PGS002033_Obesity)

ml1 <- merge(OBESITY,prsOBESITY,by=0)
rownames(ml1) <- ml1$Row.names; ml1 <- ml1[,-1]
names(ml1)
dim(ml1)

iniml3 <- merge(OBESITY,bq,by=0)
rownames(iniml3) <-iniml3$Row.names; iniml3 <- iniml3[,-1]
ml3 <- merge(iniml3[],acelImp,by=0)
rownames(ml3) <- ml3$Row.names; ml3 <- ml3[,-1]
names(ml3)
dim(ml3)

ml2 <- ml3[,-c(31:37)]
rownames(ml2) <- rownames(ml3)
names(ml2)
dim(ml2)

ml4 <- merge(ml1,ml2[,-1],by=0)
rownames(ml4) <- ml4$Row.names; ml4 <- ml4[,-1]
names(ml4)
dim(ml4)

ml5 <- merge(ml1,ml3[,-1],by=0)
rownames(ml5) <- ml5$Row.names; ml5 <- ml5[,-1]
names(ml5)
dim(ml5)

ml6 <- iniml3
names(ml6)
## SAVE

write.csv(ml1,"3.4.1_OBESITY_ML_Only_GRS.csv")
write.csv(ml2,"3.4.2_OBESITY_ML_BQpre_AcePre_Imputada.csv")
write.csv(ml3,"3.4.3_OBESITY_ML_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml4,"3.4.4_OBESITY_ML_GRSs_BQpre_Acepre_Imputada.csv")
write.csv(ml5,"3.4.5_OBESITY_ML_GRSs_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml6,"3.4.6_OBESITY_ML_GRSs_BQpre.csv")

################################################################################
#                     SBP ALTERED DATASETS
################################################################################

msData <- iniData
table(msData$Code_new_T2==rownames(prs))

SBP <- data.frame(SBP_R=msData$SBP_R)
rownames(SBP) <- rownames(msData);
SBP <- na.omit(SBP)
dim(SBP)

prsSBP <- data.frame(PGS002257_SBP=prs[,10])
rownames(prsSBP) <- rownames(prs)
head(prsSBP) == head(prs$PGS002257_SBP)

ml1 <- merge(SBP,prsSBP,by=0)
rownames(ml1) <- ml1$Row.names; ml1 <- ml1[,-1]
names(ml1)
dim(ml1)

iniml3 <- merge(SBP,bq,by=0)
rownames(iniml3) <-iniml3$Row.names; iniml3 <- iniml3[,-1]
ml3 <- merge(iniml3[],acelImp,by=0)
rownames(ml3) <- ml3$Row.names; ml3 <- ml3[,-1]
names(ml3)
dim(ml3)

ml2 <- ml3[,-c(31:37)]
rownames(ml2) <- rownames(ml3)
names(ml2)
dim(ml2)

ml4 <- merge(ml1,ml2[,-1],by=0)
rownames(ml4) <- ml4$Row.names; ml4 <- ml4[,-1]
names(ml4)
dim(ml4)

ml5 <- merge(ml1,ml3[,-1],by=0)
rownames(ml5) <- ml5$Row.names; ml5 <- ml5[,-1]
names(ml5)
dim(ml5)

ml6 <- iniml3
names(ml6)
## SAVE

write.csv(ml1,"3.5.1_SBP_ML_Only_GRS.csv")
write.csv(ml2,"3.5.2_SBP_ML_BQpre_AcePre_Imputada.csv")
write.csv(ml3,"3.5.3_SBP_ML_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml4,"3.5.4_SBP_ML_GRSs_BQpre_Acepre_Imputada.csv")
write.csv(ml5,"3.5.5_SBP_ML_GRSs_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml6,"3.5.6_SBP_ML_GRSs_BQpre.csv")


################################################################################
#                     DBP ALTERED DATASETS
################################################################################

msData <- iniData
table(msData$Code_new_T2==rownames(prs))

DBP <- data.frame(DBP_R=msData$DBP_R)
rownames(DBP) <- rownames(msData);
DBP <- na.omit(DBP)
dim(DBP)

prsDBP <- data.frame(PGS002258_DBP=prs[,11])
rownames(prsDBP) <- rownames(prs)
head(prsDBP) == head(prs$PGS002258_DBP)

ml1 <- merge(DBP,prsDBP,by=0)
rownames(ml1) <- ml1$Row.names; ml1 <- ml1[,-1]
names(ml1)
dim(ml1)

iniml3 <- merge(DBP,bq,by=0)
rownames(iniml3) <-iniml3$Row.names; iniml3 <- iniml3[,-1]
ml3 <- merge(iniml3[],acelImp,by=0)
rownames(ml3) <- ml3$Row.names; ml3 <- ml3[,-1]
names(ml3)
dim(ml3)

ml2 <- ml3[,-c(31:37)]
rownames(ml2) <- rownames(ml3)
names(ml2)
dim(ml2)

ml4 <- merge(ml1,ml2[,-1],by=0)
rownames(ml4) <- ml4$Row.names; ml4 <- ml4[,-1]
names(ml4)
dim(ml4)

ml5 <- merge(ml1,ml3[,-1],by=0)
rownames(ml5) <- ml5$Row.names; ml5 <- ml5[,-1]
names(ml5)
dim(ml5)

ml6 <- iniml3
names(ml6)

## SAVE

write.csv(ml1,"3.6.1_DBP_ML_Only_GRS.csv")
write.csv(ml2,"3.6.2_DBP_ML_BQpre_AcePre_Imputada.csv")
write.csv(ml3,"3.6.3_DBP_ML_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml4,"3.6.4_DBP_ML_GRSs_BQpre_Acepre_Imputada.csv")
write.csv(ml5,"3.6.5_DBP_ML_GRSs_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml2,"3.6.6_DBP_ML_BQpre_AcePre_Imputada.csv")

################################################################################
#                     GLUCOSE ALTERED DATASETS
################################################################################

msData <- iniData
table(msData$Code_new_T2==rownames(prs))

GLUCOSE <- data.frame(GLUCOSE_R=msData$Glucose_R)
rownames(GLUCOSE) <- rownames(msData);
GLUCOSE <- na.omit(GLUCOSE)
dim(GLUCOSE)

prsGLUCOSE <- data.frame(PGS001350_Fasting_Glucose=prs[,6])
rownames(prsGLUCOSE) <- rownames(prs)
head(prsGLUCOSE) == head(prs$PGS001350_Fasting_Glucose)

ml1 <- merge(GLUCOSE,prsGLUCOSE,by=0)
rownames(ml1) <- ml1$Row.names; ml1 <- ml1[,-1]
names(ml1)
dim(ml1)

iniml3 <- merge(GLUCOSE,bq,by=0)
rownames(iniml3) <-iniml3$Row.names; iniml3 <- iniml3[,-1]
ml3 <- merge(iniml3[],acelImp,by=0)
rownames(ml3) <- ml3$Row.names; ml3 <- ml3[,-1]
names(ml3)
dim(ml3)

ml2 <- ml3[,-c(31:37)]
rownames(ml2) <- rownames(ml3)
names(ml2)
dim(ml2)

ml4 <- merge(ml1,ml2[,-1],by=0)
rownames(ml4) <- ml4$Row.names; ml4 <- ml4[,-1]
names(ml4)
dim(ml4)

ml5 <- merge(ml1,ml3[,-1],by=0)
rownames(ml5) <- ml5$Row.names; ml5 <- ml5[,-1]
names(ml5)
dim(ml5)

ml6 <- iniml3
names(ml6)
## SAVE

write.csv(ml1,"3.7.1_GLUCOSE_ML_Only_GRS.csv")
write.csv(ml2,"3.7.2_GLUCOSE_ML_BQpre_AcePre_Imputada.csv")
write.csv(ml3,"3.7.3_GLUCOSE_ML_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml4,"3.7.4_GLUCOSE_ML_GRSs_BQpre_Acepre_Imputada.csv")
write.csv(ml5,"3.7.5_GLUCOSE_ML_GRSs_BQpre_Acepre_AcePub_Imputada.csv")
write.csv(ml6,"3.7.6_GLUCOSE_ML_GRSs_BQpre.csv")



sink(file="3_PRSs_PUBMEP_3_APPROACH_Group_Composition.txt")

table(Metabolic_Syndrome_2COMP)

table(HOMA)

table(TRYGLICERIDES)

table(HDL)

table(OBESITY)

table (SBP)

table (DBP)

table(GLUCOSE)

sink(file=NULL)

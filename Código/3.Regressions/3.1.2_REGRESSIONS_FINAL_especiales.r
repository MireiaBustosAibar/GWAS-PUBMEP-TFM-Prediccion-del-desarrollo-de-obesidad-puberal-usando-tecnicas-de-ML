
#cd /mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/TFM_REGRESSIONS
#setwd("/media/mireia/HardDiskHDD/Documentos/Universidad/Master_Bioinformatica_UOC/Curso_21_22/SEMESTRE_2/TFM/PERFORMING_ANALYSIS/4_REGRESIONS/1_INPUTS")
setwd("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/TFM_REGRESSIONS/1_INPUTS")

library("readxl")
library(tidyverse)
library(broom)

################################################
# Parte 1 . Carga de datos:
################################################

#FIRST CHUNK OF CODE.

iniData <- as.data.frame(read_excel("PUBMEP_FENOTIPOS_RAW.xlsx"));
rownames(iniData) <- iniData$Code_new_T2
dim(na.omit(iniData))

prs <- as.data.frame(read_excel("ALL_PRSs.xlsx"))
rownames(prs) <- prs[,1]; prs <- prs [,-1]

ml5 <- read.csv("3.7.5_GLUCOSE_ML_GRSs_BQpre_Acepre_AcePub_Imputada.csv", header=TRUE, row.names=1)

ml5_sm <- read.csv("5_ML_GRSs_BQpre_Acepre_AcePub_Imputada.csv", header=TRUE, row.names=1);
sm <-data.frame(Metabolic_Syndrome_2COMP=ml5_sm[,1]); rownames(sm)<-rownames(ml5_sm);

table(rownames(ml5) == rownames(iniData))


altura <- read.csv2("BASE_PUBMEP_LONGI_WIDEformat_NIÑAS_Y_NIÑOS_213inds_31_01_2022_EPIC.csv")
rownames(altura) <- altura$Code_new_T2
altura_T1 <- altura[order(altura$Code_new_T2),c(1,18)]
mod1_2Data <- merge(ml5,altura_T1,by=0)
mod1_2Data$Height_T1[c(24,64)] <- c(1.42,1.56)
rownames(mod1_2Data) <- mod1_2Data[,1]; mod1_2Data <- mod1_2Data[,-c(1,2,3,39)]

mod2Data <- merge(prs,mod1_2Data,by=0)
rownames(mod2Data) <- mod2Data[,1]; mod2Data <- mod2Data[,-c(1)]

mod3Data <- merge(mod2Data,iniData,by=0)
rownames(mod3Data) <- mod3Data[,1]; mod3Data <- mod3Data[,-c(1)]

modData <- merge(mod3Data,sm,by=0,all=TRUE)
rownames(modData) <- modData[,1]; modData <- modData[,-1]


table(rownames(modData)==rownames(ml5))
dim(modData)
names(modData)
na_count_bq <- sapply(modData, function(y) sum(length(which(is.na(y)))))
#na_count_bq
na_count_bq[which(na_count_bq>0)]

bpData <- modData[,-c(49,48)]
bpData[,c(13,49,61,65:80)] <- apply(bpData[,c(13,49,61,65:80)],2,as.factor)
allData <- bpData
names(allData)[48:64] <- paste(names(allData)[48:64],"_T2",sep="")


na_count_all <- sapply(bpData, function(y) sum(length(which(is.na(y)))))
#na_count_bq
na_count_all[which(na_count_all>0)]


#####################################################
#LA SIGUIENTE FUNCION ES PARA LUEGO SACAR LOS PLOTS
ggplotLogRegression <- function (fit) {
  require(ggplot2)
  fit$model[,1] <- as.numeric(fit$model[,1])-1
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "glm", col = "red", se=FALSE,
                method.args = list(family = binomial)) +
    labs(title = paste("Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotRegression <- function (fit) {

  require(ggplot2)

  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

#LA SIGUIENTE FUNCION ES PARA VER SI ESTAMOS EN UN INTERVALO DETERMINADO
in_interval <- function(x, interval){
  stopifnot(length(interval) == 2L)
  interval[1] < x & x < interval[2]
}

####################################################################
####################################################################
####################################################################
#
#
#               LINEAR REGRESSION WITHOUT ADJUSTING
#                       WITH BMI Z-SCORE
#
#
####################################################################
####################################################################
####################################################################


# REGRESION LINEAL SIN AJUSTAR POR BMI
#allData[,c(15,60,64)]

# REGRESION LINEAL COMPONENTES CON MAS DE 2 CATEGORIAS
#allData[,c(65,74,76)]

data_T1 <- allData

#MERGED_df_ <- allData[,c(1,13:65)]
MERGED_df_ <- allData[,c(15,60,64)]
#MERGED_df_ <- allData[,c(65,74,76)]
MERGED_df_ <- as.data.frame(lapply(MERGED_df_,as.numeric))

VariablesInteres <- names(data_T1)[1:11]

#####################################################

for (z in (1:length(VariablesInteres))){
  Name <- VariablesInteres[z]
  
  lista_variables <- names(MERGED_df_)
  
  lista_variables <- c(1:length(names(MERGED_df_)))
  
  options(warn=1) #to show warning as they occur.
  
  SEX <- allData$sex_T2
  ZSCORE <- allData$bmi_zscore_T2
  ORIGEN <- allData$Origen_T1
  TALLA <- allData$height_T2
  EDAD <- allData$age_T2
  GRS <- allData[,which(names(allData) %in% VariablesInteres[z])]
  
  
  
  FileOut <- file(paste("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/TFM_REGRESSIONS/2_RESULTS/3_RESULTS_FINAL/3.1_LINEAR/BMI_LIN_REGR_",Name,".csv"), "w")
  for (i in 1:length(names(MERGED_df_))){
    
    VARIABLE_INTERES <- MERGED_df_[,i]
    
    Work <- rbind(VARIABLE_INTERES,ZSCORE,SEX,GRS,ORIGEN,EDAD,TALLA)
    Work <- t(Work)
    nWork<-na.omit(Work)
    nWork <- as.data.frame(nWork)
    nWork <- as.data.frame(lapply(nWork,as.numeric))
    dim(nWork)
    names(nWork)
    
    
    lm1 <- lm(VARIABLE_INTERES~(GRS+as.factor(SEX)+EDAD+as.factor(ORIGEN)),data=nWork)
    
    #    cat(capture.output(summary(lm1), file="Metabolyc_Syndrome_PGS000308_BMIadj_Fasting_Insulin.txt"))
    coef0.lm1 <- summary(lm1)$coeff
    beta.1  <- coef0.lm1[2,1]
    se.1    <- coef0.lm1[2,2]
    t_value.1 <- coef0.lm1[2,3]
    pvalu.1  <- coef0.lm1[2,4]
    #   odr.1   <- exp(beta.1)
    ci.lo.1 <- (beta.1-1.96*se.1)
    ci.hi.1 <- (beta.1+1.96*se.1)
    
    print(summary(lm1))
    print("VARIABLE ESTUDIADA:")
    print(names(MERGED_df_)[i])
    print("Numero de individuos para analisis:")
    print(nrow(nWork))
    print("Tabla SEX")
    print(table(nWork$SEX))
    
    if ( pvalu.1 < 0.05 ) {
      ggplotRegression(lm1)
      ggsave(paste("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/TFM_REGRESSIONS/2_RESULTS/3_RESULTS_FINAL/3.1_LINEAR/3.1.1_LINEAR_PLOTS/",names(MERGED_df_)[i],Name,"regression_plot.pdf",sep="_"))
      psi <- ggstatsplot::ggscatterstats(
        data=nWork,
        x = GRS,
        y = VARIABLE_INTERES,
        xlab = "Polygenic Risk Score",
        ylab = print(names(MERGED_df_)[i]),
        bf.message = TRUE,
        messages = FALSE
      )
      ggsave(paste("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/TFM_REGRESSIONS/2_RESULTS/3_RESULTS_FINAL/3.1_LINEAR/3.1.1_LINEAR_PLOTS/", names(MERGED_df_)[i],Name, "regression_plot_2.pdf", sep="_"), psi)
      
    }else{}
    
    cadena0<-""
    cadena1<-""
    cadena0<-paste("VARIABLE","beta", "SE","ci.lo", "ci.hi", "t-value", "pvalue", sep='\t')
    cadena1<-paste(names(MERGED_df_)[i], round(beta.1,6), round(se.1,6), round(ci.lo.1,6), round(ci.hi.1,6), round(t_value.1,6),pvalu.1, sep='\t')
    cat( cadena0, cadena1,file=FileOut, sep="\n")
  }
  close(FileOut)
}


####################################################################
####################################################################
####################################################################
#
#
#               LINEAR REGRESSION OF CATEGORICAL
#                  WITH MORE THAN 2 CATEGORIES
#
#
####################################################################
####################################################################
####################################################################


# REGRESION LINEAL SIN AJUSTAR POR BMI
#allData[,c(15,60,64)]

# REGRESION LINEAL COMPONENTES CON MAS DE 2 CATEGORIAS
#allData[,c(65,74,76)]

data_T1 <- allData

#MERGED_df_ <- allData[,c(1,13:65)]
#MERGED_df_ <- allData[,c(15,60,64)]
MERGED_df_ <- allData[,c(65,74,76)]
MERGED_df_ <- as.data.frame(lapply(MERGED_df_,as.factor))

VariablesInteres <- names(data_T1)[1:11]

#####################################################

for (z in (1:length(VariablesInteres))){
  Name <- VariablesInteres[z]

  lista_variables <- names(MERGED_df_)

  lista_variables <- c(1:length(names(MERGED_df_)))

  options(warn=1) #to show warning as they occur.

  SEX <- allData$sex
  ZSCORE <- allData$bmi_zscore
  ORIGEN <- allData$Origen_T1
  TALLA <- allData$height
  EDAD <- allData$age
  GRS <- allData[,which(names(allData) %in% VariablesInteres[z])]



FileOut <- file(paste("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/TFM_REGRESSIONS/2_RESULTS/3_RESULTS_FINAL/3.1_LINEAR/CATEG_LIN_REGR_",Name,".csv"), "w")
for (i in 1:length(names(MERGED_df_))){

  VARIABLE_INTERES <- MERGED_df_[,i]

  Work <- rbind(VARIABLE_INTERES,ZSCORE,SEX,GRS,ORIGEN,EDAD,TALLA)
  Work <- t(Work)
  nWork<-na.omit(Work)
  nWork <- as.data.frame(nWork)
  nWork[,c(1,2,4,6)] <- as.data.frame(lapply(nWork[,c(1,2,4,6)],as.numeric))
  nWork[,c(3,5)] <- as.data.frame(lapply(nWork[,c(3,5)],as.factor))
  dim(nWork)
  names(nWork)


    lm1 <- lm(VARIABLE_INTERES~(GRS+SEX+EDAD+ORIGEN+ZSCORE),data=nWork)

#    cat(capture.output(summary(lm1), file="Metabolyc_Syndrome_PGS000308_BMIadj_Fasting_Insulin.txt"))
    coef0.lm1 <- summary(lm1)$coeff
    beta.1  <- coef0.lm1[2,1]
    se.1    <- coef0.lm1[2,2]
    t_value.1 <- coef0.lm1[2,3]
    pvalu.1  <- coef0.lm1[2,4]
 #   odr.1   <- exp(beta.1)
    ci.lo.1 <- (beta.1-1.96*se.1)
    ci.hi.1 <- (beta.1+1.96*se.1)

    print(summary(lm1))
    print("VARIABLE ESTUDIADA:")
    print(names(MERGED_df_)[i])
    print("Numero de individuos para analisis:")
    print(nrow(nWork))
    print("Tabla SEX")
    print(table(nWork$SEX))

    if ( pvalu.1 < 0.05 ) {
      ggplotRegression(lm1)
      ggsave(paste("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/TFM_REGRESSIONS/2_RESULTS/3_RESULTS_FINAL/3.1_LINEAR/3.1.1_LINEAR_PLOTS/CATEG_LIN_",names(MERGED_df_)[i],Name,"regression_plot.pdf",sep="_"))
      psi <- ggstatsplot::ggscatterstats(
        data=nWork,
        x = GRS,
        y = VARIABLE_INTERES,
        xlab = "Polygenic Risk Score",
        ylab = print(names(MERGED_df_)[i]),
        bf.message = TRUE,
        messages = FALSE
      )
      ggsave(paste("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/TFM_REGRESSIONS/2_RESULTS/3_RESULTS_FINAL/3.1_LINEAR/3.1.1_LINEAR_PLOTS/CATEG_LIN_", names(MERGED_df_)[i],Name, "regression_plot_2.pdf", sep="_"), psi)

    }else{}

    cadena0<-""
    cadena1<-""
    cadena0<-paste("VARIABLE","beta", "SE","ci.lo", "ci.hi", "t-value", "pvalue", sep='\t')
    cadena1<-paste(names(MERGED_df_)[i], round(beta.1,6), round(se.1,6), round(ci.lo.1,6), round(ci.hi.1,6), round(t_value.1,6),pvalu.1, sep='\t')
    cat( cadena0, cadena1,file=FileOut, sep="\n")
}
close(FileOut)
}

####################################################################
####################################################################
####################################################################
#
#
#                     LOGISTIC REGRESSION OF 
#                           DBP and SBP
#
#
####################################################################
####################################################################
####################################################################

data_T1 <- allData

# DBP y SBP QUE SE DEBEN AJUSTAR TAMBIÉN POR ALTURA
MERGED_df_ <- allData[,c(68,69)]

#MERGED_df_ <- allData[,c(66,67,70:73,75,77:80)]

VariablesInteres <- names(data_T1)[1:11]

#####################################################

for (z in (1:length(VariablesInteres))){
  Name <- VariablesInteres[z]

  lista_variables <- names(MERGED_df_)

  lista_variables <- c(1:length(names(MERGED_df_)))

  options(warn=1) #to show warning as they occur.

  SEX <- allData$sex_T2
  ZSCORE <- allData$bmi_zscore_T2
  ORIGEN <- allData$Origen_T1
  TALLA <- allData$height_T2
  EDAD <- allData$age_T2
  GRS <- allData[,which(names(allData) %in% VariablesInteres[z])]

  FileOut <- file(paste("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/TFM_REGRESSIONS/2_RESULTS/3_RESULTS_FINAL/3.2_LOGISTIC/BLOOD_PRESS_LOG_REGR_",Name,".csv"), "w")
  for (i in 1:length(names(MERGED_df_))){

    VARIABLE_INTERES <- MERGED_df_[,i]
    Work <- rbind(VARIABLE_INTERES,ZSCORE,SEX,GRS,ORIGEN,EDAD,TALLA)
    Work <- t(Work)
    nWork<-na.omit(Work)
    nWork <- as.data.frame(nWork)
    nWork[,c(2,4,6,7)] <- as.data.frame(lapply(nWork[,c(2,4,6,7)],as.numeric))
    nWork[,c(1,3,5)] <- as.data.frame(lapply(nWork[,c(1,3,5)],as.factor))
    dim(nWork)
    names(nWork)

    lm1 <- glm(formula= VARIABLE_INTERES~(GRS+SEX+EDAD+ORIGEN+TALLA),data=nWork, family=binomial (link=logit))

  #  cat(capture.output(summary(lm1), file="Metabolyc_Syndrome_PGS000308_BMIadj_Fasting_Insulin.txt"))
    coef0.lm1 <- summary(lm1)$coeff
    beta.1  <- coef0.lm1[2,1]
    se.1    <- coef0.lm1[2,2]
    t_value.1 <- coef0.lm1[2,3]
    pvalu.1  <- coef0.lm1[2,4]
    odr.1   <- exp(beta.1)
    ci.lo.1 <- (beta.1-1.96*se.1)
    ci.hi.1 <- (beta.1+1.96*se.1)

    print(summary(lm1))
    print("VARIABLE ESTUDIADA:")
    print(names(MERGED_df_)[i])
    print("Numero de individuos para analisis:")
    print(nrow(nWork))
    print("Tabla SEX")
    print(table(nWork$SEX))


    if ( pvalu.1 < 0.1 ) {
      ggplotLogRegression(lm1)
      ggsave(paste("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/TFM_REGRESSIONS/2_RESULTS/3_RESULTS_FINAL/3.2_LOGISTIC/3.2.1_LOGISTIC_PLOTS/", names(MERGED_df_)[i],Name,"regression_plot.pdf",sep="_"))

    }else{}


    cadena0<-""
    cadena1<-""
    cadena0<-paste("VARIABLE","beta", "SE","ci.lo", "ci.hi", "z-value", "pvalue", "OR", sep='\t')
    cadena1<-paste(names(MERGED_df_)[i], round(beta.1,6), round(se.1,6), round(ci.lo.1,6), round(ci.hi.1,6), round(t_value.1,6), pvalu.1,odr.1, sep='\t')
    cat( cadena0, cadena1,file=FileOut, sep="\n")
  }
  close(FileOut)
}

################################################################################
######## Proseek CVD in T2D ####################################################
## Christoph Nowak script 2 ####################################################
## risk prediction vs NDR risk model ###########################################
################################################################################

library(foreign)
library(xlsx)
library(data.table)
library(survival)
library(caret)
library(gbm)
library(dplyr)
library(OIsurv) 
library(ggplot2)
library(pROC)

setwd("~/Desktop/JÄ")

  # combined cohorts
COMBINE <- read.csv("COMBINE.csv")

  # the NDR calculator has total cholesterol in it (for TC:HDL ratio) and macroalbuminuria, old script didn't, add to cohorts:
cd1 <- read.dta("CARDIPP20161005STATA12.dta")
card <- cd1
card$MicroAlb <- ifelse(is.na(card$UAlb_Krea_index), 0, ifelse((card$UAlb_Krea_index > 0.2 & card$UAlb_Krea_index < 2), 1, 0))
card$MacroAlb <- ifelse((card$UAlb_Krea_index >= 2), 1, 0)
card$kol <- card$R_kol_l
card$kol[which(is.na(card$kol))] <- mean(card$R_kol_l, na.rm = T)

U77 <- read.dta("endostatin ULSAM77 outcome 2010 new alb.dta")
uls <- U77
uls$MicroAlb <- ifelse((U77$uACR > 0.2 & U77$uACR < 2), 1, 0)
uls$MacroAlb <- ifelse((U77$uACR >= 2), 1, 0)
uls$kol <- uls$v972
uls$kol[which(is.na(uls$kol))] <- mean(uls$kol, na.rm = T)

DtPv <- read.dta("Endostatin PIVUS 10 year follow-up 150129.dta")
piv <- DtPv
piv$MicroAlb <- 0 # NA
piv$MacroAlb <- 0
piv$kol <- DtPv$kolesterol

MvCov <- read.csv("MIVC DATABASE CENSORE JUNE 15.csv", dec = ".", na.string = "NA")
miv <- MvCov
miv$MicroAlb <-  ifelse(miv$Alb_Creat_Ratio > 20 & miv$Alb_Creat_Ratio  < 200, 1, 0)
miv$MacroAlb <-  ifelse(miv$Alb_Creat_Ratio >= 200, 1, 0)
miv$kol <- miv$CT / 38.67

dtsv <- read.dta("SAVA CVD outcomes STATA12.dta")
dtsv2 <- read.dta("SAVA kontroller STATA12.dta")
dtsv <- dtsv[(dtsv$serial < 1000), ]
dtsv3 <- merge(x = dtsv, y = dtsv2, by.x = "serial", by.y = "serial", all.x = T)
sav <- dtsv3
sav$MicroAlb <- 0
sav$MacroAlb <- 0
sav$kol <- sav$k_cholesterolTotal

t1 <- merge(COMBINE[which(COMBINE$Coh == "CDF"),], card[,c("loep", "kol", "MicroAlb", "MacroAlb")], by = "loep")
t2 <- merge(COMBINE[which(COMBINE$Coh == "PIV"),], piv[,c("lpnr", "kol", "MicroAlb", "MacroAlb")], by.x = "loep", by.y = "lpnr")
t3 <- merge(COMBINE[which(COMBINE$Coh == "ULS"),], uls[,c("pat", "kol", "MicroAlb", "MacroAlb")], by.x = "loep", by.y = "pat")
t4 <- merge(COMBINE[which(COMBINE$Coh == "Sv_C"),], sav[,c("serial", "kol", "MicroAlb", "MacroAlb")], by.x = "loep", by.y = "serial")
t5 <- merge(COMBINE[which(COMBINE$Coh == "MIVC"),], miv[,c("N", "kol", "MicroAlb", "MacroAlb")], by.x = "loep", by.y = "N")
t6 <- rbind(t1,t2,t3,t4,t5)

COMBINE <- t6

Outcome <- Surv(time = COMBINE$Days_to_MACE_2014, event = COMBINE$Status_MACE_tom_2014) 

COMBINE$Coh_1 <- ifelse(COMBINE$Coh == "PIV", 1, 0)
COMBINE$Coh_2 <- ifelse(COMBINE$Coh == "MIVC", 1, 0)
COMBINE$Coh_3 <- ifelse(COMBINE$Coh == "Sv_C", 1, 0)
COMBINE$Coh_4 <- ifelse(COMBINE$Coh == "ULS", 1, 0)
COMBINE$koen_a <- ifelse(COMBINE$koen_a == "Kvinna", 0, 1)

# save(COMBINE, file = "COMBINE_wholesome.Rdata")


###   About the NDR risk calculator: 5 årsrisken för hjärt-kärlsjukdom; https://www.ndr.nu/IFrameRisk/
###   För T2D: 30-75 yr, 12 variabler, rökning = pågående, mikro/makroalbuminuri är 20-300mg/ml och >300mg/ml kreat,  
###   förmaksflimmer tidigare/pågående, previous CVD = MI/stroke history
###   Original study - https://www.ncbi.nlm.nih.gov/pubmed/21719139 
###   Discovery n = 24,288 T2D, 30-74 yr, 2002 baseline, 5yr FU, retrospective register
###   Validation n = 4,906, 4yr FU, 2003-07, retrospective register
###   Events discovery 2,488/24,288 (10.24%). Events validation 522/4,906 (10.64%)
###   Model: Cox PH. All continuous variables log-transformed, except diabetes-age-of-onset + diabetes-duration, Mean-0-SD-1 scaled
###   12 variables: T2D onset age, T2D duration, log(TC/HDL), log(Hb1Ac_mmol), log(SBP), log(BMI), male-1, smoker-1, microAlb-1, macroAlb-1, AF-1, prev-CVD-1
###   Discovery:  C = 0.71, Sens/Spec upper-50%-predicted-risk  0.740 / 0.556
###                                   top quartile              0.492 / 0.783
###   Validation: C = 0.72, Sens/Spec upper-50%-predicted-risk  0.762 / 0.529
###                                   top quartile              0.512 / 0.779

load("COMBINE_wholesome.Rdata")

  # prep variables as done for NDR calculator
COMBINE$Ageonset <- COMBINE$Aalder - COMBINE$diab_dur
length(which(COMBINE$hba1c_IFCC == 0)) / nrow(COMBINE)  
COMBINE$hba1c_IFCC[which(COMBINE$hba1c_IFCC == 0)] <- mean(COMBINE$hba1c_IFCC[which(COMBINE$hba1c_IFCC != 0)]) 
COMBINE$hdl_kol[COMBINE$Coh == "MIVC"] <- COMBINE$hdl_kol[COMBINE$Coh == "MIVC"] / 38.67 # total chol to mmol/l
COMBINE$MicroAlb.y[is.na(COMBINE$MicroAlb.y)] <- 0
COMBINE$MacroAlb[is.na(COMBINE$MacroAlb)] <- 0
COMBINE$logTCHDL <- log(COMBINE$kol / COMBINE$hdl_kol) # log(TC/HDL)

for(i in grep("A_", names(COMBINE))){
  COMBINE[,i] <- scale(as.numeric(COMBINE[,i]))
  COMBINE[,i] <- (as.numeric(COMBINE[,i]))
}

prot_names <- names(COMBINE)[grep("A_", names(COMBINE))]
pr <- paste(prot_names, collapse = " + ")
COMBINE$Survived2 <- ifelse(COMBINE$Status_MACE_tom_2014 == 1, 'MACE', 'noMACE')
COMBINE$Survived2 <- as.factor(COMBINE$Survived2)

outcomeName <- "Survived2"

formula_base <- formula(paste("Outcome ~ koen_a +  Ageonset + diab_dur + scale(log(hba1c_IFCC)) + smoke +
                              MicroAlb.y + MacroAlb + scale(log(bt_syst_mean_a_B)) + prev_cvd + scale(log(BMI_VC_B)) +
                              scale(logTCHDL) +  flimmer_a + Coh_1 + Coh_2 + Coh_3 + Coh_4 "))

formula_prot <- formula(paste("Outcome ~ koen_a +  Ageonset + diab_dur + scale(log(hba1c_IFCC)) + smoke +
                              MicroAlb.y + MacroAlb + scale(log(bt_syst_mean_a_B)) + prev_cvd + scale(log(BMI_VC_B)) +
                              scale(logTCHDL) +  flimmer_a + Coh_1 + Coh_2 + Coh_3 + Coh_4 +", pr))

  # split for training and validation
set.seed(12345)
splitIndex <- createDataPartition(COMBINE$Status_MACE_tom_2014, p = .75, list = FALSE, times = 1)
trainDF <- COMBINE[ splitIndex,]
testDF  <- COMBINE[-splitIndex,]

Outcome <- Surv(time = trainDF$Days_to_MACE_2014, event = trainDF$Status_MACE_tom_2014) 

set.seed(12345)
gbm_base <- gbm(formula_base,
                data = trainDF,                   
                n.trees = 1000,               
                shrinkage = 0.01,            
                interaction.depth = 2, # to leave biological interpretation & limit complexity         
                bag.fraction = 0.5,          
                train.fraction = 0.7,       
                n.minobsinnode = 5)      
set.seed(12345)
gbm_prot <- gbm(formula_prot,
                data = trainDF,                   
                n.trees = 1000,               
                shrinkage = 0.01,            
                interaction.depth = 2,         
                bag.fraction = 0.5,          
                train.fraction = 0.7,       
                n.minobsinnode = 5)     

best.iter1 <- gbm.perf(gbm_base, method = "test") # choose  best iteration
best.iter2 <- gbm.perf(gbm_prot, method = "test")

  # custom function to convert prediction back to hazards
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

  # Test set validation
predict1 <- logit2prob(predict(gbm_base, testDF, best.iter1))
predict2 <- logit2prob(predict(gbm_prot, testDF, best.iter2))
s1 <- data.frame(cbind((predict1), 1-(predict1)))
s2 <- data.frame(cbind((predict2), 1-(predict2)))
names(s1) <- c("MACE", "noMACE")
names(s2) <- c("MACE", "noMACE")
auc1 <- roc(ifelse(testDF$Survived2 == "MACE", 1, 0), s1[[2]], ci = T)
auc2 <- roc(ifelse(testDF$Survived2 == "MACE", 1, 0), s2[[2]], ci = T)
# plot(auc1)   
# plot(auc2) 
print(auc1) 
print(auc2)

  # Training set
predict1 <- logit2prob(predict(gbm_base, trainDF, best.iter1))
predict2 <- logit2prob(predict(gbm_prot, trainDF, best.iter2))
s1 <- data.frame(cbind((predict1), 1-(predict1)))
s2 <- data.frame(cbind((predict2), 1-(predict2)))
names(s1) <- c("MACE", "noMACE")
names(s2) <- c("MACE", "noMACE")
auc1train <- roc(ifelse(trainDF$Survived2 == "MACE", 1, 0), s1[[2]], ci = T)
auc2train <- roc(ifelse(trainDF$Survived2 == "MACE", 1, 0), s2[[2]], ci = T)
# plot(auc1train) 
# plot(auc2train) 
print(auc1train) 
print(auc2train)

  # (both combined)
predict1 <- logit2prob(predict(gbm_base, COMBINE, best.iter1))
predict2 <- logit2prob(predict(gbm_prot, COMBINE, best.iter2))
s1 <- data.frame(cbind((predict1), 1-(predict1)))
s2 <- data.frame(cbind((predict2), 1-(predict2)))
names(s1) <- c("MACE", "noMACE")
names(s2) <- c("MACE", "noMACE")
auc1all <- roc(ifelse(COMBINE$Survived2 == "MACE", 1, 0), s1[[2]], ci = T)
auc2all <- roc(ifelse(COMBINE$Survived2 == "MACE", 1, 0), s2[[2]], ci = T)
# plot(auc1all) 
# plot(auc2all)  
print(auc1all) 
print(auc2all)

  # Table for export & calculcate sensitivity/specificity for 50% and 25% as in NDR paper
t1_50 <- data.frame(cbind(auc1$sensitivities, auc1$specificities, auc1$thresholds))
t2_50 <- data.frame(cbind(auc2$sensitivities, auc2$specificities, auc2$thresholds))
f1 <- t1_50[which(abs(t1_50[,3] - quantile(auc1$original.predictor, 0.5)) == 
        min(abs(t1_50[,3] - quantile(auc1$original.predictor, 0.5)))), 1:2]
f2 <- t2_50[which(abs(t2_50[,3] - quantile(auc2$original.predictor, 0.5)) == 
        min(abs(t2_50[,3] - quantile(auc2$original.predictor, 0.5)))), 1:2]
f3 <- t1_50[which(abs(t1_50[,3] - quantile(auc1$original.predictor, 0.25)) == 
              min(abs(t1_50[,3] - quantile(auc1$original.predictor, 0.25)))), 1:2]
f4 <- t2_50[which(abs(t2_50[,3] - quantile(auc2$original.predictor, 0.25)) == 
              min(abs(t2_50[,3] - quantile(auc2$original.predictor, 0.25)))), 1:2]

train1_50 <- data.frame(cbind(auc1train$sensitivities, auc1train$specificities, auc1train$thresholds))
train2_50 <- data.frame(cbind(auc2train$sensitivities, auc2train$specificities, auc2train$thresholds))
g1 <- train1_50[which(abs(train1_50[,3] - quantile(auc1train$original.predictor, 0.5)) == 
              min(abs(train1_50[,3] - quantile(auc1train$original.predictor, 0.5)))), 1:2]
g2 <- train2_50[which(abs(train2_50[,3] - quantile(auc2train$original.predictor, 0.5)) == 
              min(abs(train2_50[,3] - quantile(auc2train$original.predictor, 0.5)))), 1:2]
g3 <- train1_50[which(abs(train1_50[,3] - quantile(auc1train$original.predictor, 0.25)) == 
              min(abs(train1_50[,3] - quantile(auc1train$original.predictor, 0.25)))), 1:2]
g4 <- train2_50[which(abs(train2_50[,3] - quantile(auc2train$original.predictor, 0.25)) == 
              min(abs(train2_50[,3] - quantile(auc2train$original.predictor, 0.25)))), 1:2]

MODELING_RES <- data.frame(matrix(ncol = 10, nrow = 4))
names(MODELING_RES) <- c("Model","C_index", "LCI", "UCI", "50%_sens", "50%_spec","75%_sens", "75%_spec", "cases", "controls")

MODELING_RES[1, ] <- c("test 25% NDK", auc1$auc[1], auc1$ci[1], auc1$ci[3], f1[1,1], f1[1,2], f3[1,1], f3[1,2], length(auc1$cases), length(auc1$control))
MODELING_RES[2, ] <- c("test 25% NDK + proteins", auc2$auc[1], auc2$ci[1], auc2$ci[3], f2[1,1], f2[1,2], f4[1,1], f4[1,2], length(auc2$cases), length(auc2$control))
MODELING_RES[3, ] <- c("train 75% NDK", auc1train$auc[1], auc1train$ci[1], auc1train$ci[3], g1[1,1], g1[1,2], g3[1,1], g3[1,2], length(auc1train$cases), length(auc1train$control))
MODELING_RES[4, ] <- c("train 75% NDK + proteins", auc2train$auc[1], auc2train$ci[1], auc2train$ci[3], g2[1,1], g2[1,2], g4[1,1], g4[1,2], length(auc2train$cases), length(auc2train$control))

# View(MODELING_RES)

#  write.csv(MODELING_RES, "GBM_Results.csv")

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

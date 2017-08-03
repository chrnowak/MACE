################################################################################
######## Proseek CVD in T2D ####################################################
## Christoph Nowak script 1 ####################################################
################################################################################

library(foreign)
library(xlsx)
library(data.table)
library(nephro)
library(mice)
library(lme4)
library(survival)
library(caret)
library(rms)
library(powerSurvEpi)

setwd("~/Desktop/JÄ")

##########################################################################
######## CARDIPP #########################################################
##########################################################################

cd1 <- read.dta("CARDIPP20161005STATA12.dta")
  # female = 1; NAs either -9 or NA 

  # missing protein values = NA but 2 proteins have a few "Excluded" values. To avoid
  # later confusion in LOD/2-replacement, replace immediately by LOD/2
cd1$A_112_IL1ra[which(cd1$A_112_IL1ra == "Excluded")] <- 0.6 / 2 
cd1$A_121_hK11[which(cd1$A_121_hK11 == "Excluded")] <- 1.13 / 2
  
cd1$datum_d <- as.numeric(strtrim(as.character(cd1$datum_a), 4))  
cd1$diab_dur <- ifelse(is.na(cd1$diab_debut_a), NA, ifelse(cd1$diab_debut_a == -9, NA, as.numeric(cd1$datum_d - cd1$diab_debut_a)))
  # 48 NAs for diabetes duration 
cd1$diab_dur <- ifelse(is.na(cd1$diab_dur), mean(na.omit(cd1$diab_dur)), cd1$diab_dur)
  # prevalent CVD = MI or stroke (as available in all cohorts)
cd1$prev_cvd <- ifelse((cd1$stroke_a == 1 | cd1$hjaertinfarkt_a == 1), 1, 0)
cd1$prev_cvd <- ifelse(is.na(cd1$prev_cvd) , 0, cd1$prev_cvd) 

  # Microalbuminuria =   NICE guideline ACR >2.5  men, >3.5  women
cd1$MicroAlb <- ifelse(is.na(cd1$UAlb_Krea_index), 0, ifelse(cd1$UAlb_Krea_index < 3, 0, 1))
cd1$MicroAlb <- ifelse(cd1$koen_a == 1, 
                           ifelse(is.na(cd1$UAlb_Krea_index), 0, ifelse(cd1$UAlb_Krea_index < 2.5, 0, 1)),
                           ifelse(is.na(cd1$UAlb_Krea_index), 0, ifelse(cd1$UAlb_Krea_index < 3.5, 0, 1)))

cd1$smoke <- ifelse(is.na(cd1$roekning_p), 0, ifelse(cd1$roekning_p < 3, 0, 1)) 
  # 29 NAs for HDL 
cd1$hdl_kol_l <- ifelse((is.na(cd1$hdl_kol_l) | cd1$hdl_kol_l == -9), mean(na.omit(cd1$hdl_kol_l)), cd1$hdl_kol_l)
  # 32 NAs for LDL
cd1$ldl_kol <- ifelse(is.na(cd1$ldl_kol), mean(na.omit(cd1$ldl_kol)), cd1$ldl_kol)
  # 13  NAs HbA1c
cd1$hba1c_IFCC <- ifelse(is.na(cd1$hba1c_IFCC),mean(na.omit(cd1$hba1c_IFCC)), cd1$hba1c_IFCC)
  # AF
cd1$flimmer_a <- ifelse(is.na(cd1$flimmer_a), 0, ifelse(cd1$flimmer_a == -9, 0, cd1$flimmer_a))
  # antidiabetic drugs
cd1$ADMed <- 0
cd1$ADMed <- ifelse(cd1$insulin_a > 0, 1, cd1$ADMed)
cd1$ADMed <- ifelse(cd1$metformin_a == "Ja", 1, cd1$ADMed)
cd1$ADMed <- ifelse(cd1$glibenklamid_a == "Ja", 1, cd1$ADMed)
cd1$ADMed <- ifelse(cd1$glipizid_a == "Ja", 1, cd1$ADMed)
cd1$ADMed <- ifelse(cd1$glimeperid_a == "Ja", 1, cd1$ADMed)
cd1$ADMed <- ifelse(cd1$akarbos_a == "Ja", 1, cd1$ADMed)
cd1$ADMed <- ifelse(cd1$rosiglitazon_a == "Ja", 1, cd1$ADMed)
cd1$ADMed <- ifelse(cd1$repaglinid_a == "Ja", 1, cd1$ADMed)
cd1$ADMed <- ifelse(cd1$nateglinid_a == "Ja", 1, cd1$ADMed)
cd1$ADMed <- ifelse(cd1$sitagliptin_a == "Ja", 1, cd1$ADMed)
cd1$ADMed <- ifelse(cd1$pioglitazon_a == "Ja", 1, cd1$ADMed)
cd1$ADMed[which(is.na(cd1$ADMed))] <- 0 

  # antihypertensive drugs
cd1$AHptMed <- 0
cd1$AHptMed <- ifelse(cd1$ca_a == 1, 1, cd1$AHptMed)
cd1$AHptMed <- ifelse(cd1$betablock_a == "Ja", 1, cd1$AHptMed)
cd1$AHptMed <- ifelse(cd1$thiazid_a == 1, 1, cd1$AHptMed)
cd1$AHptMed <- ifelse(cd1$loop_a ==1, 1, cd1$AHptMed)
cd1$AHptMed[(which(is.na(cd1$AHptMed)))] <- 0

  # proteins: in original script all copy/pasted by hand - liable to error. Done by importing values from original files now.
prot <- names(cd1)[grep("^A_1", names(cd1))]
  # Covariates (names from old script)
covt <- c("loep", "Status_MACE_tom_2014", "Days_to_MACE_2014", "Aalder" ,"koen_a", "hba1c_IFCC", "BMI_VC_B", "smoke", "MicroAlb", 
          "prev_cvd", "bt_syst_mean_a_B", "ldl_kol", "hdl_kol_l", "diab_dur", "flimmer_a", "CKD_EPI_GFR_B", "ADMed", "AHptMed", "statin_a")
dtcd1 <- cd1[, c(covt, prot)]  
  
  # replace "-9" with NA
for(i in which(names(dtcd1) %in% covt)){
  dtcd1[,i] <- ifelse((is.factor(dtcd1[,i]) & (dtcd1[,i] == "-9")), NA, dtcd1[,i])
  dtcd1[,i] <- ifelse((is.numeric(dtcd1[,i]) & (dtcd1[,i] == -9)), NA, dtcd1[,i])
}

dtcd1$koen_a <- ifelse(dtcd1$koen_a == 1, "Man", "Kvinna")

  # 53 people have NAs in all proteins, exclude
dtcd1 <- dtcd1[!is.na(dtcd1$A_101_IL8), ]

  # LOD/2 replacement: 92 LOD values written in by hand in old script, now imported 
proseek_card <- data.frame(fread("CARDIPP resultat final.csv"))
load <- as.numeric(proseek_card[1, 2:93])

  # New: investigate % missing proteins
temp <- c()
pos <- grep("^A_", names(dtcd1))
for(i in pos){
  temp[i] <- length(which(dtcd1[,i] == "NAN" | is.na(dtcd1[,i])))
}
temp <- data.frame(cbind(names(dtcd1)[pos], temp[20:length(temp)] / nrow(dtcd1)))
temp2 <- temp[which(as.numeric(as.character(temp$X2)) > 0), ]

 # Histogram of % missing
# pdf("missing_prot_cardipp.pdf", 8, 5 )
barplot(as.numeric(as.character(temp2$X2)), names.arg = gsub("A_", "",temp2$X1),
        ylim = c(0,1), cex.names = 0.5, main = "missing protein in CARDIPP n=710\n exclude IL4,BNP,ITG...", cex =0.7, cex.main = 1)
abline(0.15, 0, col = "red")
# dev.off()

  ## exclude >15% missing proteins: IL4, BNP, ITG
excl_prot <- temp2$X1[which(as.numeric(as.character(temp2$X2)) > 0.15)]
dtcd1 <- dtcd1[, -which(names(dtcd1) %in% excl_prot)]

  ## replace other protein NAs with LOD / 2
  ## works, as in all files, proteins are in the same order, verifiable for 101_xxx, 102_xxx etc.
load2 <- data.frame(load[-which(gsub("-", "", proseek_card[2, ]) %in% gsub("A_", "", excl_prot))])
pos <- grep("^A_", names(dtcd1))
load2$pos <- pos
for (i in pos){
  dtcd1[which(dtcd1[,i] == "NAN"), i] <- load2[which(load2[,2] == i), 1] / 2
}

  # old script does dtcd1<-na.omit(dtcd1), reduces n=710 to 656
  # better: impute through MICE regression on all other covariates but don't imput outcome - 2 NAs in MACE - exclude so n=710 -> 708
dtcd1 <- dtcd1[-which(is.na(dtcd1$Status_MACE_tom_2014)), ]
  
pos <- grep("^A_", names(dtcd1))
cov_names <- names(dtcd1)[-pos]
d <- dtcd1[, cov_names]
  
  # Multivariate Imputation by Chained Equations (MICE) -  predictive mean matching, 50 iterations
imp <- mice(d, m = 1, method = "pmm", maxit = 50, seed = 500)
Datimp <- complete(imp, "long", include = F)
  # check:
# View(cbind(d$ldl_kol, Datimp$ldl_kol)) # looks good!

dtcd1[, which(names(dtcd1)  %in% names(Datimp))] <- Datimp[, 3:21]

dtcd1$Coh <- "CDF"

save(dtcd1, file = "CARDIPP_CN.Rdata")

################################################################################################
### DISCOVERY IN CARDIPP #######################################################################
################################################################################################

  # sex indicator 1 == Male
  # Null model: Age + Sex
CARDIPP <- dtcd1
CARDIPP$koen_a <- ifelse(CARDIPP$koen_a == "Kvinna", 0, 1) 
Outcome <- Surv(time = CARDIPP$Days_to_MACE_2014, event = CARDIPP$Status_MACE_tom_2014) 
Cox_result <- coxph(Outcome ~ koen_a + Aalder, data = CARDIPP) 

summary(Cox_result)   
cox.zph(Cox_result)   # PH assumptions: no violation, p = 0.1208
plot(cox.zph(Cox_result))

Kaplanplot = survfit(Outcome ~ koen_a , data = CARDIPP)
plot(Kaplanplot) 

  # Full model for each protein: scale(protein) + age + sex
index <- grep("^A_", names(CARDIPP))
discovery_age_sex <- data.frame(matrix(nrow = length(index), ncol = 6,
                                       dimnames = list(names(CARDIPP)[index], c("b", "se", "p", "HR", "LCI", "UCI"))))

for(i in index){
  CARDIPP[,i] <- scale(as.numeric(CARDIPP[,i]))
  CARDIPP[,i] <- (as.numeric(CARDIPP[,i]))
}

j <- 1
for (i in index){
  discovery_age_sex$b[j] <- summary(coxph(Outcome ~ CARDIPP[, i] + koen_a + Aalder, data = CARDIPP))$coeff[1,1]
  discovery_age_sex$se[j] <- summary(coxph(Outcome ~ CARDIPP[, i] + koen_a + Aalder, data = CARDIPP))$coeff[1,3]
  discovery_age_sex$p[j] <- summary(coxph(Outcome ~ CARDIPP[, i] + koen_a + Aalder, data = CARDIPP))$coeff[1,5]
  discovery_age_sex$HR[j] <- exp(discovery_age_sex$b[j])
  discovery_age_sex$LCI[j] <- exp(discovery_age_sex$b[j] - 1.96*discovery_age_sex$se[j])
  discovery_age_sex$UCI[j] <- exp(discovery_age_sex$b[j] + 1.96*discovery_age_sex$se[j])
  j <- j + 1
}

discovery_age_sex <- discovery_age_sex[order(discovery_age_sex$p, decreasing = F),]
discovery_age_sex$p_fdr <- p.adjust(discovery_age_sex$p, method = "fdr")  # False Discovery Rate
discovery_age_sex$PASS <- "fail"
discovery_age_sex$PASS[which(discovery_age_sex$p_fdr < 0.05 & discovery_age_sex$b > 0)] <- "increase"
discovery_age_sex$PASS[which(discovery_age_sex$p_fdr < 0.05 & discovery_age_sex$b < 0)] <- "decrease"

# View(discovery_age_sex)

# write.csv(discovery_age_sex, "CARDIPP_discovery_result.csv")

  # interprotein correlation - needed for power
cor1 <- cor(CARDIPP[,index])
mean(cor1[which(cor1 != 1)])  # mean correlation r = 0.299 (sd = 0.1538)
sd(cor1[which(cor1 != 1)])

  # Power calculation posthoc:
powerEpiCont(formula = CARDIPP[,20] ~ koen_a + Aalder, dat = CARDIPP, X1 = CARDIPP[,20], failureFlag = CARDIPP$Status_MACE_tom_2014,
             n = 107, theta = exp(1), alpha = 0.05)
  # Given the observed data, power to detect the observed HR is 90.551%

  # A prior power for mean-0-sd-1-protein, 80% power to ...
powerEpiCont.default(n = 708, theta = 1.42, sigma2 = 1,
                     psi = 71 / 708, rho2 = 0.299^2, alpha = 0.05) 
  # ... detect HR of 1.42

#############################################################################################
##### ULSAM 77 ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
#############################################################################################

U77 <- read.dta("endostatin ULSAM77 outcome 2010 new alb.dta")
UCVD <- read.dta("ULSAM77 cvdprev.dta")
UKr <- read.xlsx("Ulsam 77 Krea och Tim 1kim1 skickat.xlsx", sheetName = "Blad2")
U77 <- merge(U77, UCVD, by = "pat" , all.x = T)
U77 <- merge(U77, UKr, by.x = "pat", by.y = "Provid", all.x = T)
U77$DateMACE <- "20101231"
U77$DateMACE <- ifelse(U77$p001 == 1, U77$p013, U77$DateMACE)
U77$DateMACE <- ifelse(U77$p056 > 0, U77$p064, U77$DateMACE)
U77$DateMACE <- ifelse(U77$p074 > 0, U77$p077, U77$DateMACE)
U77$dat77 <- ifelse(U77$dat77 == "", NA, U77$dat77)
U77$Days_to_MACE_2014 <- difftime(as.Date(U77$DateMACE, "%Y%m%d"), as.Date(U77$dat77, "%Y%m%d"), unit = "days")
U77$Days_to_MACE_2014 <- ifelse(U77$Days_to_MACE_2014 < 0, NA, U77$Days_to_MACE_2014)
U77$Status_MACE_tom_2014 <- ifelse((U77$p042 == 1 | U77$p056 == 1 | U77$p074 == 1), 1, 0) 
U77$koen_a <- rep("Man", nrow(U77))
U77$hba1c_IFCC <- 10.93 * U77$hba1c - 23.5
  
  # added: NA replacement through MICE, not simple mean
U77$smoke <- ifelse(U77$a16 == 1, 1, 0)
U77$smoke <- ifelse(is.na(U77$smoke), 0, U77$smoke)

  # new ACR threshold by gender
U77$MicroAlb <- ifelse(is.na(U77$uACR), 0, ifelse(U77$uACR < 2.5, 0, 1))
U77$ldl_kol <- U77$v324
U77$hdl_kol_l<-U77$v302
U77$diab_dur <- 0 # impute later
U77$flimmer_a <- ifelse(is.na(U77$v730), 0, as.numeric(U77$v730 == 1)) 
U77$CKD_EPI_GFR_B <- CKDEpi.creat(creatinine = (U77$Krea_umol_L / 76.26), sex = rep(1, nrow(U77)), 
                                age = U77$age77, ethnicity = rep(0, nrow(U77)))

  # old script does not code for aHT/antiT2D (put "0"). But as variable included, it's added now
U77$ADMed <- ifelse(is.na(U77$v405), 0, ifelse(U77$v405 == 0, 0, 1))
U77$AHptMed <- ifelse(is.na(U77$v102), 0, ifelse(U77$v102 == 0, 0, 1))
U77$statin_a <- ifelse(is.na(U77$v105), 0, ifelse(U77$v105 == 0, 0, 1))

  # old script defines diabetes by the questionnaire item only, now defined also by HbA1c and FG and Medx
diab_id <- (which(U77$v378 == 1 | U77$hba1c > 6.5 | U77$pglukos >= 7 | U77$ADMed == 1)) 

covtU <- c("pat", "Status_MACE_tom_2014", "Days_to_MACE_2014", "age77", "koen_a", "hba1c_IFCC", "v290","smoke",
           "MicroAlb","cvdprev2","v013","ldl_kol","hdl_kol_l","diab_dur","flimmer_a","CKD_EPI_GFR_B","ADMed","AHptMed","statin_a")
dtU <- U77[diab_id, covtU]

  # remove those not in 77yr assessment - those with NA in pat77/mace days
dtU <- dtU[-which(is.na(dtU$Days_to_MACE_2014)), ]

  # added: impute missing by MICE 
d_u <- dtU[,covtU]
imp_u <- mice(d_u, m = 1, method = "pmm", maxit = 50, seed = 500)
Datimp_u <- complete(imp_u, "long", include = F)
dtU[, which(names(dtU)  %in% names(Datimp_u))] <- Datimp_u[, 3:21]
row.names(dtU) <- NULL
names(dtU) <- covt

  # Read proteomics data
dtUP <- read.csv2("ULSAM_140617.csv", header = T, na.string = "NAN", dec = ".") 
  # Added: investigate % missing proteins
temp <- c(rep(0,92))
for(i in 2:93){
  temp[i-1] <- length(which(dtUP[,i] == "NAN" | is.na(dtUP[,i])))
}
temp <- data.frame(cbind(names(dtUP)[2:93], temp / nrow(dtUP)))
temp2 <- temp[which(as.numeric(as.character(temp$X2)) > 0), ]

# pdf("missing_prot_ulsam.pdf", 15, 5 )
barplot(as.numeric(as.character(temp2$X2)), names.arg = gsub("A_", "",temp2$X1),
        ylim = c(0,1), cex.names = 0.5, main = "missing protein in ULSAM n=768", cex =0.7, cex.main = 1)
abline(0.15, 0, col = "red")
# dev.off()

excl_u <- temp[which(as.numeric(as.character(temp$X2)) >= 0.15), ]$X1
  # 15% missing: BNP, ITGBxx, IL4 (as in cardip) plus: CTSB (22%), HSP27 (72%), PTX3 (15.4%), mAmp (26%)

  # LOD / 2 replace. old scipt: written in by hand, now imported
proseek_ulsam <- data.frame(fread("Resultat ULSAM 140617.csv"))
  # messy: protein differently labeled in different files: IL4, IL.4, IL_4 etc., so:
proseek_ulsam_excl <- proseek_ulsam[,-which(substr(proseek_ulsam[2,],1,3) %in% substr(excl_u,2,4))]
dtUP_excl <- dtUP[,-which(substr(names(dtUP),2,4) %in% substr(excl_u,2,4))]
load_u <- (proseek_ulsam_excl[1, 2:ncol(proseek_ulsam_excl)])
load_u <- as.numeric(gsub(",", ".", load_u)) # comma = dot
for (i in 2:ncol(dtUP_excl)){
  dtUP_excl[which(is.na(dtUP_excl[,i])), i] <- load_u[i-1] / 2
}
  
dtU <- merge(x = dtU, y = dtUP_excl, by.x = "loep", by.y = "NPX", all = F, all.x = T, all.y = F)

prot_excl_u <- prot[-which(substr(prot, 3,5) %in% substr(excl_u,2,4))] # exclude names from excluded proteins
names(dtU) <- c(covt, prot_excl_u)
dtU$Coh <- "ULS"

dtU <- na.omit(dtU) ## ok to exlude missing MACE people, 90 -> 86

save(dtU, file = "ULSAM_CN.Rdata")

########################## ########################## ########################## ##########################
#### PIVUS 70 ####### ########################## ########################## ########################## ####
########################## ########################## ########################## ##########################

DtPv <- read.dta("Endostatin PIVUS 10 year follow-up 150129.dta")
DtPv$Days_to_MACE_2014 <- difftime(as.Date(DtPv$allCVDdatestata, "%Y-%m-%d"), as.Date(DtPv$beskdat70mdy, "%Y-%m-%d"), unit = "days")
DtPv$Days_to_MACE_2014 <- ifelse(DtPv$Days_to_MACE_2014 < 0, NA, DtPv$Days_to_MACE_2014)

  # T2D definition, blodsocker + 11% for plasma values
DtPv <- DtPv[which(DtPv$diabetesmellitus == 1 | DtPv$blodsocker*1.11 >=7 | DtPv$oralaantidiabetika == 1 | DtPv$insulin == 1), ]
DtPv$koen <- ifelse(DtPv$kn == 1, "Kvinna", "Man")
DtPv$Status_MACE_tom_2014 <- as.numeric((DtPv$allCVDdatestata) < (DtPv$exitdatestata))
  
DtPv$hba1c_IFCC <- NA  # old script put "0" - no good for later analysis, better impute it later, otherwise 0 distorts picture
DtPv$MicroAlb <- NA
DtPv$diab_dur <- ifelse(is.na(DtPv$rmeddiabetes), 0, DtPv$rmeddiabetes) # impute later
DtPv$CKD_EPI_GFR_B <- CKDEpi.creat(creatinine = (DtPv$kre / 76.26), sex = as.numeric(DtPv$kn == 0),
                                   age = DtPv$exaktlder70, ethnicity = rep(0, nrow(DtPv)))
DtPv$CKD_EPI_GFR_B <- ifelse(is.na(DtPv$CKD_EPI_GFR_B), mean(na.omit(DtPv$CKD_EPI_GFR_B)), DtPv$CKD_EPI_GFR_B) # leave mean repl. as only 2 NAs

DtPv$ADMed <- ifelse(((DtPv$insulin == 1 )|(DtPv$oralaantidiabetika == 1)), 1, 0)
DtPv$AHptMed <- 0
DtPv$AHptMed <- ifelse(DtPv$acehmmare == 1, 1, DtPv$AHptMed)
DtPv$AHptMed <- ifelse(DtPv$atr1atagonist == 1, 1, DtPv$AHptMed)
DtPv$AHptMed <- ifelse(DtPv$caantagonist == 1, 1, DtPv$AHptMed)
DtPv$AHptMed <- ifelse(DtPv$betablock == 1, 1, DtPv$AHptMed)
DtPv$AHptMed <- ifelse(DtPv$diuretika == 1, 1, DtPv$AHptMed)
DtPv$statin_a <- ifelse((DtPv$statiner == 1 | DtPv$andralipidmed == 1), 1, 0)

covtP <- c("lpnr", "Status_MACE_tom_2014", "Days_to_MACE_2014", "exaktlder70", "koen", "hba1c_IFCC", "bmi", "rkarenu", "MicroAlb",
           "cvdiagnos", "manuelltsbp", "ldl", "hdl", "diab_dur", "afib", "CKD_EPI_GFR_B", "ADMed", "AHptMed", "statin_a")
Pcov <- DtPv[ , covtP]

  # proteins
PrPv <- read.csv2("PivusMDAFinal.csv", header = T, na.string = ".", dec = ".")
  # Added: investigate % missing proteins
temp <- c(rep(0,92))
for(i in 2:93){
  temp[i-1] <- length(which(PrPv[,i] == "NA" | is.na(PrPv[,i])))
}
temp <- data.frame(cbind(names(PrPv)[2:93], temp / nrow(PrPv)))
temp2 <- temp[which(as.numeric(as.character(temp$X2)) > 0), ]
# pdf("missing_prot_pivus.pdf", 20, 5 )
barplot(as.numeric(as.character(temp2$X2)), names.arg = gsub("A_", "",temp2$X1),
        ylim = c(0,1), cex.names = 0.4, main = "missing protein in pivus n=1,003", cex =0.7, cex.main = 1)
abline(0.15, 0, col = "red")
# dev.off()

excl_p <- temp[which(as.numeric(as.character(temp$X2)) >= 0.15), ]$X1
  # 15% missing: BNP, ITGBxx, IL4 plus: PTX3 (19%), betaNGF (92%), mAmp (24%), MMP7 (16.3), SIRT2 (69%), NEMO (58%), NTproBNP (17%)

  # LOD / 2 replacement, again old script had LODs written in by hand
proseek_pivus <- data.frame(fread("Resultat Pivus MDA - Final.csv", header = F))
proseek_pivus_excl <- proseek_pivus[,-which(substr(proseek_pivus[2, ], 1, 3) %in% substr(excl_p, 2, 4))]
PrPv_excl <- PrPv[,-which(substr(names(PrPv), 2, 4) %in% substr(excl_p, 2, 4))]
load_p <- (proseek_pivus_excl[1, 2:ncol(proseek_pivus_excl)])
load_p <- as.numeric(gsub(",", ".", load_p)) # comma = dot
for (i in 2:ncol(PrPv_excl)){
  PrPv_excl[which(is.na(PrPv_excl[,i])), i] <- load_p[i-1] / 2 # mistake in old script, had for P70 only LOD/sqrt(LOD)
}

dtPv <- merge(x = Pcov, y = PrPv_excl, by.x = "lpnr", by.y = "NPX", all = F, all.x = T, all.y = F)
  
  # for now: put NAs to 0, impute later from whole dataset
dtPv$hba1c_IFCC <- 0
dtPv$MicroAlb <- 0
dtPv <- na.omit(dtPv)
prot_excl_p <- prot[-which(substr(prot, 3, 5) %in% substr(excl_p, 2, 4))] # exclude names from excluded proteins
names(dtPv) <- c(covt, prot_excl_p) 

dtPv$Coh <- "PIV" 

save(dtPv, file = "PIVUS_CN.Rdata")

######################################### ######################################### 
###### MIVC ###### ###### ###### ###### ###### ###### ###### ###### ###### #######
######################################### ######################################### 

MvCov <- read.csv("MIVC DATABASE CENSORE JUNE 15.csv", dec = ".", na.string = "NA")
MvCov$Days_to_MACE_2014 <- MvCov$Time_B*30
MvCov$Status_MACE_tom_2014 <- as.numeric(MvCov$CVDEvents)
MvCov$hba1c_IFCC <- 10.93*MvCov$HbA1C-23.5
MvCov$smoke <- ifelse(MvCov$Smoking == 0, 0, 1)
MvCov$prev_cvd <- ifelse((MvCov$IHD == 1 | MvCov$Stroke == 1), 1, 0)
  # new ACR by gender
MvCov$MicroAlb <- ifelse(MvCov$Gender0female == 0, 
                       ifelse(is.na(MvCov$Alb_Creat_Ratio), 0, ifelse(MvCov$Alb_Creat_Ratio < 350, 0, 1)),
                       ifelse(is.na(MvCov$Alb_Creat_Ratio), 0, ifelse(MvCov$Alb_Creat_Ratio < 250, 0, 1)))
MvCov$diab_dur <- 0   # impute later
MvCov$afib <- 0       
MvCov$CKD_EPI_GFR_B <- CKDEpi.creat(creatinine = MvCov$Creatinine, sex = MvCov$Gender0female, 
                                    age = MvCov$Age, ethnicity = rep(0, nrow(MvCov)))
MvCov$CKD_EPI_GFR_B <- ifelse(is.na(MvCov$CKD_EPI_GFR_B), mean(na.omit(MvCov$CKD_EPI_GFR_B)), MvCov$CKD_EPI_GFR_B)
MvCov$Gender0female <- ifelse(MvCov$Gender0female == 0, "Kvinna", "Man")
MvCov$AHptMed <- 0
MvCov$AHptMed <- ifelse(MvCov$ACEARB == 1, 1, MvCov$AHptMed)
MvCov$AHptMed <- ifelse(MvCov$Calc_Block == 1, 1, MvCov$AHptMed)
MvCov$AHptMed <- ifelse(MvCov$Beta_Blok == 1, 1, MvCov$AHptMed)
MvCov$AHptMed <- ifelse(MvCov$Diuretics == 1, 1, MvCov$AHptMed)
MvCov$statin_a <- ifelse((MvCov$Statin_Cat == 1 | MvCov$Ezetemibe == 1), 1, 0)

CovM <- c( "N", "Status_MACE_tom_2014", "Days_to_MACE_2014", "Age", "Gender0female", "hba1c_IFCC", "BMI", "smoke", "MicroAlb",
           "prev_cvd", "SBP", "LDL", "HDL", "diab_dur", "afib", "CKD_EPI_GFR_B", "Insulin", "AHptMed", "statin_a")
MvCov <- MvCov[MvCov$DM == 1, CovM ]
names(MvCov) <- covt

  # proteins
MvPr <- read.xlsx("Njurkohort resultat final.xlsx", sheetName = "NPX (Page 2)", startRow = 2, 
                  colClasses = c("character", rep("numeric", 92), "character"))
MvPr <- MvPr[substr(MvPr$NPX, 1, 4) == "MIVC", ]
MvPr$loep <- as.numeric(substr(MvPr$NPX, 5, 8))
  # old protein had LOD written in by hand, now import
loadMv <- data.frame(fread("LOD_MvPr.csv"))
  # Added: investigate % missing proteins
temp <- c(rep(0, 92))
for(i in 2:93){
  temp[i-1] <- length(which(MvPr[,i] == "NaN" | is.na(MvPr[,i] | is.nan(MvPr[,i]))))
}
temp <- data.frame(cbind(names(MvPr)[2:93], temp / nrow(MvPr)))
temp2 <- temp[which(as.numeric(as.character(temp$X2)) > 0), ]
# pdf("missing_prot_MIVIC.pdf", 15, 5 )
barplot(as.numeric(as.character(temp2$X2)), names.arg = gsub("A_", "",temp2$X1),
        ylim = c(0,1), cex.names = 0.5, main = "missing protein in MIVIC n=285", cex =0.7, cex.main = 1)
abline(0.15, 0, col = "red")
# dev.off()

excl_mivic <- temp[which(as.numeric(as.character(temp$X2)) >= 0.15), ]$X1
  # 15% missing: X154_mAmP  22.1%,   X169_ITGB1BP2 93.333%,   X180_IL.4  97.9%,   X188_BNP, 24.9%

  ## easier scripting - replace LOD/2 first, than exclude 15%-missing
for(i in 2:93){
  MvPr[, i] <- ifelse(is.nan(MvPr[,i]), loadMv[1, i]/2, MvPr[, i])
  MvPr[, i] <- ifelse(is.na(MvPr[, i]), loadMv[1, i]/2, MvPr[, i])
} 
MvPr <- MvPr[, -which(names(MvPr) %in% c("X154_mAmP", "X169_ITGB1BP2", "X180_IL.4", "X188_BNP"))]
Mvdt <- MvPr[, -which(names(MvPr) %in% c("X.Chip.ID.", "NPX"))]
Mvdt <- merge(MvCov, Mvdt, by.x = "loep", by.y = "loep", all = F, all.x = T, all.y = F)
Mvdt <- na.omit(Mvdt)
prot_miv <- prot[-which(prot %in% c("A_154_mAmP", "A_169_ITGB1BP2", "A_180_IL4", "A_188_BNP"))]
names(Mvdt) <- c(covt, prot_miv)

Mvdt$Coh <- "MIVC"
Mvdt2 <- na.omit(Mvdt)

dtU2 <- dtU[, which(names(dtU) %in% names(dtPv))]
dtPv2 <- dtPv[, which(names(dtPv) %in% names(dtU2))]
  # check:
# rbind(names(dtU2), names(dtPv2))

dtAll <- rbind(dtU2, dtPv2)
    ### in old script: mg/ml not converted to mmol/l :
Mvdt2$ldl_kol <- Mvdt2$ldl_kol/38.67
Mvdt2$hdl_kol <- Mvdt2$hdl_kol/38.67
Mvdt2 <-  Mvdt[, which(names(Mvdt) %in% names(dtAll))]
  # check:
# rbind(names(dtAll), names(Mvdt2))

dtAll <- rbind(dtAll, Mvdt2)

save(Mvdt2, file = "MIVC_CN.Rdata")

################################################################################################
###  SAVa  #####################################################################################
################################################################################################

dtsv <- read.dta("SAVA CVD outcomes STATA12.dta")
dtsv2 <- read.dta("SAVA kontroller STATA12.dta")
dtsv <- dtsv[(dtsv$serial < 1000), ]
dtsv3 <- merge(x = dtsv, y = dtsv2, by.x = "serial", by.y = "serial", all.x = T)
dtsv3$Status_MACE_tom_2014 <- ifelse(dtsv3$faileventCV_PH_01 == "Yes", 1, 0)
dtsv3$koen_a <- ifelse(dtsv3$k_sex == "Man", "Man", "Kvinna")
dtsv3$Diabetes <- ifelse(dtsv3$k_diabetes_Q == "yes", 1, 0)
dtsv3 <- dtsv3[dtsv3$Diabetes == 1, ]
dtsv3$Days_to_MACE_2014 <- as.numeric(difftime(dtsv3$date_faileventCV_PH_01, dtsv3$k_dateInclusion, unit = "days")) 
dtsv3$smoke <- ifelse(dtsv3$k_smokehabit == "Current smoker", 1, 0)
dtsv3$MicroAlb <- 0 # missing whole cohort
dtsv3$prev_cvd <- 0
dtsv3$diab_dur <- 0
dtsv3$afib <- ifelse((dtsv3$k_rhythmECG == "Atrial fibrillation" | dtsv3$k_rhythmECG == "Atrial flutter"), 1, 0)
dtsv3$afib <- ifelse(is.na(dtsv3$afib), 0, 1)
dtsv3$CKD_EPI_GFR_B <- CKDEpi.creat(creatinine = (dtsv3$k_creatinine / 76.26), sex = as.numeric(dtsv3$k_sex == "Man"),
                                  age = dtsv3$k_age_calculated, ethnicity = rep(0, nrow(dtsv3)))
dtsv3$CKD_EPI_GFR_B <- ifelse(is.na(dtsv3$CKD_EPI_GFR_B), mean(na.omit(dtsv3$CKD_EPI_GFR_B)), dtsv3$CKD_EPI_GFR_B)
dtsv3$ADMed <- ifelse(dtsv3$k_antidiabetic == "yes", 1, 0)
dtsv3$AHptMed <- 0
dtsv3$AHptMed <- ifelse(dtsv3$k_ace == "yes", 1, dtsv3$AHptMed)
dtsv3$AHptMed <- ifelse(dtsv3$k_arbi =="yes", 1, dtsv3$AHptMed)
dtsv3$AHptMed <- ifelse(dtsv3$k_caInhibitor == "yes", 1, dtsv3$AHptMed)
dtsv3$AHptMed <- ifelse(dtsv3$k_betablocker == "yes", 1, dtsv3$AHptMed)
dtsv3$AHptMed <- ifelse(dtsv3$k_diuretic == "yes", 1, dtsv3$AHptMed)
dtsv3$statin_a <- ifelse(dtsv3$k_statin == "yes", 1, 0)

  # LOD replacement already done in file we have - not great, how many missing not clear
  # but OK, given the small size of SAVa in relation to the whole sample
CovS <- c("serial", "Status_MACE_tom_2014", "Days_to_MACE_2014", "k_age_calculated", "koen_a", "k_HbA1c_IFCC", "k_bmi1", 
          "smoke", "MicroAlb", "prev_cvd", "k_sbpRightMean", "k_ldlCholesterol", "k_hdlCholesterol", "diab_dur", "afib", "CKD_EPI_GFR_B",
          "ADMed", "AHptMed", "statin_a")
dtsv4 <- dtsv3[, CovS]
dtsv4 <- cbind(dtsv4, dtsv3[, 599:690]) 
names(dtsv4) <- c(covt, prot)
dtsv4$Coh <- "Sv_C"
dtsv4_2 <- dtsv4[, which(names(dtsv4) %in% names(dtAll))]
dtsv4_2 <- na.omit(dtsv4_2)

save(dtsv4_2, file = "SAVa_CN.Rdata")

dtAll <- rbind(dtAll, dtsv4_2)

save(dtAll, file = "MergedReplicationSamples.Rdata")

# write.csv(dtAll, "dtAll_merged_repl_samples.csv")

############################################################################################################################
############################################################################################################################
#### REPLICATION ANALYSIS ############################################################################################# ####
############################################################################################################################
############################################################################################################################

load("MergedReplicationSamples.Rdata")
disc <- read.csv("CARDIPP_discovery_result.csv")
disc_pr <- disc$X[which(disc$PASS != "fail")]

  # There are proteins that passed QC in discovery, but failed in the replication samples
disc_pr[-which(disc_pr %in% names(dtAll))]
  # exclude them - 35 out of 39 FDR-"discovered" protein can be tested for replication

REPL <- dtAll # cohort dummies, ULSAM = 0
REPL$Coh_1 <- ifelse(REPL$Coh == "PIV", 1, 0)
REPL$Coh_2 <- ifelse(REPL$Coh == "MIVC", 1, 0)
REPL$Coh_3 <- ifelse(REPL$Coh == "Sv_C", 1, 0)
REPL$koen_a <- ifelse(REPL$koen_a == "Kvinna", 0, 1)

 # Null model: Age + Sex + cohort
Outcome <- Surv(time = REPL$Days_to_MACE_2014, event = REPL$Status_MACE_tom_2014) 
Cox_result <- coxph(Outcome ~ koen_a + Aalder + Coh_1 + Coh_2 + Coh_3, data = REPL) 

summary(Cox_result)   
cox.zph(Cox_result)         # PH assumptions not violated p = 0.779
plot(cox.zph(Cox_result))

Kaplan = survfit(Outcome ~ koen_a , data = REPL)
Kaplan = survfit(Outcome ~ Coh_1 + Coh_2 + Coh_3 , data = REPL)
plot(Kaplan)

  # select discovered FDR significant proteins
index <- which(names(REPL) %in% disc_pr)
# disc_pr[which(disc_pr %in% names(REPL)[index] == F)]
replication_age_sex <- data.frame(matrix(nrow = length(index), ncol = 6,
                                       dimnames = list(names(REPL)[index], c("b", "se", "p", "HR", "LCI", "UCI"))))

for(i in index){
  REPL[,i] <- scale(as.numeric(REPL[,i]))
  REPL[,i] <- (as.numeric(REPL[,i]))
}

j <- 1
for (i in index){
  replication_age_sex$b[j] <- summary(coxph(Outcome ~ REPL[, i] + koen_a + Aalder + Coh_1 + Coh_2 + Coh_3, data = REPL))$coeff[1,1]
  replication_age_sex$se[j] <- summary(coxph(Outcome ~ REPL[, i] + koen_a + Aalder + Coh_1 + Coh_2 + Coh_3, data = REPL))$coeff[1,3]
  replication_age_sex$p[j] <- summary(coxph(Outcome ~ REPL[, i] + koen_a + Aalder + Coh_1 + Coh_2 + Coh_3, data = REPL))$coeff[1,5]
  replication_age_sex$HR[j] <- exp(replication_age_sex$b[j])
  replication_age_sex$LCI[j] <- exp(replication_age_sex$b[j] - 1.96*replication_age_sex$se[j])
  replication_age_sex$UCI[j] <- exp(replication_age_sex$b[j] + 1.96*replication_age_sex$se[j])
  j <- j + 1
}

replication_age_sex <- replication_age_sex[order(replication_age_sex$p, decreasing = F),]
replication_age_sex$p_fdr <- p.adjust(replication_age_sex$p, method = "fdr")
replication_age_sex$PASS <- "fail"
replication_age_sex$PASS[which(replication_age_sex$p_fdr < 0.05 & replication_age_sex$b > 0)] <- "increase"
replication_age_sex$PASS[which(replication_age_sex$p_fdr < 0.05 & replication_age_sex$b < 0)] <- "decrease"

# View(replication_age_sex)

# write.csv(replication_age_sex, "replication_age_sex")

#### #### #### #### #### #### #### #### #### #### #### #### 
#### prepare all-cohort-merger for prediction script # #### 
#### #### #### #### #### #### #### #### #### #### #### #### 

excl <- which(names(CARDIPP) %in% names(dtAll) == F)
CARD <- CARDIPP[, -excl]
COMBINE <- rbind(CARD, dtAll)

#  write.csv(COMBINE, "COMBINE.csv")

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### Sensitivity: Covariate adjustment #### #### #### #### #### #### #### #####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

## (if running on all cohorts combined, upload the GBM-file that has missing TC and HbA1c in):
# load("COMBINE_wholesome.Rdata") # MicroAlb.x
  
COMBINE <- REPL
# COMBINE <- REPL[-which(REPL$prev_cvd == 1), ] # to exclude all 90 with previous CVD

##  Directed Acyclic Graph: exposure is protein, outcome is mace, covariates are all available ones
##  Minimal sufficient adjustment set for estimating the total effect of Protein on MACE:
##  AF, Age, BMI, Microalbuminuria, SBP, Sex, Smoking, T2D_duration, previous_CV

Outcome <- Surv(time = COMBINE$Days_to_MACE_2014, event = COMBINE$Status_MACE_tom_2014) 
Cox_result <- coxph(Outcome ~ koen_a  + Aalder +  flimmer_a + BMI_VC_B  + MicroAlb + 
                      bt_syst_mean_a_B + smoke  + diab_dur  + prev_cvd +
                      Coh_1 + Coh_2 + Coh_3, data = COMBINE) 
summary(Cox_result)   
cox.zph(Cox_result)         # PH assumptions not violated p = 0.942
plot(cox.zph(Cox_result))

  # for FDR-replicated protein from age/sex-model
index <- which(names(COMBINE) %in% disc_pr)
ADJUSTED <- data.frame(matrix(nrow = length(index), ncol = 6,
                                         dimnames = list(names(COMBINE)[index], c("b", "se", "p", "HR", "LCI", "UCI"))))

for(i in index){
  COMBINE[,i] <- scale(as.numeric(COMBINE[,i]))
  COMBINE[,i] <- (as.numeric(COMBINE[,i]))
}

j <- 1
for (i in index){
  ADJUSTED$b[j] <- summary(coxph(Outcome ~ COMBINE[, i] + koen_a  + Aalder +  flimmer_a + BMI_VC_B  + MicroAlb +
                                   bt_syst_mean_a_B + smoke  + diab_dur  + prev_cvd +
                                   Coh_1 + Coh_2 + Coh_3, data = COMBINE))$coeff[1,1]
  ADJUSTED$se[j] <- summary(coxph(Outcome ~ COMBINE[, i] + koen_a  +  Aalder + flimmer_a + BMI_VC_B  + MicroAlb +
                                    bt_syst_mean_a_B + smoke  + diab_dur  + prev_cvd +
                                    Coh_1 + Coh_2 + Coh_3, data = COMBINE))$coeff[1,3]
  ADJUSTED$p[j] <- summary(coxph(Outcome ~ COMBINE[, i] + koen_a  + Aalder + flimmer_a + BMI_VC_B  + MicroAlb +
                                   bt_syst_mean_a_B + smoke  + diab_dur  + prev_cvd +
                                   Coh_1 + Coh_2 + Coh_3, data = COMBINE))$coeff[1,5]
  ADJUSTED$HR[j] <- exp(ADJUSTED$b[j])
  ADJUSTED$LCI[j] <- exp(ADJUSTED$b[j] - 1.96*ADJUSTED$se[j])
  ADJUSTED$UCI[j] <- exp(ADJUSTED$b[j] + 1.96*ADJUSTED$se[j])
  j <- j + 1
}

ADJUSTED <- ADJUSTED[order(ADJUSTED$p, decreasing = F),]
ADJUSTED$p_fdr <- p.adjust(ADJUSTED$p, method = "fdr")
ADJUSTED$PASS <- "fail"
ADJUSTED$PASS[which(ADJUSTED$p_fdr < 0.05 & ADJUSTED$b > 0)] <- "increase"
ADJUSTED$PASS[which(ADJUSTED$p_fdr < 0.05 & ADJUSTED$b < 0)] <- "decrease"

# View(ADJUSTED)

# write.csv(ADJUSTED, "replicated_DAG_adjusted.csv")

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

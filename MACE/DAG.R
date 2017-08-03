################################################################################
##### Directed Acyclid Graph  ##################################################
################################################################################

## For browser-based:

dag {
  "LDL-C" [pos="0.296,0.447"]
  AF [pos="0.653,-0.006"]
  Age [pos="-0.013,0.458"]
  Antidiabetic_drugs [pos="0.609,1.003"]
  Antihypertensive_drugs [pos="0.475,0.019"]
  BMI [pos="0.022,0.083"]
  HbA1c [pos="0.232,0.255"]
  MACE [outcome,pos="0.897,0.596"]
  Microalbuminuria [pos="0.424,0.772"]
  Protein [exposure,pos="0.643,0.774"]
  SBP [pos="0.173,-0.052"]
  Sex [pos="-0.019,0.259"]
  Smoking [pos="0.100,0.019"]
  Statin [pos="0.145,0.518"]
  T2D_duration [pos="0.419,1.001"]
  Total_cholesterol [pos="0.156,0.736"]
  eGFR [pos="0.292,0.591"]
  previous_CVD [pos="0.471,0.462"]
  "LDL-C" -> HbA1c
  "LDL-C" -> MACE
  AF -> MACE
  AF -> Protein
  AF -> previous_CVD
  Age -> HbA1c
  Age -> Microalbuminuria [pos="-0.006,1.030"]
  Age -> Protein
  Age -> Statin
  Age -> T2D_duration [pos="0.011,1.032"]
  Age -> eGFR [pos="0.190,0.362"]
  Antidiabetic_drugs -> MACE [pos="0.786,0.894"]
  Antihypertensive_drugs -> MACE
  Antihypertensive_drugs -> previous_CVD
  BMI -> "LDL-C"
  BMI -> HbA1c
  BMI -> Protein
  BMI -> T2D_duration [pos="-0.023,0.950"]
  BMI -> Total_cholesterol [pos="0.032,0.498"]
  BMI -> eGFR [pos="0.068,0.271"]
  HbA1c -> MACE [pos="0.351,0.233"]
  HbA1c -> previous_CVD
  Microalbuminuria -> MACE
  Microalbuminuria -> Protein
  Protein -> MACE
  SBP -> Antihypertensive_drugs
  SBP -> MACE
  SBP -> Protein [pos="0.465,0.177"]
  SBP -> eGFR [pos="0.071,0.206"]
  SBP -> previous_CVD
  Sex -> HbA1c
  Sex -> MACE [pos="0.256,-0.052"]
  Sex -> Microalbuminuria [pos="0.463,0.344"]
  Sex -> Protein
  Sex -> eGFR
  Sex -> previous_CVD
  Smoking -> MACE
  Smoking -> Protein [pos="0.514,0.173"]
  Smoking -> SBP
  Smoking -> eGFR [pos="0.084,0.299"]
  Smoking -> previous_CVD
  Statin -> "LDL-C"
  Statin -> T2D_duration
  Statin -> Total_cholesterol
  T2D_duration -> Antidiabetic_drugs
  T2D_duration -> MACE [pos="0.774,0.832"]
  T2D_duration -> Protein
  Total_cholesterol -> MACE
  eGFR -> MACE
  eGFR -> Microalbuminuria
  eGFR -> previous_CVD [pos="0.383,0.462"]
  previous_CVD -> MACE
  previous_CVD -> Protein
}

## for R:

testImplications <- function( covariance.matrix, sample.size ){
  library(ggm)
  tst <- function(i){ pcor.test( pcor(i,covariance.matrix), length(i)-2, sample.size )$pvalue }
  tos <- function(i){ paste(i,collapse=" ") }
  implications <- list(c("MACE","Age","BMI","T2D_duration","Sex","Microalbuminuria","Smoking","SBP","LDL-C","HbA1c","Protein","Statin","eGFR","previous_CVD","AF"),
                       c("MACE","Age","T2D_duration","Total_cholesterol","LDL-C","SBP","eGFR","HbA1c","Protein","Sex","Smoking","Microalbuminuria","previous_CVD","AF"),
                       c("MACE","BMI","T2D_duration","Sex","Microalbuminuria","Smoking","SBP","LDL-C","HbA1c","Protein","eGFR","previous_CVD","AF","Total_cholesterol"),
                       c("MACE","Statin","BMI","T2D_duration","Total_cholesterol","LDL-C","Age"),
                       c("MACE","Statin","T2D_duration","Sex","Microalbuminuria","Smoking","SBP","LDL-C","HbA1c","Protein","eGFR","previous_CVD","AF","Total_cholesterol"),
                       c("Protein","HbA1c","AF","SBP","Smoking","Sex","eGFR","previous_CVD","BMI","Age","LDL-C"),
                       c("Protein","HbA1c","AF","SBP","Smoking","Sex","eGFR","previous_CVD","BMI","Age","Statin"),
                       c("Protein","HbA1c","BMI","T2D_duration","Age","AF","SBP","Smoking","Sex","eGFR","previous_CVD"),
                       c("Protein","HbA1c","AF","SBP","Smoking","Sex","previous_CVD","Microalbuminuria","Age","BMI","LDL-C"),
                       c("Protein","HbA1c","AF","SBP","Smoking","Sex","previous_CVD","BMI","Age","Statin","Microalbuminuria"),
                       c("Protein","HbA1c","BMI","T2D_duration","Age","AF","SBP","Smoking","Sex","previous_CVD","Microalbuminuria"),
                       c("Protein","LDL-C","HbA1c","Statin","BMI","Sex","Age"),
                       c("Protein","LDL-C","BMI","T2D_duration","Age","HbA1c","Sex"),
                       c("Protein","LDL-C","AF","SBP","Smoking","Sex","eGFR","previous_CVD","BMI","Age","Statin"),
                       c("Protein","LDL-C","AF","SBP","Smoking","Sex","eGFR","previous_CVD","BMI","T2D_duration","Age"),
                       c("Protein","LDL-C","AF","SBP","Smoking","Sex","previous_CVD","Microalbuminuria","Age","BMI","Statin"),
                       c("Protein","LDL-C","AF","SBP","Smoking","Sex","previous_CVD","Microalbuminuria","Age","BMI","T2D_duration"),
                       c("Protein","Total_cholesterol","Statin","BMI"),
                       c("Protein","Total_cholesterol","BMI","T2D_duration","LDL-C","Age"),
                       c("Protein","Total_cholesterol","HbA1c","BMI","Sex","Age","T2D_duration"),
                       c("Protein","Total_cholesterol","AF","SBP","Smoking","Sex","eGFR","previous_CVD","BMI","Age","T2D_duration"),
                       c("Protein","Total_cholesterol","AF","SBP","Smoking","Sex","previous_CVD","Microalbuminuria","Age","BMI","T2D_duration"),
                       c("Protein","Antihypertensive_drugs","AF","Smoking","Sex","HbA1c","eGFR","previous_CVD","SBP"),
                       c("Protein","Antihypertensive_drugs","eGFR","previous_CVD","AF","SBP","Smoking","Sex","BMI","Age","LDL-C"),
                       c("Protein","Antihypertensive_drugs","BMI","Sex","Age","Statin","eGFR","previous_CVD","AF","SBP","Smoking"),
                       c("Protein","Antihypertensive_drugs","Age","BMI","T2D_duration","Sex","eGFR","previous_CVD","AF","SBP","Smoking"),
                       c("Protein","Antihypertensive_drugs","previous_CVD","Microalbuminuria","AF","SBP","Smoking","Sex","HbA1c","Age","BMI"),
                       c("Protein","Antihypertensive_drugs","previous_CVD","Microalbuminuria","AF","SBP","Smoking","Sex","Age","BMI","LDL-C"),
                       c("Protein","Antihypertensive_drugs","BMI","Sex","Age","Statin","previous_CVD","Microalbuminuria","AF","SBP","Smoking"),
                       c("Protein","Antihypertensive_drugs","Age","BMI","T2D_duration","Sex","previous_CVD","Microalbuminuria","AF","SBP","Smoking"),
                       c("Protein","Antidiabetic_drugs","T2D_duration"),
                       c("Protein","Statin","BMI","T2D_duration","LDL-C","Age"),
                       c("Protein","Statin","HbA1c","BMI","Sex","Age","T2D_duration"),
                       c("Protein","Statin","AF","SBP","Smoking","Sex","eGFR","previous_CVD","BMI","Age","T2D_duration"),
                       c("Protein","Statin","AF","SBP","Smoking","Sex","previous_CVD","Microalbuminuria","Age","BMI","T2D_duration"),
                       c("Protein","eGFR","AF","SBP","Smoking","Sex","HbA1c","previous_CVD","Microalbuminuria","Age","BMI"),
                       c("Protein","eGFR","previous_CVD","AF","SBP","Smoking","Sex","BMI","Age","LDL-C","Microalbuminuria"),
                       c("Protein","eGFR","BMI","Sex","Age","Statin","previous_CVD","AF","SBP","Smoking","Microalbuminuria"),
                       c("Protein","eGFR","Age","BMI","T2D_duration","Sex","previous_CVD","AF","SBP","Smoking","Microalbuminuria"),
                       c("Sex","Age"),
                       c("Sex","BMI"),
                       c("Sex","Smoking"),
                       c("Sex","SBP"),
                       c("Sex","T2D_duration"),
                       c("Sex","LDL-C"),
                       c("Sex","Total_cholesterol"),
                       c("Sex","Antihypertensive_drugs"),
                       c("Sex","Antidiabetic_drugs"),
                       c("Sex","Statin"),
                       c("Sex","AF"),
                       c("Age","BMI"),
                       c("Age","Smoking"),
                       c("Age","previous_CVD","SBP","Smoking","Sex","HbA1c","eGFR"),
                       c("Age","SBP"),
                       c("Age","LDL-C","Statin"),
                       c("Age","Total_cholesterol","Statin"),
                       c("Age","Antihypertensive_drugs"),
                       c("Age","Antidiabetic_drugs","T2D_duration"),
                       c("Age","AF"),
                       c("HbA1c","Smoking"),
                       c("HbA1c","Microalbuminuria","eGFR","Age","Sex"),
                       c("HbA1c","Microalbuminuria","Sex","Age","BMI"),
                       c("HbA1c","SBP"),
                       c("HbA1c","T2D_duration","Statin","Age","BMI"),
                       c("HbA1c","T2D_duration","Age","BMI","LDL-C"),
                       c("HbA1c","Total_cholesterol","Statin","BMI"),
                       c("HbA1c","Total_cholesterol","BMI","LDL-C","Age"),
                       c("HbA1c","Antihypertensive_drugs"),
                       c("HbA1c","Antidiabetic_drugs","T2D_duration"),
                       c("HbA1c","Antidiabetic_drugs","Statin","Age","BMI"),
                       c("HbA1c","Antidiabetic_drugs","Age","BMI","LDL-C"),
                       c("HbA1c","Statin","BMI","LDL-C","Age"),
                       c("HbA1c","eGFR","Sex","Age","BMI"),
                       c("HbA1c","AF"),
                       c("BMI","Smoking"),
                       c("BMI","Microalbuminuria","eGFR","Age","Sex"),
                       c("BMI","previous_CVD","SBP","Smoking","Sex","HbA1c","eGFR"),
                       c("BMI","SBP"),
                       c("BMI","Antihypertensive_drugs"),
                       c("BMI","Antidiabetic_drugs","T2D_duration"),
                       c("BMI","Statin"),
                       c("BMI","AF"),
                       c("Smoking","Microalbuminuria","eGFR","Age","Sex"),
                       c("Smoking","T2D_duration"),
                       c("Smoking","LDL-C"),
                       c("Smoking","Total_cholesterol"),
                       c("Smoking","Antihypertensive_drugs","SBP"),
                       c("Smoking","Antidiabetic_drugs"),
                       c("Smoking","Statin"),
                       c("Smoking","AF"),
                       c("Microalbuminuria","previous_CVD","SBP","Smoking","Sex","HbA1c","eGFR"),
                       c("Microalbuminuria","previous_CVD","Age","Sex","eGFR"),
                       c("Microalbuminuria","SBP","eGFR","Sex","Age"),
                       c("Microalbuminuria","T2D_duration","Age","BMI"),
                       c("Microalbuminuria","T2D_duration","eGFR","Sex","Age"),
                       c("Microalbuminuria","LDL-C","Statin","BMI"),
                       c("Microalbuminuria","LDL-C","BMI","Age"),
                       c("Microalbuminuria","LDL-C","eGFR","Sex","Age"),
                       c("Microalbuminuria","Total_cholesterol","Statin","BMI"),
                       c("Microalbuminuria","Total_cholesterol","BMI","Age"),
                       c("Microalbuminuria","Total_cholesterol","eGFR","Sex","Age"),
                       c("Microalbuminuria","Antihypertensive_drugs","SBP"),
                       c("Microalbuminuria","Antihypertensive_drugs","eGFR","Sex","Age"),
                       c("Microalbuminuria","Antidiabetic_drugs","T2D_duration"),
                       c("Microalbuminuria","Antidiabetic_drugs","Age","BMI"),
                       c("Microalbuminuria","Antidiabetic_drugs","eGFR","Sex","Age"),
                       c("Microalbuminuria","Statin","Age"),
                       c("Microalbuminuria","AF"),
                       c("previous_CVD","T2D_duration","Statin","Age","BMI"),
                       c("previous_CVD","T2D_duration","Age","BMI","LDL-C"),
                       c("previous_CVD","T2D_duration","HbA1c","BMI","Sex","Age"),
                       c("previous_CVD","T2D_duration","HbA1c","Sex","SBP","eGFR","Smoking"),
                       c("previous_CVD","LDL-C","BMI","Age","HbA1c","Sex"),
                       c("previous_CVD","LDL-C","Smoking","SBP","eGFR","HbA1c","Sex"),
                       c("previous_CVD","Total_cholesterol","Statin","BMI"),
                       c("previous_CVD","Total_cholesterol","BMI","LDL-C","Age"),
                       c("previous_CVD","Total_cholesterol","HbA1c","BMI","Sex","Age"),
                       c("previous_CVD","Total_cholesterol","Smoking","SBP","eGFR","HbA1c","Sex"),
                       c("previous_CVD","Antidiabetic_drugs","T2D_duration"),
                       c("previous_CVD","Antidiabetic_drugs","Statin","Age","BMI"),
                       c("previous_CVD","Antidiabetic_drugs","Age","BMI","LDL-C"),
                       c("previous_CVD","Antidiabetic_drugs","HbA1c","BMI","Sex","Age"),
                       c("previous_CVD","Antidiabetic_drugs","HbA1c","Sex","SBP","eGFR","Smoking"),
                       c("previous_CVD","Statin","BMI","LDL-C","Age"),
                       c("previous_CVD","Statin","HbA1c","BMI","Sex","Age"),
                       c("previous_CVD","Statin","Smoking","SBP","eGFR","HbA1c","Sex"),
                       c("SBP","T2D_duration"),
                       c("SBP","LDL-C"),
                       c("SBP","Total_cholesterol"),
                       c("SBP","Antidiabetic_drugs"),
                       c("SBP","Statin"),
                       c("SBP","AF"),
                       c("T2D_duration","LDL-C","Statin","BMI"),
                       c("T2D_duration","Total_cholesterol","Statin","BMI"),
                       c("T2D_duration","Antihypertensive_drugs"),
                       c("T2D_duration","eGFR","Age","BMI"),
                       c("T2D_duration","AF"),
                       c("LDL-C","Total_cholesterol","Statin","BMI"),
                       c("LDL-C","Antihypertensive_drugs"),
                       c("LDL-C","Antidiabetic_drugs","T2D_duration"),
                       c("LDL-C","Antidiabetic_drugs","BMI","Statin"),
                       c("LDL-C","eGFR","Age","BMI"),
                       c("LDL-C","eGFR","BMI","Statin"),
                       c("LDL-C","AF"),
                       c("Total_cholesterol","Antihypertensive_drugs"),
                       c("Total_cholesterol","Antidiabetic_drugs","T2D_duration"),
                       c("Total_cholesterol","Antidiabetic_drugs","BMI","Statin"),
                       c("Total_cholesterol","eGFR","Age","BMI"),
                       c("Total_cholesterol","eGFR","BMI","Statin"),
                       c("Total_cholesterol","AF"),
                       c("Antihypertensive_drugs","Antidiabetic_drugs"),
                       c("Antihypertensive_drugs","Statin"),
                       c("Antihypertensive_drugs","eGFR","SBP"),
                       c("Antihypertensive_drugs","AF"),
                       c("Antidiabetic_drugs","Statin","T2D_duration"),
                       c("Antidiabetic_drugs","eGFR","Age","BMI"),
                       c("Antidiabetic_drugs","eGFR","T2D_duration"),
                       c("Antidiabetic_drugs","AF"),
                       c("Statin","eGFR","Age"),
                       c("Statin","AF"),
                       c("eGFR","AF"))
  data.frame( implication=unlist(lapply(implications,tos)),
              pvalue=unlist( lapply( implications, tst ) ) )
  
}
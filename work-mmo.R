library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(edgeR)
library(tidyverse)


BRCA_patient_clinical = readRDS("~/Desktop/TCGA-BRCA-metadata/TCGA-BRCA-metadata/BRCA_patient_clinical.RDS")
BRCA_subtype = readRDS("~/Desktop/TCGA-BRCA-metadata/TCGA-BRCA-metadata/BRCA_subtype.RDS")
BRCA_TranscriptomeProfiling = readRDS("~/Desktop/TCGA-BRCA-metadata/TCGA-BRCA-metadata/BRCA_TranscriptomeProfiling.RDS")
BRCA_miRNASeqraw = readRDS("~/Desktop/TCGA-BRCA-metadata/TCGA-BRCA-metadata/BRCA_miRNASeqraw.RDS")
BRCA_RNASeqraw = readRDS("~/Desktop/TCGA-BRCA-metadata/TCGA-BRCA-metadata/BRCA_RNASeqraw.RDS")


colnames(BRCA_patient_clinical)
colnames(BRCA_subtype)

table(BRCA_subtype$vital_status, useNA = "always")
table(BRCA_patient_clinical$vital_status, useNA = "always")

p = BRCA_patient_clinical
s = BRCA_subtype

colnames(BRCA_patient_clinical)
colnames(BRCA_subtype)
colnames(BRCA_patient_clinical[,62:85])


table(p$paper_Tumor_Grade, useNA = "always")
table(p$paper_CNV.Clusters, useNA = "always")
table(p$paper_miRNA.Clusters, useNA = "always")




na.analyze = function(p) {
  NA_ascharacter = c()
  NA_asmissing = c()
  NA_analysis = as.data.frame(colnames(p))
  for (c in 1:nrow(NA_analysis)) {
    missing_counter = 0
    character_counter = 0 
    for (e in 1:nrow(p[c])) {
      missing_stat = is.na(p[e,c])
      missing_aschar_stat= p[e,c]
      if (missing_stat == TRUE) {
        next
      } else if (missing_aschar_stat == "NA") {
        character_counter = character_counter + 1
      }
    }
    NA_ascharacter = append(NA_ascharacter,character_counter)
    missing_counter = sum(is.na(p[c]))
    NA_asmissing = append(NA_asmissing,missing_counter)
  }
  
  NA_analysis = NA_analysis %>%  mutate(NA_ascharacter_count = NA_ascharacter) %>% mutate(NA_asmissing_count = NA_asmissing)
  return(NA_analysis)
}




na.analyzed_patient = na.analyze(patient)
na.analyzed_subtype = na.analyze(subtype)


na.analyzed_patient = na.analyzed_patient %>% mutate(total_NA = NA_ascharacter_count + NA_asmissing_count) %>% 
  mutate(Available_from1226 = 1226- total_NA)
  
na.analyzed_subtype = na.analyzed_subtype %>% mutate(total_NA = NA_ascharacter_count + NA_asmissing_count) %>% 
  mutate(Available_from1087 = 1087- total_NA)

withinpatient = na.analyzed_patient[62:85,]
common_comparison_table = cbind(withinpatient,na.analyzed_subtype)
View(common_comparison_table)








sum(duplicated(patient$paper_patient, fromLast = TRUE) | duplicated(patient$paper_patient)) #149
sum(duplicated(patient$paper_patient)) #142
sum(duplicated(patient$patient)) #131

sum(duplicated(subtype$patient, fromLast = TRUE) | duplicated(subtype$patient)) #0
sum(duplicated(subtype$patient)) #0
View(na.analyzed_patient[62:85,])

# CD8A,CD3E,GZMB,PRF1 immune signature; survival 
install.packages("survminer")
install.packages("survival")

library(survminer)
library(survival)
brca_surv = brca_surv[1:1548,]
brca_surv = select(BRCA_metadata,c(patient,days_to_last_follow_up,days_to_death,vital_status,
                                   immune_signature))
brca_surv = brca_surv[1:1548,]
median_immune_signature = median(brca_surv$immune_signature, na.rm = TRUE)
brca_surv = brca_surv %>% 
  mutate(immune_category = case_when(
    immune_signature > median_immune_signature ~ "high" , 
    immune_signature < median_immune_signature ~ "low"
  )) %>% 
  mutate(days_to_event = case_when(
    is.na(days_to_death) == FALSE ~ days_to_death,
    TRUE ~ days_to_last_follow_up
  )) %>% 
  mutate(status = case_when(
    vital_status == "Alive" ~ 1,
    vital_status =="Dead" ~ 0
  ))


fit <- surv_fit(Surv(days_to_event, status) ~ immune_category, data = brca_surv)
ggsurvplot(fit,
           data=brca_surv,
           risk.table = TRUE,
           pval = TRUE,
           conf.int = TRUE,
           xlab = "Time in days",
           ggtheme = theme_bw() )



























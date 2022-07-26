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
#
colnames(BRCA_patient_clinical)
colnames(BRCA_subtype)
colnames(BRCA_patient_clinical[,62:85])

#To check whether there is any NA values that do not take place as a character in the BRCA subtype data. 
for (i in 1:24) {
  if (sum(is.na(p[i])) == 0 ) {
    print(i)
  }
}




counter = 0 
for (x in 1:nrow(s[6])) {
  if( s[x,6] == "NA") {
    counter = counter + 1
  } 
}
print(counter)











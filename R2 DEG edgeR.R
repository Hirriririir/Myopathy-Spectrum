library(vroom)
library(tidyverse)
library(readxl)
library(edgeR)

merge_L <- vroom("Integration data/All_muscle_L_combatseq_tmm.csv")
Meta_all <- read_excel("Meta/adata_adj.obs.xlsx")
merge_data <- merge_L[,Meta_all$Sample_id]


# edgeR #-----
unique(Meta_all$Phenotype_4)
mygroups <- Meta_all$Phenotype_4

y <- DGEList(counts=merge_data, genes=merge_L$Symbol, group = mygroups) #
y <- calcNormFactors(y, method ='TMM')
y <- estimateDisp(y)

# export edgeR tmm count to scanoy (layer = 'edgeR') #----

tmmcount <- cpm(y)
rownames(tmmcount) <- merge_L$Symbol
tmmcount
write.csv(tmmcount, "Integration data/All_muscle_L_combatseq_tmm_edgeR.csv")


## DEG common pathway Phenotype_5#----
unique(Meta_all$Phenotype_5)

et_Myopathy_Healthy <- exactTest(y, pair=c('Healthy', 'Myopathy') ) 
write.csv(topTags(et_Myopathy_Healthy, n=10000), "DEG/EdgeR Myopathy vs. Healthy (fast).csv",row.names = FALSE)

et_Myopathy_Wasting <- exactTest(y, pair=c('Wasting', 'Myopathy') ) 
write.csv(topTags(et_Myopathy_Wasting, n=10000), "DEG/EdgeR Myopathy vs. Wasting (slow).csv",row.names = FALSE)


## DEG unique pathway Phenotype_4#----
unique(Meta_all$Phenotype_2)
mygroups <- Meta_all$Phenotype_2

y <- DGEList(counts=merge_data, genes=merge_L$Symbol, group = mygroups) #
y <- calcNormFactors(y, method ='TMM')
y <- estimateDisp(y)

# CDM
et_CDM_fast <- exactTest(y, pair=c('Control (fast death)', 'CDM') ) 
write.csv(topTags(et_CDM_fast, n=10000), "DEG/EdgeR CDM vs. fast.csv",row.names = FALSE)

# DM1
et_DM1_fast <- exactTest(y, pair=c('Control (fast death)', 'DM1') ) 
write.csv(topTags(et_DM1_fast, n=10000), "DEG/EdgeR DM1 vs. fast.csv",row.names = FALSE)

# FSHD
et_FSHD_fast <- exactTest(y, pair=c('Control (fast death)', 'FSHD') ) 
write.csv(topTags(et_FSHD_fast, n=10000), "DEG/EdgeR FSHD vs. fast.csv",row.names = FALSE)

# IBM
et_IBM_fast <- exactTest(y, pair=c('Control (fast death)', 'IBM') ) 
write.csv(topTags(et_IBM_fast, n=10000), "DEG/EdgeR IBM vs. fast.csv",row.names = FALSE)

# LGMD R12
et_LGMDR12_fast <- exactTest(y, pair=c('Control (fast death)', 'LGMD R12') ) 
write.csv(topTags(et_LGMDR12_fast, n=10000), "DEG/EdgeR LGMDR12 vs. fast.csv",row.names = FALSE)

# Titinopathy
et_Titinopathy_fast <- exactTest(y, pair=c('Control (fast death)', 'Titinopathy') ) 
write.csv(topTags(et_Titinopathy_fast, n=10000), "DEG/EdgeR Titinopathy vs. fast.csv",row.names = FALSE)

# Control (slow death)
et_slow_fast <- exactTest(y, pair=c('Control (fast death)', 'Control (slow death)') ) 
write.csv(topTags(et_slow_fast, n=10000), "DEG/EdgeR slow vs. fast.csv",row.names = FALSE)

#  FSHD Vs. Control (accident death)
et_FSHD_accident <- exactTest(y, pair=c('Control (accident death)', 'FSHD') ) 
write.csv(topTags(et_FSHD_accident, n=10000), "DEG/EdgeR FSHD vs. accident death.csv",row.names = FALSE)


# DEG batch validation #----
unique(Meta_all$Phenotype_6)
mygroups <- Meta_all$Phenotype_6

y <- DGEList(counts=merge_data, genes=merge_L$Symbol, group = mygroups) #
y <- calcNormFactors(y, method ='TMM')
y <- estimateDisp(y)



#  IBM Vs. Control (amputee)
et_IBM_Amputee <- exactTest(y, pair=c('Control (amputee)', 'IBM') ) 
write.csv(topTags(et_IBM_Amputee, n=10000), "Validation/[Batch-combatseq] EdgeR IBM vs. Amputee.csv",row.names = FALSE)

# Helsinki: Titinopathy vs Control (amputee)
et_Titinopathy_Amputee <- exactTest(y, pair=c('Control (amputee)', 'Titinopathy') ) 
write.csv(topTags(et_Titinopathy_Amputee, n=10000), "Validation/[Batch-combatseq] EdgeR Titinopathy vs. Amputee.csv",row.names = FALSE)

# GEO: CDM vs.Control (pediatric)
et_CDM_Pediatric <- exactTest(y, pair=c('Control (pediatric)', 'CDM') ) 
write.csv(topTags(et_CDM_Pediatric, n=10000), "Validation/[Batch-combatseq] EdgeR CDM vs. Control (pediatric).csv",row.names = FALSE)

# GTEx:  Control (slow death) vs. Control (unexpected death)
et_slow_unexpected <- exactTest(y, pair=c('Control (unexpected death)', 'Control (slow death)') ) 
write.csv(topTags(et_slow_unexpected, n=10000), "Validation/[Batch-combatseq] EdgeR slow vs. unexpected.csv",row.names = FALSE)

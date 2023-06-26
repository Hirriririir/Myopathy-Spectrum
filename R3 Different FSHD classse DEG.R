library(vroom)
library(tidyverse)
library(readxl)
library(edgeR)

merge_L <- vroom("Integration data/All_muscle_L_combatseq_tmm.csv")
Meta_all <- read_excel("Meta/adata_adj.obs.xlsx")
merge_data <- merge_L[,Meta_all$Sample_id]

# FSHD DEG #-----
unique(Meta_all$Phenotype_6)
mygroups <- Meta_all$Phenotype_6

y <- DGEList(counts=merge_data, genes=merge_L$Symbol, group = mygroups) #
y <- calcNormFactors(y, method ='TMM')
y <- estimateDisp(y)


# FSHD vs. fast death
et_FSHD_fast <- exactTest(y, pair=c('Control (fast death)', 'FSHD') ) 
write.csv(topTags(et_FSHD_fast, n=10000), "DEG/EdgeR FSHD vs. fast.csv",row.names = FALSE)

#  FSHD Vs. Control (accident death)
et_FSHD_accident <- exactTest(y, pair=c('Control (accident death)', 'FSHD') ) 
write.csv(topTags(et_FSHD_accident, n=10000), "DEG/EdgeR FSHD vs. accident death.csv",row.names = FALSE)

#  FSHD Vs. Control (slow death)
et_FSHD_slow <- exactTest(y, pair=c('Control (slow death)', 'FSHD') ) 
write.csv(topTags(et_FSHD_slow, n=10000), "DEG/EdgeR FSHD vs. slow death.csv",row.names = FALSE)

#  FSHD Vs. Control (diseased)
et_FSHD_diseased <- exactTest(y, pair=c('Control (diseased)', 'FSHD') ) 
write.csv(topTags(et_FSHD_diseased, n=10000), "DEG/EdgeR FSHD vs. diseased.csv",row.names = FALSE)

#  FSHD Vs. Control (amputee)
et_FSHD_Amputee <- exactTest(y, pair=c('Control (amputee)', 'FSHD') ) 
write.csv(topTags(et_FSHD_Amputee, n=10000), "DEG/EdgeR FSHD vs. Amputee.csv",row.names = FALSE)

# Random selection #-----

Fast_index <- which(Meta_all$Phenotype_4=='Control (fast death)')
Fast_index <- which(Meta_all$Phenotype_4=='Control (fast death)')
Slow_index <- which(Meta_all$Phenotype_4=='Control (slow death)')

## randomly select 8 fast death as controls
Meta_all$Phenotype_4[sample(Fast_index, size=24, replace =F)] <- 'fast_8'
Meta_all$Phenotype_4[sample(Fast_index, size=24, replace =F)] <- 'fast_24'
#Meta_all$Phenotype_4[sample(FSHD_index, size=30, replace =F)] <- 'FSHD_30'
Meta_all$Phenotype_4[sample(Slow_index, size=24, replace =F)] <- 'slow_24'

Meta_all[Meta_all$Phenotype_4=='fast_8', 'Sample_id']


unique(Meta_all$Phenotype_4)
mygroups <- Meta_all$Phenotype_4
y <- DGEList(counts=merge_data, genes=merge_L$Symbol, group = mygroups) #
y <- calcNormFactors(y, method ='TMM')
y <- estimateDisp(y)

## FSHD vs. fast_8 
et <- exactTest(y, pair=c('fast_8', 'FSHD') ) 
write.csv(topTags(et, n=10000), "DEG/EdgeR FSHD vs. fast_8.csv",row.names = FALSE)

## FSHD vs. Control (amputee) 
et <- exactTest(y, pair=c('Control (amputee)', 'FSHD') ) 
write.csv(topTags(et, n=10000), "DEG/EdgeR FSHD vs. Amputee.csv",row.names = FALSE)

## FSHD vs. Control (slow death) 
et <- exactTest(y, pair=c('Control (slow death)', 'FSHD') ) 
write.csv(topTags(et, n=10000), "DEG/EdgeR FSHD vs. slow death.csv",row.names = FALSE)

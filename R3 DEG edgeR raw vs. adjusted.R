library(vroom)
library(tidyverse)
library(readxl)
library(edgeR)
library(ggvenn)

merge_L <- vroom("Integration data/All_muscle_L_raw.csv")
Meta_all <- read_excel("Meta/adata_adj.obs.xlsx")
merge_data <- merge_L[,Meta_all$Sample_id]


# DEG batch validation #----
unique(Meta_all$Phenotype_2)
mygroups <- Meta_all$Phenotype_2

y <- DGEList(counts=merge_data, genes=merge_L$Symbol, group = mygroups) #
y <- calcNormFactors(y, method ='TMM')
y <- estimateDisp(y)

----# IBM #----

## IBM Vs. 9 Control (amputee) #----
GSE151757 <- vroom("Integration data/GSE151757_9199.csv")
GSE151757.meta <- read_excel("GEO data/GSE151757.meta.xlsx")
IBM_data <- GSE151757[,GSE151757.meta$Sample_id]
IBM_groups <- GSE151757.meta$Phenotype

y_IBM <- DGEList(counts=IBM_data, genes=GSE151757$Symbol, group = IBM_groups) #
y_IBM <- calcNormFactors(y_IBM, method ='TMM')
y_IBM <- estimateDisp(y_IBM)

et_IBM_Amputee <- exactTest(y_IBM, pair=c('AMP', 'IBM') ) 
write.csv(topTags(et_IBM_Amputee, n=10000), "Validation/[Batch-raw] EdgeR IBM vs. 9 Amputee.csv",row.names = FALSE)

IBM_integration <- read_csv('DEG/EdgeR IBM vs. fast.csv')
IBM_original <- read_csv('Validation/[Batch-raw] EdgeR IBM vs. 9 Amputee.csv')


## IBM Vs. 24 Control (amputee) #----
et_IBM_Amputee <- exactTest(y, pair=c('Control (amputee)', 'IBM') ) 
write.csv(topTags(et_IBM_Amputee, n=10000), "Validation/[Batch-raw] EdgeR IBM vs. Amputee.csv",row.names = FALSE)

## IBM veen #----
IBM_integration <- read_csv('DEG/EdgeR IBM vs. fast.csv')
IBM_original <- read_csv('Validation/[Batch-raw] EdgeR IBM vs. Amputee.csv')

x <- list(
  IBM_integration_DEG = IBM_integration %>% filter(IBM_integration$FDR < 0.05, abs(IBM_integration$logFC) > 0.5) %>% pull(genes),
  
  IBM_original_DEG = IBM_original %>% filter(IBM_original$FDR < 0.05,abs(IBM_original$logFC) > 0.5) %>% pull(genes)
)

ggvenn(
  x, 
  fill_color = c("#440154", '#35b779'),
  stroke_size = 1, set_name_size = 6, text_size = 8,
)

# CDM #----

## CDM vs.Control (pediatric)#----
et_CDM_Pediatric <- exactTest(y, pair=c('Control (pediatric)', 'CDM') ) 
write.csv(topTags(et_CDM_Pediatric, n=10000), "Validation/[Batch-raw] EdgeR CDM vs. Control (pediatric).csv",row.names = FALSE)

## CDM veen #----
CDM_integration <- read_csv('DEG/EdgeR CDM vs. fast.csv')
CDM_original <- read_csv('Validation/[Batch-raw] EdgeR CDM vs. Control (pediatric).csv')

x <- list(
  CDM_integration_DEG = CDM_integration %>% filter(CDM_integration$FDR < 0.05, abs(CDM_integration$logFC) > 0.5) %>% pull(genes),
  
  CDM_original_DEG = CDM_original %>% filter(CDM_original$FDR < 0.05,abs(CDM_integration$logFC) > 0.5) %>% pull(genes)
)

ggvenn(
  x, 
  fill_color = c("#440154", '#35b779'),
  stroke_size = 1, set_name_size = 6, text_size = 8,
)

CDM_integration_DEG = CDM_integration %>% filter(CDM_integration$FDR < 0.05, abs(IBM_integration$logFC) > 0.5)

# FSHD #----

## FSHD vs.Control (FSHD)#----
et_FSHD_FSHDCONTROL <- exactTest(y, pair=c('Control (FSHD)', 'FSHD') ) 
write.csv(topTags(et_FSHD_FSHDCONTROL, n=10000), "Validation/[Batch-raw] EdgeR FSHD vs. Control (FSHD).csv",row.names = FALSE)

## FSHD_2019 34 vs.Control (FSHD_2019) 9#----
et_FSHD_FSHDCONTROL <- exactTest(y, pair=c('Control (FSHD_2019)', 'FSHD_2019') ) 
write.csv(topTags(et_FSHD_FSHDCONTROL, n=10000), "Validation/[Batch-raw] EdgeR FSHD_2019 vs. Control (FSHD_2019).csv",row.names = FALSE)

## FSHD_2020 27 vs.Control (FSHD_2020) 8#----
et_FSHD_FSHDCONTROL <- exactTest(y, pair=c('Control (FSHD_2020)', 'FSHD_2020') ) 
write.csv(topTags(et_FSHD_FSHDCONTROL, n=10000), "Validation/[Batch-raw] EdgeR FSHD_2020 vs. Control (FSHD_2020).csv",row.names = FALSE)


## FSHD veen #----
FSHD_integration <- read_csv('DEG/EdgeR FSHD vs. fast.csv')
FSHD_original <- read_csv('Validation/[Batch-raw] EdgeR FSHD vs. Control (FSHD).csv')

x <- list(
  FSHD_integration_DEG = FSHD_integration %>% filter(FSHD_integration$FDR < 0.05, abs(FSHD_integration$logFC) > 0.5) %>% pull(genes),
  
  FSHD_original_DEG = FSHD_original %>% filter(FSHD_original$FDR < 0.05,abs(FSHD_original$logFC) > 0.5) %>% pull(genes)
)

ggvenn(
  x, 
  fill_color = c("#440154", '#35b779'),
  stroke_size = 1, set_name_size = 6, text_size = 8,
)

# LGMDR12 #----

## LGMDR12 vs.Control (LGMDR12)#----
et_LGMDR12_LGMDR12CONTROL <- exactTest(y, pair=c('Control (LGMD R12)', 'LGMD R12') ) 
write.csv(topTags(et_LGMDR12_LGMDR12CONTROL, n=10000), "Validation/[Batch-raw] EdgeR LGMDR12 vs. Control (LGMD R12).csv",row.names = FALSE)

## LGMDR12 veen #----
LGMDR12_integration <- read_csv('DEG/EdgeR LGMDR12 vs. fast.csv')
LGMDR12_original <- read_csv('Validation/[Batch-raw] EdgeR LGMDR12 vs. Control (LGMD R12).csv')

x <- list(
  LGMDR12_integration_DEG = LGMDR12_integration %>% filter(LGMDR12_integration$FDR < 0.05, abs(LGMDR12_integration$logFC) > 0.5) %>% pull(genes),
  
  LGMDR12_original_DEG = LGMDR12_original %>% filter(LGMDR12_original$FDR < 0.05,abs(LGMDR12_original$logFC) > 0.5) %>% pull(genes)
)

ggvenn(
  x, 
  fill_color = c("#440154", '#35b779'),
  stroke_size = 1, set_name_size = 6, text_size = 8,
)




# GEO: CDM vs.Control (pediatric)
et_CDM_Pediatric <- exactTest(y, pair=c('Control (pediatric)', 'CDM') ) 
write.csv(topTags(et_CDM_Pediatric, n=10000), "Validation/[Batch-raw] EdgeR CDM vs. Control (pediatric).csv",row.names = FALSE)

# GTEx:  Control (slow death) vs. Control (unexpected death)
et_slow_unexpected <- exactTest(y, pair=c('Control (unexpected death)', 'Control (slow death)') ) 
write.csv(topTags(et_slow_unexpected, n=10000), "Validation/[Batch-raw] EdgeR slow vs. unexpected.csv",row.names = FALSE)



library(vroom)
library(tidyverse)
library(readxl)
library(edgeR)

merge_L <- vroom("Integration data/All_muscle_L_raw.csv")
Meta_all <- read_excel("Meta/adata_adj.obs.xlsx")
merge_data <- merge_L[,Meta_all$Sample_id]


# DEG batch validation #----
unique(Meta_all$Phenotype_6)
mygroups <- Meta_all$Phenotype_6

y <- DGEList(counts=merge_data, genes=merge_L$Symbol, group = mygroups) #
y <- calcNormFactors(y, method ='TMM')
y <- estimateDisp(y)

# Helsinki: Titinopathy vs Control (amputee)
et_Titinopathy_Amputee <- exactTest(y, pair=c('Control (amputee)', 'Titinopathy') ) 
write.csv(topTags(et_Titinopathy_Amputee, n=10000), "Validation/[Batch-raw] EdgeR Titinopathy vs. Amputee.csv",row.names = FALSE)

# GEO: CDM vs.Control (pediatric)
et_CDM_Pediatric <- exactTest(y, pair=c('Control (pediatric)', 'CDM') ) 
write.csv(topTags(et_CDM_Pediatric, n=10000), "Validation/[Batch-raw] EdgeR CDM vs. Control (pediatric).csv",row.names = FALSE)

# GTEx:  Control (slow death) vs. Control (unexpected death)
et_slow_unexpected <- exactTest(y, pair=c('Control (unexpected death)', 'Control (slow death)') ) 
write.csv(topTags(et_slow_unexpected, n=10000), "Validation/[Batch-raw] EdgeR slow vs. unexpected.csv",row.names = FALSE)


# Batch validation veen
library(ggvenn)

## Helsinki
Helsinki_combatseq <- read_csv('Validation/[Batch-combatseq] EdgeR Titinopathy vs. Amputee.csv')
Helsinki_raw <- read_csv('Validation/[Batch-raw] EdgeR Titinopathy vs. Amputee.csv')

x <- list(
  H_combatseq = Helsinki_combatseq %>% filter(Helsinki_combatseq$FDR < 0.05,Helsinki_combatseq$logFC > 1) %>% pull(genes),
  
  H_raw = Helsinki_raw %>% filter(Helsinki_raw$FDR < 0.05,Helsinki_raw$logFC > 1) %>% pull(genes)
)

ggvenn(
  x, 
  fill_color = c('#35b779', "#fde725"),
  stroke_size = 0.5, set_name_size = 4
)

## GEO
GEO_combatseq <- read_csv('Validation/[Batch-combatseq] EdgeR CDM vs. Control (pediatric).csv')
GEO_raw <- read_csv('Validation/[Batch-raw] EdgeR CDM vs. Control (pediatric).csv')

x <- list(
  GEO_combatseq = GEO_combatseq %>% filter(GEO_combatseq$FDR < 0.05,GEO_combatseq$logFC > 1) %>% pull(genes),
  
  GEO_raw = GEO_raw %>% filter(GEO_raw$FDR < 0.05,GEO_raw$logFC > 1) %>% pull(genes)
)

ggvenn(
  x, 
  fill_color = c('#35b779', "#fde725"),
  stroke_size = 0.5, set_name_size = 4
)

## GTEx
GTEx_combatseq <- read_csv('Validation/[Batch-combatseq] EdgeR slow vs. unexpected.csv')
GTEx_raw <- read_csv('Validation/[Batch-raw] EdgeR slow vs. unexpected.csv')

x <- list(
  GTEx_combatseq = GTEx_combatseq %>% filter(GTEx_combatseq$FDR < 0.05,GTEx_combatseq$logFC > 1) %>% pull(genes),
  
  GTEx_raw = GTEx_raw %>% filter(GTEx_raw$FDR < 0.05,GTEx_raw$logFC > 1) %>% pull(genes)
)

ggvenn(
  x, 
  fill_color = c('#35b779', "#fde725"),
  stroke_size = 0.5, set_name_size = 4
)

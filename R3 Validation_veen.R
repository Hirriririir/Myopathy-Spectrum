library(vroom)
library(tidyverse)
library(readxl)
library(ggvenn)


# FSHD GSE140261 #----
## Supplemental Table 5a. The robustly up-regulated genes in FSHD sample categorized in the High Class. The DUX4-induced,  inflammatory, immune, extraceullular matrix and stress biomarkers are tagged in column J-N.

GSE140261_FSHD_high_original <- read_excel('Validation/FSHD PMID32083293 suppl_table_5_candidate_biomarkers_and_enriched_go_ddaa031.xlsx', sheet='FSHD high class DEG')

GSE140261_FSHD_IGhigh_original <- read_excel('Validation/FSHD PMID32083293 suppl_table_5_candidate_biomarkers_and_enriched_go_ddaa031.xlsx', sheet='FSHD IG-High class DEG')

GSE140261_FSHD_moderate_original <- read_excel('Validation/FSHD PMID32083293 suppl_table_5_candidate_biomarkers_and_enriched_go_ddaa031.xlsx', sheet='FSHD moderate class DEG')

GSE140261_FSHD_mild_original <- read_excel('Validation/FSHD PMID32083293 suppl_table_7_mild_fshd_164_potential_markers_ddaa031.xlsx', sheet='FSHD mild class DEG')


EdgeR_result <- read_csv('DEG/EdgeR FSHD vs. fast.csv')


x <- list(
  FSHD_integration = EdgeR_result %>% filter(EdgeR_result$FDR < 0.05,EdgeR_result$logFC > 1) %>% pull(genes),
  
  FSHD_high = GSE140261_FSHD_high_original %>% filter(GSE140261_FSHD_high_original$padj < 0.05,GSE140261_FSHD_high_original$log2FoldChange > 1) %>% pull(gene_name),
  
  FSHD_IGhigh = GSE140261_FSHD_IGhigh_original %>% filter(GSE140261_FSHD_IGhigh_original$padj < 0.05,GSE140261_FSHD_IGhigh_original$log2FoldChange > 1) %>% pull(gene_name),
  
  FSHD_moderate = GSE140261_FSHD_moderate_original %>% filter(GSE140261_FSHD_moderate_original$padj < 0.05,GSE140261_FSHD_moderate_original$log2FoldChange > 1) %>% pull(gene_name), 
  
  FSHD_mild = GSE140261_FSHD_mild_original %>% filter(GSE140261_FSHD_mild_original$padj < 0.05,GSE140261_FSHD_mild_original$log2FoldChange > 1) %>% pull(gene_name)
  )

ggvenn(
  x, 
  fill_color = c('#252526', "#440154", '#3d6682','#35b779', "#fde725"),
  stroke_size = 0.5, set_name_size = 4
)


intersect(x$FSHD_integration, x$FSHD_moderate)


# bar plot

df <- data.frame(class=c("Integration", "FSHD (high)", "FSHD (IG-high)", "FSHD (moderate)", 'FSHD (mild)'),
                 count=c(22, 5, 4, 1, 0))

ggplot(data=df, aes(x=class, y=count)) +
  geom_bar(stat="identity")+
  theme_minimal()

ggplot(data=df, aes(x='class', y='count')) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=len), vjust=-0.3, size=3.5)+
  theme_minimal()
# load data #----
library(vroom)
library(tidyverse)
library(readxl)

# Data import (all muscles) #-----

merge_L <- vroom("Integration data/All_muscle_L_raw.csv")
merge_data <- merge_L[2:1222]
merge_data[is.na(merge_data)] <- 0
merge_data <- as.matrix(merge_data)

colnames(merge_data)


# ComBat-seq #-----

Meta_GEO <- read_excel("Meta/RNA-seq integration meta.xlsx", sheet = "GEO")
Meta_GTEx = read_excel('Meta/RNA-seq integration meta.xlsx', sheet='GTEx')
Meta_Helsinki <- read_excel("Meta/RNA-seq integration meta.xlsx", sheet = "Helsinki")

Meta_all <- bind_rows(Meta_GEO, Meta_GTEx, Meta_Helsinki)
batch <- Meta_all$Geo_accession #  'Geo_accession' represents 'batch' here

merge_data_adjusted <- sva::ComBat_seq(merge_data, batch=batch, group=NULL)
merge_L[2:1222] <- merge_data_adjusted

vroom_write(merge_L, "Integration data/All_muscle_L_combatseq.csv")




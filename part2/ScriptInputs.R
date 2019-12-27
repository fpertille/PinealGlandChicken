################################
# Author: F?bio P?rtille
# Script description: PINEAL MEDIP ANALYSIS
# Date: OCTOBER 10, 2019
# Files: Fabio P?rtille
#        input: ??
#         output: "???"
################################
#to get a node
#salloc -A b2014219 -t 200:00:00 -p core -n 1 &
#going to appropriate folder within the node
#  cd /home/fibio/Project/MeDIP/RJF/Stacks
#start R in that folder
#Brazilian server:
#qsub -I -q long
#qsub -I -q short
#################################
#Brazil
#################################
. /home/program2/bin/.miniconda3/etc/profile.d/conda.sh
conda activate R-3.5_env

module load HPC/R-3.5.1

R

#installing Modifier###

#devtools::install_git(url = "https://gitlab.com/Gustafsson-lab/MODifieR.git")
#devtools::install_git(url = "https://gitlab.com/Gustafsson-lab/MODifieR.git", args ="--no-multiarch")
#devtools::install_git(url = "https://gitlab.com/Gustafsson-lab/MODifieR.git", args ="--no-lock", args ="--no-multiarch")
#unlink("H:/Documents/R/win-library/3.6/00LOCK-MODifieR", recursive = TRUE)

library(MODifieR)
library(plyr)
library(data.table)
library(limma)
library(Glimma)
library(edgeR)
library(org.Gg.eg.db)
library(readr)
library ("ChIPseeker")
library(RVenn)
library(purrr)
library(ggplot2)
library (ChIPpeakAnno)
library (GenomicFeatures)
#library("mason")

###########
#first saving all the files and than, I will need to run the deep analysis.
setwd("H:/Medicina Veterinária/2017.POSDOC-ESALQ-USP/ExperimentoPineal/INTEGRATION/")
#First saving objetcts from the last scripts
#save(overlaps, file="overlaps.rda")
#save(overlaps_5k, file="overlaps_5k.rda")
#save(overlaps_10k, file="overlaps_10k.rda")
#save(annot_exprs, file= "annot_exprs.rda")
#save(res_comp1 , file="res_comp1.rda")
#save(res_comp2 , file="res_comp2.rda")
#save(res_comp3 , file="res_comp3.rda")
#save(res_comp4 , file="res_comp4.rda")
#save(res_comp5 , file="res_comp5.rda")
#save(res_comp1_medips , file="res_comp1_medips.rda")
#save(res_comp2_medips , file="res_comp2_medips.rda")
#save(res_comp3_medips , file="res_comp3_medips.rda")
#save(res_comp4_medips , file="res_comp4_medips.rda")
#save(res_comp5_medips , file="res_comp5_medips.rda")
#save(res_comp1_medips_base , file="res_comp1_medips_base.rda")
#save(res_comp2_medips_base , file="res_comp2_medips_base.rda")
#save(res_comp3_medips_base , file="res_comp3_medips_base.rda")
#save(res_comp4_medips_base , file="res_comp4_medips_base.rda")
#save(res_comp5_medips_base , file="res_comp5_medips_base.rda")
#save(res_comp1_mirs , file="res_comp1_mirs.rda")
#save(res_comp2_mirs , file="res_comp2_mirs.rda")
#save(res_comp3_mirs , file="res_comp3_mirs.rda")
#save(res_comp4_mirs , file="res_comp4_mirs.rda")
#save(res_comp5_mirs , file="res_comp5_mirs.rda")

#saving frames
#write.table(overlaps , "overlaps.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
#write.table(res_comp1 , "res_comp1.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp2 , "res_comp2.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp3 , "res_comp3.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp4 , "res_comp4.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp5 , "res_comp5.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp1_medips , "res_comp1_medips.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp2_medips , "res_comp2_medips.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp3_medips , "res_comp3_medips.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp4_medips , "res_comp4_medips.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp5_medips , "res_comp5_medips.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp1_medips_base , "res_comp1_medips_base.txt" , sep="\t" , row.names = T , quote = F, col.names = NA)
#write.table(res_comp2_medips_base , "res_comp2_medips_base.txt" , sep="\t" , row.names = T , quote = F, col.names = NA)
#write.table(res_comp3_medips_base , "res_comp3_medips_base.txt" , sep="\t" , row.names = T , quote = F, col.names = NA)
#write.table(res_comp4_medips_base , "res_comp4_medips_base.txt" , sep="\t" , row.names = T , quote = F, col.names = NA)
#write.table(res_comp5_medips_base , "res_comp5_medips_base.txt" , sep="\t" , row.names = T , quote = F, col.names = NA)
#write.table(res_comp1_mirs , "res_comp1_mirs.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp2_mirs , "res_comp2_mirs.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp3_mirs , "res_comp3_mirs.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp4_mirs , "res_comp4_mirs.txt" , sep="\t" , row.names = F , quote = F)
#write.table(res_comp5_mirs , "res_comp5_mirs.txt" , sep="\t" , row.names = F , quote = F)

#Reading objects - run in the server FROM HERE
#setwd("/home/fabio/postdoc/PINEAL/Ranalysis/InputsForDeep")
GFM_mod <- readRDS(file="GFM_mod.rds", refhook = NULL)

load(file="overlaps.rda")
load(file="overlaps_5k.rda")
load(file="overlaps_10k.rda")
load(file="annot_exprs.rda")
load(file="res_comp1.rda")
load(file="res_comp2.rda")
load(file="res_comp3.rda")
load(file="res_comp4.rda")
load(file="res_comp5.rda")

load(file="res_comp1_medips.rda")
load(file="res_comp2_medips.rda")
load(file="res_comp3_medips.rda")
load(file="res_comp4_medips.rda")
load(file="res_comp5_medips.rda")

load(file="res_comp1_medips_base.rda")
load(file="res_comp2_medips_base.rda")
load(file="res_comp3_medips_base.rda")
load(file="res_comp4_medips_base.rda")
load(file="res_comp5_medips_base.rda")

load(file="res_comp1_mirs.rda")
load(file="res_comp2_mirs.rda")
load(file="res_comp3_mirs.rda")
load(file="res_comp4_mirs.rda")
load(file="res_comp5_mirs.rda")

#reading frames
#library(readr)
#overlaps <- read.table("overlaps.txt", header = T)


##########################################################################
#Merging the overlaps with the genes from gene expression
##########################################################################
###merging with gene###comp5####
overlaps_main_info <- overlaps[,c("peaks1" , "peaks2" , "overlapFeature" , "shortestDistance")]
res_comp5_medips_base$peaks1 <- rownames(res_comp5_medips_base)
res_comp5_medips_g <- merge(res_comp5_medips_base , overlaps_main_info , all.x=T)

###merging with gene_5k###
overlaps_5k_main_info <- overlaps_5k[,c("peaks1" , "peaks2" , "overlapFeature" , "shortestDistance")]
res_comp5_medips_g5 <- merge(res_comp5_medips_base , overlaps_5k_main_info , all.x=T)

###merging with gene_10k###
overlaps_10k_main_info <- overlaps_10k[,c("peaks1" , "peaks2" , "overlapFeature" , "shortestDistance")]
res_comp5_medips_g10 <- merge(res_comp5_medips_base , overlaps_10k_main_info , all.x=T)

###merging with gene###comp1####
res_comp1_medips_base$peaks1 <- rownames(res_comp1_medips_base)
res_comp1_medips_g <- merge(res_comp1_medips_base , overlaps_main_info , all.x=T)
###merging with gene_5k###
res_comp1_medips_g5 <- merge(res_comp1_medips_base , overlaps_5k_main_info , all.x=T)
###merging with gene_10k###
res_comp1_medips_g10 <- merge(res_comp1_medips_base , overlaps_10k_main_info , all.x=T)
###merging with gene###comp2####
res_comp2_medips_base$peaks1 <- rownames(res_comp2_medips_base)
res_comp2_medips_g <- merge(res_comp2_medips_base , overlaps_main_info , all.x=T)
###merging with gene_5k###
res_comp2_medips_g5 <- merge(res_comp2_medips_base , overlaps_5k_main_info , all.x=T)
###merging with gene_10k###
res_comp2_medips_g10 <- merge(res_comp2_medips_base , overlaps_10k_main_info , all.x=T)
###merging with gene###comp3####
res_comp3_medips_base$peaks1 <- rownames(res_comp3_medips_base)
res_comp3_medips_g <- merge(res_comp3_medips_base , overlaps_main_info , all.x=T)
###merging with gene_5k###
res_comp3_medips_g5 <- merge(res_comp3_medips_base , overlaps_5k_main_info , all.x=T)
###merging with gene_10k###
res_comp3_medips_g10 <- merge(res_comp3_medips_base , overlaps_10k_main_info , all.x=T)
###merging with gene###comp4####
res_comp4_medips_base$peaks1 <- rownames(res_comp4_medips_base)
res_comp4_medips_g <- merge(res_comp4_medips_base , overlaps_main_info , all.x=T)
###merging with gene_5k###
res_comp4_medips_g5 <- merge(res_comp4_medips_base , overlaps_5k_main_info , all.x=T)
###merging with gene_10k###
res_comp4_medips_g10 <- merge(res_comp4_medips_base , overlaps_10k_main_info , all.x=T)

###merging with gene###comp5####
res_comp5_medips_base$peaks1 <- rownames(res_comp5_medips_base)
res_comp5_medips_g <- merge(res_comp5_medips_base , overlaps_main_info , all.x=T)
###merging with gene_5k###
res_comp5_medips_g5 <- merge(res_comp5_medips_base , overlaps_5k_main_info , all.x=T)
###merging with gene_10k###
res_comp5_medips_g10 <- merge(res_comp5_medips_base , overlaps_10k_main_info , all.x=T)

#################################
Gallus_network_PPI_700 <- read_table2("Gallus_network_PPI_700.txt")
#################################
#Sys.setenv(PATH = paste("/Users/fabio/miniconda3/bin", Sys.getenv("PATH"), sep=":"))
#set.seed(143)
#if(.Platform$OS.type == "windows") Sys.setenv(PATH= paste("C:/Users/fabio/Miniconda3/Lib/bin",Sys.getenv()["PATH"],sep=";"))
#library(reticulate)


#getting the number of each test we did!
length(unique(res_comp1_medips_g$peaks2[res_comp1_medips_g$P.Value<0.05]))
length(unique(res_comp2_medips_g$peaks2[res_comp2_medips_g$P.Value<0.05]))
length(unique(res_comp3_medips_g$peaks2[res_comp3_medips_g$P.Value<0.05]))
length(unique(res_comp4_medips_g$peaks2[res_comp4_medips_g$P.Value<0.05]))
length(unique(res_comp5_medips_g$peaks2[res_comp5_medips_g$P.Value<0.05]))


length(unique(res_comp1_medips_g5$peaks2[res_comp1_medips_g5$P.Value<0.05]))
length(unique(res_comp2_medips_g5$peaks2[res_comp2_medips_g5$P.Value<0.05]))
length(unique(res_comp3_medips_g5$peaks2[res_comp3_medips_g5$P.Value<0.05]))
length(unique(res_comp4_medips_g5$peaks2[res_comp4_medips_g5$P.Value<0.05]))
length(unique(res_comp5_medips_g5$peaks2[res_comp5_medips_g5$P.Value<0.05]))


length(unique(res_comp1_medips_g10$peaks2[res_comp1_medips_g10$P.Value<0.05]))
length(unique(res_comp2_medips_g10$peaks2[res_comp2_medips_g10$P.Value<0.05]))
length(unique(res_comp3_medips_g10$peaks2[res_comp3_medips_g10$P.Value<0.05]))
length(unique(res_comp4_medips_g10$peaks2[res_comp4_medips_g10$P.Value<0.05]))
length(unique(res_comp5_medips_g10$peaks2[res_comp5_medips_g10$P.Value<0.05]))

length(unique(res_comp1_mirs$mirbase_id[res_comp1_mirs$P.Val<0.05]))
length(unique(res_comp2_mirs$mirbase_id[res_comp2_mirs$P.Val<0.05]))
length(unique(res_comp3_mirs$mirbase_id[res_comp3_mirs$P.Val<0.05]))
length(unique(res_comp4_mirs$mirbase_id[res_comp4_mirs$P.Val<0.05]))
length(unique(res_comp5_mirs$mirbase_id[res_comp5_mirs$P.Val<0.05]))

length(unique(res_comp1$ensembl_id[res_comp1$adj.P.Val<0.2]))
length(unique(res_comp2$ensembl_id[res_comp2$adj.P.Val<0.2]))
length(unique(res_comp3$ensembl_id[res_comp3$adj.P.Val<0.2]))
length(unique(res_comp4$ensembl_id[res_comp4$adj.P.Val<0.2]))
length(unique(res_comp5$ensembl_id[res_comp5$adj.P.Val<0.2]))

length(unique(res_comp1_medips_base$loc_id[res_comp1_medips_base$adj.P.Val<0.2]))
length(unique(res_comp2_medips_base$loc_id[res_comp2_medips_base$adj.P.Val<0.2]))
length(unique(res_comp3_medips_base$loc_id[res_comp3_medips_base$adj.P.Val<0.2]))
length(unique(res_comp4_medips_base$loc_id[res_comp4_medips_base$adj.P.Val<0.2]))
length(unique(res_comp5_medips_base$loc_id[res_comp5_medips_base$adj.P.Val<0.2]))

length(unique(res_comp1_medips_g$peaks2[res_comp1_medips_g$adj.P.Val<0.2]))
length(unique(res_comp2_medips_g$peaks2[res_comp2_medips_g$adj.P.Val<0.2]))
length(unique(res_comp3_medips_g$peaks2[res_comp3_medips_g$adj.P.Val<0.2]))
length(unique(res_comp4_medips_g$peaks2[res_comp4_medips_g$adj.P.Val<0.2]))
length(unique(res_comp5_medips_g$peaks2[res_comp5_medips_g$adj.P.Val<0.2]))

length(unique(res_comp1_medips_g5$peaks2[res_comp1_medips_g5$adj.P.Val<0.2]))
length(unique(res_comp2_medips_g5$peaks2[res_comp2_medips_g5$adj.P.Val<0.2]))
length(unique(res_comp3_medips_g5$peaks2[res_comp3_medips_g5$adj.P.Val<0.2]))
length(unique(res_comp4_medips_g5$peaks2[res_comp4_medips_g5$adj.P.Val<0.2]))
length(unique(res_comp5_medips_g5$peaks2[res_comp5_medips_g5$adj.P.Val<0.2]))


length(unique(res_comp1_medips_g10$peaks2[res_comp1_medips_g10$adj.P.Val<0.2]))
length(unique(res_comp2_medips_g10$peaks2[res_comp2_medips_g10$adj.P.Val<0.2]))
length(unique(res_comp3_medips_g10$peaks2[res_comp3_medips_g10$adj.P.Val<0.2]))
length(unique(res_comp4_medips_g10$peaks2[res_comp4_medips_g10$adj.P.Val<0.2]))
length(unique(res_comp5_medips_g10$peaks2[res_comp5_medips_g10$adj.P.Val<0.2]))

length(unique(res_comp1_mirs$mirbase_id[res_comp1_mirs$adj.P.Val<0.2]))
length(unique(res_comp2_mirs$mirbase_id[res_comp2_mirs$adj.P.Val<0.2]))
length(unique(res_comp3_mirs$mirbase_id[res_comp3_mirs$adj.P.Val<0.2]))
length(unique(res_comp4_mirs$mirbase_id[res_comp4_mirs$adj.P.Val<0.2]))
length(unique(res_comp5_mirs$mirbase_id[res_comp5_mirs$adj.P.Val<0.2]))

length(unique(res_comp1_mirs$ensembl_id[res_comp1_mirs$P.Val<0.05 & res_comp1_mirs$pred_score>90]))
length(unique(res_comp2_mirs$ensembl_id[res_comp2_mirs$P.Val<0.05 & res_comp2_mirs$pred_score>90]))
length(unique(res_comp3_mirs$ensembl_id[res_comp3_mirs$P.Val<0.05 & res_comp3_mirs$pred_score>90]))
length(unique(res_comp4_mirs$ensembl_id[res_comp4_mirs$P.Val<0.05 & res_comp4_mirs$pred_score>90]))
length(unique(res_comp5_mirs$ensembl_id[res_comp5_mirs$P.Val<0.05 & res_comp5_mirs$pred_score>90]))

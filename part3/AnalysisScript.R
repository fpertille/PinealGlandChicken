################################
# Author: Fabio Pertille
# Script description: PINEAL MEDIP ANALYSIS
# Date: OCTOBER 18, 2019
# Files: Fabio Partille
#        input: from ScriptInput.R 
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
#. /home/program2/bin/.miniconda3/etc/profile.d/conda.sh
#conda activate R-3.5_env

#module load HPC/R-3.5.1
#R
setwd("H:/Medicina Veterinária/2017.POSDOC-ESALQ-USP/ExperimentoPineal/INTEGRATION/")
library (BSgenome.Ggallus.UCSC.galGal6)
library(org.Gg.eg.db)
##############################################################################
##############################################################################
##############################################################################
#########################SEX DIFFERENCES######################################
##############################################################################

library(stringr)
DMRsig_SMvF<- subset(res_comp1_medips, adj.P.Val<0.2)
collumn <- as.data.frame(str_split_fixed(DMRsig_SMvF$loc_id, "_", 3))
names(collumn) <- c("chr", "start", "end")
DMRsig_SMvF <- cbind(collumn, DMRsig_SMvF)
DMRsig_SMvF$strand <- "*"
DMRsig_SMvF$Test <- "SMvF";
DMRsig_SMvF$start <- as.numeric(as.character(DMRsig_SMvF$start))
DMRsig_SMvF$end <- as.numeric(as.character((DMRsig_SMvF$end)))
GR_DMRsig_SMvF <- with(DMRsig_SMvF, GRanges(chr, IRanges(start, end), strand, loc_id, logFC, AveExpr, t, P.Value, adj.P.Val, Test))

DMRsig_CMvF<- subset(res_comp2_medips, adj.P.Val<0.2)
collumn <- as.data.frame(str_split_fixed(DMRsig_CMvF$loc_id, "_", 3))
names(collumn) <- c("chr", "start", "end")
DMRsig_CMvF <- cbind(collumn, DMRsig_CMvF)
DMRsig_CMvF$strand <- "*"
DMRsig_CMvF$Test <- "CMvF";
DMRsig_CMvF$start <- as.numeric(as.character(DMRsig_CMvF$start))
DMRsig_CMvF$end <- as.numeric(as.character((DMRsig_CMvF$end)))
GR_DMRsig_CMvF <- with(DMRsig_CMvF, GRanges(chr, IRanges(start, end), strand, loc_id, logFC, AveExpr, t, P.Value, adj.P.Val, Test))
CSMvF <- rbind(DMRsig_CMvF, DMRsig_SMvF)
GR_CSMvF <- with(CSMvF, GRanges(chr, IRanges(start, end), strand, loc_id, logFC, AveExpr, t, P.Value, adj.P.Val, Test))

ol <- findOverlapsOfPeaks(GR_DMRsig_SMvF, GR_DMRsig_CMvF)
GR_SexLap1 <- as.data.frame (subsetByOverlaps(GR_DMRsig_SMvF, GR_DMRsig_CMvF))
GR_SexLap2 <- as.data.frame (subsetByOverlaps(GR_DMRsig_CMvF, GR_DMRsig_SMvF))

write.table(GR_SexLap1 , "SexLap1.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
write.table(GR_SexLap2 , "SexLap2.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
#write.table(as.data.frame(ol$overlappingPeaks$`GR_DMRsig_SMvF///GR_DMRsig_CMvF`) , "SexLap.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
grl <- splitAsList(GR_CSMvF, GR_CSMvF$Test)
tiff("Venn_OverlapedCSMvF.tiff", units="cm", width=30.00, height=30.00, res=600, compression = "lzw")
#CairoPDF("venn.pdf")
res <- makeVennDiagram(Peaks=grl, NameOfPeaks=names(grl), main = "DMRs between CMvF and SMvF contrasts ", col = "black", alpha = 0.7, cex= 5.5, ignore.strand=T, scaled=F, fill = c("darkorange1", "deepskyblue2"))
#fill = c("Yellow", "cyan"), 
dev.off()

#merged peaklist
GR_SexLap <- as.data.frame(ol$peaklist[["GR_DMRsig_SMvF///GR_DMRsig_CMvF"]])
GR_SexUniqueC <- as.data.frame(ol$peaklist[["GR_DMRsig_CMvF"]])
GR_SexUniqueS <- as.data.frame(ol$peaklist[["GR_DMRsig_SMvF"]])
#all the peaks
SexLap <- as.data.frame(ol$overlappingPeaks[["GR_DMRsig_SMvF///GR_DMRsig_CMvF"]])
SexLap <- as.data.frame(GR_SexLap[!duplicated(GR_SexLap[,c('loc_id')]),])
SexUnique <- as.data.frame(ol$uniquePeaks)
SexUnique <- as.data.frame(GR_SexUnique[!duplicated(GR_SexUnique[,c('loc_id')]),])

#####annotating that peaks#######
library ("ChIPseeker")
library (GenomicFeatures)
supportedUCSCtables(genome = "galGal6")
gg_txdb <- makeTxDbFromUCSC(genome="galGal6", tablename="ensGene")

peak_file <- as.data.frame(SexLap)
peak_file <- as.data.frame(SexUnique)
####################################
peak_file <- subset(peak_file, select=c(seqnames, start,end))
names(peak_file) <- c("CHR", "BP", "BP2")

peak_anno_function <- function(input_data){
  b2 <- input_data[,c("CHR" , "BP" , "BP2")]
  write.table(b2 , "annotate_file.txt" , sep = "\t" , col.names = TRUE , row.names = FALSE)
  peakAnno <- annotatePeak("H:/Medicina Veterinária/2017.POSDOC-ESALQ-USP/ExperimentoPineal/INTEGRATION/annotate_file.txt",
                           TxDb=gg_txdb, annoDb="org.Gg.eg.db")
  return(peakAnno)
}

peak_chicken <- peak_anno_function(peak_file)
plotAnnoPie(peak_chicken)
info_chicken <- peak_chicken@anno@elementMetadata

SexLap <- cbind(SexLap, info_chicken , all=T)
SexUnique <- cbind(SexUnique, info_chicken , all=T)

write.table(SexLap, "SexLap.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
write.table(SexUnique, "SexUnique.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)

#END###########END####################END###################END###############
##############################################################################
##########################DMRbyLOCATION USING GFM#############################
##############################################################################
###DMR regard to genes
peak_file <- as.data.frame(res_comp5_medips_g10)
#peak_file <- as.data.frame(res_comp5_medips_g5)
#peak_file <- as.data.frame(res_comp5_medips_g)

library(dplyr)
library(tidyr)
library(stringr)
peak_file <- as.data.frame(str_split_fixed(peak_file$peaks1, "_", 3))
names(peak_file) <- c("CHR", "BP", "BP2")

library ("ChIPseeker")
supportedUCSCtables(genome = "galGal6")
gg_txdb <- makeTxDbFromUCSC(genome="galGal6", tablename="ensGene")
####################################

peak_anno_function <- function(input_data){
  b2 <- input_data[,c("CHR" , "BP" , "BP2")]
  write.table(b2 , "annotate_file.txt" , sep = "\t" , col.names = TRUE , row.names = FALSE)
  peakAnno <- annotatePeak("H:/Medicina Veterinária/2017.POSDOC-ESALQ-USP/ExperimentoPineal/INTEGRATION/annotate_file.txt",
                           TxDb=gg_txdb, annoDb="org.Gg.eg.db")
  return(peakAnno)
}

peak_chicken <- peak_anno_function(peak_file)
plotAnnoPie(peak_chicken)
info_chicken <- peak_chicken@anno@elementMetadata

DMRall_GFM10_annotation <- cbind(res_comp5_medips_g10, info_chicken , all=T)
#DMRall_GFM5_annotation <- cbind(res_comp5_medips_g5, info_chicken , all=T)
#DMRall_GFMg_annotation <- cbind(res_comp5_medips_g, info_chicken , all=T)

DMRall_GFM_annotation <- as.data.frame(DMRall_GFM_annotation)
#DMRall_GFM5_annotation <- as.data.frame(DMRall_GFM5_annotation)
#DMRall_GFMg_annotation <- as.data.frame(DMRall_GFMg_annotation)
#to filter inside genes use filter to exon, intron, 3 and 3´UTR , and <1kb dowstream

#write.table(DMRall_GFM5_annotation, "DMR_gene_GFM.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
#write.table(DMRall_GFMg_annotation, "DMR_gene5_GFM.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
#write.table(DMRall_GFM10_annotation, "DMR_gene10_GFM.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)

#I SEPARETED THE WITHIN GENES FROM TSS-5K AND 5-10KG USING EXCEL 
###DMR gene LAP####
x <- scan("DMRgeneLAP.txt", what="", sep="\n")
y <- strsplit(x, "[[:space:]]+")
names(y) <- sapply(y, `[[`, 1)
y <- lapply(y, `[`, -1)
DMR_gene_LAP<-y

DMR_gene_LAPs <- c()
DMR_gene_LAPs$DMRWithinGENE <- annot_exprs$gene_id[annot_exprs$ensembl_id %in% DMR_gene_LAP$DMRWithinGENE]
DMR_gene_LAPs$DMRTSS_5K<- annot_exprs$gene_id[annot_exprs$ensembl_id %in% DMR_gene_LAP$DMRTSS_5K]
DMR_gene_LAPs$DMR5K_10K<- annot_exprs$gene_id[annot_exprs$ensembl_id %in% DMR_gene_LAP$DMR5K_10K]

go_enrich <- compareCluster(geneClusters = DMR_gene_LAPs, 
                                fun = "enrichGO" , 
                                OrgDb = org.Gg.eg.db , 
                                ont="BP" , 
                                pvalueCutoff=1, 
                                qvalueCutoff=1, 
                                readable=T)
tiff("PathwayGO_GFM_DMR_Genelap.tiff", units="cm", width=40.00, height=40.00, res=600, compression = "lzw")
dotplot(go_mod_enrich , showCategory=20)
dev.off()

kegg_enrich <- compareCluster(geneClusters = DMR_gene_LAPs, 
                                  fun = "enrichKEGG" , 
                                  org="gga",
                                  pvalueCutoff=1, 
                                  qvalueCutoff=1)
tiff("PathwaygGFM_DMR_Genelap.tiff", units="cm", width=25.00, height=15.00, res=600, compression = "lzw")
dotplot(kegg_mod_enrich)
dev.off()

setwd("H:/Medicina Veterinária/2017.POSDOC-ESALQ-USP/ExperimentoPineal/INTEGRATION/Pathway/")
write.table(as.data.frame(go_enrich@compareClusterResult) , "PathwayGO_Enrich_DMRgene.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
write.table(as.data.frame(kegg_enrich@compareClusterResult) , "PathwayKEGG_Enrich_DMRgene.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
setwd("H:/Medicina Veterinária/2017.POSDOC-ESALQ-USP/ExperimentoPineal/INTEGRATION/")

##############################################################################
##########################DMRbyLOCATION USING GFM#############################
##############################################################################
#END###########END####################END###################END##############

library(UpSetR)
library(RColorBrewer)
### expression ###
tiff("VennBar_RNA.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
RNA_seq_list <- c()
RNA_seq_list$SMvF <- na.omit(unique(res_comp1$ensembl_id[res_comp1$P.Value<0.05]))
RNA_seq_list$CMvF <- na.omit(unique(res_comp2$ensembl_id[res_comp2$P.Value<0.05]))
RNA_seq_list$MSvC <- na.omit(unique(res_comp3$ensembl_id[res_comp3$P.Value<0.05]))
RNA_seq_list$FSvC <- na.omit(unique(res_comp4$ensembl_id[res_comp4$P.Value<0.05]))
RNA_seq_list$GFM <- na.omit(unique(res_comp5$ensembl_id[res_comp5$P.Value<0.05]))
upset(data = fromList(RNA_seq_list) , keep.order = T, sets= c("GFM", "FSvC", "MSvC", "CMvF","SMvF"),
      matrix.color = "blue", point.size = 5, text.scale= c(1,1, 1, 1, 1.5, 1.5),
      sets.bar.color = brewer.pal(5, "Set3"))
dev.off() 

tiff("VennBar_DMR.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
DMR_medips <- c()
DMR_medips$SMvF <- na.omit(unique(res_comp1_medips_base$loc_id[res_comp1_medips_base$P.Value<0.05]))
DMR_medips$CMvF <- na.omit(unique(res_comp2_medips_base$loc_id[res_comp2_medips_base$P.Value<0.05]))
DMR_medips$MSvC <- na.omit(unique(res_comp3_medips_base$loc_id[res_comp3_medips_base$P.Value<0.05]))
DMR_medips$FSvC <- na.omit(unique(res_comp4_medips_base$loc_id[res_comp4_medips_base$P.Value<0.05]))
DMR_medips$GFM <- na.omit(unique(res_comp5_medips_base$loc_id[res_comp5_medips_base$P.Value<0.05]))
upset(data = fromList(DMR_medips) , keep.order = T, sets= c("GFM", "FSvC", "MSvC", "CMvF","SMvF"),
      matrix.color = "blue", point.size = 5, text.scale= c(1,1, 1, 1, 1.5, 1.3),
      sets.bar.color = brewer.pal(5, "Set3"))
dev.off() 


tiff("VennBar_MIR.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
mir_seq <- c()
mir_seq$SMvF <- na.omit(unique(res_comp1_mirs$mirbase_id[res_comp1_mirs$P.Value<0.05]))
mir_seq$CMvF <- na.omit(unique(res_comp2_mirs$mirbase_id[res_comp2_mirs$P.Value<0.05]))
mir_seq$MSvC <- na.omit(unique(res_comp3_mirs$mirbase_id[res_comp3_mirs$P.Value<0.05]))
mir_seq$FSvC <- na.omit(unique(res_comp4_mirs$mirbase_id[res_comp4_mirs$P.Value<0.05]))
mir_seq$GFM <- na.omit(unique(res_comp5_mirs$mirbase_id[res_comp5_mirs$P.Value<0.05]))
upset(data = fromList(mir_seq) , keep.order = T, sets= c("GFM", "FSvC", "MSvC", "CMvF","SMvF"),
      matrix.color = "blue", point.size = 5, text.scale= c(1,1, 1, 1, 1.5, 1.5),
      sets.bar.color = brewer.pal(5, "Set3"))
dev.off() 
########################################################################################################
tiff("VennBar_SMvF.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
SMvF <- c()
SMvF$DEG <- na.omit(unique(res_comp1$ensembl_id[res_comp1$P.Value<0.05]))
SMvF$DMR <- as.character(na.omit(unique(res_comp1_medips_g$peaks2[res_comp1_medips_g$P.Value<0.05])))
SMvF$DMR_5k <- as.character(na.omit(unique(res_comp1_medips_g5$peaks2[res_comp1_medips_g5$P.Value<0.05])))
SMvF$DMR_10k <- as.character(na.omit(unique(res_comp1_medips_g10$peaks2[res_comp1_medips_g10$P.Value<0.05])))
SMvF$DMiR <- na.omit(unique(res_comp1_mirs$ensembl_id[res_comp1_mirs$P.Value<0.05 & res_comp1_mirs$pred_score > 90]))
upset(data = fromList(SMvF) , keep.order = T, sets= c("DMR", "DMR_5k", "DMR_10k", "DEG","DMiR"),
      matrix.color = "blue", point.size = 5, text.scale= c(1,1, 1, 1, 1.5, 1.5),
      sets.bar.color = brewer.pal(5, "Set3"))
dev.off() 

tiff("VennBar_CMvF.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
CMvF <- c()
CMvF$DEG <- na.omit(unique(res_comp2$ensembl_id[res_comp2$P.Value<0.05]))
CMvF$DMR <- as.character(na.omit(unique(res_comp2_medips_g$peaks2[res_comp2_medips_g$P.Value<0.05])))
CMvF$DMR_5k <- as.character(na.omit(unique(res_comp2_medips_g5$peaks2[res_comp2_medips_g5$P.Value<0.05])))
CMvF$DMR_10k <- as.character(na.omit(unique(res_comp2_medips_g10$peaks2[res_comp2_medips_g10$P.Value<0.05])))
CMvF$DMiR <- na.omit(unique(res_comp2_mirs$ensembl_id[res_comp2_mirs$P.Value<0.05 & res_comp2_mirs$pred_score > 90]))
upset(data = fromList(CMvF) , keep.order = TRUE, sets= c("DMR", "DMR_5k", "DMR_10k", "DEG","DMiR"),
      matrix.color = "blue", point.size = 5, text.scale= c(1,1, 1, 1, 1.5, 1.5),
      sets.bar.color = brewer.pal(5, "Set3"))
dev.off()

tiff("VennBar_MSvC.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
MSvC <- c()
MSvC$DEG <- na.omit(unique(res_comp3$ensembl_id[res_comp3$P.Value<0.05]))
MSvC$DMR <- as.character(na.omit(unique(res_comp3_medips_g$peaks2[res_comp3_medips_g$P.Value<0.05])))
MSvC$DMR_5k <- as.character(na.omit(unique(res_comp3_medips_g5$peaks2[res_comp3_medips_g5$P.Value<0.05])))
MSvC$DMR_10k <- as.character(na.omit(unique(res_comp3_medips_g10$peaks2[res_comp3_medips_g10$P.Value<0.05])))
MSvC$DMiR <- na.omit(unique(res_comp3_mirs$ensembl_id[res_comp3_mirs$P.Value<0.05 & res_comp3_mirs$pred_score > 90]))
upset(data = fromList(MSvC) , keep.order = T, sets= c("DMR", "DMR_5k", "DMR_10k", "DEG","DMiR"),
      matrix.color = "blue", point.size = 5, text.scale= c(1,1, 1, 1, 1.5, 1.5),
      sets.bar.color = brewer.pal(5, "Set3"))
dev.off()

tiff("VennBar_FSvC.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
FSvC <- c()
FSvC$DEG <- na.omit(unique(res_comp4$ensembl_id[res_comp4$P.Value<0.05]))
FSvC$DMR <- as.character(na.omit(unique(res_comp4_medips_g$peaks2[res_comp4_medips_g$P.Value<0.05])))
FSvC$DMR_5k <- as.character(na.omit(unique(res_comp4_medips_g5$peaks2[res_comp4_medips_g5$P.Value<0.05])))
FSvC$DMR_10k <- as.character(na.omit(unique(res_comp4_medips_g10$peaks2[res_comp4_medips_g10$P.Value<0.05])))
FSvC$DMiR <- na.omit(unique(res_comp4_mirs$ensembl_id[res_comp4_mirs$P.Value<0.05 & res_comp4_mirs$pred_score > 90]))
upset(data = fromList(FSvC) , keep.order = T, sets= c("DMR", "DMR_5k", "DMR_10k", "DEG","DMiR"),
      matrix.color = "blue", point.size = 5, text.scale= c(1,1, 1, 1, 1.5, 1.5),
      sets.bar.color = brewer.pal(5, "Set3"))
dev.off()

tiff("VennBar_GFM.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
GFM <- c()
GFM$DEG <- na.omit(unique(res_comp5$ensembl_id[res_comp5$P.Value<0.05]))
GFM$DMR <- as.character(na.omit(unique(res_comp5_medips_g$peaks2[res_comp5_medips_g$P.Value<0.05])))
GFM$DMR_5k <- as.character(na.omit(unique(res_comp5_medips_g5$peaks2[res_comp5_medips_g5$P.Value<0.05])))
GFM$DMR_10k <- as.character(na.omit(unique(res_comp5_medips_g10$peaks2[res_comp5_medips_g10$P.Value<0.05])))
GFM$DMiR <- na.omit(unique(res_comp5_mirs$ensembl_id[res_comp5_mirs$P.Value<0.05 & res_comp5_mirs$pred_score > 90]))
upset(data = fromList(GFM) , keep.order = T, sets= c("DMR", "DMR_5k", "DMR_10k", "DEG","DMiR"),
              matrix.color = "blue", point.size = 5, text.scale= c(1,1, 1, 1, 1.5, 1.5),
              sets.bar.color = brewer.pal(5, "Set3"))
dev.off()

####################################################################
##make a comparison list
library(RVenn)
library(purrr)
library(ggplot2)
library (ggvenn)

SMvF <- c()
SMvF$DEG <- na.omit(unique(res_comp1$ensembl_id[res_comp1$P.Value<0.05]))
SMvF$DMR <- as.character(na.omit(unique(res_comp1_medips_g10$peaks2[res_comp1_medips_g10$P.Value<0.05])))
SMvF$DMiR <- na.omit(unique(res_comp1_mirs$ensembl_id[res_comp1_mirs$P.Value<0.05 & res_comp1_mirs$pred_score > 90]))
tiff("Venn_SMvF.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
ggvenn(SMvF)+ scale_fill_manual(values = c("red", "blue", "yellow"))
dev.off()

CMvF <- c()
CMvF$DEG <- na.omit(unique(res_comp2$ensembl_id[res_comp2$P.Value<0.05]))
CMvF$DMR <- as.character(na.omit(unique(res_comp2_medips_g10$peaks2[res_comp2_medips_g10$P.Value<0.05])))
CMvF$DMiR <- na.omit(unique(res_comp2_mirs$ensembl_id[res_comp2_mirs$P.Value<0.05 & res_comp2_mirs$pred_score > 90]))
tiff("Venn_CMvF.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
ggvenn(CMvF)+ scale_fill_manual(values = c("red", "blue", "yellow"))
dev.off()

MSvC <- c()
MSvC$DEG <- na.omit(unique(res_comp3$ensembl_id[res_comp3$P.Value<0.05]))
MSvC$DMR <- as.character(na.omit(unique(res_comp3_medips_g10$peaks2[res_comp3_medips_g10$P.Value<0.05])))
MSvC$DMiR <- na.omit(unique(res_comp3_mirs$ensembl_id[res_comp3_mirs$P.Value<0.05 & res_comp3_mirs$pred_score > 90]))
tiff("Venn_MSvC.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
ggvenn(MSvC) + scale_fill_manual(values = c("red", "blue", "yellow"))
dev.off()

FSvC <- c()
FSvC$DEG <- na.omit(unique(res_comp4$ensembl_id[res_comp4$P.Value<0.05]))
FSvC$DMR <- as.character(na.omit(unique(res_comp4_medips_g10$peaks2[res_comp4_medips_g10$P.Value<0.05])))
FSvC$DMiR <- na.omit(unique(res_comp4_mirs$ensembl_id[res_comp4_mirs$P.Value<0.05 & res_comp4_mirs$pred_score > 90]))
tiff("Venn_FSvC.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
ggvenn(FSvC)+ scale_fill_manual(values = c("red", "blue", "yellow"))
dev.off()

GFM <- c()
GFM$DEG <- na.omit(unique(res_comp5$ensembl_id[res_comp5$P.Value<0.05]))
GFM$DMR <- as.character(na.omit(unique(res_comp5_medips_g10$peaks2[res_comp5_medips_g10$P.Value<0.05])))
GFM$DMiR <- na.omit(unique(res_comp5_mirs$ensembl_id[res_comp5_mirs$P.Value<0.05 & res_comp5_mirs$pred_score > 90]))
tiff("Venn_GFM.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
ggvenn(GFM) + scale_fill_manual(values = c("red", "blue", "yellow"))
dev.off()


##############################with intersection##########################
pathway_list_overlaps <- c()
pathway_list_overlaps$RNA_MiR <- annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DEG[GFM$DEG %in% GFM$DMiR]]
pathway_list_overlaps$RNA_DMR <- annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DEG[GFM$DEG %in% GFM$DMR]]
pathway_list_overlaps$DMR_MiR <- annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DMR[GFM$DMR %in% GFM$DMiR]]


library(clusterProfiler)
library(org.Gg.eg.db)
go_enrich <- compareCluster(geneClusters = pathway_list_overlaps , 
                            fun = "enrichGO" , 
                            OrgDb = org.Gg.eg.db , 
                            ont="BP" , 
                            pvalueCutoff=1, 
                            qvalueCutoff=0.5, 
                            readable=T)
tiff("PathwayGO_IntersectionGFM.tiff", units="cm", width=35.00, height=25.00, res=600, compression = "lzw")
dotplot(go_enrich , showCategory=15)
dev.off()

go_g_enrich <- compareCluster(geneClusters = pathway_list_overlaps , 
                              fun = "groupGO" , 
                              OrgDb = org.Gg.eg.db , 
                              ont="BP" , 
                              level = 10 , 
                              readable=T)
tiff("PathwaygroupsGO_IntersectionGFM.tiff", units="cm", width=35.00, height=25.00, res=600, compression = "lzw")
dotplot(go_g_enrich)
dev.off()

kegg_enrich <- compareCluster(geneClusters = pathway_list_overlaps , 
                              fun = "enrichKEGG" , 
                              org="gga",
                              pvalueCutoff=1, 
                              qvalueCutoff=1)
tiff("PathwaygKEGG_IntersectionGFM.tiff", units="cm", width=35.00, height=25.00, res=600, compression = "lzw")
dotplot(kegg_enrich)
dev.off()

##############################UNIQUE##########################
`%notin%` <- Negate(`%in%`)

pathway_list_overlaps <- c()
pathway_list_overlaps$RNA <- na.omit(unique(annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DEG[GFM$DEG %notin% GFM$DMiR & GFM$DEG %notin% GFM$DMR]]))
pathway_list_overlaps$DMR <- na.omit(unique(annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DMR[GFM$DMR %notin% GFM$DMiR & GFM$DMR %notin% GFM$DEG]]))
pathway_list_overlaps$MiR <- na.omit(unique(annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DMiR[GFM$DMiR %notin% GFM$DEG & GFM$DMiR %notin% GFM$DMR]]))

library(clusterProfiler)
library(org.Gg.eg.db)
go_enrich <- compareCluster(geneClusters = pathway_list_overlaps , 
                            fun = "enrichGO" , 
                            OrgDb = org.Gg.eg.db , 
                            ont="BP" , 
                            pvalueCutoff=1, 
                            qvalueCutoff=0.5, 
                            readable=T)
tiff("PathwayGO_UniqueGFM.tiff", units="cm", width=35.00, height=25.00, res=600, compression = "lzw")
dotplot(go_enrich , showCategory=15)
dev.off()

go_g_enrich <- compareCluster(geneClusters = pathway_list_overlaps , 
                              fun = "groupGO" , 
                              OrgDb = org.Gg.eg.db , 
                              ont="BP" , 
                              level = 10 , 
                              readable=T)
tiff("PathwaygroupsGO_UniqueGFM.tiff", units="cm", width=25.00, height=15.00, res=600, compression = "lzw")
dotplot(go_g_enrich)
dev.off()

kegg_enrich <- compareCluster(geneClusters = pathway_list_overlaps , 
                              fun = "enrichKEGG" , 
                              org="gga",
                              pvalueCutoff=1, 
                              qvalueCutoff=0.8)
tiff("PathwaygKEGG_UniqueGFM.tiff", units="cm", width=25.00, height=15.00, res=600, compression = "lzw")
dotplot(kegg_enrich)
dev.off()

##############################BOTH##########################
pathway_list_overlaps <- c()
pathway_list_overlaps$RNA <- na.omit(unique(annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DEG[GFM$DEG %notin% GFM$DMiR & GFM$DEG %notin% GFM$DMR]]))
pathway_list_overlaps$DMR <- na.omit(unique(annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DMR[GFM$DMR %notin% GFM$DMiR & GFM$DMR %notin% GFM$DEG]]))
pathway_list_overlaps$MiR <- na.omit(unique(annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DMiR[GFM$DMiR %notin% GFM$DEG & GFM$DMiR %notin% GFM$DMR]]))
pathway_list_overlaps$RNA_MiR <- annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DEG[GFM$DEG %in% GFM$DMiR]]
pathway_list_overlaps$RNA_DMR <- annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DEG[GFM$DEG %in% GFM$DMR]]
pathway_list_overlaps$DMR_MiR <- annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DMR[GFM$DMR %in% GFM$DMiR]]

go_enrich <- compareCluster(geneClusters = pathway_list_overlaps , 
                            fun = "enrichGO" , 
                            OrgDb = org.Gg.eg.db , 
                            ont="BP" , 
                            pvalueCutoff=1, 
                            qvalueCutoff=0.9, 
                            readable=T)
tiff("PathwayGO_GFM.tiff", units="cm", width=35.00, height=40.00, res=600, compression = "lzw")
dotplot(go_enrich , showCategory=15)
dev.off()

go_g_enrich <- compareCluster(geneClusters = pathway_list_overlaps , 
                              fun = "groupGO" , 
                              OrgDb = org.Gg.eg.db , 
                              ont="BP" , 
                              level = 10 , 
                              readable=T)
tiff("PathwaygroupsGO_GFM.tiff", units="cm", width=30.00, height=15.00, res=600, compression = "lzw")
dotplot(go_g_enrich)
dev.off()

kegg_enrich <- compareCluster(geneClusters = pathway_list_overlaps , 
                              fun = "enrichKEGG" , 
                              org="gga",
                              pvalueCutoff=1, 
                              qvalueCutoff=0.8)
tiff("PathwaygKEGG_GFM.tiff", units="cm", width=30.00, height=20.00, res=600, compression = "lzw")
dotplot(kegg_enrich)
dev.off()


#################################################################################################
#################################################################################################
#################################################################################################
#######################TO EXTRACT THE INFORMATION FROM THE INTERSECTIONS#########################
#################################################################################################
#################################################################################################
#################################################################################################
pathway_list_overlaps <- c()
pathway_list_overlaps$RNA_MiR <- annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DEG[GFM$DEG %in% GFM$DMiR]]
pathway_list_overlaps$RNA_DMR <- annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DEG[GFM$DEG %in% GFM$DMR]]
pathway_list_overlaps$DMR_MiR <- annot_exprs$gene_id[annot_exprs$ensembl_id %in% GFM$DMR[GFM$DMR %in% GFM$DMiR]]

d1 <- res_comp5[res_comp5$P.Value<0.05 & res_comp5$gene_id %in% pathway_list_overlaps$RNA_MiR, ]
d2 <- res_comp5[res_comp5$P.Value<0.05 & res_comp5$gene_id %in% pathway_list_overlaps$DMR_MiR, ]
d3 <- res_comp5[res_comp5$P.Value<0.05 & res_comp5$gene_id %in% pathway_list_overlaps$RNA_DMR, ]
d1$Intersection <- "RNA_MiR"
d2$Intersection <- "DMR_MiR"
d3$Intersection <- "RNA_DMR"
RNAinter<- rbind(d1, d3)

d1 <- res_comp5_mirs[res_comp5_mirs$P.Value<0.05 & res_comp5_mirs$gene_id %in% pathway_list_overlaps$RNA_MiR, ]
d2 <- res_comp5_mirs[res_comp5_mirs$P.Value<0.05 & res_comp5_mirs$gene_id %in% pathway_list_overlaps$DMR_MiR, ]
#d3 <- res_comp5_mirs[res_comp5_mirs$P.Value<0.05 & res_comp5_mirs$gene_id %in% pathway_list_overlaps$RNA_DMR, ]
d1$Intersection <- "RNA_MiR"
d2$Intersection <- "DMR_MiR"
#d3$Intersection <- "RNA_DMR"
DMiRinter<- rbind(d1, d2)
DMiRinter<- na.omit(DMiRinter)
DMiRinter <- merge(DMiRinter, annot_exprs[, c("symbol", "ensembl_id")], by="ensembl_id")
DMiRinter <- distinct(DMiRinter, symbol,  .keep_all = TRUE)

d1 <- res_comp5[res_comp5$P.Value<0.05 & res_comp5$gene_id %in% pathway_list_overlaps$RNA_MiR, ]
d2 <- res_comp5[res_comp5$P.Value<0.05 & res_comp5$gene_id %in% pathway_list_overlaps$DMR_MiR, ]
d3 <- res_comp5[res_comp5$P.Value<0.05 & res_comp5$gene_id %in% pathway_list_overlaps$RNA_DMR, ]
d1$Intersection <- "RNA_MiR"
#d2$Intersection <- "DMR_MiR"
d3$Intersection <- "RNA_DMR"
DEGinter<- rbind(d1, d3)
DEGinter<- na.omit(DEGinter)

pathway_list_overlaps <- c()
pathway_list_overlaps$RNA_MiR <- annot_exprs$ensembl_id[annot_exprs$ensembl_id %in% GFM$DEG[GFM$DEG %in% GFM$DMiR]]
pathway_list_overlaps$RNA_DMR <- annot_exprs$ensembl_id[annot_exprs$ensembl_id %in% GFM$DEG[GFM$DEG %in% GFM$DMR]]
pathway_list_overlaps$DMR_MiR <- annot_exprs$ensembl_id[annot_exprs$ensembl_id %in% GFM$DMR[GFM$DMR %in% GFM$DMiR]]

#d1 <- res_comp5_medips_g10[res_comp5_medips_g10$P.Value<0.05 & res_comp5_medips_g10$peaks2 %in% pathway_list_overlaps$RNA_MiR, ]
d2 <- res_comp5_medips_g10[res_comp5_medips_g10$P.Value<0.05 & res_comp5_medips_g10$peaks2 %in% pathway_list_overlaps$DMR_MiR, ]
d3 <- res_comp5_medips_g10[res_comp5_medips_g10$P.Value<0.05 & res_comp5_medips_g10$peaks2 %in% pathway_list_overlaps$RNA_DMR, ]
#d1$Intersection <- "RNA_MiR"
d2$Intersection <- "DMR_MiR"
d3$Intersection <- "RNA_DMR"
DMRinter<- rbind(d2, d3)
DMRinter<- na.omit(DMRinter)
DMRinter$ensembl_id <- DMRinter$peaks2
DMRinter <- merge(DMRinter, annot_exprs[, c("symbol", "ensembl_id")], by="ensembl_id")



write.table(RNAinter , "DEGinter.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
write.table(DMiRinter , "DMiRinter.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
write.table(DMRinter , "DMRinter.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)

#################################################################################################
#################################################################################################
#################################################################################################
##########################################MODULE#################################################
#################################################################################################
#################################################################################################
#################################################################################################

#ppi_net2 <- Gallus_network_PPI_700[Gallus_network_PPI_700$Score >= 900 ,]
##############make_modules#################
#Sys.setenv(PATH = paste("/Users/badt/anaconda3/bin", Sys.getenv("PATH"), sep=":"))
library(MODifieR)
library(plyr)
library(data.table)
#set.seed(143)

####RNA#######
input_tab <- as.data.frame(res_comp5)
input_tab <- input_tab[, c("gene_id" , "P.Value")]
colnames(input_tab) <- c("gene" , "pvalue")
input_tab <- na.omit(input_tab)
input_tab$gene <- as.character(input_tab$gene)
input_tab <- ddply(input_tab , .(gene) , numcolwise(min))
input_tab <- input_tab[order(input_tab$pvalue , decreasing = F) , ]
MODifieR_input <- MODifieR:::create_custom_input_object(diff_genes = input_tab )
diamond_module_rna_gfm <- diamond(MODifieR_input = MODifieR_input, 
                                  ppi_network = Gallus_network_PPI_700 , include_seed = T)

####MEDIPS#######
input_tab <- as.data.frame(res_comp5_medips_g10)
input_tab <- input_tab[, c("peaks2" , "P.Value")]
colnames(input_tab) <- c("ensembl_id" , "P.Value")
input_tab$ensembl_id <- as.character(input_tab$ensembl_id)
input_tab <- na.omit(input_tab)
input_tab <- merge(input_tab , annot_exprs , all.x=T)
input_tab <- input_tab[, c("gene_id" , "P.Value")]
colnames(input_tab) <- c("gene" , "pvalue")
input_tab <- na.omit(input_tab)
input_tab$gene <- as.character(input_tab$gene)
input_tab <- ddply(input_tab , .(gene) , numcolwise(min))
input_tab <- input_tab[order(input_tab$pvalue , decreasing = F) , ]
MODifieR_input <- MODifieR:::create_custom_input_object(diff_genes = input_tab )
diamond_module_dna_gfm <- diamond(MODifieR_input = MODifieR_input, 
                                  ppi_network = Gallus_network_PPI_700 , include_seed = T)


####MIRS#######
input_tab <- as.data.frame(res_comp5_mirs)
input_tab <- input_tab[, c("gene_id" , "P.Value")]
colnames(input_tab) <- c("gene" , "pvalue")
input_tab <- na.omit(input_tab)
input_tab$gene <- as.character(input_tab$gene)
input_tab <- ddply(input_tab , .(gene) , numcolwise(min))
input_tab <- input_tab[order(input_tab$pvalue , decreasing = F) , ]
MODifieR_input <- MODifieR:::create_custom_input_object(diff_genes = input_tab )
diamond_module_mir_gfm <- diamond(MODifieR_input = MODifieR_input, 
                                  ppi_network = Gallus_network_PPI_700 , include_seed = T)

#############################################ATENTION#############################################
#####for GFM_mod obeject we will use the loaded object in ScriptInputsR because something is happening different across windows and Mac ###
GFM_mod <- readRDS(file="GFM_mod.rds", refhook = NULL)######################################
names(GFM_mod ) <- c("DEG", "DMR", "DMiR")
##################################################################################################
#GFM_mod <- c()
#GFM_mod$DEG <- diamond_module_rna_gfm$module_genes
#GFM_mod$Medips <- diamond_module_dna_gfm$module_genes
#GFM_mod$mirs <- diamond_module_mir_gfm$module_genes

tiff("Venn_GFM_mod.tiff", units="cm", width=15.00, height=15.00, res=600, compression = "lzw")
ggvenn(GFM_mod) + scale_fill_manual(values = c("red", "blue", "yellow"))
dev.off()

##############################with intersection##########################
pathway_list_overlaps_mod <- c()
pathway_list_overlaps_mod$RNA_MiR <- GFM_mod$DEG[GFM_mod$DEG %in% GFM_mod$DMiR]
pathway_list_overlaps_mod$RNA_DMR <- GFM_mod$DEG[GFM_mod$DEG %in% GFM_mod$DMR ]
pathway_list_overlaps_mod$DMR_DMiR <- GFM_mod$DMR[GFM_mod$DMR %in% GFM_mod$DMiR ]

library(clusterProfiler)
library(org.Gg.eg.db)
go_mod_enrich <- compareCluster(geneClusters = pathway_list_overlaps_mod , 
                                fun = "enrichGO" , 
                                OrgDb = org.Gg.eg.db , 
                                ont="BP" , 
                                pvalueCutoff=1, 
                                qvalueCutoff=0.5, 
                                readable=T)
tiff("PathwayGO_IntersectionGFM_mod.tiff", units="cm", width=30.00, height=30.00, res=600, compression = "lzw")
dotplot(go_mod_enrich , showCategory=20)
dev.off()



go_g_mod_enrich <- compareCluster(geneClusters = pathway_list_overlaps_mod , 
                                  fun = "groupGO" , 
                                  OrgDb = org.Gg.eg.db , 
                                  ont="BP" , 
                                  level = 10 , 
                                  readable=T)
tiff("PathwaygroupsGO_IntersectionGFM_mod.tiff", units="cm", width=30.00, height=20.00, res=600, compression = "lzw")
dotplot(go_g_mod_enrich)
dev.off()

kegg_mod_enrich <- compareCluster(geneClusters = pathway_list_overlaps_mod , 
                                  fun = "enrichKEGG" , 
                                  org="gga",
                                  pvalueCutoff=1, 
                                  qvalueCutoff=1)
tiff("PathwaygKEGG_IntersectionGFM_mod.tiff", units="cm", width=25.00, height=15.00, res=600, compression = "lzw")
dotplot(kegg_mod_enrich)
dev.off()

##############################UNIQUE##########################

pathway_list_overlaps_mod <- c()
pathway_list_overlaps_mod$RNA <- GFM_mod$DEG[GFM_mod$DEG %notin% GFM_mod$DMiR & GFM_mod$DEG %notin% GFM_mod$DMR]
pathway_list_overlaps_mod$DMR <- GFM_mod$DMR[GFM_mod$DMR %notin% GFM_mod$DMiR & GFM_mod$DMR %notin% GFM_mod$DEG]
pathway_list_overlaps_mod$DMiR <- GFM_mod$DMiR[GFM_mod$DMiR %notin% GFM_mod$DEG & GFM_mod$DMiR %notin% GFM_mod$DMR]

go_mod_enrich <- compareCluster(geneClusters = pathway_list_overlaps_mod , 
                                fun = "enrichGO" , 
                                OrgDb = org.Gg.eg.db , 
                                ont="BP" , 
                                pvalueCutoff=1, 
                                qvalueCutoff=0.5, 
                                readable=T)
tiff("PathwayGO_UniqueGFM_mod.tiff", units="cm", width=30.00, height=30.00, res=600, compression = "lzw")
dotplot(go_mod_enrich , showCategory=20)
dev.off()



go_g_mod_enrich <- compareCluster(geneClusters = pathway_list_overlaps_mod , 
                                  fun = "groupGO" , 
                                  OrgDb = org.Gg.eg.db , 
                                  ont="BP" , 
                                  level = 10 , 
                                  readable=T)
tiff("PathwaygroupsGO_UniqueGFM_mod.tiff", units="cm", width=30.00, height=20.00, res=600, compression = "lzw")
dotplot(go_g_mod_enrich)
dev.off()

kegg_mod_enrich <- compareCluster(geneClusters = pathway_list_overlaps_mod , 
                                  fun = "enrichKEGG" , 
                                  org="gga",
                                  pvalueCutoff=1, 
                                  qvalueCutoff=1)
tiff("PathwaygKEGG_UniqueGFM_mod.tiff", units="cm", width=25.00, height=15.00, res=600, compression = "lzw")
dotplot(kegg_mod_enrich)
dev.off()

##############################BOTH##########################
pathway_list_overlaps_mod <- c()
pathway_list_overlaps_mod$RNA <- GFM_mod$DEG[GFM_mod$DEG %notin% GFM_mod$DMiR & GFM_mod$DEG %notin% GFM_mod$DMR]
pathway_list_overlaps_mod$DMR <- GFM_mod$DMR[GFM_mod$DMR %notin% GFM_mod$DMiR & GFM_mod$DMR %notin% GFM_mod$DEG]
pathway_list_overlaps_mod$DMiR <- GFM_mod$DMiR[GFM_mod$DMiR %notin% GFM_mod$DEG & GFM_mod$DMiR %notin% GFM_mod$DMR]
pathway_list_overlaps_mod$RNA_DMiR <- GFM_mod$DEG[GFM_mod$DEG %in% GFM_mod$DMiR ]
pathway_list_overlaps_mod$RNA_DMR <- GFM_mod$DEG[GFM_mod$DEG %in% GFM_mod$DMR ]
pathway_list_overlaps_mod$DMR_DMiR <- GFM_mod$DMR[GFM_mod$DMR %in% GFM_mod$DMiR ]

go_mod_enrich <- compareCluster(geneClusters = pathway_list_overlaps_mod , 
                                fun = "enrichGO" , 
                                OrgDb = org.Gg.eg.db , 
                                ont="BP" , 
                                pvalueCutoff=1, 
                                qvalueCutoff=0.9, 
                                readable=T)
tiff("PathwayGO_GFM_mod.tiff", units="cm", width=40.00, height=50.00, res=600, compression = "lzw")
dotplot(go_mod_enrich , showCategory=20)
dev.off()



go_g_mod_enrich <- compareCluster(geneClusters = pathway_list_overlaps_mod , 
                                  fun = "groupGO" , 
                                  OrgDb = org.Gg.eg.db , 
                                  ont="BP" , 
                                  level = 10 , 
                                  readable=T)
tiff("PathwaygroupsGO_GFM_mod.tiff", units="cm", width=35.00, height=20.00, res=600, compression = "lzw")
dotplot(go_g_mod_enrich)
dev.off()

kegg_mod_enrich <- compareCluster(geneClusters = pathway_list_overlaps_mod , 
                                  fun = "enrichKEGG" , 
                                  org="gga",
                                  pvalueCutoff=1, 
                                  qvalueCutoff=1)
tiff("PathwaygKEGG_GFM_mod.tiff", units="cm", width=30.00, height=25.00, res=600, compression = "lzw")
dotplot(kegg_mod_enrich)
dev.off()


#####Writing pathway data.frames#####
write.table(as.data.frame(go_enrich@compareClusterResult) , "PathwayGO_Enrich.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
write.table(as.data.frame(kegg_enrich@compareClusterResult) , "PathwayKegg_Enrich.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
write.table(as.data.frame(go_mod_enrich@compareClusterResult) , "PathwayGO_Mod_Enrich.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)
write.table(as.data.frame(kegg_mod_enrich@compareClusterResult) , "PathwayKegg_Mod_Enrich.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)

#transform KEGG annotation geneID in gene simboll
KEGG_LIST <- as.data.frame(kegg_mod_enrich@compareClusterResult)
collnames <- c("Cluster", "geneID", "p.adjust")
KEGG_LIST <- KEGG_LIST[collnames]
KEGG_LIST <- KEGG_LIST %>% filter(p.adjust <= 0.2)
KEGG_LIST <- separate_rows(KEGG_LIST, geneID)
colnames(KEGG_LIST) <- c("cluster", "gene_id")
KEGG_LIST <- merge(KEGG_LIST , xs , all.x=T)
write.table(as.data.frame(KEGG_LIST) , "PathwayKegg_Mod_Enrich_gene_ID.txt" , sep="\t" , row.names = T, quote = F, col.names = NA)

###extracting the last gene list to use as input in CytoScape

all_genes_list <- as.data.frame(table(c(GFM_mod$DEG , GFM_mod$DMR , GFM_mod$DMiR)))
#all_genes_list <- as.data.frame(table(c(diamond_module_dna_gfm$module_genes , diamond_module_mir_gfm$module_genes , diamond_module_rna_gfm$module_genes)))
colnames(all_genes_list) <- c("gene_id" , "frequency")
xs <- as.data.frame(org.Gg.egSYMBOL2EG)
all_genes_list <- merge(all_genes_list , xs , all.x=T)
all_genes_list$DMR <- c(0)
all_genes_list$DMR[all_genes_list$gene_id %in% diamond_module_dna_gfm$seed_genes] <- 1
all_genes_list$DEG <- c(0)
all_genes_list$DEG[all_genes_list$gene_id %in% diamond_module_rna_gfm$seed_genes] <- 1
all_genes_list$DMiR <- c(0)
all_genes_list$DMiR[all_genes_list$gene_id %in% diamond_module_mir_gfm$seed_genes] <- 1
all_genes_list[,4:6] <- all_genes_list[,4:6]/rowSums(all_genes_list[,4:6])
all_genes_list[is.na(all_genes_list)] <- 0
write.table(all_genes_list , "gene_list_annotation_cytoscape.txt" , sep="\t" , row.names = F , quote = F)


#############https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
cnetplot(res_comp5, showCategory = 5, foldChange = res_comp5$logFC)

cnetplot(ego, foldChange=geneList)




#################################################################################################
#################################################################################################
#################################################################################################
##########################################FOLDCHANGE#############################################
#################################################################################################
#################################################################################################
#################################################################################################
DMRinter
DMiRinter
DEGinter

na.omit
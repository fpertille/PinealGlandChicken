rm(list=ls())
##Medips_analysis##
library(readr)
library(RMariaDB)
library(GenomicFeatures)
library(org.Gg.eg.db)
library(ChIPseeker)
library(biomaRt)

xf <- as.data.frame(org.Gg.egENSEMBL2EG)
mart <- useMart(biomart = "ensembl", dataset = "ggallus_gene_ensembl")
gene_ranges <- getBM(attributes = c( "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
                     filters    = "ensembl_gene_id",
                     values     = xf$ensembl_id, 
                     mart       = mart)

Medips_windows <- read_delim("mr.edgeR_PINEAL_allind.txt","\t", escape_double = FALSE, trim_ws = TRUE)
count_medips <- unique(Medips_windows[, 1:16])
count_medips <- as.data.frame(count_medips)
rownames(count_medips) <- paste(count_medips$chr , count_medips$start , count_medips$stop , sep="_")
count_medips <- count_medips[,-1:-4]


#plot(density(count_medips))
any(is.na(count_medips))
plot(density(as.matrix(count_medips)))
plot(density(log(as.matrix(count_medips))))


count_medips <- as.matrix(count_medips)
threshold <- count_medips > 1
table(rowSums(threshold))
keep_medips <- rowSums(threshold) >= 4
count_medips_keep <- count_medips[keep_medips ,]


plot(density(log(as.matrix(count_medips_keep))))
plot(density(as.matrix(count_medips_keep)))



transfor_lcpm <- voom(count_medips_keep)
colnames(transfor_lcpm) <- substring(colnames(transfor_lcpm) , 1,10)
plot(density(transfor_lcpm$E))

medips_pheno <- read_delim("sample_info copy 2.txt","\t", escape_double = FALSE, trim_ws = TRUE)
medips_pheno <- medips_pheno[order(match(medips_pheno$Sample_ID,colnames(transfor_lcpm))),]

group <- as.factor(medips_pheno$Gender)
group2 <- as.factor(medips_pheno$Stress)
group3 <- as.factor(paste(medips_pheno$Gender , medips_pheno$Stress , sep = "_"))

design <- model.matrix(~0+group3)
colnames(design) <- gsub("group3", "", colnames(design))
contr.matrix <- makeContrasts(
comp1= M_Stress-F_Stress ,
comp2= M_Control-F_Control ,
comp3= M_Stress-M_Control ,
comp4= F_Stress-F_Control ,
levels = colnames(design))
vfit_medips <- lmFit(transfor_lcpm, design)
vfit_medips <- contrasts.fit(vfit_medips, contrasts=contr.matrix)
efit_medips <- eBayes(vfit_medips)
res_comp1_medips_base <- topTable(efit_medips , coef = 1 , number = Inf)
res_comp1_medips_base$loc_id <- rownames(res_comp1_medips_base)
res_comp2_medips_base <- topTable(efit_medips , coef = 2 , number = Inf)
res_comp2_medips_base$loc_id <- rownames(res_comp2_medips_base)
res_comp3_medips_base <- topTable(efit_medips , coef = 3 , number = Inf)
res_comp3_medips_base$loc_id <- rownames(res_comp3_medips_base)
res_comp4_medips_base <- topTable(efit_medips , coef = 4 , number = Inf)
res_comp4_medips_base$loc_id <- rownames(res_comp4_medips_base)

design <- model.matrix(~0+group2)
colnames(design) <- gsub("group2", "", colnames(design))
design
contr.matrix <- makeContrasts(
comp1= Stress-Control ,
levels = colnames(design))
contr.matrix
vfit_medips <- lmFit(transfor_lcpm, design)
vfit_medips <- contrasts.fit(vfit_medips, contrasts=contr.matrix)
efit_medips <- eBayes(vfit_medips)
res_comp5_medips_base <- topTable(efit_medips , coef = 1 , number = Inf)
res_comp5_medips_base$loc_id <- rownames(res_comp5_medips_base)

var_id_func <- function(gse_info , column_name , splitchar){
df <- gse_info # split the strings
temp <- strsplit(df[,paste(column_name)], split=paste(splitchar))# maximum length of the list items
maxL <- max(sapply(temp, length))
temp <- data.frame(do.call(rbind, lapply(temp, function(i) c(i, rep(NA, maxL-length(i))))))# contstruct data.frame with NAs as fills
df <- cbind(x=df, temp)
gse_info <- df
rm(temp)
rm(df)
colnames(gse_info) <- gsub("x.", "", colnames(gse_info))
gse_info <- data.frame(lapply(gse_info, as.character), stringsAsFactors=FALSE)
return(gse_info)
}


peak_file <- as.data.frame(count_medips_keep)
peak_file$loc <- rownames(peak_file)
peak_file <- var_id_func(gse_info = peak_file , column_name = "loc" , splitchar = "_")
peak_file <- peak_file[,c("X1" , "X2" , "X3")]
colnames(peak_file) <- c("CHR" , "BP" , "BP2")
supportedUCSCtables(genome = "galGal6")
gg_txdb <- makeTxDbFromUCSC(genome="galGal6", tablename="ensGene")
gg_txdb <- readRDS("txdb_object.RDS")

peak_file_meth <- peak_file
rownames(peak_file_meth) <- paste(peak_file_meth$CHR , peak_file_meth$BP , peak_file_meth$BP2 , sep = "_")
library(GenomicRanges)
library(ChIPpeakAnno)
bed1 <- makeGRangesFromDataFrame(peak_file_meth , start.field = "BP" , end.field = "BP2")
bed2 <- makeGRangesFromDataFrame(gene_ranges , start.field = "start_position" , end.field = "end_position")
ovp = findOverlapsOfPeaks(bed1 , bed2)
overlaps <- ovp$overlappingPeaks$`bed1///bed2`

gene_ranges_5k <- gene_ranges
gene_ranges_5k$start_position <- gene_ranges_5k$start_position-5000
gene_ranges_5k$end_position <- gene_ranges_5k$end_position+5000
bed3_5k <- makeGRangesFromDataFrame(gene_ranges_5k , start.field = "start_position" , end.field = "end_position")
ovp_5k <- findOverlapsOfPeaks(bed1 , bed3_5k)
overlaps_5k <- ovp_5k$overlappingPeaks$`bed1///bed3_5k`


gene_ranges_10k <- gene_ranges
gene_ranges_10k$end_position <- gene_ranges_10k$end_position+10000
gene_ranges_10k$start_position <- gene_ranges_10k$start_position-10000
bed4_10k <- makeGRangesFromDataFrame(gene_ranges_10k , start.field = "start_position" , end.field = "end_position")
ovp_10k <- findOverlapsOfPeaks(bed1 , bed4_10k)
overlaps_10k <- ovp_10k$overlappingPeaks$`bed1///bed4_10k`

###merging with gene###
overlaps_main_info <- overlaps[,c("peaks1" , "peaks2" , "overlapFeature" , "shortestDistance")]
res_comp1_medips_base$peaks1 <- rownames(res_comp1_medips_base)
res_comp1_medips_g <- merge(res_comp1_medips_base , overlaps_main_info , all.x=T)
res_comp2_medips_base$peaks1 <- rownames(res_comp2_medips_base)
res_comp2_medips_g <- merge(res_comp2_medips_base , overlaps_main_info , all.x=T)
res_comp3_medips_base$peaks1 <- rownames(res_comp3_medips_base)
res_comp3_medips_g <- merge(res_comp3_medips_base , overlaps_main_info , all.x=T)
res_comp4_medips_base$peaks1 <- rownames(res_comp4_medips_base)
res_comp4_medips_g <- merge(res_comp4_medips_base , overlaps_main_info , all.x=T)
res_comp5_medips_base$peaks1 <- rownames(res_comp5_medips_base)
res_comp5_medips_g <- merge(res_comp5_medips_base , overlaps_main_info , all.x=T)

###merging with gene_5k###
overlaps_5k_main_info <- overlaps_5k[,c("peaks1" , "peaks2" , "overlapFeature" , "shortestDistance")]
res_comp1_medips_g5 <- merge(res_comp1_medips_base , overlaps_5k_main_info , all.x=T)
res_comp2_medips_g5 <- merge(res_comp2_medips_base , overlaps_5k_main_info , all.x=T)
res_comp3_medips_g5 <- merge(res_comp3_medips_base , overlaps_5k_main_info , all.x=T)
res_comp4_medips_g5 <- merge(res_comp4_medips_base , overlaps_5k_main_info , all.x=T)
res_comp5_medips_g5 <- merge(res_comp5_medips_base , overlaps_5k_main_info , all.x=T)

###merging with gene_10k###
overlaps_10k_main_info <- overlaps_10k[,c("peaks1" , "peaks2" , "overlapFeature" , "shortestDistance")]
res_comp1_medips_g10 <- merge(res_comp1_medips_base , overlaps_10k_main_info , all.x=T)
res_comp2_medips_g10 <- merge(res_comp2_medips_base , overlaps_10k_main_info , all.x=T)
res_comp3_medips_g10 <- merge(res_comp3_medips_base , overlaps_10k_main_info , all.x=T)
res_comp4_medips_g10 <- merge(res_comp4_medips_base , overlaps_10k_main_info , all.x=T)
res_comp5_medips_g10 <- merge(res_comp5_medips_base , overlaps_10k_main_info , all.x=T)












txdb <- gg_txdb
peak_anno_function <- function(input_data){
b2 <- input_data[,c("CHR" , "BP" , "BP2")]
write.table(b2 , "annotate_file.txt" , sep = "\t" , col.names = TRUE , row.names = FALSE)
peakAnno <- annotatePeak("/Users/badt/Documents/documents/Phd_thesis/Chicken/annotate_file.txt",
TxDb=txdb, annoDb="org.Gg.eg.db", tssRegion = c(-5000, 5000) , overlap = "all")
return(peakAnno)
}

peak_anno_function <- function(input_data){
  b2 <- input_data[,c("CHR" , "BP" , "BP2")]
  write.table(b2 , "annotate_file.txt" , sep = "\t" , col.names = TRUE , row.names = FALSE)
  peakAnno <- annotatePeak("annotate_file.txt",
                           TxDb=txdb, annoDb="org.Gg.eg.db", tssRegion = c(-10000, 10000) , overlap = "all")
  return(peakAnno)
}
#peak_chicken <- peak_anno_function(peak_file)


peak_chicken_5k <- peak_anno_function(peak_file)

peak_chicken_10k <- peak_anno_function(peak_file)





plotAnnoPie(peak_chicken)
#info_chicken <- peak_chicken@anno@elementMetadata
info_chicken <- readRDS("info_chicken.rds")


info_chicken$loc_id <- rownames(count_medips_keep)
res_comp1_medips$loc_id <- rownames(res_comp1_medips)
res_comp1_medips <- merge(res_comp1_medips , info_chicken , all=T)
res_comp2_medips$loc_id <- rownames(res_comp2_medips)
res_comp2_medips <- merge(res_comp2_medips , info_chicken , all=T)
res_comp3_medips$loc_id <- rownames(res_comp3_medips)
res_comp3_medips <- merge(res_comp3_medips , info_chicken , all=T)
res_comp4_medips$loc_id <- rownames(res_comp4_medips)
res_comp4_medips <- merge(res_comp4_medips , info_chicken , all=T)
res_comp5_medips$loc_id <- rownames(res_comp5_medips)
res_comp5_medips <- merge(res_comp5_medips , info_chicken , all=T)
write.table(res_comp1_medips , "Medips_StressMvsStressF.txt" , sep="\t" , row.names = F , quote = F)
write.table(res_comp2_medips , "Medips_ControlMvsControlF.txt" , sep="\t" , row.names = F , quote = F)
write.table(res_comp3_medips , "Medips_StressMvsControlM.txt" , sep="\t" , row.names = F , quote = F)
write.table(res_comp4_medips , "Medips_StressFvsControlF.txt" , sep="\t" , row.names = F , quote = F)
write.table(res_comp5_medips , "Medips_SVC_no_gender.txt" , sep="\t" , row.names = F , quote = F)



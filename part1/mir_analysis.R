library(xlsx)
library(readxl)
mir_counts <- read_excel("count_table_miR_1527_all.xlsx")
mir_counts <- as.data.frame(mir_counts)
rownames(mir_counts) <- mir_counts[,1]
mir_counts <- mir_counts[,-1]
mir_counts <- as.matrix(mir_counts)
threshold <- mir_counts > 1
table(rowSums(threshold))
keep_mirs <- rowSums(threshold) >= 4
mir_counts_keep <- mir_counts[keep_mirs ,]



plot(density(log(as.matrix(mir_counts_keep))))
plot(density(as.matrix(mir_counts_keep)))



transfor_lcpm_mirs <- voom(mir_counts_keep)
plot(density(transfor_lcpm_mirs$E))

mirs_pheno <- read_delim("targets.txt","\t", escape_double = FALSE, trim_ws = TRUE)
mirs_pheno <- mirs_pheno[order(match(mirs_pheno$Sample,colnames(transfor_lcpm_mirs))),]

group <- as.factor(mirs_pheno$sex)
group2 <- as.factor(mirs_pheno$treatment)
group3 <- as.factor(paste(mirs_pheno$sex , mirs_pheno$treatment , sep = "_"))

design <- model.matrix(~0+group3)
colnames(design) <- gsub("group3", "", colnames(design))
contr.matrix <- makeContrasts(
  comp1= male_stress-female_stress ,
  comp2= male_ctrl-female_ctrl ,
  comp3= male_stress-male_ctrl ,
  comp4= female_stress-female_ctrl ,
  levels = colnames(design))
vfit_mirs <- lmFit(transfor_lcpm_mirs, design)
vfit_mirs <- contrasts.fit(vfit_mirs, contrasts=contr.matrix)
efit_mirs <- eBayes(vfit_mirs)
res_comp1_mirs <- topTable(efit_mirs , coef = 1 , number = Inf)
res_comp2_mirs <- topTable(efit_mirs , coef = 2 , number = Inf)
res_comp3_mirs <- topTable(efit_mirs , coef = 3 , number = Inf)
res_comp4_mirs <- topTable(efit_mirs , coef = 4 , number = Inf)

design <- model.matrix(~0+group2)
colnames(design) <- gsub("group2", "", colnames(design))
design
contr.matrix <- makeContrasts(
  comp1= stress-ctrl ,
  levels = colnames(design))
contr.matrix
vfit_mirs <- lmFit(transfor_lcpm_mirs, design)
vfit_mirs <- contrasts.fit(vfit_mirs, contrasts=contr.matrix)
efit_mirs <- eBayes(vfit_mirs)
res_comp5_mirs <- topTable(efit_mirs , coef = 1 , number = Inf)

all_mir_test <- rownames(mir_counts_keep)
miRDB_v6_chicken_real_targets_df <- miRDB_v6_chicken_real_targets[miRDB_v6_chicken_real_targets$X1 %in% all_mir_test ,]
colnames(miRDB_v6_chicken_real_targets_df) <- c("mirbase_id" , "accession" , "pred_score")
miRDB_v6_chicken_real_targets_df <- merge(miRDB_v6_chicken_real_targets_df , annot_mirs)


res_comp1_mirs$mirbase_id <- rownames(res_comp1_mirs)
res_comp1_mirs <- merge(res_comp1_mirs , miRDB_v6_chicken_real_targets_df , all.x=T)


res_comp2_mirs$mirbase_id <- rownames(res_comp2_mirs)
res_comp2_mirs <- merge(res_comp2_mirs , miRDB_v6_chicken_real_targets_df , all.x=T)

res_comp3_mirs$mirbase_id <- rownames(res_comp3_mirs)
res_comp3_mirs <- merge(res_comp3_mirs , miRDB_v6_chicken_real_targets_df , all.x=T)

res_comp4_mirs$mirbase_id <- rownames(res_comp4_mirs)
res_comp4_mirs <- merge(res_comp4_mirs , miRDB_v6_chicken_real_targets_df , all.x=T)

res_comp5_mirs$mirbase_id <- rownames(res_comp5_mirs)
res_comp5_mirs <- merge(res_comp5_mirs , miRDB_v6_chicken_real_targets_df , all.x=T)


length(na.omit(unique(res_comp1_mirs$ensembl_id[res_comp1_mirs$P.Value<0.05 & res_comp1_mirs$pred_score>90])))
length(na.omit(unique(res_comp2_mirs$ensembl_id[res_comp2_mirs$P.Value<0.05 & res_comp2_mirs$pred_score>90])))
length(na.omit(unique(res_comp3_mirs$ensembl_id[res_comp3_mirs$P.Value<0.05 & res_comp3_mirs$pred_score>90])))
length(na.omit(unique(res_comp4_mirs$ensembl_id[res_comp4_mirs$P.Value<0.05 & res_comp4_mirs$pred_score>90])))
length(na.omit(unique(res_comp5_mirs$ensembl_id[res_comp5_mirs$P.Value<0.05 & res_comp5_mirs$pred_score>90])))




multimir_results <- get_multimir(org     = 'gga',
                                 mirna   = rownames(res_comp1_mirs),
                                 table   = 'validated',
                                 summary = TRUE)

par(mfrow=c(3,3))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.group2 <- group2
levels(col.group2) <-  brewer.pal(nlevels(col.group2), "Set2")
col.group2 <- as.character(col.group2)
col.group3 <- group3
levels(col.group3) <-  brewer.pal(nlevels(col.group3), "Set3")
col.group3 <- as.character(col.group3)

plotMDS(lcpm, labels=group, col=col.group)
title(main="A. gender")
plotMDS(lcpm, labels=group2, col=col.group2)
title(main="B. stress status")
plotMDS(lcpm, labels=group3, col=col.group3)
title(main="B. combined")

plotMDS(transfor_lcpm, labels=group, col=col.group)
title(main="A. gender")
plotMDS(transfor_lcpm, labels=group2, col=col.group2)
title(main="B. stress status")
plotMDS(transfor_lcpm, labels=group3, col=col.group3)
title(main="B. combined")

plotMDS(transfor_lcpm_mirs, labels=group, col=col.group)
title(main="A. gender")
plotMDS(transfor_lcpm_mirs, labels=group2, col=col.group2)
title(main="B. stress status")
plotMDS(transfor_lcpm_mirs, labels=group3, col=col.group3)
title(main="B. combined")

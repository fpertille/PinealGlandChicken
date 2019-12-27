##pathway_enrichment##
enrich_list <- c()
enrich_list$rnaseq_comp1 <- na.omit(res_comp1$ENTREZID[res_comp1$P.Value < 0.05])
enrich_list$rnaseq_comp2 <- na.omit(res_comp2$ENTREZID[res_comp2$P.Value < 0.05])
enrich_list$rnaseq_comp3 <- na.omit(res_comp3$ENTREZID[res_comp3$P.Value < 0.05])
enrich_list$rnaseq_comp4 <- na.omit(res_comp4$ENTREZID[res_comp4$P.Value < 0.05])
enrich_list$rnaseq_comp5 <- na.omit(res_comp5$ENTREZID[res_comp5$P.Value < 0.05])
enrich_list$medips_comp1 <- na.omit(res_comp1_medips$ENTREZID[res_comp1_medips$P.Value < 0.05])
enrich_list$medips_comp2 <- na.omit(res_comp2_medips$ENTREZID[res_comp2_medips$P.Value < 0.05])
enrich_list$medips_comp3 <- na.omit(res_comp3_medips$ENTREZID[res_comp3_medips$P.Value < 0.05])
enrich_list$medips_comp4 <- na.omit(res_comp4_medips$ENTREZID[res_comp4_medips$P.Value < 0.05])
enrich_list$medips_comp5 <- na.omit(res_comp5_medips$ENTREZID[res_comp5_medips$P.Value < 0.05])

library(clusterProfiler)
go_enrich <- compareCluster(geneClusters = enrich_list , 
                            fun = "enrichGO" , 
                            OrgDb = org.Gg.eg.db , 
                            ont="BP" , 
                            pvalueCutoff=1, 
                            qvalueCutoff=0.5, 
                            readable=T)
dotplot(go_enrich)



go_g_enrich <- compareCluster(geneClusters = enrich_list , 
                            fun = "groupGO" , 
                            OrgDb = org.Gg.eg.db , 
                            ont="BP" , 
                            level = 10 , 
                            readable=T)
dotplot(go_g_enrich)

kegg_enrich <- compareCluster(geneClusters = enrich_list , 
                            fun = "enrichKEGG" , 
                            org="gga",
                            pvalueCutoff=1, 
                            qvalueCutoff=1)
dotplot(kegg_enrich)

common_medips_RNA <- c()

##heat_try##
library(plyr)
library(data.table)
##comp1##
comp1_rna <- data.frame("geneId"= res_comp1$SYMBOL , "RNA_pval"=res_comp1$P.Value , stringsAsFactors = F)
comp1_rna <- ddply(comp1_rna , .(geneId) , numcolwise(min))
comp1_medips <- data.frame("geneId"= res_comp1_medips$SYMBOL , "medips_pval"=res_comp1_medips$P.Value , stringsAsFactors = F)
comp1_medips <- ddply(comp1_medips , .(geneId) , numcolwise(min))
comp_heat <- merge(comp1_rna , comp1_medips , all=T)
comp_heat <- na.omit(comp_heat)
comp_heat <- unique(comp_heat)

comp_heat_rna_side <- comp_heat[comp_heat$RNA_pval < 0.05 , ]
rownames(comp_heat_rna_side) <- comp_heat_rna_side$geneId
comp_heat_rna_side <- comp_heat_rna_side[,-1]
comp_heat_medips_side <- comp_heat[comp_heat$medips_pval < 0.05 , ]
rownames(comp_heat_medips_side) <- comp_heat_medips_side$geneId
comp_heat_medips_side <- comp_heat_medips_side[,-1]

comp_common <- comp_heat[comp_heat$RNA_pval < 0.05 & comp_heat$medips_pval < 0.05 , ]
rownames(comp_common) <- comp_common$geneId
comp_common <- comp_common[,-1]

common_medips_RNA$comp1 <- rownames(comp_common)

library(gplots)

heatmap.2(-log10(as.matrix(comp_heat_rna_side)) , trace = "none" , margins = c(15,8))
heatmap.2(-log10(as.matrix(comp_heat_medips_side)) , trace = "none" , margins = c(15,8))
heatmap.2(-log10(as.matrix(comp_common)) , trace = "none" , margins = c(15,8))

##comp2##
comp2_rna <- data.frame("geneId"= res_comp2$SYMBOL , "RNA_pval"=res_comp2$P.Value , stringsAsFactors = F)
comp2_rna <- ddply(comp2_rna , .(geneId) , numcolwise(min))
comp2_medips <- data.frame("geneId"= res_comp2_medips$SYMBOL , "medips_pval"=res_comp2_medips$P.Value , stringsAsFactors = F)
comp2_medips <- ddply(comp2_medips , .(geneId) , numcolwise(min))
comp_heat <- merge(comp2_rna , comp2_medips , all=T)
comp_heat <- na.omit(comp_heat)
comp_heat <- unique(comp_heat)

comp_heat_rna_side <- comp_heat[comp_heat$RNA_pval < 0.05 , ]
rownames(comp_heat_rna_side) <- comp_heat_rna_side$geneId
comp_heat_rna_side <- comp_heat_rna_side[,-1]
comp_heat_medips_side <- comp_heat[comp_heat$medips_pval < 0.05 , ]
rownames(comp_heat_medips_side) <- comp_heat_medips_side$geneId
comp_heat_medips_side <- comp_heat_medips_side[,-1]

comp_common <- comp_heat[comp_heat$RNA_pval < 0.05 & comp_heat$medips_pval < 0.05 , ]
rownames(comp_common) <- comp_common$geneId
comp_common <- comp_common[,-1]

common_medips_RNA$comp2 <- rownames(comp_common)

library(gplots)

heatmap.2(-log10(as.matrix(comp_heat_rna_side)) , trace = "none" , margins = c(15,8))
heatmap.2(-log10(as.matrix(comp_heat_medips_side)) , trace = "none" , margins = c(15,8))
heatmap.2(-log10(as.matrix(comp_common)) , trace = "none" , margins = c(15,8))

##comp3##
comp3_rna <- data.frame("geneId"= res_comp3$SYMBOL , "RNA_pval"=res_comp3$P.Value , stringsAsFactors = F)
comp3_rna <- ddply(comp3_rna , .(geneId) , numcolwise(min))
comp3_medips <- data.frame("geneId"= res_comp3_medips$SYMBOL , "medips_pval"=res_comp3_medips$P.Value , stringsAsFactors = F)
comp3_medips <- ddply(comp3_medips , .(geneId) , numcolwise(min))
comp_heat <- merge(comp3_rna , comp3_medips , all=T)
comp_heat <- na.omit(comp_heat)
comp_heat <- unique(comp_heat)

comp_heat_rna_side <- comp_heat[comp_heat$RNA_pval < 0.05 , ]
rownames(comp_heat_rna_side) <- comp_heat_rna_side$geneId
comp_heat_rna_side <- comp_heat_rna_side[,-1]
comp_heat_medips_side <- comp_heat[comp_heat$medips_pval < 0.05 , ]
rownames(comp_heat_medips_side) <- comp_heat_medips_side$geneId
comp_heat_medips_side <- comp_heat_medips_side[,-1]

comp_common <- comp_heat[comp_heat$RNA_pval < 0.05 & comp_heat$medips_pval < 0.05 , ]
rownames(comp_common) <- comp_common$geneId
comp_common <- comp_common[,-1]


common_medips_RNA$comp3 <- rownames(comp_common)
library(gplots)

heatmap.2(-log10(as.matrix(comp_heat_rna_side)) , trace = "none" , margins = c(15,8))
heatmap.2(-log10(as.matrix(comp_heat_medips_side)) , trace = "none" , margins = c(15,8))
heatmap.2(-log10(as.matrix(comp_common)) , trace = "none" , margins = c(15,8))

##comp4##
comp4_rna <- data.frame("geneId"= res_comp4$SYMBOL , "RNA_pval"=res_comp4$P.Value , stringsAsFactors = F)
comp4_rna <- ddply(comp4_rna , .(geneId) , numcolwise(min))
comp4_medips <- data.frame("geneId"= res_comp4_medips$SYMBOL , "medips_pval"=res_comp4_medips$P.Value , stringsAsFactors = F)
comp4_medips <- ddply(comp4_medips , .(geneId) , numcolwise(min))
comp_heat <- merge(comp4_rna , comp4_medips , all=T)
comp_heat <- na.omit(comp_heat)
comp_heat <- unique(comp_heat)

comp_heat_rna_side <- comp_heat[comp_heat$RNA_pval < 0.05 , ]
rownames(comp_heat_rna_side) <- comp_heat_rna_side$geneId
comp_heat_rna_side <- comp_heat_rna_side[,-1]
comp_heat_medips_side <- comp_heat[comp_heat$medips_pval < 0.05 , ]
rownames(comp_heat_medips_side) <- comp_heat_medips_side$geneId
comp_heat_medips_side <- comp_heat_medips_side[,-1]

comp_common <- comp_heat[comp_heat$RNA_pval < 0.05 & comp_heat$medips_pval < 0.05 , ]
rownames(comp_common) <- comp_common$geneId
comp_common <- comp_common[,-1]

common_medips_RNA$comp4 <- rownames(comp_common)

library(gplots)

heatmap.2(-log10(as.matrix(comp_heat_rna_side)) , trace = "none" , margins = c(15,8))
heatmap.2(-log10(as.matrix(comp_heat_medips_side)) , trace = "none" , margins = c(15,8))
heatmap.2(-log10(as.matrix(comp_common)) , trace = "none" , margins = c(15,8))

##comp5##
comp5_rna <- data.frame("geneId"= res_comp5$SYMBOL , "RNA_pval"=res_comp5$P.Value , stringsAsFactors = F)
comp5_rna <- ddply(comp5_rna , .(geneId) , numcolwise(min))
comp5_medips <- data.frame("geneId"= res_comp5_medips$SYMBOL , "medips_pval"=res_comp5_medips$P.Value , stringsAsFactors = F)
comp5_medips <- ddply(comp5_medips , .(geneId) , numcolwise(min))
comp_heat <- merge(comp5_rna , comp5_medips , all=T)
comp_heat <- na.omit(comp_heat)
comp_heat <- unique(comp_heat)

comp_heat_rna_side <- comp_heat[comp_heat$RNA_pval < 0.05 , ]
rownames(comp_heat_rna_side) <- comp_heat_rna_side$geneId
comp_heat_rna_side <- comp_heat_rna_side[,-1]
comp_heat_medips_side <- comp_heat[comp_heat$medips_pval < 0.05 , ]
rownames(comp_heat_medips_side) <- comp_heat_medips_side$geneId
comp_heat_medips_side <- comp_heat_medips_side[,-1]

comp_common <- comp_heat[comp_heat$RNA_pval < 0.05 & comp_heat$medips_pval < 0.05 , ]
rownames(comp_common) <- comp_common$geneId
comp_common <- comp_common[,-1]

common_medips_RNA$comp5 <- rownames(comp_common)

library(gplots)

heatmap.2(-log10(as.matrix(comp_heat_rna_side)) , trace = "none" , margins = c(15,8))
heatmap.2(-log10(as.matrix(comp_heat_medips_side)) , trace = "none" , margins = c(15,8))
heatmap.2(-log10(as.matrix(comp_common)) , trace = "none" , margins = c(15,8))



###methylation_change###
req_genes <- unique(unlist(common_medips_RNA))
comp1_meth <- res_comp1_medips[res_comp1_medips$SYMBOL %in% req_genes , ]
comp1_meth <- data.frame("ID"= paste(comp1_meth$SYMBOL , comp1_meth$annotation , sep="_") , "comp1_Meth_change"= comp1_meth$logFC , stringsAsFactors = F)
comp2_meth <- res_comp2_medips[res_comp2_medips$SYMBOL %in% req_genes , ]
comp2_meth <- data.frame("ID"= paste(comp2_meth$SYMBOL , comp2_meth$annotation , sep="_") , "comp2_Meth_change"= comp2_meth$logFC , stringsAsFactors = F)
comp3_meth <- res_comp3_medips[res_comp3_medips$SYMBOL %in% req_genes , ]
comp3_meth <- data.frame("ID"= paste(comp3_meth$SYMBOL , comp3_meth$annotation , sep="_") , "comp3_Meth_change"= comp3_meth$logFC , stringsAsFactors = F)
comp4_meth <- res_comp4_medips[res_comp4_medips$SYMBOL %in% req_genes , ]
comp4_meth <- data.frame("ID"= paste(comp4_meth$SYMBOL , comp4_meth$annotation , sep="_") , "comp4_Meth_change"= comp4_meth$logFC , stringsAsFactors = F)
comp5_meth <- res_comp5_medips[res_comp5_medips$SYMBOL %in% req_genes , ]
comp5_meth <- data.frame("ID"= paste(comp5_meth$SYMBOL , comp5_meth$annotation , sep="_") , "comp5_Meth_change"= comp5_meth$logFC , stringsAsFactors = F)


comp2_meth <- comp2_meth[order(match(comp2_meth$ID,comp1_meth$ID)),]
comp3_meth <- comp3_meth[order(match(comp3_meth$ID,comp1_meth$ID)),]
comp4_meth <- comp4_meth[order(match(comp4_meth$ID,comp1_meth$ID)),]
comp5_meth <- comp5_meth[order(match(comp5_meth$ID,comp1_meth$ID)),]



   


comp1_rna <- na.omit(res_comp1[res_comp1_medips$SYMBOL %in% req_genes , ])
comp1_rna <- data.frame("ID"= paste(comp1_rna$SYMBOL , comp1_rna$annotation , sep="_") , "comp1_rna_change"= comp1_rna$logFC , stringsAsFactors = F)
comp2_rna <- res_comp2[res_comp2$SYMBOL %in% req_genes , ]
comp2_rna <- data.frame("ID"= paste(comp2_rna$SYMBOL , comp2_rna$annotation , sep="_") , "comp2_rna_change"= comp2_rna$logFC , stringsAsFactors = F)
comp3_rna <- res_comp3[res_comp3$SYMBOL %in% req_genes , ]
comp3_rna <- data.frame("ID"= paste(comp3_rna$SYMBOL , comp3_rna$annotation , sep="_") , "comp3_rna_change"= comp3_rna$logFC , stringsAsFactors = F)
comp4_rna <- res_comp4[res_comp4$SYMBOL %in% req_genes , ]
comp4_rna <- data.frame("ID"= paste(comp4_rna$SYMBOL , comp4_rna$annotation , sep="_") , "comp4_rna_change"= comp4_rna$logFC , stringsAsFactors = F)
comp5_rna <- res_comp5[res_comp5$SYMBOL %in% req_genes , ]
comp5_rna <- data.frame("ID"= paste(comp5_rna$SYMBOL , comp5_rna$annotation , sep="_") , "comp5_rna_change"= comp5_rna$logFC , stringsAsFactors = F)


comp2_rna <- comp2_rna[order(match(comp2_rna$ID,comp1_rna$ID)),]
comp3_rna <- comp3_rna[order(match(comp3_rna$ID,comp1_rna$ID)),]
comp4_rna <- comp4_rna[order(match(comp4_rna$ID,comp1_rna$ID)),]
comp5_rna <- comp5_rna[order(match(comp5_rna$ID,comp1_rna$ID)),]


for_meth_matrix <- data.frame("Stress_MvsF_meth"=comp1_meth$comp1_Meth_change , 
                              "Stress_MvsF_rna"=comp1_rna$comp1_rna_change ,
                              "Control_MvsF_meth"=comp2_meth$comp2_Meth_change , 
                              "Control_MvsF_rna"=comp2_rna$comp2_rna_change ,
                              "Male_SVC_meth"=comp3_meth$comp3_Meth_change , 
                              "Male_SVC_rna"=comp3_rna$comp3_rna_change , 
                              "Female_SVC_meth"=comp4_meth$comp4_Meth_change , 
                              "Female_SVC_rna"=comp4_rna$comp4_rna_change , 
                              "SVC_meth" =comp5_meth$comp5_Meth_change, 
                              "SVC_rna" =comp5_rna$comp5_rna_change,
                              stringsAsFactors = F)
rownames(for_meth_matrix) <- comp1_meth$ID
heatmap.2(as.matrix(for_meth_matrix), trace = "none" , margins = c(15,8) , dendrogram = "row" , Colv = F)   

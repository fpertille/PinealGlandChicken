rm(list = ls())
library(limma)
library(Glimma)
library(edgeR)
library(org.Gg.eg.db)
library(readr)
entrez_gene_matrix <- read_csv("gene_count_matrix_trimmed.csv")
entrez_gene_matrix <- as.data.frame(entrez_gene_matrix)
rownames(entrez_gene_matrix) <- entrez_gene_matrix$gene_id
entrez_gene_matrix <- entrez_gene_matrix[,-1]
entrez_gene_matrix <- as.matrix(entrez_gene_matrix)
#colnames(entrez_gene_matrix) <- substr(colnames(entrez_gene_matrix),1,nchar(colnames(entrez_gene_matrix))-39)
colnames(entrez_gene_matrix) <- substr(colnames(entrez_gene_matrix),1,nchar(colnames(entrez_gene_matrix))-48)
colnames(entrez_gene_matrix)


library(readr)
sample_info <- read_delim("sample_info.txt","\t", escape_double = FALSE, trim_ws = TRUE)
sample_info <- sample_info[sample_info$Sample_ID %in% colnames(entrez_gene_matrix) ,]
sample_info <- sample_info[order(match(sample_info$Sample_ID,colnames(entrez_gene_matrix))),]
sample_info$combined <- paste(sample_info$Gender , sample_info$Stress , sep = "_")

#medips_pheno <- medips_pheno[order(match(medips_pheno$Sample_ID,colnames(transfor_lcpm))),]
# Obtain CPMs
myCPM <- cpm(entrez_gene_matrix)
# Have a look at the output
head(myCPM)
# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 1
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 4
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- entrez_gene_matrix[keep,]
summary(keep)


##DGELIST object##
x <- DGEList(counts = counts.keep)
samplenames <-sample_info$Sample_ID
x$samples$group <- as.factor(sample_info$Gender)
x$samples$group2 <- as.factor(sample_info$Stress)
x$samples$group3 <- as.factor(sample_info$combined)
group <- as.factor(sample_info$Gender)
group2 <- as.factor(sample_info$Stress)
group3 <- as.factor(sample_info$combined)
#geneid <- rownames(x)
#genes <- select(org.Gg.eg.db, keys=geneid, columns=c("SYMBOL" , "ENTREZID"), keytype="ENSEMBL")
#genes <- genes[!duplicated(genes$ENTREZID),]
#x$genes <- genes


cpm <-cpm(x)
lcpm <- cpm(x, log=TRUE)

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)


summary(lcpm)


lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")



x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,3))
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



design <- model.matrix(~0+group3)
colnames(design) <- gsub("group3", "", colnames(design))
design

contr.matrix <- makeContrasts(
  comp1= M_Stress-F_Stress , 
  comp2= M_Control-F_Control , 
  comp3= M_Stress-M_Control , 
  comp4= F_Stress-F_Control ,
  levels = colnames(design))
contr.matrix


par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
#plotSA(efit, main="Final model: Mean-variance trend")


#summary(decideTests(efit , p.value = 0.9))

res_comp1 <- topTable(efit , coef = 1 , number = Inf)
res_comp2 <- topTable(efit , coef = 2 , number = Inf)
res_comp3 <- topTable(efit , coef = 3 , number = Inf)
res_comp4 <- topTable(efit , coef = 4 , number = Inf)


design <- model.matrix(~0+group2)
colnames(design) <- gsub("group2", "", colnames(design))
design
contr.matrix <- makeContrasts(
  comp1= Stress-Control ,
  levels = colnames(design))
contr.matrix
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
res_comp5 <- topTable(efit , coef = 1 , number = Inf)

xg <- as.data.frame(org.Gg.egENSEMBL2EG)
xs <- as.data.frame(org.Gg.egSYMBOL2EG)
annot_exprs <- merge(xg, xs , all=T)

res_comp1$ensembl_id <- rownames(res_comp1)
res_comp1 <- merge(res_comp1 , annot_exprs , all.x=T)

res_comp2$ensembl_id <- rownames(res_comp2)
res_comp2 <- merge(res_comp2 , annot_exprs , all.x=T)

res_comp3$ensembl_id <- rownames(res_comp3)
res_comp3 <- merge(res_comp3 , annot_exprs , all.x=T)

res_comp4$ensembl_id <- rownames(res_comp4)
res_comp4 <- merge(res_comp4 , annot_exprs , all.x=T)

res_comp5$ensembl_id <- rownames(res_comp5)
res_comp5 <- merge(res_comp5 , annot_exprs , all.x=T)










info_chicken_rna <- info_chicken
info_chicken_rna$geneId <-  substr(info_chicken_rna$geneId ,1,nchar(info_chicken_rna$geneId)-2)
info_chicken_rna <- info_chicken_rna[,-8]
res_comp1$geneId <- rownames(res_comp1)
res_comp1 <- merge(res_comp1 , info_chicken_rna , all=T)

res_comp2$geneId <- rownames(res_comp2)
res_comp2 <- merge(res_comp2 , info_chicken_rna , all=T)

res_comp3$geneId <- rownames(res_comp3)
res_comp3 <- merge(res_comp3 , info_chicken_rna , all=T)

res_comp4$geneId <- rownames(res_comp4)
res_comp4 <- merge(res_comp4 , info_chicken_rna , all=T)

res_comp5$geneId <- rownames(res_comp5)
res_comp5 <- merge(res_comp5 , info_chicken_rna , all=T)


###entrez_new_genome##
res_comp1$ensembl_id <- rownames(res_comp1)
res_comp1 <- merge(res_comp1 , xx , all.x=T)





write.table(res_comp1 , "RNAseq_trim_StressMvsStressF.txt" , sep="\t" , row.names = F , quote = F)
write.table(res_comp2 , "RNAseq_trim_ControlMvsControlF.txt" , sep="\t" , row.names = F , quote = F)
write.table(res_comp3 , "RNAseq_trim_StressMvsControlM.txt" , sep="\t" , row.names = F , quote = F)
write.table(res_comp4 , "RNAseq_trim_StressFvsControlF.txt" , sep="\t" , row.names = F , quote = F)
write.table(res_comp5 , "RNAseq_trim_SVC_no_gender.txt" , sep="\t" , row.names = F , quote = F)



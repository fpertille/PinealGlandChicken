load("mr.edgeR_PINEAL.rda")
Meth_pd <- read_delim("sample_info copy.txt","\t", escape_double = FALSE, trim_ws = TRUE)
meth_matrix <- mr.edgeR[,1:16]
meth_matrix <- unique(meth_matrix)
rownames(meth_matrix) <- paste(meth_matrix$chr , meth_matrix$start , meth_matrix$stop , meth_matrix$CF , sep = "_")
meth_matrix <- meth_matrix[,-1:-4]
colnames(meth_matrix) <- substring(colnames(meth_matrix) ,1, 10)
threshold <- meth_matrix>1
head(threshold)
table(rowSums(threshold))
keep <- rowSums(threshold) >= 10
counts.keep <- meth_matrix[keep,]
counts.keep <- cpm(counts.keep)

plot(density(log2(as.matrix(counts.keep))))
plot(density(as.matrix(counts.keep)))

Meth_pd <- Meth_pd[Meth_pd$Sample_ID %in% colnames(meth_matrix) ,]
Meth_pd <- Meth_pd[order(match(Meth_pd$Sample_ID,colnames(meth_matrix))),]
Meth_pd$combined <- paste(Meth_pd$Gender , Meth_pd$Stress , sep = "_")


##model without log transformation
design <- model.matrix(~0+combined , Meth_pd )
colnames(design) <- gsub("combined", "", colnames(design))
fit <- lmFit(as.matrix(counts.keep), design)

contrast.matrix <- makeContrasts(test1=F_Stress-F_Control, 
                                 test2=M_Stress-M_Control , 
                                 test3 =F_Stress-M_Stress , 
                                 test4 =F_Control-M_Control , 
                                  levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
test_1_sig_cov <- topTable(fit2, n=Inf, coef=1, adjust="BH")
test_2_sig_cov <- topTable(fit2, n=Inf, coef=2, adjust="BH")
test_3_sig_cov <- topTable(fit2, n=Inf, coef=3, adjust="BH")
test_4_sig_cov <- topTable(fit2, n=Inf, coef=4, adjust="BH")

##model with log transformation
design <- model.matrix(~0+combined , Meth_pd )
colnames(design) <- gsub("combined", "", colnames(design))
fit <- lmFit(log2(as.matrix(counts.keep)), design)

contrast.matrix <- makeContrasts(test1=F_Stress-F_Control, 
                                 test2=M_Stress-M_Control , 
                                 test3 =F_Stress-M_Stress , 
                                 test4 =F_Control-M_Control , 
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
test_5_sig_cov <- topTable(fit2, n=Inf, coef=1, adjust="BH")
test_6_sig_cov <- topTable(fit2, n=Inf, coef=2, adjust="BH")
test_7_sig_cov <- topTable(fit2, n=Inf, coef=3, adjust="BH")
test_8_sig_cov <- topTable(fit2, n=Inf, coef=4, adjust="BH")



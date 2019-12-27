network_ppi <- read_table2("9031.protein.links.detailed.v11.0.txt")
network_ppi$protein1 <- substring(network_ppi$protein1 , 6)
network_ppi$protein2 <- substring(network_ppi$protein2 , 6)
network_ppi <- network_ppi[network_ppi$combined_score >= 700 , ]
library(org.Gg.eg.db)
xx <- as.data.frame(org.Gg.egENSEMBLPROT2EG)


temp <- network_ppi[, c(1,2, 10)]
colnames(temp) <- c("prot_id" , "prot_2" , "Score")
first_merge <- merge(temp , xx , all.x=T)
first_merge <- first_merge[,-1]
first_merge <- first_merge[,c(3,1,2)]
colnames(first_merge) <- c("Node1" , "prot_id" , "Score")
second_merge <- merge(first_merge , xx , all.x=T)
second_merge <- second_merge[,-1]
second_merge <- second_merge[,c(1,3,2)]
colnames(second_merge) <- c("Node1" , "Node2" , "Score")
network_ppi <- na.omit(second_merge)
rm(first_merge)
rm(temp)
write.table(network_ppi , "Gallus_network_PPI_700.txt" , sep="\t" , row.names = F , quote = F)



##test_run
diff_gen_obj <- res_comp1[, c("ENTREZID" , "P.Value" , "adj.P.Val")]
diff_gen_obj <- diff_gen_obj[diff_gen_obj$adj.P.Val < 0.1 , ]
diff_gen_obj <- na.omit(diff_gen_obj[,-2])
rownames(diff_gen_obj) <- NULL
colnames(diff_gen_obj) <- c("gene" , "pvalue")
##MODIFIER

Sys.setenv(PATH = paste("/Users/badt/anaconda3/bin", Sys.getenv("PATH"), sep=":"))##location of python on your system
set.seed(143)
MOD_inp_obj <- MODifieRDev:::create_custom_input_object(diff_genes = diff_gen_obj)
diamond_mod_comp1 <- MODifieRDev::diamond(MODifieR_input = MOD_inp_obj , ppi_network = network_ppi , deg_cutoff = 0.1 , n_output_genes = 200 , seed_weight = 10 , include_seed =T)


##test_run
diff_gen_obj <- res_comp2[, c("ENTREZID" , "P.Value" , "adj.P.Val")]
diff_gen_obj <- diff_gen_obj[diff_gen_obj$adj.P.Val < 0.1 , ]
diff_gen_obj <- na.omit(diff_gen_obj[,-2])
rownames(diff_gen_obj) <- NULL
colnames(diff_gen_obj) <- c("gene" , "pvalue")
##MODIFIER

Sys.setenv(PATH = paste("/Users/badt/anaconda3/bin", Sys.getenv("PATH"), sep=":"))
set.seed(143)
MOD_inp_obj <- MODifieRDev:::create_custom_input_object(diff_genes = diff_gen_obj)
diamond_mod_comp2 <- MODifieR::diamond(MODifieR_input = MOD_inp_obj , ppi_network = network_ppi , deg_cutoff = 0.1 , n_output_genes = 200 , seed_weight = 10 , include_seed =T)


xx <- as.data.frame(org.Gg.egSYMBOL2EG)
write.table(xx$symbol[xx$gene_id %in% diamond_mod_comp3$module_genes] , "Male_svc_module.txt" , sep="\t" , row.names=F , quote=F)

##test_run
diff_gen_obj <- res_comp3[, c("ENTREZID" , "P.Value" , "adj.P.Val")]
diff_gen_obj <- diff_gen_obj[diff_gen_obj$P.Value< 0.05 , ]
diff_gen_obj <- na.omit(diff_gen_obj[,-3])
rownames(diff_gen_obj) <- NULL
colnames(diff_gen_obj) <- c("gene" , "pvalue")
##MODIFIER

Sys.setenv(PATH = paste("/Users/badt/anaconda3/bin", Sys.getenv("PATH"), sep=":"))
set.seed(143)
MOD_inp_obj <- MODifieRDev:::create_custom_input_object(diff_genes = diff_gen_obj)
diamond_mod_comp3 <- MODifieR::diamond(MODifieR_input = MOD_inp_obj , ppi_network = network_ppi , deg_cutoff = 0.05 , n_output_genes = 200 , seed_weight = 10 , include_seed =T)


##test_run
diff_gen_obj <- res_comp4[, c("ENTREZID" , "P.Value" , "adj.P.Val")]
diff_gen_obj <- diff_gen_obj[diff_gen_obj$P.Value< 0.05 , ]
diff_gen_obj <- na.omit(diff_gen_obj[,-3])
rownames(diff_gen_obj) <- NULL
colnames(diff_gen_obj) <- c("gene" , "pvalue")
##MODIFIER

Sys.setenv(PATH = paste("/Users/badt/anaconda3/bin", Sys.getenv("PATH"), sep=":"))
set.seed(143)
MOD_inp_obj <- MODifieRDev:::create_custom_input_object(diff_genes = diff_gen_obj)
diamond_mod_comp4 <- MODifieR::diamond(MODifieR_input = MOD_inp_obj , ppi_network = network_ppi , deg_cutoff = 0.05 , n_output_genes = 200 , seed_weight = 10 , include_seed =T)

##test_run
diff_gen_obj <- res_comp_main[, c("ENTREZID" , "P.Value" , "adj.P.Val")]
diff_gen_obj <- diff_gen_obj[diff_gen_obj$P.Value< 0.05 , ]
diff_gen_obj <- na.omit(diff_gen_obj[,-3])
rownames(diff_gen_obj) <- NULL
colnames(diff_gen_obj) <- c("gene" , "pvalue")
##MODIFIER

Sys.setenv(PATH = paste("/Users/badt/anaconda3/bin", Sys.getenv("PATH"), sep=":"))
set.seed(143)
MOD_inp_obj <- MODifieRDev:::create_custom_input_object(diff_genes = diff_gen_obj)
diamond_mod_comp_main <- MODifieR::diamond(MODifieR_input = MOD_inp_obj , ppi_network = network_ppi , deg_cutoff = 0.05 , n_output_genes = 200 , seed_weight = 10 , include_seed =T)



library(clusterProfiler)

##ONLY MODULES
for_path_modules <- list("Gender_fixed_module"= diamond_mod_comp_main$module_genes , 
                         "Male_female_stress_module"= diamond_mod_comp1$module_genes , 
                         "Male_female_Control_module"=  diamond_mod_comp2$module_genes, 
                         "Male_stress_v_control_module"=  diamond_mod_comp3$module_genes, 
                         "Female_stress_v_control_module"= diamond_mod_comp4$module_genes 
)


GO_analysis <- compareCluster(geneClusters = for_path_modules , fun = "enrichGO" , ont="BP" ,OrgDb=org.Gg.eg.db)
jpeg("GO_ANALYSIS_MODULES_ONLY.jpeg" , width = 1650 , height = 850)
dotplot(GO_analysis)
dev.off()
kegg_analysis <- compareCluster(geneClusters = for_path_modules , fun = "enrichKEGG" ,organism="gga" )
jpeg("KEGG_ANALYSIS_MODULES_ONLY.jpeg" , width = 1650 , height = 850)
dotplot(kegg_analysis , showCategory=50)
dev.off()

##ONLY DEG


for_path_DEG <- list( "Gender_fixed_genes"= res_comp_main$ENTREZID[res_comp_main$P.Value < 0.05] , 
                      "Male_female_stress_genes"= res_comp1$ENTREZID[res_comp1$adj.P.Val < 0.1] , 
                      "Male_female_Control_genes"=  res_comp2$ENTREZID[res_comp2$adj.P.Val < 0.1 ], 
                      "Male_stress_v_control_genes"=  res_comp3$ENTREZID[res_comp3$P.Value < 0.05 ], 
                      "Female_stress_v_control_genes"= res_comp4$ENTREZID[res_comp4$P.Value < 0.05 ])



GO_analysis <- compareCluster(geneClusters = for_path_DEG , fun = "enrichGO" , ont="BP" ,OrgDb=org.Gg.eg.db)
jpeg("GO_ANALYSIS_genes_ONLY.jpeg" , width = 1650 , height = 850)
dotplot(GO_analysis)
dev.off()
kegg_analysis <- compareCluster(geneClusters = for_path_DEG , fun = "enrichKEGG" ,organism="gga" )
jpeg("KEGG_ANALYSIS_genes_ONLY.jpeg" , width = 1650 , height = 850)
dotplot(kegg_analysis , showCategory=50)
dev.off()

##combined
for_path_comb <- c(for_path_modules , for_path_DEG)

GO_analysis <- compareCluster(geneClusters = for_path_comb , fun = "enrichGO" , ont="BP" ,OrgDb=org.Gg.eg.db)
jpeg("GO_ANALYSIS_combined.jpeg" , width = 1650 , height = 850)
dotplot(GO_analysis)
dev.off()
kegg_analysis <- compareCluster(geneClusters = for_path_comb , fun = "enrichKEGG" ,organism="gga" )
jpeg("KEGG_ANALYSIS_combined.jpeg" , width = 1650 , height = 850)
dotplot(kegg_analysis , showCategory=50)
dev.off()





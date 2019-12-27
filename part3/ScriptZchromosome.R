
#use this line to filter chrZ and normalize based only on Z chr
count_medips <-  dplyr::filter(count_medips, grepl("chrZ",chr))



keepZcomp1 <- res_comp1_medips_base
keepZcomp2 <- res_comp2_medips_base

keepcomp1_voom <-res_comp1_medips_base
keepcomp2_voom <-res_comp2_medips_base

keepcomp1_z <-dplyr::filter(keepcomp1_voom, grepl("chrZ",loc_id))
keepcomp2_z <-dplyr::filter(keepcomp2_voom, grepl("chrZ",loc_id))

head(keepcomp1_z)
cor.test(keepZcomp1$P.Value, keepcomp1_z$P.Value)
cor.test(keepZcomp2$P.Value, keepcomp2_z$P.Value)
plot(-log10(keepZcomp2$P.Value), -log10(keepcomp2_z$P.Value)): abline(a = 0,b = 1)
plot(-log10(keepZcomp1$P.Value), -log10(keepcomp1_z$P.Value)): abline(a = 0,b = 1)


###
#Pearson's product-moment correlation
#data:  keepZcomp1$P.Value and keepcomp1_z$P.Value
#t = 122.71, df = 7768, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.8045044 0.8196436
#sample estimates:
#  cor 
#0.8122107 


#Pearson's product-moment correlation
#data:  keepZcomp2$P.Value and keepcomp2_z$P.Value
#t = 163.24, df = 7768, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8748156 0.8848573
#sample estimates:
#      cor 
#0.8799347 


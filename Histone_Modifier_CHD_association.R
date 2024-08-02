#This code is to identify Histon_modifier genes that might show an association with CHD development


# Used Data: https://www.medrxiv.org/content/10.1101/2023.12.23.23300495v1.supplementary-material
# "media-2"
#Total cases 3876
#n_total_syndromic 1471
#n_total_nonsyndromic 2405
#Total control 45082

# Used Panel
# Cell signaling : https://www.cellsignal.com/learn-and-support/reference-tables/histone-modification-table#acetylation
# Goterms: Searched for Histone and modifier Term 
# OMIM : searchd only for Histone Term
# Reaktom: https://reactome.org/download-data Reaktom Pathway gene set 


#Columns description:

#csq_group: Variant consequence (hcLOF: high-confidence loss-of-function; missC: missense constrained; syn: Synonymous).
#SYMBOL: Gene symbol.
#n_het_cases: Number of heterozygous variants observed in all CHD probands (aCHD).
#n_het_syndromic: Number of heterozygous variants observed in syndromic CHD probands (sCHD).
#n_het_nonsyndromic: Number of heterozygous variants observed in non-syndromic CHD probands (nsCHD).
#n_het_controls: Number of heterozygous variants observed in controls.
#n_total_cases: Total number of CHD cases.
#n_total_syndromic: Total number of syndromic CHD cases.
#n_total_nonsyndromic: Total number of non-syndromic CHD cases.
#n_total_controls: Total number of controls.
#fet.p_value: P-value from Fisher exact test.
#fet.odds_ratio: Odd ratio from Fisher exact test.
#fet.ci_95_lower: Lower confidence interval from Fisher exact test.
#fet.ci_95_upper: Upper confidence interval from Fisher exact test.
#analysis: Type of proband compared to controls.
#maf: Minor variant allelic frequency.



## All list were merged in Histone_mod_Gene_list.csv
# total is 942 genes 

histone_panel<-read.csv2("Histone_mod_Gene_list.csv", header= TRUE,na.strings=c("","NA"))




aCHD_vs_controls<-read.csv2("aCHD_vs_controls.csv", header= TRUE,na.strings=c("","NA"))
sCHD_vs_controls<-read.csv2("sCHD_vs_controls.csv", header= TRUE,na.strings=c("","NA"))
nsCHD_vs_controls<-read.csv2("nsCHD_vs_controls.csv", header= TRUE,na.strings=c("","NA"))


#############       For all CHD cases       ###################

# calculate the corrected P values bonferoi and FDR 
#Adjustment via Bonferroni #
adjusted_p_value_bon <- NULL # create value for corrected p values for later 
bon<-p.adjust(aCHD_vs_controls$fet.p_value,method = "bonferroni") # apply the bonferoni correction on each p value 
p.value_bon<-c(adjusted_p_value_bon,bon) # save the adjusted p values in a vector 
aCHD_vs_controls$p.value_bon<- p.value_bon# put the saved adjusted p values in a new column 

#Adjustment for FDR #
adjusted_p_value_FDR <- NULL # create value for corrected p values for later 
FDR<-p.adjust(aCHD_vs_controls$fet.p_value,method = "BH") # apply the bonferoni correction on each p value 
p.value_FDR<-c(adjusted_p_value_FDR,FDR) # save the adjusted p values in a vector 
aCHD_vs_controls$p.value_FDR<- p.value_FDR# put the saved adjusted p values in a new column 



# take only the genes present in the Histone Modifier panel 
Hist_aCHD<-subset(aCHD_vs_controls, SYMBOL %in% histone_panel$Gene)


#filter out all with a p-value above 0,05
gene_of_intrest_bon<-filter(Hist_aCHD,p.value_bon<=0.05)
gene_of_intrest_FDR<-filter(Hist_aCHD,p.value_FDR<=0.05)

write.csv2(gene_of_intrest_FDR, file ="Hist_aCHD_pValue.csv", row.names= FALSE, quote = FALSE) # save list of genes with their corrected p-values





#filter out synonymous 
Hist_aCHD2<- Hist_aCHD[!grepl("syn",Hist_aCHD$csq_group),] #  !!!!!!!!!!!   not sure if i should do it / migth be interesting 



install.packages("qqman")
library("qqman")
# visualize reults for Hist_aCHD
qq(Hist_aCHD$fet.p_value, 
   main = "Q-Q Plot Hist_aCHD (qqman)",
   pch = 19,
   col = "black",
   cex = 1, las = 1,
   cex.axis=1.7,
   xlim=c(0, 3.5),
   ylim=c(0, 5.5)) 



##############  For Syndromic only      ###############

# calculate the corrected P values bonferoi and FDR 
#Adjustment via Bonferroni #
adjusted_p_value_bon <- NULL # create value for corrected p values for later 
bon<-p.adjust(sCHD_vs_controls$fet.p_value,method = "bonferroni") # apply the bonferoni correction on each p value 
p.value_bon<-c(adjusted_p_value_bon,bon) # save the adjusted p values in a vector 
sCHD_vs_controls$p.value_bon<- p.value_bon# put the saved adjusted p values in a new column 

#Adjustment for FDR #
adjusted_p_value_FDR <- NULL # create value for corrected p values for later 
FDR<-p.adjust(sCHD_vs_controls$fet.p_value,method = "BH") # apply the bonferoni correction on each p value 
p.value_FDR<-c(adjusted_p_value_FDR,FDR) # save the adjusted p values in a vector 
sCHD_vs_controls$p.value_FDR<- p.value_FDR# put the saved adjusted p values in a new column 



# take only the genes present in the Histone Modifier panel 
Hist_sCHD<-subset(sCHD_vs_controls, SYMBOL %in% histone_panel$Gene)


#filter out all with a p-value above 0,05
gene_of_intrest_bon<-filter(Hist_sCHD,p.value_bon<=0.05)
gene_of_intrest_FDR<-filter(Hist_sCHD,p.value_FDR<=0.05)

write.csv2(gene_of_intrest_FDR, file ="Hist_sCHD_pValue.csv", row.names= FALSE, quote = FALSE) # save list of genes with their corrected p-values





#filter out synonymous 
Hist_sCHD2<- Hist_sCHD[!grepl("syn",Hist_sCHD$csq_group),] #  !!!!!!!!!!!   not sure if i should do it / migth be interesting 


# visualize reults for Hist_aCHD
qq(Hist_sCHD$fet.p_value, 
   main = "Q-Q Plot Hist_sCHD (qqman)",
   pch = 19,
   col = "black",
   cex = 1, las = 1,
   cex.axis=1.7,
   xlim=c(0, 3.5),
   ylim=c(0, 5.5)) 



###########  For nonsyndromic CHD         ##############

# calculate the corrected P values bonferoi and FDR 
#Adjustment via Bonferroni #
adjusted_p_value_bon <- NULL # create value for corrected p values for later 
bon<-p.adjust(nsCHD_vs_controls$fet.p_value,method = "bonferroni") # apply the bonferoni correction on each p value 
p.value_bon<-c(adjusted_p_value_bon,bon) # save the adjusted p values in a vector 
nsCHD_vs_controls$p.value_bon<- p.value_bon# put the saved adjusted p values in a new column 

#Adjustment for FDR #
adjusted_p_value_FDR <- NULL # create value for corrected p values for later 
FDR<-p.adjust(nsCHD_vs_controls$fet.p_value,method = "BH") # apply the bonferoni correction on each p value 
p.value_FDR<-c(adjusted_p_value_FDR,FDR) # save the adjusted p values in a vector 
nsCHD_vs_controls$p.value_FDR<- p.value_FDR# put the saved adjusted p values in a new column 



# take only the genes present in the Histone Modifier panel 
Hist_nsCHD<-subset(nsCHD_vs_controls, SYMBOL %in% histone_panel$Gene)


#filter out all with a p-value above 0,05
gene_of_intrest_bon<-filter(Hist_nsCHD,p.value_bon<=0.05)
gene_of_intrest_FDR<-filter(Hist_nsCHD,p.value_FDR<=0.05)

# ------------>>>>   No genes of interest by only looking at non-syndromic cases 














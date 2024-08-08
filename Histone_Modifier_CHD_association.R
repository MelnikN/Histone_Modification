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

rm(list = ls())

## All list were merged in Histone_mod_Gene_list.csv
# total is 912 genes (Tier1 = 265 / Tier2 = 647) 

histone_panel<-read.csv2("raw_data/Histone_mod_Gene_list.csv", header= TRUE,na.strings=c("","NA"))


###########################       First for CHD case control data set       ###############

aCHD_vs_controls<-read.csv2("raw_data/aCHD_vs_controls.csv", header= TRUE,na.strings=c("","NA"))
sCHD_vs_controls<-read.csv2("raw_data/sCHD_vs_controls.csv", header= TRUE,na.strings=c("","NA"))
nsCHD_vs_controls<-read.csv2("raw_data/nsCHD_vs_controls.csv", header= TRUE,na.strings=c("","NA"))

#seperate the table according to their variant type (syn/LOF/missense)
# and calculate afterwards the corrected p_values 



#############       For all CHD cases       ###################


####     SYN 
aCHD_syn<- aCHD_vs_controls[grepl("syn",aCHD_vs_controls$csq_group),] 

adjusted_p_value_bon <- NULL # create value for corrected p values for later 
bon<-p.adjust(aCHD_syn$fet.p_value,method = "bonferroni") # apply the bonferoni correction on each p value 
p.value_bon<-c(adjusted_p_value_bon,bon) # save the adjusted p values in a vector 
aCHD_syn$p.value_bon<- p.value_bon# put the saved adjusted p values in a new column 

#Adjustment for FDR #
adjusted_p_value_FDR <- NULL # create value for corrected p values for later 
FDR<-p.adjust(aCHD_syn$fet.p_value,method = "BH") # apply the bonferoni correction on each p value 
p.value_FDR<-c(adjusted_p_value_FDR,FDR) # save the adjusted p values in a vector 
aCHD_syn$p.value_FDR<- p.value_FDR# put the saved adjusted p values in a new column 

# add another column to describe from which analysis and data set the gene come 
aCHD_syn$Situation <- with(aCHD_syn, paste(csq_group,analysis,"Case_Control", sep=" / "))




####  LOF
aCHD_LOF<-aCHD_vs_controls[grepl("hcLOF",aCHD_vs_controls$csq_group),]

#Adjustment via Bonferroni #
adjusted_p_value_bon <- NULL # create value for corrected p values for later 
bon<-p.adjust(aCHD_LOF$fet.p_value,method = "bonferroni") # apply the bonferoni correction on each p value 
p.value_bon<-c(adjusted_p_value_bon,bon) # save the adjusted p values in a vector 
aCHD_LOF$p.value_bon<- p.value_bon# put the saved adjusted p values in a new column 

#Adjustment for FDR #
adjusted_p_value_FDR <- NULL # create value for corrected p values for later 
FDR<-p.adjust(aCHD_LOF$fet.p_value,method = "BH") # apply the bonferoni correction on each p value 
p.value_FDR<-c(adjusted_p_value_FDR,FDR) # save the adjusted p values in a vector 
aCHD_LOF$p.value_FDR<- p.value_FDR# put the saved adjusted p values in a new column 

# add another column to describe from which analysis and data set the gene come 
aCHD_LOF$Situation <- with(aCHD_LOF, paste(csq_group,analysis,"Case_Control", sep=" / "))



####      Missense variants
aCHD_missC<-aCHD_vs_controls[grepl("missC",aCHD_vs_controls$csq_group),]

#Adjustment via Bonferroni #
adjusted_p_value_bon <- NULL # create value for corrected p values for later 
bon<-p.adjust(aCHD_missC$fet.p_value,method = "bonferroni") # apply the bonferoni correction on each p value 
p.value_bon<-c(adjusted_p_value_bon,bon) # save the adjusted p values in a vector 
aCHD_missC$p.value_bon<- p.value_bon# put the saved adjusted p values in a new column 

#Adjustment for FDR #
adjusted_p_value_FDR <- NULL # create value for corrected p values for later 
FDR<-p.adjust(aCHD_missC$fet.p_value,method = "BH") # apply the bonferoni correction on each p value 
p.value_FDR<-c(adjusted_p_value_FDR,FDR) # save the adjusted p values in a vector 
aCHD_missC$p.value_FDR<- p.value_FDR# put the saved adjusted p values in a new column 

# add another column to describe from which analysis and data set the gene come 
aCHD_missC$Situation <- with(aCHD_missC, paste(csq_group,analysis,"Case_Control", sep=" / "))





#############################################


##############  For Syndromic only      ###############

####     SYN 
sCHD_syn<- sCHD_vs_controls[grepl("syn",sCHD_vs_controls$csq_group),] 

adjusted_p_value_bon <- NULL # create value for corrected p values for later 
bon<-p.adjust(sCHD_syn$fet.p_value,method = "bonferroni") # apply the bonferoni correction on each p value 
p.value_bon<-c(adjusted_p_value_bon,bon) # save the adjusted p values in a vector 
sCHD_syn$p.value_bon<- p.value_bon# put the saved adjusted p values in a new column 

#Adjustment for FDR #
adjusted_p_value_FDR <- NULL # create value for corrected p values for later 
FDR<-p.adjust(sCHD_syn$fet.p_value,method = "BH") # apply the bonferoni correction on each p value 
p.value_FDR<-c(adjusted_p_value_FDR,FDR) # save the adjusted p values in a vector 
sCHD_syn$p.value_FDR<- p.value_FDR# put the saved adjusted p values in a new column 

# add another column to describe from which analysis and data set the gene come 
sCHD_syn$Situation <- with(sCHD_syn, paste(csq_group,analysis,"Case_Control", sep=" / "))




####  LOF
sCHD_LOF<-sCHD_vs_controls[grepl("hcLOF",sCHD_vs_controls$csq_group),]

#Adjustment via Bonferroni #
adjusted_p_value_bon <- NULL # create value for corrected p values for later 
bon<-p.adjust(sCHD_LOF$fet.p_value,method = "bonferroni") # apply the bonferoni correction on each p value 
p.value_bon<-c(adjusted_p_value_bon,bon) # save the adjusted p values in a vector 
sCHD_LOF$p.value_bon<- p.value_bon# put the saved adjusted p values in a new column 

#Adjustment for FDR #
adjusted_p_value_FDR <- NULL # create value for corrected p values for later 
FDR<-p.adjust(sCHD_LOF$fet.p_value,method = "BH") # apply the bonferoni correction on each p value 
p.value_FDR<-c(adjusted_p_value_FDR,FDR) # save the adjusted p values in a vector 
sCHD_LOF$p.value_FDR<- p.value_FDR# put the saved adjusted p values in a new column 

# add another column to describe from which analysis and data set the gene come 
sCHD_LOF$Situation <- with(sCHD_LOF, paste(csq_group,analysis,"Case_Control", sep=" / "))



####      Missense variants
sCHD_missC<-sCHD_vs_controls[grepl("missC",sCHD_vs_controls$csq_group),]

#Adjustment via Bonferroni #
adjusted_p_value_bon <- NULL # create value for corrected p values for later 
bon<-p.adjust(sCHD_missC$fet.p_value,method = "bonferroni") # apply the bonferoni correction on each p value 
p.value_bon<-c(adjusted_p_value_bon,bon) # save the adjusted p values in a vector 
sCHD_missC$p.value_bon<- p.value_bon# put the saved adjusted p values in a new column 

#Adjustment for FDR #
adjusted_p_value_FDR <- NULL # create value for corrected p values for later 
FDR<-p.adjust(sCHD_missC$fet.p_value,method = "BH") # apply the bonferoni correction on each p value 
p.value_FDR<-c(adjusted_p_value_FDR,FDR) # save the adjusted p values in a vector 
sCHD_missC$p.value_FDR<- p.value_FDR# put the saved adjusted p values in a new column 

# add another column to describe from which analysis and data set the gene come 
sCHD_missC$Situation <- with(sCHD_missC, paste(csq_group,analysis,"Case_Control", sep=" / "))




##############  For Non-Syndromic only      ###############

####     SYN 
nsCHD_syn<- nsCHD_vs_controls[grepl("syn",nsCHD_vs_controls$csq_group),] 

adjusted_p_value_bon <- NULL # create value for corrected p values for later 
bon<-p.adjust(nsCHD_syn$fet.p_value,method = "bonferroni") # apply the bonferoni correction on each p value 
p.value_bon<-c(adjusted_p_value_bon,bon) # save the adjusted p values in a vector 
nsCHD_syn$p.value_bon<- p.value_bon# put the saved adjusted p values in a new column 

#Adjustment for FDR #
adjusted_p_value_FDR <- NULL # create value for corrected p values for later 
FDR<-p.adjust(nsCHD_syn$fet.p_value,method = "BH") # apply the bonferoni correction on each p value 
p.value_FDR<-c(adjusted_p_value_FDR,FDR) # save the adjusted p values in a vector 
nsCHD_syn$p.value_FDR<- p.value_FDR# put the saved adjusted p values in a new column 

# add another column to describe from which analysis and data set the gene come 
nsCHD_syn$Situation <- with(nsCHD_syn, paste(csq_group,analysis,"Case_Control", sep=" / "))




####  LOF
nsCHD_LOF<-nsCHD_vs_controls[grepl("hcLOF",nsCHD_vs_controls$csq_group),]

#Adjustment via Bonferroni #
adjusted_p_value_bon <- NULL # create value for corrected p values for later 
bon<-p.adjust(nsCHD_LOF$fet.p_value,method = "bonferroni") # apply the bonferoni correction on each p value 
p.value_bon<-c(adjusted_p_value_bon,bon) # save the adjusted p values in a vector 
nsCHD_LOF$p.value_bon<- p.value_bon# put the saved adjusted p values in a new column 

#Adjustment for FDR #
adjusted_p_value_FDR <- NULL # create value for corrected p values for later 
FDR<-p.adjust(nsCHD_LOF$fet.p_value,method = "BH") # apply the bonferoni correction on each p value 
p.value_FDR<-c(adjusted_p_value_FDR,FDR) # save the adjusted p values in a vector 
nsCHD_LOF$p.value_FDR<- p.value_FDR# put the saved adjusted p values in a new column 

# add another column to describe from which analysis and data set the gene come 
nsCHD_LOF$Situation <- with(nsCHD_LOF, paste(csq_group,analysis,"Case_Control", sep=" / "))



####      Missense variants
nsCHD_missC<-nsCHD_vs_controls[grepl("missC",nsCHD_vs_controls$csq_group),]

#Adjustment via Bonferroni #
adjusted_p_value_bon <- NULL # create value for corrected p values for later 
bon<-p.adjust(nsCHD_missC$fet.p_value,method = "bonferroni") # apply the bonferoni correction on each p value 
p.value_bon<-c(adjusted_p_value_bon,bon) # save the adjusted p values in a vector 
nsCHD_missC$p.value_bon<- p.value_bon# put the saved adjusted p values in a new column 

#Adjustment for FDR #
adjusted_p_value_FDR <- NULL # create value for corrected p values for later 
FDR<-p.adjust(nsCHD_missC$fet.p_value,method = "BH") # apply the bonferoni correction on each p value 
p.value_FDR<-c(adjusted_p_value_FDR,FDR) # save the adjusted p values in a vector 
nsCHD_missC$p.value_FDR<- p.value_FDR# put the saved adjusted p values in a new column 

# add another column to describe from which analysis and data set the gene come 
nsCHD_missC$Situation <- with(nsCHD_missC, paste(csq_group,analysis,"Case_Control", sep=" / "))






# next step is to merge everything, filter for histone panel and add the tier 
#bind the table with each other 
sum_table_CHD<- bind_rows(nsCHD_LOF,nsCHD_syn,nsCHD_missC,
                          sCHD_LOF,sCHD_syn,sCHD_missC,
                          aCHD_LOF,aCHD_syn,aCHD_missC)
# number is correct nothing got lost

Hist_sum_table_CHD<-subset(sum_table_CHD, SYMBOL %in% histone_panel$Gene)
# add the genepanel columns 
Hist_sum_table_CHD<- merge(x=Hist_sum_table_CHD, y=histone_panel, by.x = "SYMBOL", by.y = "Gene")


Hist_sum_table_CHD <- Hist_sum_table_CHD %>% arrange(Hist_sum_table_CHD$p.value_FDR) # order list from lowest p-value to highest. 

write.csv2(Hist_sum_table_CHD, file ="significant_gene_tables/all_histone_genes_CaseControl_merged.csv", row.names= FALSE, quote = FALSE) # save list of genes with their corrected p-values



# filter for significant genes 
gene_of_intrest_bon<-filter(Hist_sum_table_CHD,p.value_bon<=0.05)
gene_of_intrest_FDR<-filter(Hist_sum_table_CHD,p.value_FDR<=0.05)





# safe significant_gene_tables
write.csv2(gene_of_intrest_FDR, file ="significant_gene_tables/significant_histone_genes_CaseControl_merged.csv", row.names= FALSE, quote = FALSE) # save list of genes with their corrected p-values






####################  try to make qq plot work 

# entweder neues column machen wo nur die ersten x rows einen namen haben 
# oder ich mache eine loop zum beschreiben der punkte die x mal läft von oben runter 

# test first only on 
hist_aCHD_LOF<-Hist_sum_table_CHD[grepl("hcLOF / all_cases / Case_Control",Hist_sum_table_CHD$Situation),]
hist_aCHD_syn<-Hist_sum_table_CHD[grepl("syn / all_cases / Case_Control",Hist_sum_table_CHD$Situation),]



#YEAH !!!!!!!    labeling the points works know correctly 
# create qq Plot with fastqq
fastqq (hist_aCHD_LOF$fet.p_value,
        main="Title",
        p2=NULL,
        logtransform=TRUE,
        speedup=TRUE,
        lambda=TRUE,
        maxP=NULL, # truncation of p values visible if nedd be 
        cex=1,   # punkt größe
        cex.axis=1# achsen beschriftungs größe
)
#calculate observed pvalue / to get the xy coordinates   
hist_aCHD_LOF$obs.p<- -log10(hist_aCHD_LOF$fet.p_value)
#calculate expected p-value
hist_aCHD_LOF$exp.p<- -log10((rank(hist_aCHD_LOF$fet.p_value, ties.method="first")-.5)/(length(hist_aCHD_LOF$fet.p_value)-1))
#add the gene Symbol via coordinates but only for teh first two / defined in head()
text(head(hist_aCHD_LOF$exp.p,n=2),head(hist_aCHD_LOF$obs.p,n=2),labels=head(hist_aCHD_LOF$SYMBOL,n=2), cex= 0.7, pos=3)



#text(hist_aCHD_LOF$exp.p,hist_aCHD_LOF$obs.p,labels=hist_aCHD_LOF$SYMBOL, cex= 0.7, pos=3)





##############    Function test ---> Funktioniert 

create_qq <- function(df, x_num){
  fastqq (df$fet.p_value,
          main=sprintf("QQ-Plot \n %s \n", head(hist_aCHD_LOF$Situation,n=1)),
          sub=head(df$Situation,n=1),
          p2=NULL,
          logtransform=TRUE,
          speedup=TRUE,
          lambda=TRUE,
          maxP=NULL, # truncation of p values visible if nedd be 
          cex=1,   # punkt größe
          cex.axis=1# achsen beschriftungs größe
  )
  #calculate observed pvalue / to get the xy coordinates   
  df$obs.p<- -log10(df$fet.p_value)
  #calculate expected p-value
  df$exp.p<- -log10((rank(df$fet.p_value, ties.method="first")-.5)/(length(df$fet.p_value)-1))
  #add the gene Symbol via coordinates but only for teh first two / defined in head()
  text(head(df$exp.p,n=x_num),head(df$obs.p,n=x_num),labels=head(df$SYMBOL,n=x_num), cex= 0.7, pos=3)
}

#test run 
create_qq(hist_aCHD_syn,4)
#it works perfectly 

############





















x<- -log10(hist_aCHD_syn$fet.p_value) # das ist die y achse 
y<- -log10(exp.pvalues) # das ist die x achse 
y
#Calculate expectations
exp.pvalues<-(rank(hist_aCHD_syn$fet.p_value, ties.method="first")+.5)/(length(hist_aCHD_syn$fet.p_value)+1)
x
qqPlot(hist_aCHD_syn$fet.p_value)
2.6937269489 #y
2.75202673 #x
text(2.6937269489,2.75202673,labels="X", cex= 0.7, pos=3)
text(2.4718781993,2.71444269,labels="Y", cex= 0.7, pos=3)
text(2.3257501636,2.47108330,labels="Z", cex= 0.7, pos=3)

hist_aCHD_syn$exp<- -log10(hist_aCHD_syn$fet.p_value)

 #x2
 #y2
qqplot()
plot(y,x,aps=1,ylim = c(0,3))
abline(0,1)
# generating a qq plot with ggplot2
ggplot(hist_aCHD_syn, aes(sample = -log(fet.p_value))) +
  geom_qq() +
  geom_qq_line()


ggqqplot(hist_aCHD_syn$fet.p_value,color="black",add=c("qqline","none"),)

qqnorm(hist_aCHD_syn$fet.p_value)

fastqq (hist_aCHD_LOF$fet.p_value,
        p2=NULL,
        logtransform=TRUE,
        pairwisecompare=TRUE,
        speedup=TRUE,
        lambda=TRUE,
        maxP=14,
        cex=1,   # 
        cex.axis=2) # achsen beschriftungs größe 

fastqq (hist_aCHD_syn$fet.p_value,
        p2=NULL,
        logtransform=TRUE,
        speedup=TRUE,
        lambda=TRUE,
        maxP=NULL, # truncation of p values visible if nedd be 
        cex=1,   # punkt größe
        cex.axis=1# achsen beschriftungs größe
        )  






text(2,2,labels="here", cex= 0.7, pos=3)
inst


library("qqman")
# visualize reults for Hist_aCHD
qq(Hist_aCHD_LOF$fet.p_value, 
   label=Hist_aCHD_LOF$SYMBOL,
   main = "Q-Q Plot Hist_aCHD (qqman)",
   pch = 19,
   col = "black",
   cex = 1, las = 1,
   cex.axis=1.7,
   xlim=c(0, 3.5),
   ylim=c(0, 5.5)) 


sCHD_syn


gg<-qq(hist_aCHD_LOF$fet.p_value, 
   main = "Q-Q Plot sCHD_syn (qqman)",
   pch = 19,
   col = "black",
   cex = 1, las = 1,
   cex.axis=1.7) 



Top_Hits = head(arrange(sCHD_syn,fet.p_value),5) 

# Add column label, containing the gene name for the top hits or nothing for all others
sCHD_syn$label = if_else(sCHD_syn$SYMBOL %in% Top_Hits$SYMBOL,  
                         sCHD_syn$SYMBOL, NULL)

ggrepel::geom_text_repel(aes(label = label),
                         size = 3, show.legend = FALSE) 




install.packages("devtools")
devtools::install_github('kaustubhad/fastman',build_vignettes = TRUE)

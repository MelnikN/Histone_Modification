#This code is to identify Histon_modifier genes that might show an association with CHD development


# Used Data: https://www.medrxiv.org/content/10.1101/2023.12.23.23300495v1.supplementary-material
# "media-2"
#Total cases 3876
#n_total_syndromic 1471
#n_total_nonsyndromic 2405
#Total control 45082

# Used Panel
# wiki Pathway: https://www.wikipathways.org/pathways/WP2369.html
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






####################  QQ-Plot generation 

# separate for the situations for aCHD 
aCHD_LOF<-sum_table_CHD[grepl("hcLOF / all_cases / Case_Control",sum_table_CHD$Situation),]
aCHD_syn<-sum_table_CHD[grepl("syn / all_cases / Case_Control",sum_table_CHD$Situation),]
aCHD_missC<- sum_table_CHD[grepl("missC / all_cases / Case_Control",sum_table_CHD$Situation),]

# separate for the situations for sCHD 
sCHD_LOF<-sum_table_CHD[grepl("hcLOF / syndromic / Case_Control",sum_table_CHD$Situation),]
sCHD_syn<-sum_table_CHD[grepl("syn / syndromic / Case_Control",sum_table_CHD$Situation),]
sCHD_missC <-sum_table_CHD[grepl("missC / syndromic / Case_Control",sum_table_CHD$Situation),]

# separate for the situations for nsCHD 
nsCHD_LOF<-sum_table_CHD[grepl("hcLOF / nonsyndromic / Case_Control",sum_table_CHD$Situation),]
nsCHD_syn<-sum_table_CHD[grepl("syn / nonsyndromic / Case_Control",sum_table_CHD$Situation),]
nsCHD_missC <-sum_table_CHD[grepl("missC / nonsyndromic / Case_Control",sum_table_CHD$Situation),]





##############    Function test ---> Funktioniert 
# the following function is a alternative version from me to optimize xachsis for labeling gene points 
fastqqx<-function (p1, p2 = NULL, colour, logtransform = TRUE, pairwisecompare = TRUE, 
          speedup = TRUE, lambda = TRUE, maxP = 14, fix_zero = TRUE, 
          cex = 0.6, cex.axis = 0.9, xlab, ylab, ...)
  {
  if (!is.numeric(p1)) {
    stop("p1 must be numeric.")
  }
  if (missing(colour)) {
    colour = "black"
  }
  f = !is.na(p1) & !is.na(p1) & !is.nan(p1) & !is.na(p1) & 
    !is.null(p1) & is.finite(p1)
  p1 = p1[f]
  if (logtransform) {
    if (any(p1 < 0)) {
      warning("negative p-values found in p1. will be truncated to 0 first and then converted to minimum p-value/2.")
    }
    if (any(p1 == 0)) {
      warning("some p-values in p1 equal to zero. check fix_zero behavior.")
    }
    if (any(p1 > 1)) {
      warning("some p-values > 1 in p1. will be truncated to 1.")
    }
    f = p1 < 0
    p1[f] = 0
    f = p1 > 1
    p1[f] = 1
    f = (p1 == 0)
    if (any(f)) {
      if (fix_zero) {
        mx = min(p1[!f])
        p1[f] = mx
      }
      else {
        p1 = p1[!f]
      }
    }
  }
  lmb = qchisq(median(p1), 1, lower.tail = F)/0.4549364
  if (length(p2) > 1) {
    if (length(p2) != length(p1)) {
      stop("p1 and p2 must have same length.")
    }
    else {
      if (!is.numeric(p2)) {
        stop("p2 must be numeric.")
      }
      f = !is.na(p2) & !is.na(p2) & !is.nan(p2) & !is.na(p2) & 
        !is.null(p2) & is.finite(p2)
      p2 = p2[f]
      if (logtransform) {
        if (any(p2 < 0)) {
          warning("negative p-values found in p1. will be truncated to 0 first and then converted to minimum p-value/2.")
        }
        if (any(p2 == 0)) {
          warning("some p-values in p2 equal to zero. check fix_zero behavior.")
        }
        if (any(p2 > 1)) {
          warning("some p-values > 1 in p1. will be truncated to 1.")
        }
        f = (p2 < 0)
        p2[f] = 0
        f = (p2 > 1)
        p2[f] = 1
        f = (p2 == 0)
        if (any(f)) {
          if (fix_zero) {
            mx = min(p2[!f])
            p2[f] = mx
          }
          else {
            p2 = p2[!f]
          }
        }
      }
      lmb = qchisq(median(p1), 1, lower.tail = F)/qchisq(median(p2), 
                                                         1, lower.tail = F)
      if (pairwisecompare) {
        if (logtransform) {
          p2 = -log10(p2)
          p1 = -log10(p1)
          xlbl = expression(-log[10](italic(p2)))
          ylbl = expression(-log[10](italic(p1)))
        }
        else {
          xlbl = expression(italic(p2))
          ylbl = expression(italic(p1))
        }
      }
      else {
        if (logtransform) {
          p2 = -log10(sort(p2))
          p1 = -log10(sort(p1))
          xlbl = expression(-log[10](italic(p2)))
          ylbl = expression(-log[10](italic(p1)))
        }
        else {
          p2 = sort(p2)
          p1 = sort(p1)
          xlbl = expression(italic(p2))
          ylbl = expression(italic(p1))
        }
      }
    }
  }
  else {
    if (logtransform) {
      p2 = -log10(ppoints(length(p1)))
      p1 = -log10(sort(p1))
      xlbl = expression(Expected ~ ~-log[10](italic(p1)))
      ylbl = expression(Observed ~ ~-log[10](italic(p1)))
    }
    else {
      p2 = ppoints(length(p1))
      p1 = sort(p1)
      xlbl = expression(Expected ~ ~italic(p1))
      ylbl = expression(Observed ~ ~italic(p1))
    }
  }
  if (!is.null(maxP)) {
    f = p1 <= -maxP
    p1[f] = -maxP
    f = p1 >= maxP
    p1[f] = maxP
    f = p2 <= -maxP
    p2[f] = -maxP
    f = p2 >= maxP
    p2[f] = maxP
  }
  else {
    f = p1 == Inf
    p1[f] = sort(p1)[length(p1) - 1]
    f = p1 == -Inf
    p1[f] = sort(p1)[2]
    f = p2 == Inf
    p2[f] = sort(p2)[length(p2) - 1]
    f = p2 == -Inf
    p2[f] = sort(p2)[2]
  }
  if (length(colour) > 1) {
    if (length(colour) != length(p1)) {
      stop("p-value and colour must have same length.")
    }
    else {
      m = data.frame(p1, p2, colour)
    }
    if (speedup) {
      digs = 3
      m$p1 = round(m$p1, digits = digs)
      m$p2 = round(m$p2, digits = digs)
      f = duplicated(m)
      m = m[!f, ]
    }
    rm(p1, p2, colour)
    if (!missing(xlab)) {
      xlbl = xlab
    }
    if (!missing(ylab)) {
      ylbl = ylab
    }
    fac1 = 1.2
    fac2 = 0.4
    xbnd = c(0, max(m$p2) + fac2)
    ybnd = c(0, max(m$p1) + fac1)
    plot(m$p2, m$p1, col = m$colour, pch = 20, cex = cex, 
         cex.axis = cex.axis, las = 1, xaxs = "i", yaxs = "i", 
         xlim = xbnd, ylim = ybnd, xlab = xlbl, ylab = ylbl, 
         ...)
    abline(0, 1, col = "red")
    if (lambda) {
      text(x = max(m$p2) * 0.05, y = max(m$p1), labels = bquote(lambda == 
                                                                  .(round(lmb, digits = 4))), adj = c(0, 1))
    }
    return(lmb)
  }
  else {
    m = as.data.frame(cbind(p1, p2))
    if (speedup) {
      digs = 3
      m$p1 = round(m$p1, digits = digs)
      m$p2 = round(m$p2, digits = digs)
      f = duplicated(m)
      m = m[!f, ]
    }
    rm(p1, p2)
    if (!missing(xlab)) {
      xlbl = xlab
    }
    if (!missing(ylab)) {
      ylbl = ylab
    }
    fac1 = 1
    fac2 = 0.4
    xbnd = c(0, max(m$p2) + fac2)
    ybnd = c(0, max(m$p1) + fac1)
    plot(m$p2, m$p1, col = colour, pch = 20, cex = cex, cex.axis = cex.axis, 
         las = 1, xaxs = "i", yaxs = "i", xlim = xbnd, ylim = ybnd, 
         xlab = xlbl, ylab = ylbl, ...)
    abline(0, 1, col = "red")
    if (lambda) {
      text(x = max(m$p2) * 0.05, y = max(m$p1), labels = bquote(lambda == 
                                                                  .(round(lmb, digits = 4))), adj = c(0, 1))
    }
    return(lmb)
  }
}

# the following function is used to generate QQ-Plots and label gene points 
create_qq <- function(df, x_num){
  df <-subset(df, fet.p_value >= 0.0000000000000000000000000000000000000000000000000001 & fet.p_value <= 0.9999999999) # remove 1 and zeros from the df 
  hist_vec<-histone_panel$Gene # create vector for hist panel 
  df<-df %>%        # check if gene is in hist panel 
    mutate(Hist_Panel= 
             case_when(
               SYMBOL %in% hist_vec ~"green",
               .default =  "black"
             ))
  df<- merge(x=df, y=histone_panel, by.x = "SYMBOL", by.y = "Gene",all=TRUE)
  df<-subset(df, !is.na(csq_group))
  
  df$counts[is.na(df$counts)] <- 0
  df<-df %>% 
    mutate(farbe= 
             case_when(
               counts <=1 ~"salmon",
               counts >=2 ~"black"
             ))
  
  #calculate observed pvalue / to get the xy coordinates   
  df$obs.p<- -log10(df$fet.p_value)
  #calculate expected p-value
  df$exp.p<- -log10((rank(df$fet.p_value, ties.method="first")-.5)/(length(df$fet.p_value)-1))
  #add the gene Symbol via coordinates but only for teh first x rows / defined in head()
  df <- df %>% arrange(df$fet.p_value)
  
  fastqqx (df$fet.p_value,
          main=sprintf("QQ-Plot \n %s \n", head(df$Situation,n=1)),
          sub=head(df$Situation,n=1),
          p2=NULL,
          logtransform=TRUE,
          speedup=TRUE,
          lambda=TRUE,
          fix_zero = TRUE, # exclude the 0 from the Data to make the lamda better 
          maxP=NULL, # truncation of p values visible if nedd be 
          cex=1,   # punkt größe
          col=df$farbe,
          cex.axis=1# achsen beschriftungs größe
  )
  
  text(head(df$exp.p,n=x_num),head(df$obs.p,n=x_num),labels=head(df$SYMBOL,n=x_num), cex= 0.7, pos=3,col=df$Hist_Panel)
  }


#Create all QQ Plots for vizualising 

png(filename="figures/QQ_PLOT_aCHD_syn.png",width =9 ,height = 5,units="in",res = 300)
create_qq(aCHD_syn,4)
dev.off()
png(filename="figures/QQ_PLOT_aCHD_LOF.png",width =9 ,height = 5,units="in",res = 300)
create_qq(aCHD_LOF,4)
dev.off()
png(filename="figures/QQ_PLOT_aCHD_missC.png",width =9 ,height = 5,units="in",res = 300)
create_qq(aCHD_missC,4)
dev.off()

png(filename="figures/QQ_PLOT_sCHD_syn.png",width =9 ,height = 5,units="in",res = 300)
create_qq(sCHD_syn,4)
dev.off()

png(filename="figures/QQ_PLOT_sCHD_LOF.png",width =9 ,height = 5,units="in",res = 300)
create_qq(sCHD_LOF,4)###
dev.off()

png(filename="figures/QQ_PLOT_sCHD_missC.png",width =9 ,height = 5,units="in",res = 300)
create_qq(sCHD_missC,4)
dev.off()


png(filename="figures/QQ_PLOT_nsCHD_syn.png",width =9 ,height = 5,units="in",res = 300)
create_qq(nsCHD_syn,4)
dev.off()

png(filename="figures/QQ_PLOT_nsCHD_LOF.png",width =9 ,height = 5,units="in",res = 300)
create_qq(nsCHD_LOF,4)
dev.off()

png(filename="figures/QQ_PLOT_nsCHD_missC.png",width =9 ,height = 5,units="in",res = 300)
create_qq(nsCHD_missC,4)
dev.off()
############
























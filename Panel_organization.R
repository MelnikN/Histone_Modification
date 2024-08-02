#This code is to identify Histon_modifier genes that might show an association with CHD development


# Used Data:

# Cell signaling : https://www.cellsignal.com/learn-and-support/reference-tables/histone-modification-table#acetylation
# Goterms: Searched for Histone and modifier Term 
# OMIM : searchd only for Histone Term
# Reaktom: https://reactome.org/download-data Reaktom Pathway gene set 



#######       Data List organization  ######

getwd()

#load data
cellsignal <- read.csv2("CellSignal_Histonmod_correct.csv", header= TRUE)
OMIM <- read.csv2("OMIM_Histone_mod_human.csv", header= TRUE,na.strings=c("","NA"))
GoTerms <- data.frame(read.table("GoTerms_Histon_Modifier2.txt",fill=T,header=F))
Reaktom<-read.csv2("Reaktom_Pathway_Histone.csv", header= TRUE)

#The following codes doesn't create a data frame as i want. How can I achive it ? 
#Reaktom <-GSA.read.gmt(paste("ReactomePathways.gmt",sep = "/"))






names(GoTerms) <-c("UniProtKB","Gene","Synonym","Type")
GoTerms<-GoTerms[grepl("Uni",GoTerms$UniProtKB),]

HistPenel<- GoTerms[2]
HistPenel<-HistPenel %>%
  rowwise() %>%
  mutate(Source = "GoTerms")


########      cell signal  #######


#str_replace_all(cellsignal, fixed(" "), "")
cellsignal<-separate_rows(cellsignal,Histone.modifying.Enzymes,sep=",") # seperate by , and make a new row for each word
cellsignal<-separate_rows(cellsignal,Histone.modifying.Enzymes,sep="/") # seperate by , and make a new row for each word
cellsignal <- cellsignal %>% mutate( Histone.modifying.Enzymes= str_squish(Histone.modifying.Enzymes)) # remove the space at the beginning



#cellsignal$Histone.modifying.Enzymes<- gsub('\\s+',"",cellsignal$Histone.modifying.Enzymes)

#as.data.frame(apply(cellsignal, 2, function(x) gsub('\\s+', '', x)))

cellsignal_panel<-distinct(cellsignal,Histone.modifying.Enzymes) # get the unique names from all genes

cellsignal_panel<-cellsignal_panel %>%
  rowwise() %>%
  mutate(Source = "CellSignal")

names(cellsignal_panel) <-c("Gene","Source")



#######      OMIM     ########
# ist ansich jetzt fertig mit approved symbol 

#Mouse.Gene..from.MGI.


#OMIM<-subset(OMIM, is.na(Mouse.Gene..from.MGI.))
#OMIM<-subset(OMIM, !is.na(Approved.Symbol))
#write.csv2(OMIM, file ="OMIM_Histone_mod_human.csv")

#OMIM$Mouse.Gene..from.MGI. <- NULL
OMIM_panel<-distinct(OMIM,Approved.Symbol) # get the unique names from all genes

OMIM_panel<-OMIM_panel %>%
  rowwise() %>%
  mutate(Source = "OMIM")

names(OMIM_panel) <-c("Gene","Source")


#########       Reaktom       ##############



Reaktom<-separate_rows(Reaktom,Gene_Name,sep=" ") # seperate by , and make a new row for each word
Reaktom_panel<-distinct(Reaktom,Gene_Name) # get the unique names from all genes
Reaktom_panel<-Reaktom_panel %>%
  rowwise() %>%
  mutate(Source = "Reaktom")
names(Reaktom_panel) <-c("Gene","Source")






############        Create the Gene List with Ranks       ##############


#merge all panels 
Panel<- bind_rows(HistPenel,OMIM_panel,Reaktom_panel,cellsignal_panel)
# number is correct. i dndt lose any gene names 




sum_table<-Panel %>%
  group_by(Gene) %>%  
  summarise(counts = n())
# this code groups all variant that show the same symbol and consequence 
# it counts only how often a unique ID/ Term occurs in the 





# summarize all sources in one column and get distinct genes only 
x<-Panel %>%
  group_by(Gene) %>%
  summarise(vlist = paste0(Source, collapse = ",")) %>% 
  distinct(Gene, .keep_all = TRUE)
# merge the counts and all sources 
my_list<-merge(sum_table,x,by="Gene", all=T)

write.csv2(my_list, file ="Histone_mod_Gene_list.csv",, row.names= FALSE)













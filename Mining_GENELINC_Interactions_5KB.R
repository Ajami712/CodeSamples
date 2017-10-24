#==========================================
# In this script, I take pre-processed public Hi-C data
# and look for all chromosomal interactions with a particular gene of interest.
# The genomic region of the gene is defined by the gene body itself as well as spans of the genome that had features of interest in acetylomic and methylomic data.
# Once I identify the regions of the genome that interact with this gene, I cross-reference those regions with gene annotations.
# This helps us define possible regulatory relationships between other genes and this gene.
#==========================================



# 1. Define genomic region of interest

ChromsBefore = c("1","2","3","4","5","6","7","8","9","10","11","12","13")
ChromsAfter = c("15","16","17","18","19","20","21","22","X")
GENELINC_Windows = c(56240000,56245000,56250000,56255000,56260000,56265000,56270000,56275000,56280000,56285000)

# 2. Load gene annotations

setwd("C:\\Users\\Sabah\\Documents\\Nassim\\Grad School Year 3\\EC Project\\Project Route 3 - Basic Analysis of New 10pt data\\Specific Analysis 6 - GENELINC\\4C connections 5KBres")

GeneName_Reference = read.csv("ListOfAllGenes.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
GeneName_Reference = GeneName_Reference[c("Associated.Gene.Name","EntrezGene.ID")]
GeneName_Reference = GeneName_Reference[complete.cases(GeneName_Reference$EntrezGene.ID),]
rownames(GeneName_Reference) = GeneName_Reference[,"EntrezGene.ID"]
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
tx_hg19 <- as.data.frame(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
tx_hg19["Associated.Gene.Name"] = GeneName_Reference[tx_hg19$gene_id,"Associated.Gene.Name"]

colnames(tx_hg19) = c("Chromosome.Name","Gene.Start..bp.","Gene.End..bp.","width","Strand","EntrezGene.ID","Associated.Gene.Name")
tx_hg19$width = NULL
tx_hg19["Chromosome.Name"] = substring(tx_hg19[,"Chromosome.Name"],4)
GeneLocations = tx_hg19


# 3. Define the annotate function, which assigns annotations to gene regions

Annotate <- function(regionBegin,regionEnd,Annotations = ChromosomalGeneLocations){
  NumRegions = length(regionBegin)
  FinalAnnotations = vector()
  for(i in 1:NumRegions){
    InvalidGenes_1 = Annotations[Annotations$Gene.End..bp. < regionBegin[i],]
    InvalidGenes_2 = Annotations[Annotations$Gene.Start..bp. > regionEnd[i],]
    ValidGenes = Annotations[!(Annotations$EntrezGene.ID %in% 
                                 c(InvalidGenes_1[,"EntrezGene.ID"],InvalidGenes_2[,"EntrezGene.ID"])),]
    if(dim(ValidGenes)[1] > 0){
      FinalAnnotations = c(FinalAnnotations,paste(ValidGenes[,"Associated.Gene.Name"],collapse=","))
    }else{
      InvalidGenes_1["Difference"] = regionBegin[i]-InvalidGenes_1["Gene.End..bp."]
      InvalidGenes_2["Difference"] = InvalidGenes_2["Gene.Start..bp."] - regionEnd[i]
      InvalidGenes = rbind(InvalidGenes_1,InvalidGenes_2)
      closest = head(InvalidGenes[InvalidGenes["Difference"] == min(InvalidGenes["Difference"]),],1)
      AnnoToAdd = paste(closest[,"Difference"],"bases away from",closest[,"Associated.Gene.Name"])
      FinalAnnotations = c(FinalAnnotations,AnnoToAdd)
    }
  }
  
  return(FinalAnnotations)
}


# 4. Search all chromosomes outside of chromosome 14 for interacting regions and then annotate them.

Found_Genes = vector()
for(i in ChromsBefore){
  Interactions = read.csv(paste("chr",i,"_14_5kb.RAWobserved",sep=""),sep="\t",stringsAsFactors=FALSE,header=FALSE)
  colnames(Interactions) = c(paste("chr",i,sep=""),"chr14","Signal")
  Interactions = Interactions[Interactions$chr14 %in% GENELINC_Windows,]
  Interactions = Interactions[order(Interactions$Signal,decreasing=TRUE),]
  
  ChromosomalGeneLocations = GeneLocations[GeneLocations$Chromosome.Name == as.numeric(i),]
  Interactions["WindowEnd"] = Interactions[paste("chr",i,sep="")] + 5000
  Interactions["Annotations"] = Annotate(Interactions[,paste("chr",i,sep="")],Interactions[,"WindowEnd"])
  
  Found_Genes = c(Found_Genes, Interactions[!grepl("bases away from",Interactions[,"Annotations"]),"Annotations"])
  Name = paste("annotated_chr",i,"_14.csv",sep="")
  write.csv(Interactions,Name)
}



for(i in ChromsAfter){
  Interactions = read.csv(paste("chr14_",i,"_5kb.RAWobserved",sep=""),sep="\t",stringsAsFactors=FALSE,header=FALSE)
  colnames(Interactions) = c("chr14",paste("chr",i,sep=""),"Signal")
  Interactions = Interactions[Interactions$chr14 %in% GENELINC_Windows,]
  Interactions = Interactions[order(Interactions$Signal,decreasing=TRUE),]
  
  ChromosomalGeneLocations = GeneLocations[GeneLocations$Chromosome.Name == i,]
  Interactions["WindowEnd"] = Interactions[paste("chr",i,sep="")] + 5000
  Interactions["Annotations"] = Annotate(Interactions[,paste("chr",i,sep="")],Interactions[,"WindowEnd"])
  
  Found_Genes = c(Found_Genes, Interactions[!grepl("bases away from",Interactions[,"Annotations"]),"Annotations"])
  Name = paste("annotated_chr14_",i,".csv",sep="")
  write.csv(Interactions,Name)
}




Found_Genes = unique(Found_Genes[!is.na(Found_Genes)])
write.csv(Found_Genes,"GENELINC Targets.csv")

setwd("C:\\Users\\Sabah\\Documents\\Nassim\\Grad School Year 3\\EC Project\\Project Route 3 - Basic Analysis of New 10pt data\\Part 2 - DE Analysis\\Gene Level Analysis")

LFC_Combined = read.csv("AllGenes_LFC.csv",header=TRUE,sep=",",row.names=1,stringsAsFactors=FALSE)
P_Combined = read.csv("AllGenes_P.csv",header=TRUE,sep=",",row.names=1,stringsAsFactors=FALSE)
PAdj_Combined = read.csv("AllGenes_PAdj.csv",header=TRUE,sep=",",row.names=1,stringsAsFactors=FALSE)

IDs = LFC_Combined[c("RefSeq","BaseMean")]

vsPS_Data_P = P_Combined[,colnames(P_Combined)[grepl("to.PS",colnames(P_Combined))]]
OSvsPS_Data_P = vsPS_Data_P[,colnames(vsPS_Data_P)[grepl("OS",colnames(vsPS_Data_P))]]


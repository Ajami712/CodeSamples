#============================
# This script runs basic DESEQ2 analysis on ten time point RNA-seq data.
# The experimental conditions are shown as "OS", "PS", and "ST".
# The ten time points are 1, 2, 3, 4, 6, 9, 12, 16, 20, and 24 hours.
#============================

# 1. Load Data
setwd("C:\\Users\\Sabah\\Documents\\Nassim\\Grad School Year 3\\EC Project\\Project Route 3 - Basic Analysis of New 10pt data\\Part 1 - Quality Checks")
Data = read.csv("RawCounts.csv",header=TRUE,sep=",",row.names=1,stringsAsFactors=FALSE)
setwd("C:\\Users\\Sabah\\Documents\\Nassim\\Grad School Year 3\\EC Project\\Project Route 3 - Basic Analysis of New 10pt data\\Part 2 - DE Analysis\\Gene Level Analysis")
GeneNames = Data["GeneName"]

# 2. The data is already aligned and annotated, but counts have been found on the transcript level.
# Here, we combine counts of transcripts from the same gene, in order to get gene-level counts.
Data = Data[order(Data["GeneName"]),]
Gene_Level_Data = as.data.frame(matrix(0.0,nrow=0,ncol=ncol(Data)),stringsAsFactors=FALSE)
colnames(Gene_Level_Data) = colnames(Data)
Genes = unique(Data[,"GeneName"])

Gene_Unique_Identifiers = data.frame(row.names=Genes,RefSeq = Genes,stringsAsFactors=FALSE)
for(gene in Genes){
  Subset = Data[Data["GeneName"] == gene,]
  Subset$GeneName = NULL
  Subset[gene,] = colSums(Subset)
  
  Gene_Unique_Identifiers[gene,"RefSeq"] = rownames(Subset)[1]
  Gene_Level_Data = rbind(Gene_Level_Data,Subset[gene,])
}


Data = Gene_Level_Data


# 3. DESeq2 Pipeline

library(DESeq2)
SumEx_Data = SummarizedExperiment(round(as.matrix(Data)))
colData(SumEx_Data)$Category = c("ST1","ST1","OS1","OS2",
                                 "OS3","OS4","OS6","OS9",
                                 "OS12","OS16","OS20","OS24",
                                 "PS1","PS2","PS3","PS4","PS6","PS9","PS12","PS16","PS20","PS24",
                                 "OS1","OS2","OS3","OS4","OS6","OS9","OS12","OS16","OS20","OS24",
                                 "PS1","PS2","PS3","PS4","PS6","PS9","PS12","PS16","PS20","PS24")
ddsSE = DESeqDataSet(SumEx_Data,design=~Category)
ddsSE = DESeq(ddsSE)

Pairs2Compare = matrix(c(c("OS1","PS1"),
                         c("OS2","PS2"),
                         c("OS3","PS3"),
                         c("OS4","PS4"),
                         c("OS6","PS6"),
                         c("OS9","PS9"),
                         c("OS12","PS12"),
                         c("OS16","PS16"),
                         c("OS20","PS20"),
                         c("OS24","PS24"),
                         c("PS2","PS1"),
                         c("PS3","PS1"),
                         c("PS4","PS1"),
                         c("PS6","PS1"),
                         c("PS9","PS1"),
                         c("PS12","PS1"),
                         c("PS16","PS1"),
                         c("PS20","PS1"),
                         c("PS24","PS1"),
                         c("OS2","OS1"),
                         c("OS3","OS1"),
                         c("OS4","OS1"),
                         c("OS6","OS1"),
                         c("OS9","OS1"),
                         c("OS12","OS1"),
                         c("OS16","OS1"),
                         c("OS20","OS1"),
                         c("OS24","OS1"),
                         c("OS1","ST1"),
                         c("OS2","ST1"),
                         c("OS3","ST1"),
                         c("OS4","ST1"),
                         c("OS6","ST1"),
                         c("OS9","ST1"),
                         c("OS12","ST1"),
                         c("OS16","ST1"),
                         c("OS20","ST1"),
                         c("OS24","ST1"),
                         c("PS1","ST1"),
                         c("PS2","ST1"),
                         c("PS3","ST1"),
                         c("PS4","ST1"),
                         c("PS6","ST1"),
                         c("PS9","ST1"),
                         c("PS12","ST1"),
                         c("PS16","ST1"),
                         c("PS20","ST1"),
                         c("PS24","ST1")),ncol=48)

AllSigs = list()
for(i in 1:ncol(Pairs2Compare)){
  res = results(ddsSE,contrast=c("Category",Pairs2Compare[,i]))
  res = as.data.frame(res)
  res[is.na(res)] = 1
  res_sigs = rownames(res[res$pvalue < 0.05,])
  ID = paste(Pairs2Compare[,i][1]," to ",Pairs2Compare[,i][2],sep="")
  AllSigs[[ID]] = res_sigs
  
  name = paste(Pairs2Compare[,i][1]," to ",Pairs2Compare[,i][2],".csv",sep="")
  write.csv(res,name)
}


# 4. Organize DESeq2 results into spreadsheets for easy access and downstream analysis.

FinalSigs = unique(unlist(AllSigs))

LFC_Combined = Gene_Unique_Identifiers
LFC_Combined$Filler = 0
LFC_Combined = LFC_Combined[FinalSigs,]
P_Combined = LFC_Combined
PAdj_Combined = LFC_Combined



LFC_Combined_Master = Gene_Unique_Identifiers
LFC_Combined_Master$Filler = 0
P_Combined_Master = LFC_Combined_Master
PAdj_Combined_Master = LFC_Combined_Master

for(i in 1:ncol(Pairs2Compare)){
  res = results(ddsSE,contrast=c("Category",Pairs2Compare[,i]))
  res = as.data.frame(res)
  res[is.na(res)] = 1
  res_DE = res[FinalSigs,]
  name = paste(Pairs2Compare[,i][1]," to ",Pairs2Compare[,i][2],sep="")
  if(i == 1){
    LFC_Combined["BaseMean"] = res_DE["baseMean"]
    P_Combined["BaseMean"] = res_DE["baseMean"]
    PAdj_Combined["BaseMean"] = res_DE["baseMean"]
    LFC_Combined_Master["BaseMean"] = res["baseMean"]
    P_Combined_Master["BaseMean"] = res["baseMean"]
    PAdj_Combined_Master["BaseMean"] = res["baseMean"]
  }
  LFC_Combined[name] = res_DE["log2FoldChange"]
  P_Combined[name] = res_DE["pvalue"]
  PAdj_Combined[name] = res_DE["padj"]
  LFC_Combined_Master[name] = res["log2FoldChange"]
  P_Combined_Master[name] = res["pvalue"]
  PAdj_Combined_Master[name] = res["padj"]
}
LFC_Combined$Filler = NULL
P_Combined$Filler = NULL
PAdj_Combined$Filler = NULL

LFC_Combined_Master$Filler = NULL
P_Combined_Master$Filler = NULL
PAdj_Combined_Master$Filler = NULL


# BEGIN MANUAL FIXES - Some genes get assigned to a refseqID that aren't in the databases (usually NRs)
# This happens because of splice variants.
# We can manually fix these by replacing their RefSeqIDs with alternatively spliced RefseqIDs that
# we know are in the databases (e.g. KEGG, CPDB).
Genes2Fix = c('ADIPOR1',	'ANAPC13',	'ATG13',	'BPHL',	'BTF3L4',	'C10orf58',	'C12orf32',	'C7orf49',	'CDKN1A',	'CLU',	'CRAT',
  'CSAG2',	'CYB5D2',	'DAXX',	'EFHD1',	'EXOC7',	'FAM178A',	'FBXO31',	'FRG2',	'GATSL2',	'GHRL',	'GPATCH8',
  'IL9R',	'INPP5F',	'MAPK8',	'MGC34034',	'MTMR2',	'NARS2',	'NCRNA00032',	'NDUFB11',	'NOL8',	'NT5C3',	'NUDT7',
  'OBFC2A',	'PARVG',	'PCYT2',	'PDCD6IP',	'PDIK1L',	'PGAP2',	'PLSCR4',	'POMGNT1',	'RHBG',	'RPL10',	'SDK1',
  'SPANXB1',	'SPANXD',	'TPI1',	'WWP2',	'XAGE1C',	'ZFP161',	'ZMYND11',	'ZNF630')	

RefSeq4Fix = c('NM_015999',	'NM_015391',	'NM_001142673',	'NM_004332',	'NM_001136497',	'NM_032333',	'NR_027365',	'NM_024033',	'NM_078467',	'NM_001831',	'NM_000755',
               'NM_001080848',	'NM_144611',	'NM_001350',	'NM_025202',	'NM_001013839',	'NM_001136123',	'NM_024735',	'NM_001005217',	'NM_001145064',	'NM_001134941',	'NM_001002909',
               'NM_002186',	'NM_014937',	'NM_002750',	'NR_027030',	'NM_201281',	'NM_024678',	'NR_026679',	'NM_001135998',	'NM_017948',	'NM_016489',	'NM_001105663',
               'NM_001031716',	'NM_022141',	'NM_001184917',	'NM_001162429',	'NM_152835',	'NR_027018',	'NM_001177304',	'NM_017739',	'NM_020407',	'NM_006013',	'NM_152744',
               'NM_032461',	'NM_032417',	'NM_001159287',	'NM_199424',	'NM_001097597',	'NM_003409',	'NM_006624',	'NM_001190255'
)

for(i in 1:length(Genes2Fix)){
  LFC_Combined[Genes2Fix[i],"RefSeq"] = paste(RefSeq4Fix[i],"_chrNA",sep="")
  P_Combined[Genes2Fix[i],"RefSeq"] = paste(RefSeq4Fix[i],"_chrNA",sep="")
  PAdj_Combined[Genes2Fix[i],"RefSeq"] = paste(RefSeq4Fix[i],"_chrNA",sep="")
  LFC_Combined_Master[Genes2Fix[i],"RefSeq"] = paste(RefSeq4Fix[i],"_chrNA",sep="")
  P_Combined_Master[Genes2Fix[i],"RefSeq"] = paste(RefSeq4Fix[i],"_chrNA",sep="")
  PAdj_Combined_Master[Genes2Fix[i],"RefSeq"] = paste(RefSeq4Fix[i],"_chrNA",sep="")
}

# END MANUAL FIXES




write.csv(LFC_Combined,"CombinedDE_LFC.csv")
write.csv(P_Combined,"CombinedDE_P.csv")
write.csv(PAdj_Combined,"CombinedDE_PAdj.csv")
write.csv(Gene_Unique_Identifiers,"Combined_Unique IDs for Genes.csv")


write.csv(LFC_Combined_Master,"AllGenes_LFC.csv")
write.csv(P_Combined_Master,"AllGenes_P.csv")
write.csv(PAdj_Combined_Master,"AllGenes_PAdj.csv")


# number of total DE genes
AllSigGenes = list()
for(i in 1:10){
  AllSigGenes[[i]] = unique(GeneNames[AllSigs[[i]],"GeneName"])
}

FinalSigGenes = unique(unlist(AllSigGenes))

par(mar = c(5,4,4,2)+0.1)
Number_Of_DE_Genes = as.numeric(summary(AllSigGenes)[,"Length"])
Time_Points = c(1,2,3,4,6,9,12,16,20,24)
plot(Time_Points,Number_Of_DE_Genes,ylim=c(0,4000),col='red',type="b")

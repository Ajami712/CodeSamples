#============================
# This script takes the output from DE Analysis_GeneLevel_Annotated.R and plots genes of interest.
# Genes to plot are specified in the Genes_I_Like vector by gene name.
#============================

# 1. Load data.

setwd("C:\\Users\\Sabah\\Documents\\Nassim\\Grad School Year 3\\EC Project\\Project Route 3 - Basic Analysis of New 10pt data\\Part 2 - DE Analysis\\Gene Level Analysis")

LFC_Combined = read.csv("AllGenes_LFC.csv",header=TRUE,sep=",",row.names=1,stringsAsFactors=FALSE)
P_Combined = read.csv("AllGenes_P.csv",header=TRUE,sep=",",row.names=1,stringsAsFactors=FALSE)
PAdj_Combined = read.csv("AllGenes_PAdj.csv",header=TRUE,sep=",",row.names=1,stringsAsFactors=FALSE)

IDs = LFC_Combined[c("RefSeq","BaseMean")]

vsST_Data_LFC = LFC_Combined[,colnames(LFC_Combined)[grepl("to.ST1",colnames(LFC_Combined))]]
vsST_Data_P = P_Combined[,colnames(P_Combined)[grepl("to.ST1",colnames(P_Combined))]]
vsST_Data_PAdj = PAdj_Combined[,colnames(PAdj_Combined)[grepl("to.ST1",colnames(PAdj_Combined))]]

PSvsST_Data_LFC = vsST_Data_LFC[,colnames(vsST_Data_LFC)[grepl("PS",colnames(vsST_Data_LFC))]]
PSvsST_Data_P = vsST_Data_P[,colnames(vsST_Data_P)[grepl("PS",colnames(vsST_Data_P))]]
PSvsST_Data_PAdj = vsST_Data_PAdj[,colnames(vsST_Data_PAdj)[grepl("PS",colnames(vsST_Data_PAdj))]]
OSvsST_Data_LFC = vsST_Data_LFC[,colnames(vsST_Data_LFC)[grepl("OS",colnames(vsST_Data_LFC))]]
OSvsST_Data_P = vsST_Data_P[,colnames(vsST_Data_P)[grepl("OS",colnames(vsST_Data_P))]]
OSvsST_Data_PAdj = vsST_Data_PAdj[,colnames(vsST_Data_PAdj)[grepl("OS",colnames(vsST_Data_PAdj))]]

Time = c(1,2,3,4,6,9,12,16,20,24)


# 2. Select genes.

Genes_I_Like = c('AP1M1',
                 'ARL8B',
                 'ASCL1',
                 'ATG4C',
                 'ATG5',
                 'ATOH1',
                 'BMP4',
                 'BMP7',
                 'BNIP3',
                 'CD83',
                 'CDKN1A',
                 'CNIH1',
                 'CSNK2A1',
                 'CYLD',
                 'DENND2A',
                 'DLGAP1',
                 'DMRT1',
                 'DNAJC19',
                 'DNMT3A',
                 'DPYSL4',
                 'DTNB',
                 'FOS',
                 'FOXD3',
                 'FOXJ2',
                 'GATA2',
                 'GJB1',
                 'GPR101',
                 'GSC',
                 'HOOK2',
                 'HOXB13',
                 'HOXB4',
                 'HOXB5',
                 'HOXB9',
                 'HOXC6',
                 'HOXC8',
                 'JARID2',
                 'JUNB',
                 'KLF2',
                 'KLF4',
                 'KLF9',
                 'LEFTY1',
                 'LEFTY2',
                 'MKRN1',
                 'MSRA',
                 'MYCL',
                 'NANOG',
                 'NOTCH3',
                 'NRP2',
                 'OTX2',
                 'PAWR',
                 'PDGFRA',
                 'PGS1',
                 'PHC1',
                 'PLA2G1B',
                 'PPP1R12A',
                 'PRDM1',
                 'PXN',
                 'PYCR2',
                 'RABIF',
                 'RAD23B',
                 'RAD51C',
                 'RPN1',
                 'RPS6KA1',
                 'RYBP',
                 'SALL1',
                 'SDE2',
                 'SERPINA3',
                 'SIX1',
                 'SIX4',
                 'SLC2A3',
                 'SMARCAD1',
                 'SNAI1',
                 'SNCG',
                 'SOCS3',
                 'SOX18',
                 'SOX2',
                 'SP1',
                 'SP7',
                 'SRSF3',
                 'STAB2',
                 'TCEA2',
                 'TCF15',
                 'TEX14',
                 'TFAP2C',
                 'TMED10',
                 'TMED10',
                 'TPT1',
                 'TRPM3',
                 'UBE2V1',
                 'ZIC3',
                 'ZMYM3')

IDs_I_Like = Genes_I_Like



# 3. plot genes.

for(elem in IDs_I_Like){
  Significant = (PSvsST_Data_P[elem,] < 0.05)
  Very_Significant = (PSvsST_Data_PAdj[elem,] < 0.05)
  Significant_OS = (OSvsST_Data_P[elem,] < 0.05)
  Very_Significant_OS = (OSvsST_Data_PAdj[elem,] < 0.05)
  Sizes = c(1,1,1,1,1,1,1,1,1,1)
  Sizes[Significant] = 2
  Sizes[Very_Significant] = 4
  
  Sizes_OS = c(1,1,1,1,1,1,1,1,1,1)
  Sizes_OS[Significant_OS] = 2
  Sizes_OS[Very_Significant_OS] = 4
  
  Title = paste(elem," - ",IDs[elem,"RefSeq"]," - Basemean: ",round(IDs[elem,"BaseMean"],2),sep="")
    
  plot(Time,PSvsST_Data_LFC[elem,],type='b',cex=Sizes,col='blue',ylim=c(-3,3),ylab="LFC",main=Title)
  par(new=T)
  plot(Time,OSvsST_Data_LFC[elem,],type='b',cex=Sizes_OS,col='red',ylim=c(-3,3),ylab="LFC")
  legend(20,3,c("PS vs ST","OS vs ST"),lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red"))
  abline(0,0,col="gray")
  
}

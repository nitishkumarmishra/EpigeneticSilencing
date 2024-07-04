## This code developed by methdology proposed by using Hui Shen
## Creating problem when there is no CpG in genomica range of given gene
library('org.Hs.eg.db')
library(Homo.sapiens)
library(biomaRt)
library(FDb.InfiniumMethylation.hg19)
#######################################################
genes <- read.csv("TSS-HGNC-name.csv") ### TSS-HGNC.CSV is the list of TCGA genes in RNASeqV2 level-3 data
yy<-as.character(genes$hgnc, na.rm=FALSE)
keytypes(Homo.sapiens)
#x <- org.Hs.egSYMBOL2EG
#mapped_genes <- mappedkeys(x)
#mapped_genes <- keys(Homo.sapiens, keytype="SYMBOL") #Both are same
hm450 <- get450k()
TSS.nearest <- getNearestTSS(hm450)##  Make file for nearest TSS for each CpGs
gene <- intersect(yy, TSS.nearest.uniq.gene)


################################################
##List of HGNC symbol whose corresponding ENTREZ-ID available
txs <- transcriptsBy(Homo.sapiens, 'gene', col='GENEID')
#use Entrez GeneID as a string, not numeric or factor
getProbes<-function(geneID){
  tmp <- org.Hs.egSYMBOL2EG[[geneID]]
  temp<-txs[[tmp]]
  if(is.null(temp)){
    probes<-NULL
  }else{
    upstream.probes<-names(subsetByOverlaps(hm450,flank(temp,1500, start=TRUE)))
    downstream.probes<-names(subsetByOverlaps(hm450,flank(temp,-1500,start=TRUE)))
    probes<-unique(c(upstream.probes,downstream.probes))
  }
  return(probes)
}
# gene1 <- gene[1:10]
# sapply(gene1, getProbes)

############################ Final line to get list of probes for all genes ##########
##listprobe is character vector of all probes for yy, while lisprobes1 is list of CpGs for yy
listprobes <- sapply(gene, getProbes)
#listprobes[['A1BG']]## Print CpG ID's for A1BG
#listprobes <- unlist(sapply(yy, getProbes))
#listprobes1 <-lapply(yy, getProbes)

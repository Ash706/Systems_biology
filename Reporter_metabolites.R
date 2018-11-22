---
  title: "Reporter metbolites for differentially expressed genes under palmitate treatment based on human metabolic network"
author: "Ashfaq Ali"
output:
  pdf_document:
  fig_caption: yes
html_notebook: default
html_document: default
word_document: default
---
  ## The differental expression analyses was performed using Fisher exact test in pairwise samples and fold changes and p-values were produced using meta analyses of fold chnages and p-values.
  
  ```{r, include=FALSE}
library(piano)
library(Hmisc)
library(limma)
library(impute)
library(readr)
library(stringr)
library(biomaRt)
```

```{r, echo=FALSE}
Differentially_expressed_HOMA <- read_csv("C:/Users/AALI0054/Documents/Data/Naba/GE_Bioinformatics/HOMA_associations_24wk_2.csv")
Differentially_expressed_HOMA<-Differentially_expressed_HOMA[!is.na(Differentially_expressed_HOMA$ensembl_gene_id),]
rownames(Differentially_expressed_HOMA)<-make.unique(Differentially_expressed_HOMA$ensembl_gene_id)

## Load the the gene metabolite association file
HMR2_gm_GSC <- read_delim("Z:/Ashfaq_Ali/GSM/Analyses_pipeline/Base_files/HMR2_Generic/HMR2_gm_GSC.sif",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
HMR2_gm_GSC<-HMR2_gm_GSC[,c(3,1)] ## Extract the gene names( 3rd columns)  and associated metabolites (column 1)

HMR2_gp_GSC <- read_delim("Z:/Ashfaq_Ali/GSM/Analyses_pipeline/Base_files/HMR2_Generic/Gene_pathway_HMR2_generic.csv",  ",", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
HMR2_gp_GSC<-HMR2_gp_GSC[, c(2,1)]

myPval<- Differentially_expressed_HOMA$P.Value
names(myPval)<-Differentially_expressed_HOMA$ensembl_gene_id
myFC<-Differentially_expressed_HOMA$logFC
names(myFC)<-Differentially_expressed_HOMA$ensembl_gene_id
direction<-ifelse(myFC>0, 1,-1)
myGsc<-loadGSC(HMR2_gm_GSC)
myGsc_gp<-loadGSC(HMR2_gp_GSC)

gene_stats<-cbind(myPval,myFC)
rownames(gene_stats)<-Differentially_expressed_HOMA$ensembl_gene_id
reporter_HOMA<-runGSA(geneLevelStats = myPval, directions = direction,gsc=myGsc, geneSetStat="reporter", signifMethod="nullDist", nPerm=1000, gsSizeLim=c(5,100))
pathways_HOMA<-runGSA(geneLevelStats= myPval,directions = direction, gsc=myGsc_gp, geneSetStat = "reporter", signifMethod="nullDist", nPerm=1000, gsSizeLim=c(5,100) )
```

## Metabolic pathways enriched under palmitate treatmesnt in isolated human islets

### List of differentally activated pathways in human islets after palmitate treatment

```{r, echo=FALSE}
gsaResTab_palm_path <- GSAsummaryTable(pathways_palm)
ii <- which(gsaResTab_palm_path[,11]<0.05)
gsaResTab_palm_path[ii,]
write.csv(gsaResTab_palm_path, "Z:/Ashfaq_Ali/GSM/Analyses_pipeline/Palmitat_vs_control/Pathway_list_palmitate.csv", row.names = FALSE)
```

```{r, echo=FALSE}
gsaResTab_palm_metabo <- GSAsummaryTable(reporter_palm)
ii <- which(gsaResTab_palm_metabo[,11]<0.05)
gsaResTab_palm_metabo[ii,]
write.csv(gsaResTab_palm_metabo, "Z:/Ashfaq_Ali/GSM/Analyses_pipeline/Palmitat_vs_control/Pathway_list_palmitate_2.csv", row.names = FALSE)
```
### A network representation of metabolic pathways enriched and interconnectivity between the pathways based on shared number of genes

```{r, echo=FALSE}
pdf("Z:/Ashfaq_Ali/GSM/Analyses_pipeline/Palmitat_vs_control/Palm_pathways.pdf", width = 8)
networkPlot(reporter_HOMA, class="distinct", direction="up", significance=0.01, label="names",adjusted = TRUE, lay = 4 , cexLabel = 1  , ncharLabel = 100,edgeColor = "yellow4", nodeSize = c(5,30), overlap = 1 , main = "Metabolic pathways over-represented in the palmitate condition")
networkPlot(pathways_HOMA, class="distinct", direction="down", significance=0.1, label="names",adjusted = FALSE, lay = 4, cexLabel = 0.9, edgeColor = "darkgreen", main = "Reporter Pathways in Palmitate treatment", edgeWidth = c(1,15), ncharLabel = 70, cexLegend = 1.2)
dev.off()
```

GSAheatmap(pathways_HOMA, cutoff = 4, ncharLabel = 100, columnnames = "abbr", adjusted = TRUE, cellnote = "pvalue", colorkey =TRUE,cex = 1)

GSAheatmap(reporter_HOMA, cutoff = 4, ncharLabel = 100, columnnames = "abbr", adjusted = TRUE, cellnote = "pvalue", colorkey =TRUE,cex = 1)


biocLite("biomaRt")
library(biomaRt)
library(stringr)
Differentially_expressed_Palmitate <- read_csv("Z:/AQAI - Ashfaq Ali/GSM/Palmitate_control/Differentially_expressed_Palmitate.csv")
bm <- useMart("ensembl")
bm <- useDataset("hsapiens_gene_ensembl", mart=bm)

ids_Palm<-Differentially_expressed_Palmitate$Gene_ID
Etrez_Palm <- getBM(attributes=c("hgnc_symbol", "entrezgene"), filters=c("hgnc_symbol"), ids_Palm, mart=bm)
head(merge.data.frame(x=Differentially_expressed_Palmitate,y=Etrez_Palm, by.x="Gene_ID", by.y="hgnc_symbol", all.x=TRUE))

Entrez_genes<-merge.data.frame(x=Differentially_expressed_Palmitate,y=Etrez_Palm, by.x="Gene_ID", by.y="hgnc_symbol", all.x=TRUE)
Entrez_genes$expression<-ifelse(Entrez_genes$Fold_change>0, 1,0)
Entrez_genes$entrezgene<-paste(Entrez_genes$entrezgene, ".","1", sep = "")
write.csv(Entrez_genes, "Z:/AQAI - Ashfaq Ali/GSM/Palmitate_control/Differentially_expressed_ENTREZ.csv", row.names = FALSE)                                

Differentially_expressed_Cytokine <- read_csv("Z:/AQAI - Ashfaq Ali/GSM/Cytokine_control/Differentially_expressed_cytokine.csv")
ids_cyt<-Differentially_expressed_Cytokine$`Gene name`
Etrez_cyt <- getBM(attributes=c("hgnc_symbol", "entrezgene"), filters=c("hgnc_symbol"), ids_cyt, mart=bm)

head(Etrez_cyt)

Entrez_genes_cyt<-merge.data.frame(x=Differentially_expressed_Cytokine,y=Etrez_cyt, by.x="Gene name", by.y="hgnc_symbol", all.x=TRUE)
Entrez_genes_cyt$expression<-ifelse(Entrez_genes_cyt$`Log2 fold change`>0, 1,-1)
Entrez_genes_cyt$entrezgene<-paste(Entrez_genes_cyt$entrezgene, ".","1", sep = "")
write.csv(Entrez_genes_cyt, "Z:/AQAI - Ashfaq Ali/GSM/Cytokine_control/Differentially_expressed_ENTREZ_cyt.csv", row.names = FALSE)                                

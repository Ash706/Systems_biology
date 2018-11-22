###Import the files with significant probes with fold changes
NBW_foldchanges_adj_p_2 <- read.delim("~/Desktop/Christa_ANders/Dante_calculations/NBW_foldchanges_adj_p_2.txt")
LBW_foldchages_adj_pvalue_2 <- read.delim("~/Desktop/Christa_ANders/Dante_calculations/LBW_foldchages_adj_pvalue_2.txt")

###Import the original files to extract q_values and annotations
LBW_adypocytes_expression_paired_NBW <- read.delim("~/Desktop/Christa_ANders/LBW_adypocytes_expression_paired_NBW.txt")
LBW_adypocytes_expression_paired_LBW <- read.delim("~/Desktop/Christa_ANders/LBW_adypocytes_expression_paired_LBW.txt")
####merge the two data sets and extract relevant columns
NBW_qval_fold<- merge.data.frame(x=NBW_foldchanges_adj_p_2,y=LBW_adypocytes_expression_paired_NBW, by.x="Illumina_ID", 
                                 by.y="IlluminaID")

LBW_qval_fold<- merge.data.frame(x=LBW_foldchages_adj_pvalue_2,y=LBW_adypocytes_expression_paired_LBW, by.x="Illumina_ID", 
                                 by.y="IlluminaID")
#Select significant probes
NBW_sigz<-NBW_qval_fold[NBW_qval_fold$qvalues_stemVSdiff<0.05,]
LBW_sigz<-LBW_qval_fold[LBW_qval_fold$qvalues_stemVSdiff<0.05, ]
##Remove the data files
rm(NBW_foldchanges_adj_p_2,LBW_foldchages_adj_pvalue_2,LBW_adypocytes_expression_paired_NBW, LBW_adypocytes_expression_paired_LBW)
rm(NBW_qval_fold,LBW_qval_fold)
#Write the output to a file
write.table(NBW_sigz, file="~/Desktop/Christa_ANders/R_calc/NBW_qval_foldc_Anno.txt",sep="\t", row.names=F, col.names = TRUE)
write.table(LBW_sigz, file="~/Desktop/Christa_ANders/R_calc/LBW_qval_foldc_Anno.txt",sep="\t", row.names=F, col.names = TRUE)
#Select relevant columns include illumina ID, Symbol, qvalues, Averrage expressions, fold changes etc
NBW_sig_select<- NBW_sigz[,c("Illumina_ID","Symbol","Fold.Change", "Absolute.Fold.Change","qvalues_stemVSdiff","diff_mean" , "diff_sd", "stem_mean","stem_sd", "ENSEMBL")]
LBW_sig_select<- LBW_sigz[,c("Illumina_ID","Symbol","Fold.Change", "Absolute.Fold.Change","qvalues_stemVSdiff","diff_mean" , "diff_sd", "stem_mean","stem_sd", "ENSEMBL")]
#Get gene coordinates
source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
bm <- useMart("ensembl")
bm <- useDataset("hsapiens_gene_ensembl", mart=bm)
ids_Palm<-Differentially_expressed_Palmitate$Gene_ID
Etrez_Palm <- getBM(attributes=c("hgnc_symbol", "entrezgene", "hgnc_id"), filters=c("hgnc_symbol"), ids_Palm, mart=bm)
eid2fb <- getBM(attributes=c("chromosome_name", "start_position" , "end_position", "hgnc_symbol", "entrezgene", "hgnc_id"), filters=c("hgnc_symbol"), ids_Palm, mart=bm)


NBW_gene_exp_symbol_loc_fold_qval<-merge.data.frame(x=NBW_sig_select,y=eid2fb, by.x="Symbol", 
                 by.y="hgnc_symbol", all.x=TRUE)
#extract location based on gene symbol for low birthweight data
ids_LBW<-paste(as.vector(LBW_sig_select$Symbol),collapse = ",")
eid2fb_LBW <- getBM(attributes=c("chromosome_name", "start_position" , "end_position", "hgnc_symbol", "entrezgene"), filters=c("hgnc_symbol"), ids_LBW, mart=bm)
LBW_gene_exp_symbol_loc_fold_qval<-merge.data.frame(x=LBW_sig_select,y=eid2fb, by.x="Symbol", 
                                                    by.y="hgnc_symbol", all.x=TRUE)
##write the files to the text

write.table(NBW_gene_exp_symbol_loc_fold_qval, file="~/Desktop/Christa_ANders/R_calc/To_Mandy/NBW_gene_exp_symbol_loc_fold_qval.txt",sep="\t", row.names=F, col.names = TRUE)
write.table(LBW_gene_exp_symbol_loc_fold_qval, file="~/Desktop/Christa_ANders/R_calc/To_Mandy/LBW_gene_exp_symbol_loc_fold_qval.txt",sep="\t", row.names=F, col.names = TRUE)

####Import paired samples methylation files for paired samples
LBW_adypocytes_methylation_3_pairs_removed_paired_NBW_full_dataset <- read.delim("~/Desktop/Christa_ANders/Illumina_methylation/Illumina Methylation/LBW_adypocytes_methylation_3_pairs_removed_paired_NBW_full_dataset.txt")
NBW_Meth_sigz<-LBW_adypocytes_methylation_3_pairs_removed_paired_NBW_full_dataset[LBW_adypocytes_methylation_3_pairs_removed_paired_NBW_full_dataset$qvalues_stemVSdiff<0.05,]
rm(LBW_adypocytes_methylation_3_pairs_removed_paired_NBW_full_dataset)
NBW_Meth_sig_delta<- NBW_Meth_sigz[,c("ILMNID","qvalues_stemVSdiff","diff_minus_stem_mean_beta" ,"GENOME_BUILD", "CHR","MAPINFO","UCSC_CPG_ISLANDS_NAME","UCSC_REFGENE_NAME")]
write.table(NBW_Meth_sig_delta, file="~/Desktop/Christa_ANders/R_calc/NBW_Meth_sig_delta.txt", sep="\t")
#aaa<-strsplit(as.character(NBW_Meth_sig_delta$UCSC_REFGENE_NAME),";") #split the column by separator
##to convert the list to columns do following 
#library("plyr")
#df <- ldply(aaa)
#colnames(df) <- c("REFGEN1", "REFGEN2","REFGEN3")
LBW_adypocytes_methylation_3_pairs_removed_paired_LBW_full_dataset <- read.delim("~/Desktop/Christa_ANders/Illumina_methylation/Illumina Methylation/LBW_adypocytes_methylation_3_pairs_removed_paired_LBW_full_dataset.txt")
LBW_Meth_sigz<-LBW_adypocytes_methylation_3_pairs_removed_paired_LBW_full_dataset[LBW_adypocytes_methylation_3_pairs_removed_paired_LBW_full_dataset$qvalues_stemVSdiff<0.05,]
rm(LBW_adypocytes_methylation_3_pairs_removed_paired_LBW_full_dataset)
LBW_Meth_sig_delta<- LBW_Meth_sigz[,c("ILMNID","qvalues_stemVSdiff","diff_minus_stem_mean_beta" ,"GENOME_BUILD", "CHR","MAPINFO","UCSC_CPG_ISLANDS_NAME","UCSC_REFGENE_NAME")]
write.table(LBW_Meth_sig_delta, file="~/Desktop/Christa_ANders/R_calc/LBW_Meth_sig_delta.txt", sep="\t")
######to perform gene leve methylation changes analyses
NBW_Meth_genes1 <- read.delim("~/Desktop/Christa_ANders/R_calc/Meth_genes/NBW_Meth_genes1.txt")
LBW_Meth_sig_delta <- read.delim("~/Desktop/Christa_ANders/R_calc/Meth_genes/LBW_Meth_sig_delta.txt")
###Sum delta methylation values over genes for NBW people
NBW_Meth_genes_delta_sums<-tapply(NBW_Meth_genes1$diff_minus_stem_mean_beta,NBW_Meth_genes1$UCSC_REFGENE_NAME,sum, simplify=T)
NBW_Meth_genes_delta_sums<-as.data.frame(NBW_Meth_genes_delta_sums)
NBW_Meth_genes_delta_sums<-NBW_Meth_genes_delta_sums[-1,]
NBW_Meth_genes_delta_sums<-as.data.frame(NBW_Meth_genes_delta_sums)
NBW_Meth_genes_delta_sums<-cbind(rownames(NBW_Meth_genes_delta_sums),NBW_Meth_genes_delta_sums)
###Sum delta methylation values over genes LBW people
LBW_Meth_genes_delta_sums<-tapply(LBW_Meth_sig_delta$diff_minus_stem_mean_beta,LBW_Meth_sig_delta$UCSC_REFGENE_NAME,sum, simplify=T)
LBW_Meth_genes_delta_sums<-as.data.frame(LBW_Meth_genes_delta_sums)
LBW_Meth_genes_delta_sums<-LBW_Meth_genes_delta_sums[-1,]
LBW_Meth_genes_delta_sums<-as.data.frame(LBW_Meth_genes_delta_sums)
LBW_Meth_genes_delta_sums<-cbind(rownames(LBW_Meth_genes_delta_sums),LBW_Meth_genes_delta_sums)
#####count number of significant methylation for NBW
NBW_Meth_genes_count<-count(NBW_Meth_genes1, "UCSC_REFGENE_NAME")[-1,]
NBW_Meth_cout_sum<-cbind(NBW_Meth_genes_count,NBW_Meth_genes_delta_sums)
LBW_Meth_genes_count<-count(LBW_Meth_sig_delta, "UCSC_REFGENE_NAME")[-1,]
LBW_Meth_cout_sum<-cbind(LBW_Meth_genes_count,LBW_Meth_genes_delta_sums)
##### fetch chromosolam locations for NBW genes methylated
ids_NBW_Meth<-paste(as.vector(NBW_Meth_genes_count$UCSC_REFGENE_NAME),collapse = ",")
eid2fb_NBW <- getBM(attributes=c("chromosome_name", "start_position" , "end_position", "hgnc_symbol"), filters=c("hgnc_symbol"), ids_NBW_Meth, mart=bm)
NBW_Meth_count_sum_loc<-merge.data.frame(x=NBW_Meth_cout_sum,y=eid2fb_NBW, by.x="UCSC_REFGENE_NAME", 
                                                    by.y="hgnc_symbol", all.x=TRUE, all.y=TRUE)

##### fetch chromosolam locations for LBW genes methylated
ids_LBW_Meth<-paste(as.vector(LBW_Meth_genes_count$UCSC_REFGENE_NAME),collapse = ",")
eid2fb_LBW <- getBM(attributes=c("chromosome_name", "start_position" , "end_position", "hgnc_symbol"), filters=c("hgnc_symbol"), ids_LBW_Meth, mart=bm)
LBW_Meth_count_sum_loc<-merge.data.frame(x=LBW_Meth_cout_sum,y=eid2fb_LBW, by.x="UCSC_REFGENE_NAME", 
                                         by.y="hgnc_symbol", all.x=TRUE, all.y=TRUE)
write.table(NBW_Meth_count_sum_loc, file="~/Desktop/Christa_ANders/R_calc/To_Mandy/NBW_Meth_count_sum_loc.txt", sep="\t")
write.table(LBW_Meth_count_sum_loc, file="~/Desktop/Christa_ANders/R_calc/To_Mandy/LBW_Meth_count_sum_loc.txt", sep="\t")

#write a file with count and gene expression fold chages  and q values
colnames(LBW_Meth_genes_count)<-c("Gena_Name", "LBW_Meth_sites")
colnames(NBW_Meth_genes_count)<-c("Gena_Name", "NBW_Meth_sites")
Meth_count_gene_fold<- merge.data.frame(x=LBW_Meth_genes_count,y=NBW_Meth_genes_count, by.x="Gena_Name", by.y="Gena_Name", all.x=TRUE, all.y=TRUE)
NBW_genes<-NBW_sig_select[, c(2,3,5)]
colnames(NBW_genes)<-c("Gena_Name", "NBW_stemvsdiff_fold_changes", "NBW_qvalues_stemVSdiff")
LBW_genes<-LBW_sig_select[, c(2,3,5)]
colnames(LBW_genes)<-c("Gena_Name", "LBW_stemvsdiff_fold_changes", "LBW_qvalues_stemVSdiff")
Gene_fold_NBW_LBW<- merge.data.frame(x=NBW_genes,y=LBW_genes, by.x="Gena_Name", by.y="Gena_Name", all.x=TRUE, all.y=TRUE)

Meth_count_Gene_fold_NBW_LBW<- merge.data.frame(x=Meth_count_gene_fold,y=Gene_fold_NBW_LBW, by.x="Gena_Name", by.y="Gena_Name", all.x=TRUE, all.y=TRUE)
Meth_count_Gene_fold_NBW_LBW_1_emth<-subset(Meth_count_Gene_fold_NBW_LBW, LBW_Meth_sites>0 | NBW_Meth_sites>0)
Meth_count_Gene_fold_NBW_LBW_5_emth<-subset(Meth_count_Gene_fold_NBW_LBW, LBW_Meth_sites>4 | NBW_Meth_sites>4)
Meth_count_Gene_fold_NBW_LBW_5_emth_comp_data<-na.omit(Meth_count_Gene_fold_NBW_LBW_5_emth)
write.table(Meth_count_Gene_fold_NBW_LBW_5_emth_comp_data, file="~/Desktop/Christa_ANders/R_calc/To_Mandy/Meth_count_Gene_fold_NBW_LBW_5_emth_comp_data", sep="\t", row.names = FALSE)
write.table(Meth_count_Gene_fold_NBW_LBW_5_emth, file="~/Desktop/Christa_ANders/R_calc/To_Mandy/Meth_count_Gene_fold_NBW_LBW_5_emth.txt", sep="\t", row.names = FALSE)


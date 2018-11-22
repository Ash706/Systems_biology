#### Generate files for Kiwi using Piano

DE_genes_palm_controls_EdgeR <- read_csv("Z:/Ashfaq_Ali/GSM/Saeed_pipeline/Codes/MATLAB/New_files/DE-genes-palm-controls-EdgeR.csv")
DE_genes_palm_controls_EdgeR<-DE_genes_palm_controls_EdgeR[!(is.na(DE_genes_palm_controls_EdgeR$`P-value`)),]
rownames(DE_genes_palm_controls_EdgeR)<-DE_genes_palm_controls_EdgeR$`ENSG ID`
HMR2_gm_GSC <- read_delim("Z:/Ashfaq_Ali/GSM/Saeed_pipeline/Codes/MATLAB/KIWI/kiwi/data/HMR2_gm_GSC.sif",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)


HMR2_gm_GSC<-HMR2_gm_GSC[,c(3,1)]

myPval<- DE_genes_palm_controls_EdgeR$`P-value`
names(myPval)<-DE_genes_palm_controls_EdgeR$`ENSG ID`
myFC<-DE_genes_palm_controls_EdgeR$`log2 FC`
names(myFC)<-DE_genes_palm_controls_EdgeR$`ENSG ID`

myGsc<-loadGSC(HMR2_gm_GSC)
gene_stats<-cbind(myPval,myFC)
rownames(gene_stats)<-DE_genes_palm_controls_EdgeR$`ENSG ID`
reporter<-runGSA(myPval,myFC, gsc=myGsc, geneSetStat="reporter", signifMethod="nullDist", nPerm=1000, gsSizeLim=c(5,100))

nw <- networkPlot(reporter, class="distinct", direction="both", significance=0.001, label="names")

gsaResTab <- GSAsummaryTable(reporter)
grep("p \\(",colnames(gsaResTab),value=T)
grep("p \\(",colnames(gsaResTab))
ii <- which(gsaResTab[,10]<0.0001)
gsaResTab$Name[ii]
minPval <- apply(gsaResTab[,c(4,7,10,14,18)],1,min,na.rm=TRUE)
minPval
ii <- which(minPval<0.0001)
gsaResTabSign <- gsaResTab[ii,c(1,4,7,10,14,18)]
gsaResTabSign[1:20,]

myTval <- pfc$resTable[["aerobic_Clim - anaerobic_Clim"]]$t
names(myTval) <- pfc$resTable[["aerobic_Clim - anaerobic_Clim"]]$ProbesetID

# gsaRes1 <- runGSA(myTval,geneSetStat="mean",gsc=myGsc, nPerm=1000,gsSizeLim=c(10,800))
# gsaRes2 <- runGSA(myTval,geneSetStat="median",gsc=myGsc, nPerm=1000,gsSizeLim=c(10,800))
# gsaRes3 <- runGSA(myTval,geneSetStat="sum",gsc=myGsc,nPerm=1000,gsSizeLim=c(10,800))
# gsaRes4 <- runGSA(myTval,geneSetStat="maxmean",gsc=myGsc,nPerm=1000,gsSizeLim=c(10,800))
gsaRes5 <- runGSA(myPval,myFC,geneSetStat="fisher",gsc=myGsc,nPerm=1000,gsSizeLim=c(10,800))
gsaRes6 <- runGSA(myPval,myFC,geneSetStat="stouffer",gsc=myGsc,nPerm=1000,gsSizeLim=c(10,800))
gsaRes7 <- runGSA(myPval,myFC,geneSetStat="tailStrength",gsc=myGsc,nPerm=1000,gsSizeLim=c(10,800))

resList <- list(reporter,gsaRes5,gsaRes6,gsaRes7)
names(resList) <- c("reporter","fisher", "stouffer","tailStrength")
ch <- consensusHeatmap(resList,cutoff=30,method="mean")

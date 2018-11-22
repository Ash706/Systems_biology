library(tidyr)
library(dplyr)

Gena_pathway_association_HMR2<-Gene_pathways_HMR2 %>% 
  rename(x=`GENE ASSOCIATION`) %>% # rename because separate_rows() seems to not like columns with spaces
  separate_rows(x, sep=";") %>% 
  rename(`GENE ASSOCIATION` =x) 
write.csv(Gena_pathway_association_HMR2, "Z:/Ashfaq_Ali/GSM/Analyses_pipeline/Base_files/HMR2_Generic/Gene_pathway_HMR2_generic.csv", row.names = FALSE)

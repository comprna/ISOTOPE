args = commandArgs(trailingOnly=TRUE)

#Get the gene name assoicated to the ensembl_gene_ids
library(biomaRt)
cat("Starting execution...\n")
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

#Load the list of genes
# listAttributes(ensembl)
# file <- read.table(file="/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v5/neoskipping_junctions_filtered_mut.tab",sep="\t",header = TRUE)
file <- read.table(file=args[1],sep="\t",header = TRUE)
genes <- unique(file$Gene_id)
genes_associated_names <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),filters = 'ensembl_gene_id', values = genes, mart = ensembl)

#Associate the new gene names to the previous list
file_merged <- merge(file,genes_associated_names,by.x="Gene_id",by.y="ensembl_gene_id",all.x=TRUE)

#Save the new file
# write.table(file_merged,file="/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v5/neoskipping_junctions_filtered_mut2.tab",quote=FALSE,row.names = FALSE,sep="\t")
write.table(file_merged,file=args[2],quote=FALSE,row.names = FALSE,sep="\t")
cat(paste0("Done. Created ",args[2],".\n"))

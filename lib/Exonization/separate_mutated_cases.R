args = commandArgs(trailingOnly=TRUE)

####Separate the mutated cases from the non-mutated####
#Load the file from check_mutations_nearby.py 
# exonizations <- read.table(file="/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/exonizations_by_sample_coverage_mut.tab",sep="\t",header = TRUE)
exonizations <- read.table(file=args[1],sep="\t",header = TRUE)
#Get the ones that has a mutation an significant p_value from coverageBed
exonizations_mutated <- exonizations[which(exonizations$mut_coincidence=="True" & exonizations$pvalue<=0.05),]
#Get also the lines without mutation but significant p_vaue
exonizations_not_mutated <- exonizations[which(exonizations$mut_coincidence=="False" & exonizations$pvalue<=0.05),]
#Save both files
# write.table(exonizations_mutated,file="/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/mutated_exonizations.tab",quote=FALSE,sep="\t",row.names=FALSE)
# write.table(exonizations_not_mutated,file="/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/non_mutated_exonizations.tab",quote=FALSE,sep="\t",row.names=FALSE)
write.table(exonizations_mutated,file=args[2],quote=FALSE,sep="\t",row.names=FALSE)
write.table(exonizations_not_mutated,file=args[3],quote=FALSE,sep="\t",row.names=FALSE)
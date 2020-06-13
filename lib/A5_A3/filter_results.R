args = commandArgs(trailingOnly=TRUE)
####Filter the A5_A3 junctions with peptide change####

file <- read.table(file=args[1],sep="\t",header = TRUE)
#Associate the Ensembl id
#Load the file with the associated Gene names
conversion <- read.table(file=args[2],sep="\t",header = TRUE)
merge_file <- merge(file,conversion,by.x="Gene",by.y="Gene.stable.ID",all.x=TRUE)
#How many above TPM > 1?
file_f <- merge_file[which(as.numeric(merge_file$Transcript_TPM)>1),]
#How many junctions giving rise to exons larger than 500nt?
file_f2 <- file_f[which(as.numeric(file_f$New_Exon_length)<=500),]

#Save this file
write.table(file_f2,file=args[3],sep="\t",quote = FALSE,row.names = FALSE)

#How many peptide changes without NMD?
table(file_f2$Peptide_change,file_f2$NMD)
#Get only the cases in which there is a peptide change without NMD
file_f3 <- file_f2[which(file_f2$Peptide_change=="True" & file_f2$NMD=="False"),]

#Save this file (for running run_NetMHC.py)
write.table(file_f3,file=args[4],sep="\t",quote=FALSE,row.names = FALSE)

cat("How many falling on repeated regions?")
table(file_f3$Repeats)
length(which(file_f3$Repeats!="No repeat"))

cat("How many with associated mutations?")
table(file_f3$mut_coincidence)

cat("How many gets stalled?")
table(file_f3$Stalling)

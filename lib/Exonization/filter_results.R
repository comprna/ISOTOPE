args = commandArgs(trailingOnly=TRUE)
####Filter the exonizations with peptide change ####

junctions <- read.table(file=args[1],sep="\t",header = TRUE)
#Get transcripts with a minimum expresion TPM > 1
junctions_f <- junctions[which(junctions$Transcript_TPM>1),]
#Get only the cases in which there is significance (All cases are singificant, this has been prefiltered before)
junctions_f2 <- junctions_f[which(junctions_f$pvalue<=0.05),]
#Save this file
write.table(junctions_f2,file=args[2],sep="\t",quote = FALSE,row.names = FALSE)

cat("How many exons?")
length(unique(junctions_f2$New_exon))
cat("How many samples?")
length(unique(junctions_f2$Sample_id))
cat("How many peptide changes without NMD?")
table(junctions_f2$Peptide_change,junctions_f2$NMD)
#Get only the cases in which there is a peptide change without NMD
junctions_f3 <- junctions_f2[which(junctions_f2$Peptide_change=="True" & junctions_f2$NMD=="False"),]
junctions_f3$Index <- seq(1:nrow(junctions_f3))

#Save this file (for running run_NetMHC.py)
write.table(junctions_f3,file=args[3],sep="\t",quote = FALSE,row.names = FALSE)

args = commandArgs(trailingOnly=TRUE)
####Filter the neoskippings with peptide change ####

junctions <- read.table(file=args[1],sep="\t",header = TRUE)
#Get transcripts with a minimum expresion TPM > 1
junctions_f <- junctions[which(junctions$Transcript_TPM>1),]
#Save this file
write.table(junctions_f,file=args[2],sep="\t",quote = FALSE,row.names = FALSE)


cat("How many exons?")
length(unique(junctions_f$Neoskipping_junction))
cat("How many samples?")
length(unique(junctions_f$Sample_id))
cat("How many peptide changes without NMD?")
table(junctions_f$Peptide_change,junctions_f$NMD)
#Get only the cases in which there is a peptide change without NMD
junctions_f2 <- junctions_f[which(junctions_f$Peptide_change=="True" & junctions_f$NMD=="False"),]
junctions_f2$Index <- seq(1:nrow(junctions_f2))

#Save this file (for running run_NetMHC.py)
write.table(junctions_f2,file=args[3],sep="\t",quote = FALSE,row.names = FALSE)

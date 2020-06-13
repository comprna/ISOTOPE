args = commandArgs(trailingOnly=TRUE)

#Extract the coordinates from the RI events:
input_file <- as.data.frame(read.table(file=args[1],sep="\t",header = TRUE))

#Remove all the introns == 0
mean <- apply(input_file,1,function(x)mean(x, na.rm = TRUE))
input_file_f <- input_file[which(mean>=0.1),]
input_file_f$Row.names <- rownames(input_file_f)
#Reformat the ids for the U12 introns
u12 <- grepl("U12",as.character(input_file_f$Row.names))
u12ids <- as.character(input_file_f[u12,]$Row.names)
u12ids2 <- gsub("[.]",":",u12ids)
input_file_f[u12,ncol(input_file_f)] <- u12ids2
#Reformat the ids also for the non U12 introns
nonu12ids <- as.character(input_file_f[!u12,]$Row.names)
nonu12ids2 <- sapply(nonu12ids, function(x)paste0("chr",strsplit(x,":")[[1]][1],":",strsplit(x,":")[[1]][2],"(+):",strsplit(x,":")[[1]][3]))
input_file_f[!u12,ncol(input_file_f)] <- nonu12ids2
Index <- as.character(input_file_f$Row.names)
#Extract the cordinates
chr <- unlist(lapply(as.character(input_file_f$Row.names),function(x)strsplit(x,":")[[1]][1]))
start <- as.numeric(unlist(lapply(unlist(lapply(as.character(input_file_f$Row.names),function(x)strsplit(x,":")[[1]][2])),function(x)strsplit(x,"-")[[1]][1])))
end <- as.numeric(unlist(lapply(unlist(lapply(unlist(lapply(as.character(input_file_f$Row.names),function(x)strsplit(x,":")[[1]][2])),function(x)strsplit(x,"[(]")[[1]][1])),function(x)strsplit(x,"-")[[1]][2])))
strand <- unlist(lapply(as.character(input_file_f$Row.names),function(x)strsplit(x,":")[[1]][2]))
strand <- substr(strand,nchar(strand)-1,nchar(strand)-1)
merge_table_final <- cbind(Index,chr,strand,start,end,input_file_f[,-ncol(input_file_f)])

#Save this file
write.table(merge_table_final,file=args[2],sep="\t",quote=FALSE,row.names = FALSE)
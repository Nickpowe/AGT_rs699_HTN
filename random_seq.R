setwd("/filepath/AGT/RNAduplex")
library(stringi)
library(tidyverse)

mir122_5p <- "uggagugugacaaugguguuug"
mir122_5p_rev_comp <- "caaacaccattgtcacactcca"

rs699_wt_t <- "GTAATTCTAGGAGCTCTGACAGGATGGAAGACTGGCTGCTCCCTGATGGGAGCCAGTGTGGACAGCACCCTGGCTTTCGTTCTAGAGTCGGG"
rs699_var_c <- "GTAATTCTAGGAGCTCTGACAGGATGGAAGACTGGCTGCTCCCTGACGGGAGCCAGTGTGGACAGCACCCTGGCTTTCGTTCTAGAGTCGGG"
perf_match <- "GTAATTCTAGGAGCTCTGACAGGATGGAAGACTGGCcaaacaccattgtcacactccaTGGACAGCACCCTGGCTTTCGTTCTAGAGTCGGG"

x <- data.frame()
for(n in 1:1000){
  y1 <- stringi::stri_rand_shuffle(rs699_var_c) 
  y <- data.frame(sequence=c(y1,mir122_5p))
  x <- rbind(x,y) 
}

x1 <- data.frame()
for(n in 1:1000){
  y1 <- stringi::stri_rand_shuffle(rs699_wt_t) 
  y <- data.frame(sequence=c(y1,mir122_5p))
  x1 <- rbind(x1,y) 
}

x <- rbind(x,x1)

write.table(x, file="rand_sequences.txt",row.names = F,col.names = F,quote = F)
write.table(data.frame(c(perf_match,mir122_5p)), file="perf_sequences.txt",row.names = F,col.names = F,quote = F)

system("RNAduplex <rand_sequences.txt > out.txt")

out <- read.csv("/filepath/AGT/RNAduplex/out.txt", header=FALSE)

out[4:5] <- str_split_fixed(out$V3,"\\(",2)
out[,5] <- gsub("\\)","",out[,5])

summary(as.numeric(out$V5))

system("RNAduplex <perf_sequences.txt > out_perf.txt")

-39.3--7.1
1/32.2

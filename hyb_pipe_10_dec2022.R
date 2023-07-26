arg = commandArgs(trailingOnly=TRUE) #this allows arguments from linux command to be passed to R: "Rscript hyb_pipe.R arg1 arg2 arg3"

a1 <- arg[1] #should be filepath to the hyb document
ID <- arg[2] #should be output filename identifier

setwd("/filepath/")

library(stringi)
library(httr)
library(jsonlite)
library(xml2)

library(tidyverse)
library(ensembldb)
library(AnnotationHub)

#######################################################
## get ensembl database loaded and prepared to use
#######################################################

dbs_avail <- query(AnnotationHub(), "EnsDb.Hsapiens") ## use this to find newer ensembl data
db_latest <- rownames(as.data.frame(dbs_avail@.db_uid))[nrow(as.data.frame(dbs_avail@.db_uid))]
db_to_use <- "AH98047"
if(db_to_use!=db_latest){warning("not useing latest ensb database")}

supportedFilters() ## use this to see what can be filtered on

edb <- AnnotationHub()[[db_to_use]]
rm(db_latest,db_to_use,dbs_avail)

#######################################################
## read in data
#######################################################

b <- read.delim(a1, header=FALSE) #loads the .hyb file into R and names it "b"

#filter to mi_m and m_mi only
b1 <- dplyr::filter(b, grepl("hsa-", b$V4))
b2 <- dplyr::filter(b1, grepl("ENSG", b1$V10))
b2$V16 <- "mir_first"
c1 <- dplyr::filter(b, grepl("ENSG", b$V4))
c2 <- dplyr::filter(c1, grepl("hsa-", c1$V10))
c2$V16 <- "m_first"
c3 <- c2[c(1,2,3,10,11,12,13,14,15,4,5,6,7,8,9,16)]#getting all mRNAs and microRNAs in same columns
colnames(c3)[4:15] <- c("V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15")
d <- rbind(b2,c3)
rm(b1, b2, c1, c2, c3)

#splitting gene info into columns
d1 <- str_split_fixed(d$V10, "\\|", 4) #splitting text for later naming hybrids
d2 <- merge(d, d1, by=0, all=T) #merging on R's 'column.names' aka by=0
rownames(d2) <- d2$Row.names; d2$Row.names <- NULL #puts the column numbers (which were put into a new column) back into R's 'column.names'
rm(d, d1)

#making hybrid name column
d3 <- add_column(d2, Hybrid=0)
d3$Hybrid <- stri_join(d3$V4.x, d3$V3.y, sep=" ") #joins columns to make hybrid names
names(d3) <- c("Unique Assigned Sequence Number", 
               "Sequence", "Hybridization Delta G", 
               "microRNA reference Header", "microRNA start position of sequence", 
               "microRNA end position of sequence", "reference microRNA start coordinates", 
               "reference microRNA end coordinates", "microRNA mapping score", "Gene reference Header", "Gene start position of sequence",
               "Gene end position of sequence", "reference Gene start coordinates", "reference Gene end coordinates", "Gene mapping score",
               "directionality", "ENSG", "ENST", "Gene", "Strand", "Hybrid")
rm(d2)



#remove UMI duplicates
d3$umi <- substr(d3$Sequence,1,17)
d3 <- d3[order(d3$`Gene mapping score`),] # this should result in choosing the one with the lower probability that it was misaligned to reference.
rownames(d3) <- NULL
d3 <- d3 %>% distinct(umi, .keep_all = T) # CAUTION - SKIP IF UMI not in first part of read

#counting duplicate hybrids
e <- paste(d3$Hybrid)
e1 <- data.frame(table(e)) #the table function counts duplicates
colnames(e1) <- c("Hybrid", "Total Unique Sequence Reads per Hybrid") #renaming the column headers
filename <- paste0(ID,"_hybrid_counts",".csv")
write.csv(e1, file=filename, row.names=F, quote=F)
rm(e)
d3 <- merge(d3, e1, by="Hybrid", all = F)
d3 <- d3[c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,1,22)]
rm(e1)
e1 <- read.csv(paste0(ID,"_hybrid_counts",".csv"))

#loop for transcript genome conversion into columns
 
d4 <- d3[rowSums(is.na(d3))>0,]
print(paste(nrow(d4),"rows were NA prior to genome fetch loop"),sep=" ")
print("start loop for transcript to genome conversion into columns")

d3$chrom <- NA
d3$start1 <- NA
d3$end1 <- NA
d3$start2 <- NA
d3$end2 <- NA
d3$start3 <- NA
d3$end3 <- NA

###create GRanges object ###
print("creating GRanges object")
rng_tx <- IRanges(start = d3$`reference Gene start coordinates`, end = d3$`reference Gene end coordinates`, names = d3$ENST)
rng_gnm <- transcriptToGenome(rng_tx, edb)

######## extracting genomic coords ######
print("running loop")
for (n in 1:nrow(d3)){
  print(n)
  if(length(rng_gnm@listData[[n]]@ranges)==0){next}
  
  d3$chrom[n] <- as.character(rng_gnm@listData[[n]]@seqnames[1]@values)
  d3$start1[n] <- rng_gnm@listData[[n]]@ranges[1]@start
  d3$end1[n] <- as.data.frame(rng_gnm@listData[[n]]@ranges[1])$end
  
  if(length(rng_gnm@listData[[n]]@ranges)>1){
    d3$start2[n] <- rng_gnm@listData[[n]]@ranges[2]@start
    d3$end2[n] <- as.data.frame(rng_gnm@listData[[n]]@ranges[2])$end
  } else {next}
  
  if(length(rng_gnm@listData[[n]]@ranges)>2){
    d3$start3[n] <- rng_gnm@listData[[n]]@ranges[3]@start
    d3$end3[n] <- as.data.frame(rng_gnm@listData[[n]]@ranges[3])$end
  }
}

# checkpoint
write.table(d3, paste0("checkpoint1_",ID,".csv"))

#### finding the ones that didn't match the ensdb #####
d4 <- subset(d3, is.na(d3$start1))

#assign different ENST based on FIRST ENST in the ensdb
for(n in 1:nrow(d4)){
  if(nrow(d4)==0){break}
  x <- transcriptsBy(edb, by = "gene", filter = GeneIdFilter(d4$ENSG[n]))
  d4$ENST[n] <-  as.data.frame(x@unlistData)[1,6]
}

d5 <- subset(d4, is.na(d4$ENST))
for(n in 1:nrow(d5)){
  if(nrow(d5)==0){break}
  x <- transcriptsBy(edb, by = "gene", filter = SymbolFilter(d5$Gene[n]))
  d5$ENST[n] <-  as.data.frame(x@unlistData)[1,6]
}

#subset to those that worked and find ranges
d6 <- rbind(subset(d5, !is.na(d5$ENST)),subset(d4, !is.na(d4$ENST)))
rng_tx <- IRanges(start = d6$`reference Gene start coordinates`, end = d6$`reference Gene end coordinates`, names = d6$ENST)
if(length(rng_tx)>1){rng_gnm <- transcriptToGenome(rng_tx, edb)}

for (n in 1:nrow(d6)){
  if(nrow(d6)==0){break}
 print(n)
  if(length(rng_gnm@listData[[n]]@ranges)==0){next}
  
  d6$chrom[n] <- as.character(rng_gnm@listData[[n]]@seqnames[1]@values)
  d6$start1[n] <- rng_gnm@listData[[n]]@ranges[1]@start
  d6$end1[n] <- as.data.frame(rng_gnm@listData[[n]]@ranges[1])$end
  
  if(length(rng_gnm@listData[[n]]@ranges)>1){
    d6$start2[n] <- rng_gnm@listData[[n]]@ranges[2]@start
    d6$end2[n] <- as.data.frame(rng_gnm@listData[[n]]@ranges[2])$end
  } else {next}
  
  if(length(rng_gnm@listData[[n]]@ranges)>2){
    d6$start3[n] <- rng_gnm@listData[[n]]@ranges[3]@start
    d6$end3[n] <- as.data.frame(rng_gnm@listData[[n]]@ranges[3])$end
  }
}

d7 <- rbind(subset(d3, !is.na(d3$start1)),subset(d6, !is.na(d6$start1)))

#subset to those that didnt work and repeat using SECOND ENST option in ensdb
d8 <- subset(d6, is.na(d6$start1))

for(n in 1:nrow(d8)){
  if(nrow(d8)==0){break}
  x <- transcriptsBy(edb, by = "gene", filter = GeneIdFilter(d8$ENSG[n]))
  d8$ENST[n] <-  as.data.frame(x@unlistData)[2,6]
}

d5 <- subset(d8, is.na(d8$ENST))
for(n in 1:nrow(d5)){
  if(nrow(d5)==0){break}
  x <- transcriptsBy(edb, by = "gene", filter = SymbolFilter(d5$Gene[n]))
  d5$ENST[n] <-  as.data.frame(x@unlistData)[2,6]
}

d6 <- rbind(subset(d5, !is.na(d5$ENST)),subset(d8, !is.na(d8$ENST)))
rng_tx <- IRanges(start = d6$`reference Gene start coordinates`, end = d6$`reference Gene end coordinates`, names = d6$ENST)
if(length(rng_tx)>1){rng_gnm <- transcriptToGenome(rng_tx, edb)}

for (n in 1:nrow(d6)){
  if(nrow(d6)==0){break}
  if(length(rng_gnm@listData[[n]]@ranges)==0){next}
  
  d6$chrom[n] <- as.character(rng_gnm@listData[[n]]@seqnames[1]@values)
  d6$start1[n] <- rng_gnm@listData[[n]]@ranges[1]@start
  d6$end1[n] <- as.data.frame(rng_gnm@listData[[n]]@ranges[1])$end
  
  if(length(rng_gnm@listData[[n]]@ranges)>1){
    d6$start2[n] <- rng_gnm@listData[[n]]@ranges[2]@start
    d6$end2[n] <- as.data.frame(rng_gnm@listData[[n]]@ranges[2])$end
  } else {next}
  
  if(length(rng_gnm@listData[[n]]@ranges)>2){
    d6$start3[n] <- rng_gnm@listData[[n]]@ranges[3]@start
    d6$end3[n] <- as.data.frame(rng_gnm@listData[[n]]@ranges[3])$end
  }
}
d7 <- rbind(d7,subset(d6, !is.na(d6$start1)))

#subset to those that didnt work and repeat using THIRD ENST option in ensdb
d8 <- subset(d6, is.na(d6$start1))

for(n in 1:nrow(d8)){
  if(nrow(d8)==0){break}
  x <- transcriptsBy(edb, by = "gene", filter = GeneIdFilter(d8$ENSG[n]))
  d8$ENST[n] <-  as.data.frame(x@unlistData)[3,6]
}

d5 <- subset(d8, is.na(d8$ENST))
for(n in 1:nrow(d5)){
  if(nrow(d5)==0){break}
  x <- transcriptsBy(edb, by = "gene", filter = SymbolFilter(d5$Gene[n]))
  d5$ENST[n] <-  as.data.frame(x@unlistData)[3,6]
}

d6 <- rbind(subset(d5, !is.na(d5$ENST)),subset(d8, !is.na(d8$ENST)))
rng_tx <- IRanges(start = d6$`reference Gene start coordinates`, end = d6$`reference Gene end coordinates`, names = d6$ENST)
if(length(rng_tx)>1){rng_gnm <- transcriptToGenome(rng_tx, edb)}

for (n in 1:nrow(d6)){
  if(nrow(d6)==0){break}
  if(length(rng_gnm@listData[[n]]@ranges)==0){next}
  
  d6$chrom[n] <- as.character(rng_gnm@listData[[n]]@seqnames[1]@values)
  d6$start1[n] <- rng_gnm@listData[[n]]@ranges[1]@start
  d6$end1[n] <- as.data.frame(rng_gnm@listData[[n]]@ranges[1])$end
  
  if(length(rng_gnm@listData[[n]]@ranges)>1){
    d6$start2[n] <- rng_gnm@listData[[n]]@ranges[2]@start
    d6$end2[n] <- as.data.frame(rng_gnm@listData[[n]]@ranges[2])$end
  } else {next}
  
  if(length(rng_gnm@listData[[n]]@ranges)>2){
    d6$start3[n] <- rng_gnm@listData[[n]]@ranges[3]@start
    d6$end3[n] <- as.data.frame(rng_gnm@listData[[n]]@ranges[3])$end
  }
}

d7 <- rbind(d7,subset(d6, !is.na(d6$start1)))

# ###use REST API to try and find the missings
# d4 <- subset(d3, is.na(d3$start1))
# d8 <- subset(d4, d4$`Unique Assigned Sequence Number` %in% subset(d6, is.na(d6$start1))$`Unique Assigned Sequence Number`)
# #d8 <- subset(d6, is.na(d6$start1))
# 
# for (val in 1:nrow(d8)){
#   if(nrow(d8)==0){break}
#   skip_to_next <- FALSE
#   z <- d8[val,13]
#   y <- d8[val,14]
#   x <- stringr::str_extract(string = d8[val,18], pattern = "ENST[0-99999999999]+")
#   w <- content(GET(paste("https://rest.ensembl.org/map/cdna/",x,"/",z,"..",y,"?content-type=application/json",sep = "")))
#   tryCatch(w <- as.data.frame(w), error=function(e) { skip_to_next <<- TRUE })
#   if(skip_to_next) { next }
#   d8[val,23] <- as.character(w[1,"mappings.seq_region_name"])
#   d8[val,24] <- w[1,"mappings.start"]
#   d8[val,25] <- w[1,"mappings.end"]
# }
# 
# vec1 <- c(seq(1,23),"X","Y")
# 
# d8 <- subset(d8, d8$chrom %in% vec1)
# 
# d7 <- rbind(d8, d7)

#### collapsing intron-spanners ####
f1 <- subset(d7, d7$start1>0)
f2 <- subset(d7, d7$start2>0)
f3 <- subset(d7, d7$start3>0)

f1 <- f1[-c(26:29)]

f2 <- f2[-c(24:25)]
f2 <- f2[-c(26:27)]

f3 <- f3[-c(24:27)]

colnames(f1)[24:25] <- c("start","end")
colnames(f2)[24:25] <- c("start","end")
colnames(f3)[24:25] <- c("start","end")

f4 <- rbind(f1,f2,f3)

print("writing file Hybrid Full ...")
filename <- paste(ID,"_hybrid_full", ".csv", sep="")
write.csv(f4, file=filename,row.names=F)

#make bed file
print("start make bed")

vec <- c(paste0("chr",seq(1,23)),"chrX","chrY")

bed <- f4[c(23,24,25,21)]
names(bed) <- c("chrom","chromStart","chromEnd","name")
bed$chrom <- paste('chr',bed$chrom,sep="")

t <- as.data.frame(table(bed$chrom))

bed <- subset(bed, bed$chrom %in% vec)

bed <- bed[order(bed$chrom),]

bed$name <- gsub(" ","_",bed$name)

filename <- paste0(ID,".bed")
write.table(bed, file = filename, row.names=F, quote=F, col.names=F, sep = "\t") #this is in 1-based since coming from Hyb and ensembl
#!!!!!!!!! bedtools genomecov is expecting 0-based coordinates 
#!!!!!!!!! - as of 11/3/2022 bedtools genomecov was fed 1-based coordinates
system(paste0("sh bedgraph.sh ",ID,".bed ",ID,".bdg")) 

system(paste0("sed -i '1i track type=bed name='",ID,"_bed' description='",ID,"_bed'' ",filename))
system(paste0("sed -i '1i track type=bedGraph name='",ID,"_BedGraph' description='",ID,"_bdg'' ",paste0(ID,".bdg")))

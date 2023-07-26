setwd("/filepath/GTEx/plink/agt")

####################################################################################
# run feature counts with "-J" option to count exon boundaries from bam files
####################################################################################

# run "/filepath/GTEx/scripts/bam_feature_liver.sh"
# subset to just the header info: from /filepath/GTEx/feature_counts_rna_seq_liver run "head -1 liver_read_matrix2.txt.jcounts > jcounts_header.txt"
# subset to just the agt info: from /filepath/GTEx/feature_counts_rna_seq_liver run "grep 'ENSG00000135744\|PrimaryGene' liver_read_matrix2.txt.jcounts > jcounts_agt.txt"

####################################################################################
# pull in header and agt files
####################################################################################
y1 <- read.table(file="/filepath/GTEx/feature_counts_rna_seq_liver/jcounts_agt.txt", stringsAsFactors = FALSE)
y <- data.frame(t(y1[1,]))
y$X1 <- gsub("/filepath/GTEx/bam_rna_seq_liver/","",y$X1)
y[,2:4] <- str_split_fixed(y[,1],"\\-",3)
y$colnames <- paste0(y$V2,".",y$V3)
colnames(y1) <- y$colnames

jcounts <- select_if(y1, colnames(y1) %in% c(colnames(y1)[1:8],rownames(gtex_agt_snps)))

jcounts <- jcounts[-c(1),]

jcount_exon2 <- as.data.frame(t(subset(jcounts, jcounts$Site2_location.==230709995)))
colnames(jcount_exon2) <- jcount_exon2[4,]
jcount_exon2 <- jcount_exon2[-c(1:8),]

for(n in 1:ncol(jcount_exon2)){
  jcount_exon2[n] <- as.numeric(unlist(jcount_exon2[n]))
}

jcount_exon2$sum <- rowSums(jcount_exon2)
jcount_exon2$percent_alt <- 1-(jcount_exon2$`230706200`/jcount_exon2$sum)

####################################################################################
# pull in genotype
####################################################################################
gtex_agt_snps <- read.delim("/filepath/GTEx/plink/agt/gtex_agt_snps.raw", stringsAsFactors=FALSE)
gtex_agt_snps <- gtex_agt_snps[c(2,7,8)]
gtex_agt_snps$IID <- gsub("-",".",gtex_agt_snps$IID)
rownames(gtex_agt_snps) <- gtex_agt_snps$IID
gtex_agt_snps <- gtex_agt_snps[c(2,3)]
gtex_agt_snps <- 2-gtex_agt_snps

colnames(gtex_agt_snps) <- c("rs699_G","rs5051_T")

####################################################################################
# merge
####################################################################################
jc <- merge(gtex_agt_snps, jcount_exon2, by=0)

summary(lm(percent_alt ~ rs699_G, data=jc))

#library(rtracklayer)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)
library(rstatix)
library(ggtext)
#library(VariantAnnotation)
library(tidyverse)
library(ggrepel)
library(grid)
library(gridExtra)
library(MASS)
library(leaps)
library(tidyverse)
library(caret)
library(stringr)
library(dplyr)
library(ggplot2)
setwd("/filepath/GTEx/plink/agt")

######################################################################
### Genotype
######################################################################
gtex_agt_snps <- read.delim("/filepath/GTEx/plink/agt/gtex_agt_snps.raw", stringsAsFactors=FALSE)
gtex_agt_snps <- gtex_agt_snps[c(2,7,8)]
gtex_agt_snps$IID <- gsub("-",".",gtex_agt_snps$IID)
rownames(gtex_agt_snps) <- gtex_agt_snps$IID
gtex_agt_snps <- gtex_agt_snps[c(2,3)]
gtex_agt_snps <- 2-gtex_agt_snps

colnames(gtex_agt_snps) <- c("rs699_G","rs5051_T")


######################################################################
### Liver
######################################################################
liver_read_matrix <- read.delim("/filepath/GTEx/feature_counts_rna_seq_liver/liver_read_matrix.txt", comment.char="#")
a <- as.data.frame(colnames(liver_read_matrix))
a[,1] <- gsub("X.N.project.GTEx.bam_rna_seq_liver.GTEX.","GTEX-",a[,1])
a[,2:3] <- str_split_fixed(a[,1],"\\.",2)
colnames(liver_read_matrix) <- a[,2]
row.names(liver_read_matrix) <- liver_read_matrix$Geneid
liver_read_matrix <- liver_read_matrix[-c(1)]
liver_read_matrix <- liver_read_matrix[-c(1:4)]

lengths <- liver_read_matrix[1]

liver_read_matrix <- liver_read_matrix[-c(1)]

a <- as.data.frame(colnames(liver_read_matrix))
for(n in 1:nrow(a)){
  a[n,1] <- colnames(liver_read_matrix[n])
  a[n,2] <- colSums(liver_read_matrix[n])
}

liver_read_matrix_TPM <- liver_read_matrix
for(n in 1:nrow(a)){
  x <- as.character(a[n,1])
  liver_read_matrix_TPM[x] <- (liver_read_matrix[x]/(a[n,2]))*1000000
}

#write.table(liver_read_matrix_TPM,file = "/filepath/GTEx/final_data_frames/liver_read_matrix_TPM")
liver_read_matrix_TPM <- read.csv("/filepath/GTEx/final_data_frames/liver_read_matrix_TPM", sep="")

liver_read_matrix_logTPM <- log10(liver_read_matrix_TPM)

#write.table(liver_read_matrix_logTPM,file = "/filepath/GTEx/final_data_frames/liver_read_matrix_logTPM")
liver_read_matrix_logTPM <- read.csv("/filepath/GTEx/final_data_frames/liver_read_matrix_logTPM", sep="")

############ log ################
agt <- subset(liver_read_matrix_logTPM, grepl("ENSG00000135744",rownames(liver_read_matrix_logTPM)))

agt <- as.data.frame(t(agt))
colnames(agt) <- "AGT"

a <- merge(agt,gtex_agt_snps,by=0)
a$hap <- paste0(a$rs699,"+",a$rs5051)
a$hap[a$hap=="0+0"] <- 6
a$hap[a$hap=="0+2"] <- 1
a$hap[a$hap=="0+1"] <- 3
a$hap[a$hap=="1+0"] <- 8
a$hap[a$hap=="1+1"] <- 5
a$hap[a$hap=="1+2"] <- 2
a$hap[a$hap=="2+0"] <- 9
a$hap[a$hap=="2+1"] <- 7
a$hap[a$hap=="2+2"] <- 4

aa <- a %>% group_by(hap) %>% summarize(mean=mean(AGT),n=n())

liver_log <- ggplot()+
  geom_boxplot(data=a,mapping=aes(x=hap,y=AGT),outlier.shape = NA)+
  geom_point(data=a,mapping=aes(x=hap,y=AGT),alpha=0.6,col="black",size=1,position=position_jitter(h=0,w=.1))+
  geom_point(data=aa,mapping=aes(x=hap,y=mean),alpha=0.8,col="red",size=2,position="identity")+
  theme_classic()+theme(text = element_text(size = 6))

liver_log

############ linear ################
agt <- subset(liver_read_matrix_TPM, grepl("ENSG00000135744",rownames(liver_read_matrix_TPM)))

agt <- as.data.frame(t(agt))
colnames(agt) <- "AGT"

a <- merge(agt,gtex_agt_snps,by=0)
a$hap <- paste0(a$rs699,"+",a$rs5051)
a$hap[a$hap=="0+0"] <- 6
a$hap[a$hap=="0+2"] <- 1
a$hap[a$hap=="0+1"] <- 3
a$hap[a$hap=="1+0"] <- 8
a$hap[a$hap=="1+1"] <- 5
a$hap[a$hap=="1+2"] <- 2
a$hap[a$hap=="2+0"] <- 9
a$hap[a$hap=="2+1"] <- 7
a$hap[a$hap=="2+2"] <- 4

aa <- a %>% group_by(hap) %>% summarize(mean=mean(AGT),n=n())

liver_linear <- ggplot()+
  geom_boxplot(data=a,mapping=aes(x=hap,y=AGT),outlier.shape = NA)+
  geom_point(data=a,mapping=aes(x=hap,y=AGT),alpha=0.6,col="black",size=1,position=position_jitter(h=0,w=.1))+
  geom_point(data=aa,mapping=aes(x=hap,y=mean),alpha=0.8,col="red",size=2,position="identity")+
  theme_classic()+theme(text = element_text(size = 6))

liver_linear

plot1 <- ggarrange(liver_log, liver_linear)
plot1
ggsave("gtex_agt_rs699_rs5051_liver.pdf",plot = plot1,device = "pdf",width = 90, height = 30, units = "mm")

######################################################################
### Cerebellum
######################################################################
ids <- read.delim("/filepath/GTEx/txt_analysis_supplement/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt")
ids <- subset(ids, ids$SMTSD=="Brain - Cerebellum")

cerebellum_read_matrix <- read.delim("/filepath/GTEx/feature_counts_rna_seq_cerebellum/cerebellum_read_matrix.txt", comment.char="#")
a <- as.data.frame(colnames(cerebellum_read_matrix))

a[,2] <- gsub("X.N.project.GTEx.bam_rna_seq_cerebellum.GTEX.","GTEX-",a[,1])
a[,3] <- gsub(".Aligned.sortedByCoord.out.patched.md.bam","",a[,2])
a[,4] <- gsub("\\.","\\-",a[,3])
ab <- subset(a, a$V4 %in% c("Geneid","Length",ids$SAMPID))

cerebellum_read_matrix <- select_if(cerebellum_read_matrix, colnames(cerebellum_read_matrix) %in% ab[,1])

a <- as.data.frame(colnames(cerebellum_read_matrix))
a[,1] <- gsub("X.N.project.GTEx.bam_rna_seq_cerebellum.GTEX.","GTEX-",a[,1])
a[,2:3] <- str_split_fixed(a[,1],"\\.",2)
colnames(cerebellum_read_matrix) <- a[,2]
row.names(cerebellum_read_matrix) <- cerebellum_read_matrix$Geneid


lengths <- cerebellum_read_matrix["Length"]

cerebellum_read_matrix <- cerebellum_read_matrix[-c(1:2)]

a <- as.data.frame(colnames(cerebellum_read_matrix))
for(n in 1:nrow(a)){
  a[n,1] <- colnames(cerebellum_read_matrix[n])
  a[n,2] <- colSums(cerebellum_read_matrix[n])
}

cerebellum_read_matrix_TPM <- cerebellum_read_matrix
for(n in 1:nrow(a)){
  x <- as.character(a[n,1])
  cerebellum_read_matrix_TPM[x] <- (cerebellum_read_matrix[x]/(a[n,2]))*1000000
}

write.table(cerebellum_read_matrix_TPM,file = "/filepath/GTEx/final_data_frames/cerebellum_read_matrix_TPM")
#cerebellum_read_matrix_TPM <- read.csv("/filepath/GTEx/final_data_frames/cerebellum_read_matrix_TPM", sep="")

cerebellum_read_matrix_logTPM <- log10(cerebellum_read_matrix_TPM)

write.table(cerebellum_read_matrix_logTPM,file = "/filepath/GTEx/final_data_frames/cerebellum_read_matrix_logTPM")
#cerebellum_read_matrix_logTPM <- read.csv("/filepath/GTEx/final_data_frames/cerebellum_read_matrix_logTPM", sep="")

############ log ################
agt <- subset(cerebellum_read_matrix_logTPM, grepl("ENSG00000135744",rownames(cerebellum_read_matrix_logTPM)))

agt <- as.data.frame(t(agt))
colnames(agt) <- "AGT"
rownames(agt) <- gsub("-",".",rownames(agt))

a <- merge(agt,gtex_agt_snps,by=0)
a$hap <- paste0(a$rs699,"+",a$rs5051)
a$hap[a$hap=="0+0"] <- 6
a$hap[a$hap=="0+2"] <- 1
a$hap[a$hap=="0+1"] <- 3
a$hap[a$hap=="1+0"] <- 8
a$hap[a$hap=="1+1"] <- 5
a$hap[a$hap=="1+2"] <- 2
a$hap[a$hap=="2+0"] <- 9
a$hap[a$hap=="2+1"] <- 7
a$hap[a$hap=="2+2"] <- 4

a$AGT <- as.numeric(a$AGT)

aa <- a %>% group_by(hap) %>% summarize(mean=mean(AGT),n=n())

cerebellum_log <- ggplot()+
  geom_boxplot(data=a,mapping=aes(x=hap,y=AGT),outlier.shape = NA)+
  geom_point(data=a,mapping=aes(x=hap,y=AGT),alpha=0.6,col="black",size=1,position=position_jitter(h=0,w=.1))+
  geom_point(data=aa,mapping=aes(x=hap,y=mean),alpha=0.8,col="red",size=2,position="identity")+
  theme_classic()+theme(text = element_text(size = 6))

cerebellum_log

############ linear ################
agt <- subset(cerebellum_read_matrix_TPM, grepl("ENSG00000135744",rownames(cerebellum_read_matrix_TPM)))

agt <- as.data.frame(t(agt))
colnames(agt) <- "AGT"
rownames(agt) <- gsub("-",".",rownames(agt))

a <- merge(agt,gtex_agt_snps,by=0)
a$hap <- paste0(a$rs699,"+",a$rs5051)
a$hap[a$hap=="0+0"] <- 6
a$hap[a$hap=="0+2"] <- 1
a$hap[a$hap=="0+1"] <- 3
a$hap[a$hap=="1+0"] <- 8
a$hap[a$hap=="1+1"] <- 5
a$hap[a$hap=="1+2"] <- 2
a$hap[a$hap=="2+0"] <- 9
a$hap[a$hap=="2+1"] <- 7
a$hap[a$hap=="2+2"] <- 4

a$AGT <- as.numeric(a$AGT)

aa <- a %>% group_by(hap) %>% summarize(mean=mean(AGT),n=n())

cerebellum_linear <- ggplot()+
  geom_boxplot(data=a,mapping=aes(x=hap,y=AGT),outlier.shape = NA)+
  geom_point(data=a,mapping=aes(x=hap,y=AGT),alpha=0.6,col="black",size=1,position=position_jitter(h=0,w=.1))+
  geom_point(data=aa,mapping=aes(x=hap,y=mean),alpha=0.8,col="red",size=2,position="identity")+
  theme_classic()+theme(text = element_text(size = 6))

cerebellum_linear

plot1 <- ggarrange(cerebellum_log, cerebellum_linear)
plot1
ggsave("gtex_agt_rs699_rs5051_cerebellum.pdf",plot = plot1,device = "pdf",width = 90, height = 30, units = "mm")

######################################################################
### sigmoid
######################################################################
ids <- read.delim("/filepath/GTEx/txt_analysis_supplement/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt")
ids <- subset(ids, ids$SMTSD=="Colon - Sigmoid")

sigmoid_read_matrix <- read.delim("/filepath/GTEx/feature_counts_rna_seq_sigmoid/sigmoid_read_matrix.txt", comment.char="#")
a <- as.data.frame(colnames(sigmoid_read_matrix))

a[,2] <- gsub("X.N.project.GTEx.bam_rna_seq_sigmoid.GTEX.","GTEX-",a[,1])
a[,3] <- gsub(".Aligned.sortedByCoord.out.patched.md.bam","",a[,2])
a[,4] <- gsub("\\.","\\-",a[,3])
ab <- subset(a, a$V4 %in% c("Geneid","Length",ids$SAMPID))

sigmoid_read_matrix <- select_if(sigmoid_read_matrix, colnames(sigmoid_read_matrix) %in% ab[,1])

a <- as.data.frame(colnames(sigmoid_read_matrix))
a[,1] <- gsub("X.N.project.GTEx.bam_rna_seq_sigmoid.GTEX.","GTEX-",a[,1])
a[,2:3] <- str_split_fixed(a[,1],"\\.",2)
colnames(sigmoid_read_matrix) <- a[,2]
row.names(sigmoid_read_matrix) <- sigmoid_read_matrix$Geneid


lengths <- sigmoid_read_matrix["Length"]

sigmoid_read_matrix <- sigmoid_read_matrix[-c(1:2)]

a <- as.data.frame(colnames(sigmoid_read_matrix))
for(n in 1:nrow(a)){
  a[n,1] <- colnames(sigmoid_read_matrix[n])
  a[n,2] <- colSums(sigmoid_read_matrix[n])
}

sigmoid_read_matrix_TPM <- sigmoid_read_matrix
for(n in 1:nrow(a)){
  x <- as.character(a[n,1])
  sigmoid_read_matrix_TPM[x] <- (sigmoid_read_matrix[x]/(a[n,2]))*1000000
}

write.table(sigmoid_read_matrix_TPM,file = "/filepath/GTEx/final_data_frames/sigmoid_read_matrix_TPM")
#sigmoid_read_matrix_TPM <- read.csv("/filepath/GTEx/final_data_frames/sigmoid_read_matrix_TPM", sep="")

sigmoid_read_matrix_logTPM <- log10(sigmoid_read_matrix_TPM)

write.table(sigmoid_read_matrix_logTPM,file = "/filepath/GTEx/final_data_frames/sigmoid_read_matrix_logTPM")
#sigmoid_read_matrix_logTPM <- read.csv("/filepath/GTEx/final_data_frames/sigmoid_read_matrix_logTPM", sep="")

############ log ################
agt <- subset(sigmoid_read_matrix_logTPM, grepl("ENSG00000135744",rownames(sigmoid_read_matrix_logTPM)))

agt <- as.data.frame(t(agt))
colnames(agt) <- "AGT"
rownames(agt) <- gsub("-",".",rownames(agt))

a <- merge(agt,gtex_agt_snps,by=0)
a$hap <- paste0(a$rs699,"+",a$rs5051)
a$hap[a$hap=="0+0"] <- 6
a$hap[a$hap=="0+2"] <- 1
a$hap[a$hap=="0+1"] <- 3
a$hap[a$hap=="1+0"] <- 8
a$hap[a$hap=="1+1"] <- 5
a$hap[a$hap=="1+2"] <- 2
a$hap[a$hap=="2+0"] <- 9
a$hap[a$hap=="2+1"] <- 7
a$hap[a$hap=="2+2"] <- 4

a$AGT <- as.numeric(a$AGT)

aa <- a %>% group_by(hap) %>% summarize(mean=mean(AGT),n=n())

sigmoid_log <- ggplot()+
  geom_boxplot(data=a,mapping=aes(x=hap,y=AGT),outlier.shape = NA)+
  geom_point(data=a,mapping=aes(x=hap,y=AGT),alpha=0.6,col="black",size=1,position=position_jitter(h=0,w=.1))+
  geom_point(data=aa,mapping=aes(x=hap,y=mean),alpha=0.8,col="red",size=2,position="identity")+
  theme_classic()+theme(text = element_text(size = 6))

sigmoid_log

############ linear ################
agt <- subset(sigmoid_read_matrix_TPM, grepl("ENSG00000135744",rownames(sigmoid_read_matrix_TPM)))

agt <- as.data.frame(t(agt))
colnames(agt) <- "AGT"
rownames(agt) <- gsub("-",".",rownames(agt))

a <- merge(agt,gtex_agt_snps,by=0)
a$hap <- paste0(a$rs699,"+",a$rs5051)
a$hap[a$hap=="0+0"] <- 6
a$hap[a$hap=="0+2"] <- 1
a$hap[a$hap=="0+1"] <- 3
a$hap[a$hap=="1+0"] <- 8
a$hap[a$hap=="1+1"] <- 5
a$hap[a$hap=="1+2"] <- 2
a$hap[a$hap=="2+0"] <- 9
a$hap[a$hap=="2+1"] <- 7
a$hap[a$hap=="2+2"] <- 4

a$AGT <- as.numeric(a$AGT)

aa <- a %>% group_by(hap) %>% summarize(mean=mean(AGT),n=n())

sigmoid_linear <- ggplot()+
  geom_boxplot(data=a,mapping=aes(x=hap,y=AGT),outlier.shape = NA)+
  geom_point(data=a,mapping=aes(x=hap,y=AGT),alpha=0.6,col="black",size=1,position=position_jitter(h=0,w=.1))+
  geom_point(data=aa,mapping=aes(x=hap,y=mean),alpha=0.8,col="red",size=2,position="identity")+
  theme_classic()+theme(text = element_text(size = 6))

sigmoid_linear

plot1 <- ggarrange(sigmoid_log, sigmoid_linear)
plot1
ggsave("gtex_agt_rs699_rs5051_sigmoid.pdf",plot = plot1,device = "pdf",width = 90, height = 30, units = "mm")

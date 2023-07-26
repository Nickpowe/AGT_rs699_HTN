setwd("/filepath/UKB")
library(ggplot2)
library(dplyr)
library(ggpubr)

############################################################################################################################
### Pull in and name data pieces
############################################################################################################################

#system("sh /filepath/UKB/chr1_calls/agt_plink_filter.sh")
rs699 <- read.delim("/filepath/UKB/chr1_calls/rs699.raw")

#system("sh /filepath/UKB/chr1_imputed/agt_plink_filter.sh")
rs5051 <- read.delim("/filepath/UKB/chr1_imputed/rs5051.raw")

#system("./ukbconv ukb52321.enc_ukb r -s4080 -osys_bp_aut") 
sys_bp <- read.table("/filepath/UKB/sys_bp_aut.tab", header=TRUE, sep="\t")
sys_bp$sys_bp <- rowMeans(sys_bp[2:9], na.rm = TRUE)
colnames(sys_bp)[1] <- "IID"

#system("./ukbconv ukb52321.enc_ukb r -s4079 -odia_bp_aut")
dia_bp <- read.table("/filepath/UKB/dia_bp_aut.tab", header=TRUE, sep="\t")
dia_bp$dia_bp <- rowMeans(dia_bp[2:9], na.rm = TRUE)
colnames(dia_bp)[1] <- "IID"

#system("./ukbconv ukb52321.enc_ukb r -s2966 -oage_htn_dg")
age_htn_dg <- read.table("/filepath/UKB/age_htn_dg.tab", header=TRUE, sep="\t")
colnames(age_htn_dg)[1] <- "IID"
age_htn_dg[age_htn_dg=="-1"] <- NA
age_htn_dg[age_htn_dg=="-3"] <- NA
age_htn_dg$age_htn_dg <- rowMeans(age_htn_dg[2:5], na.rm = TRUE)

#system("./ukbconv ukb52321.enc_ukb r -s21022 -oage")
age <- read.table("/filepath/UKB/age.tab", header=TRUE, sep="\t")
colnames(age) <- c("IID","age")

#system("./ukbconv ukb52321.enc_ukb r -s21000 -oethnic_background")
eth_rac <- read.table("/filepath/UKB/ethnic_background.tab", header=TRUE, sep="\t")
lvl.1001 <- c(-3,-1,1,2,3,4,5,6,1001,1002,1003,2001,2002,2003,2004,3001,3002,3003,3004,4001,4002,4003)
lbl.1001 <- c("Prefer not to answer","Do not know","White","Mixed","Asian or Asian British","Black or Black British","Chinese","Other ethnic group","British","Irish","Any other white background","White and Black Caribbean","White and Black African","White and Asian","Any other mixed background","Indian","Pakistani","Bangladeshi","Any other Asian background","Caribbean","African","Any other Black background")
eth_rac$f.21000.0.0 <- ordered(eth_rac$f.21000.0.0, levels=lvl.1001, labels=lbl.1001)
eth_rac$f.21000.1.0 <- ordered(eth_rac$f.21000.1.0, levels=lvl.1001, labels=lbl.1001)
eth_rac$f.21000.2.0 <- ordered(eth_rac$f.21000.2.0, levels=lvl.1001, labels=lbl.1001)
rm(lvl.1001,lbl.1001)
colnames(eth_rac)[1] <- "IID"

eth_rac$col5 <- paste(eth_rac$f.21000.0.0, eth_rac$f.21000.1.0, eth_rac$f.21000.2.0, sep = ",")
tab <- data.frame(table(eth_rac$col5))
#write.csv(tab, file = "eth_rac_table.csv")
#manually decided among the three instances how to code ethnicity and race into broader groups
eth_rac_table <- read.csv("/filepath/UKB/eth_rac_table.csv")[c(2,5)]
colnames(eth_rac_table) <- c("col5","eth_rac")
eth_rac <- merge(eth_rac, eth_rac_table, by="col5")

rm(tab)

#system("./ukbconv ukb52321.enc_ukb r -s6177 -omed_for_bp_chol_db")
med_for_bp_chol_db_males <- read.table("/filepath/UKB/med_for_bp_chol_db.tab", header=TRUE, sep="\t")
colnames(med_for_bp_chol_db_males)[1] <- "IID"

med_for_bp_chol_db_males[2:13][med_for_bp_chol_db_males[2:13]!=2] <- 0
med_for_bp_chol_db_males[2:13][med_for_bp_chol_db_males[2:13]==2] <- 1

med_for_bp_chol_db_males$on_bp_med <- rowSums(med_for_bp_chol_db_males[2:13], na.rm = TRUE)
med_for_bp_chol_db_males$on_bp_med[med_for_bp_chol_db_males$on_bp_med>0] <- 1

#system("./ukbconv ukb52321.enc_ukb r -s6153 -omed_for_bp_chol_db_females")
med_for_bp_chol_db_females <- read.table("/filepath/UKB/med_for_bp_chol_db_females.tab", header=TRUE, sep="\t")
colnames(med_for_bp_chol_db_females)[1] <- "IID"

med_for_bp_chol_db_females[2:17][med_for_bp_chol_db_females[2:17]!=2] <- 0
med_for_bp_chol_db_females[2:17][med_for_bp_chol_db_females[2:17]==2] <- 1

med_for_bp_chol_db_females$on_bp_med <- rowSums(med_for_bp_chol_db_females[2:17], na.rm = TRUE)
med_for_bp_chol_db_females$on_bp_med[med_for_bp_chol_db_females$on_bp_med>0] <- 1

#combine males and females
med_for_bp_chol_db <- merge(med_for_bp_chol_db_females[c("IID","on_bp_med")],med_for_bp_chol_db_males[c("IID","on_bp_med")],by="IID")
med_for_bp_chol_db$on_bp_med <- rowSums(med_for_bp_chol_db[2:3])
med_for_bp_chol_db <- med_for_bp_chol_db[c(1,4)]

#system("./ukbconv ukb52321.enc_ukb r -s131286 -odate_htn_dg")
htn_dg <- read.table("/filepath/UKB/date_htn_dg.tab", header=TRUE, sep="\t")
htn_dg[htn_dg$f.131286.0.0 %in% c(1900,1901,1902,1903,2037)] <- NA #"Code has no event date","Code has event date before participant\'s date of birth","Code has event date matching participant\'s date of birth","Code has event date after participant\'s date of birth and falls in the same calendar year as date of birth","Code has event date in the future and is presumed to be a place-holder or other system default"
htn_dg$f.131286.0.0[!is.na(htn_dg$f.131286.0.0)] <- 1
htn_dg$f.131286.0.0[is.na(htn_dg$f.131286.0.0)] <- 0
htn_dg$f.131286.0.0 <- as.numeric(htn_dg$f.131286.0.0)
nrow(subset(htn_dg, htn_dg$f.131286.0.0==1))
colnames(htn_dg)[1] <- "IID"

htn_dg <- merge(htn_dg, age_htn_dg[c("IID","age_htn_dg")],by="IID")
htn_dg$age_htn_dg[!is.na(htn_dg$age_htn_dg)] <- 1
htn_dg$age_htn_dg[is.na(htn_dg$age_htn_dg)] <- 0

nrow(subset(htn_dg, htn_dg$f.131286.0.0 == 0 & htn_dg$age_htn_dg == 1))
htn_dg$htn_dg <- rowSums(htn_dg[2:3])
htn_dg$htn_dg[htn_dg$htn_dg>0] <- 1

#system("./ukbconv ukb52321.enc_ukb r -s21001 -obmi_kg_m2")
bmi <- read.table("/filepath/UKB/bmi_kg_m2.tab", header=TRUE, sep="\t")
bmi$bmi <- rowMeans(bmi[2:5],na.rm = TRUE)
colnames(bmi)[1] <- "IID"

#system("./ukbconv ukb52321.enc_ukb r -s20116 -osmk_stat")
smk_stat <- read.table("/filepath/UKB/smk_stat.tab", header=TRUE, sep="\t")
#Prefer not to answer=-3,Never=0,Previous=1,Current=2

smk_stat$f.20116.0.0[is.na(smk_stat$f.20116.0.0)] <- smk_stat$f.20116.1.0[is.na(smk_stat$f.20116.0.0)]
smk_stat$f.20116.0.0[is.na(smk_stat$f.20116.0.0)] <- smk_stat$f.20116.2.0[is.na(smk_stat$f.20116.0.0)]
smk_stat$f.20116.0.0[is.na(smk_stat$f.20116.0.0)] <- smk_stat$f.20116.3.0[is.na(smk_stat$f.20116.0.0)]

smk_stat <- smk_stat[1:2]
smk_stat$f.20116.0.0[smk_stat$f.20116.0.0==-3] <- NA
colnames(smk_stat) <- c("IID","smk_stat")

#system("./ukbconv ukb52321.enc_ukb r -s22032 -oactivity")
act <- read.table("/filepath/UKB/activity.tab", header=TRUE, sep="\t")
#0=low, 1=moderate, 2=high
colnames(act) <- c("IID","act_more")

#system("./ukbconv ukb52321.enc_ukb r -s1558 -oalc")
alc <- read.table("/filepath/UKB/alc.tab", header=TRUE, sep="\t")
#lvl.100402 <- c(-3,1,2,3,4,5,6)
#lbl.100402 <- c("Prefer not to answer","Daily or almost daily","Three or four times a week","Once or twice a week","One to three times a month","Special occasions only","Never")

alc$f.1558.0.0[is.na(alc$f.1558.0.0)] <- alc$f.1558.1.0[is.na(alc$f.1558.0.0)]
alc$f.1558.0.0[is.na(alc$f.1558.0.0)] <- alc$f.1558.2.0[is.na(alc$f.1558.0.0)]
alc$f.1558.0.0[is.na(alc$f.1558.0.0)] <- alc$f.1558.3.0[is.na(alc$f.1558.0.0)]

alc <- alc[1:2]
alc$f.1558.0.0[alc$f.1558.0.0==-3] <- NA
colnames(alc) <- c("IID","alc_less")

#system("./ukbconv ukb52321.enc_ukb r -s1478 -osalt_add")
salt_add <- read.table("/filepath/UKB/salt_add.tab", header=TRUE, sep="\t")
#lvl.100394 <- c(-3,1,2,3,4)
#lbl.100394 <- c("Prefer not to answer","Never/rarely","Sometimes","Usually","Always")

salt_add$f.1478.0.0[is.na(salt_add$f.1478.0.0)] <- salt_add$f.1478.1.0[is.na(salt_add$f.1478.0.0)]
salt_add$f.1478.0.0[is.na(salt_add$f.1478.0.0)] <- salt_add$f.1478.2.0[is.na(salt_add$f.1478.0.0)]
salt_add$f.1478.0.0[is.na(salt_add$f.1478.0.0)] <- salt_add$f.1478.3.0[is.na(salt_add$f.1478.0.0)]

salt_add <- salt_add[1:2]
salt_add$f.1478.0.0[salt_add$f.1478.0.0==-3] <- NA
colnames(salt_add) <- c("IID","salt_add")

#system("./ukbconv ukb52321.enc_ukb r -s30710 -ocrp")
crp <- read.table("/filepath/UKB/crp.tab", header=TRUE, sep="\t")
colnames(crp) <- c("IID","crp")


#system("./ukbconv ukb52321.enc_ukb r -s30630 -oapo_a")
apo_a <- read.table("/filepath/UKB/apo_a.tab", header=TRUE, sep="\t")
colnames(apo_a) <- c("IID","apo_a")

#system("./ukbconv ukb52321.enc_ukb r -s30640 -oapo_b")
apo_b <- read.table("/filepath/UKB/apo_b.tab", header=TRUE, sep="\t")
colnames(apo_b) <- c("IID","apo_b")

#system("")

#system("./ukbconv ukb52321.enc_ukb r -s26414 -oedu_england")
#system("./ukbconv ukb52321.enc_ukb r -s26431 -oedu_scotland")
#system("./ukbconv ukb52321.enc_ukb r -s26421 -oedu_wales")
edu <- read.table("/filepath/UKB/edu_england.tab", header=TRUE, sep="\t")
colnames(edu) <- c("IID","edu")

edu2 <- read.table("/filepath/UKB/edu_scotland.tab", header=TRUE, sep="\t")
colnames(edu2) <- c("IID","edu")

edu3 <- read.table("/filepath/UKB/edu_wales.tab", header=TRUE, sep="\t")
colnames(edu3) <- c("IID","edu")

edu <- merge(edu, edu2, by="IID")
edu <- merge(edu, edu3, by="IID")

edu$edu.x[is.na(edu$edu.x)] <- edu$edu.y[is.na(edu$edu.x)]
edu$edu.x[is.na(edu$edu.x)] <- edu$edu[is.na(edu$edu.x)]
edu <- edu[1:2]
colnames(edu) <- c("IID","edu")

rm(edu2, edu3)

#system("./ukbconv ukb52321.enc_ukb r -s26411 -oinc_england")
#system("./ukbconv ukb52321.enc_ukb r -s26428 -oinc_scotland")
#system("./ukbconv ukb52321.enc_ukb r -s26418 -oinc_wales")
#system("")
#system("./ukbconv ukb52321.enc_ukb r -s20003 -oall_meds")

############################################################################################################################
### Create Master Data Frame
############################################################################################################################

df <- merge(rs699[c("IID","SEX","rs699_G")],rs5051[c("IID","rs5051_T")],by="IID")
df <- merge(df, age,by="IID")
df <- merge(df, age_htn_dg[c("IID","age_htn_dg")], by="IID")
df <- merge(df, sys_bp[c("IID","sys_bp")], by="IID")
df <- merge(df, dia_bp[c("IID","dia_bp")], by="IID")
df <- merge(df, eth_rac[c("IID","eth_rac")], by="IID")
df <- merge(df, med_for_bp_chol_db[c("IID","on_bp_med")], by="IID")
df$mbp <- (df$sys_bp+df$dia_bp)/2

df$rs5051_T_round <- round(df$rs5051_T, digits = 0)
df$group <- NA

x1 <- subset(df, df$rs699_G==0 & df$rs5051_T_round==2)
x1$group <- "group1"

x2 <- subset(df, df$rs699_G==1 & df$rs5051_T_round==2)
x2$group <- "group2"

x3 <- subset(df, df$rs699_G==0 & df$rs5051_T_round==1)
x3$group <- "group3"

x4 <- subset(df, df$rs699_G==2 & df$rs5051_T_round==2)
x4$group <- "group4"

x5 <- subset(df, df$rs699_G==1 & df$rs5051_T_round==1)
x5$group <- "group5"

x6 <- subset(df, df$rs699_G==0 & df$rs5051_T_round==0)
x6$group <- "group6"

x7 <- subset(df, df$rs699_G==2 & df$rs5051_T_round==1)
x7$group <- "group7"

x8 <- subset(df, df$rs699_G==1 & df$rs5051_T_round==0)
x8$group <- "group8"

x9 <- subset(df, df$rs699_G==2 & df$rs5051_T_round==0)
x9$group <- "group9"

df <- rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9)
rm(x1,x2,x3,x4,x5,x6,x7,x8,x9)

df <- merge(df, pul_cor_sys[c("IID","pul_cor_sys")], by="IID")
df <- merge(df, pul_cor_dia[c("IID","pul_cor_dia")], by="IID")

df$mbp_pul_cor <- (df$pul_cor_sys+df$pul_cor_dia)/2

df <- merge(df, htn_dg[c("IID","htn_dg")], by="IID")
df <- merge(df, bmi[c("IID","bmi")], by="IID")

df <- merge(df, smk_stat[c("IID","smk_stat")], by="IID")
# df <- merge(df, veg_con[c("IID","veg_con")], by="IID")
# df <- merge(df, meat_con[c("IID","meat_con")], by="IID")
df <- merge(df, alc[c("IID","alc_less")], by="IID")
# df <- merge(df, carb_est_yst[c("IID","carb_est_yst")], by="IID")
# df <- merge(df, sod_est_yst[c("IID","sod_est_yst")], by="IID")
df <- merge(df, act[c("IID","act_more")], by="IID")
df <- merge(df, salt_add[c("IID","salt_add")], by="IID")
df <- merge(df, crp[c("IID","crp")], by="IID")
df <- merge(df, apo_a[c("IID","apo_a")], by="IID")
df <- merge(df, apo_b[c("IID","apo_b")], by="IID")
df <- merge(df, edu[c("IID","edu")], by="IID")

#create haplotype group
df$hap <- gsub("group","",df$group)
df$hap <- as.numeric(df$hap)

df$hap2 <- df$hap
df$hap2[df$hap2==2 | df$hap2==3] <- 2
df$hap2[df$hap2>=4 & df$hap2<=6] <- 3
df$hap2[df$hap2==7 | df$hap2==8] <- 4
df$hap2[df$hap2==9] <- 5

############################################################################################################################
### Systolic Blood Pressure 
############################################################################################################################

### rs699 ###
summary(lm(sys_bp ~ rs699_G, df))
summary(lm(sys_bp ~ rs699_G+SEX+age+bmi, df))
summary(lm(sys_bp ~ rs699_G+SEX+age+bmi+smk_stat, df))
summary(lm(sys_bp ~ rs699_G+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, df))

### rs5051 ###
summary(lm(sys_bp ~ rs5051_T, df))
summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, df))

### haplotype group ###
x <- df %>% group_by(hap) %>% summarise(mean = mean(sys_bp,na.rm = TRUE),n = n())
x <- df %>% group_by(hap2) %>% summarise(mean = mean(sys_bp,na.rm = TRUE),n = n())
x <- df %>% group_by(rs5051_T_round) %>% summarise(mean = mean(sys_bp,na.rm = TRUE),n = n())

summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, df))
summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, df))

#test for normality
ggdensity(df$sys_bp)

#download regression estimates
y <- as.data.frame(summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, df))$coefficients)
y1 <- as.data.frame(summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, df))$coefficients)
y2 <- as.data.frame(summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, df))$coefficients)

y <- rbind(y, y1, y2)
write.csv(y, file="sys_regression.csv")

### haplotype groups stratified by on_bp_med ###
summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$on_bp_med==1)))# on_bp_med
summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$on_bp_med==0)))# not on_bp_med

summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$on_bp_med==1)))# on_bp_med
summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$on_bp_med==0)))# not on_bp_med

summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$on_bp_med==1)))# on_bp_med
summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$on_bp_med==0)))# not on_bp_med

### haplotype groups stratified by age >=50 ###
summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$age>=50)))# age>=50
summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$age<50)))# age<50

summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$age>=50)))# age>=50
summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$age<50)))# age<50

summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$age>=50)))# age>=50
summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$age<50)))# age<50

### haplotype groups stratified by race ###
summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Black")))
summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="White_and_Black")))
summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Asian")))
summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="White_and_Asian")))
summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Indian")))
summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="White")))
summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Prefer_not_to_answer")))
summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Mixed_other")))
summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Other")))                                                                                         

summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Black")))
summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="White_and_Black")))
summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Asian")))
summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="White_and_Asian")))
summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Indian")))
summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="White")))
summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Prefer_not_to_answer")))
summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Mixed_other")))
summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Other")))  

summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Black")))
summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="White_and_Black")))
summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Asian")))
summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="White_and_Asian")))
summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Indian")))
summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="White")))
summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Prefer_not_to_answer")))
summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Mixed_other")))
summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Other")))  

summary(lm(sys_bp ~ rs699_G+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Black")))
summary(lm(sys_bp ~ rs699_G+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="White_and_Black")))
summary(lm(sys_bp ~ rs699_G+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Asian")))
summary(lm(sys_bp ~ rs699_G+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="White_and_Asian")))
summary(lm(sys_bp ~ rs699_G+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Indian")))
summary(lm(sys_bp ~ rs699_G+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="White")))
summary(lm(sys_bp ~ rs699_G+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Prefer_not_to_answer")))
summary(lm(sys_bp ~ rs699_G+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Mixed_other")))
summary(lm(sys_bp ~ rs699_G+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, subset(df, df$eth_rac=="Other")))  

summary(lm(sys_bp ~ rs5051_T, subset(df, df$eth_rac=="Black")))
summary(lm(sys_bp ~ rs5051_T, subset(df, df$eth_rac=="White")))

black <- ggplot(data = subset(df, df$eth_rac=="Black"), aes(x = as.factor(rs5051_T_round), y = sys_bp))+
  geom_jitter(color="black", size=0.4, alpha=0.5)+
  geom_boxplot(outlier.shape=NA)+theme_classic()

ggsave(black, filename = "rs5051_sys_black.pdf",units = "mm",width = 60, height = 60)

x <- subset(df, df$eth_rac=="Black") %>% group_by(rs5051_T_round) %>% summarise(mean = mean(sys_bp,na.rm = TRUE),n = n())
x <- subset(df, df$eth_rac=="Black") %>% group_by(rs699_G) %>% summarise(mean = mean(sys_bp,na.rm = TRUE),n = n())
x <- subset(df, df$eth_rac=="Black") %>% group_by(hap) %>% summarise(mean = mean(sys_bp,na.rm = TRUE),n = n())
x <- subset(df, df$eth_rac=="Black") %>% group_by(hap2) %>% summarise(mean = mean(sys_bp,na.rm = TRUE),n = n())

maf_rs5051 <- sum(df$rs5051_T)/(nrow(df)*2)
maf_rs699 <- sum(df$rs699_G)/(nrow(df)*2)

maf_rs5051_black <- sum(subset(df, df$eth_rac=="Black")$rs5051_T)/(nrow(subset(df, df$eth_rac=="Black"))*2)
maf_rs699_black <- sum(subset(df, df$eth_rac=="Black")$rs699_G)/(nrow(subset(df, df$eth_rac=="Black"))*2)

maf_rs5051_white <- sum(subset(df, df$eth_rac=="White")$rs5051_T)/(nrow(subset(df, df$eth_rac=="White"))*2)
maf_rs699_white <- sum(subset(df, df$eth_rac=="White")$rs699_G)/(nrow(subset(df, df$eth_rac=="White"))*2)

ld_pair <- data.frame(hap=c("AC","AT","GC","GT","AC","AT","GC","GT"),group=c(rep("African (n=1320",4),rep("European(n=1006)",4))
                      ,count=c(108,20,3,1189,583,9,8,406))

ggplot(data=ld_pair,aes(x=hap,y=count, fill=group,label=count))+
  geom_col(position = position_dodge())+
  geom_text(size=4, position =position_dodge(1),vjust=-.5)+
  theme_classic()

ggsave(plot = last_plot(), filename = "rs5051_rs699_ldpair.pdf",units = "mm",width = 120, height = 120)

############################################################################################################################
### Age HTN diagnosed 
############################################################################################################################
x <- df %>% group_by(hap2) %>% summarise(mean = mean(age_htn_dg,na.rm = TRUE),n = n())
x <- df %>% group_by(hap) %>% summarise(mean = mean(age_htn_dg,na.rm = TRUE),n = n())
x <- df %>% group_by(rs5051_T_round) %>% summarise(mean = mean(age_htn_dg,na.rm = TRUE),n = n())

### covariates alone ###
summary(lm(age_htn_dg ~ SEX, data = df)) # it is true that women seek medical treatment earlier than men, so sex as a covariate may not be appropriate here
summary(lm(age_htn_dg ~ age, data = df)) # age doesn't make sense for this analysis

### hap analysis ###
summary(lm(age_htn_dg ~ hap2, data = df))
summary(lm(age_htn_dg ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, df))

summary(lm(age_htn_dg ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, df))
summary(lm(age_htn_dg ~ rs5051_T, df))

ggplot(data = df, aes(x = group_cont, y = age_ht_dg))+
  geom_jitter(color="black", size=0.4, alpha=0.5)+
  geom_smooth(method = lm)+
  theme_classic()+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9), labels = c(1,2,3,4,5,6,7,8,9))+
  scale_y_continuous(breaks = c(20,30,40,45,50,55,60,70,80), labels = c(20,30,40,45,50,55,60,70,80))

############################################################################################################################
### logistic regression in people on BP meds
############################################################################################################################
summary(glm(on_bp_med ~ hap2, data = df,family = "binomial"))
summary(glm(on_bp_med ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df,family = "binomial"))

x <- df %>% group_by(hap2) %>% summarise(mean = mean(on_bp_med),n = n())

############################################################################################################################
### logistic regression in people with HTN diagnosis
############################################################################################################################
summary(glm(htn_dg ~ hap2, data = df,family = "binomial"))
summary(glm(htn_dg ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df,family = "binomial"))

summary(glm(htn_dg ~ rs5051_T_round+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df,family = "binomial"))

x <- df %>% group_by(hap) %>% summarise(mean = mean(htn_dg),n = n())

############################################################################################################################
### Diastolic Blood Pressure 
############################################################################################################################

### rs699 alone ###
summary(lm(dia_bp ~ rs699_G, df))
summary(lm(dia_bp ~ rs699_G + SEX + age, df))

### covariates alone ###
summary(lm(dia_bp ~ SEX, df))
summary(lm(dia_bp ~ on_bp_med, df))
summary(lm(dia_bp ~ SEX+age, df))
summary(lm(dia_bp ~ SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, df))

### haplotype group ###
x <- df %>% group_by(group) %>% summarise(mean = mean(dia_bp,na.rm = TRUE),n = n())

summary(lm(dia_bp ~ hap, data = df))
summary(lm(dia_bp ~ hap + SEX, data = df))
summary(lm(dia_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df)) #effect size increases after adding sex and age as covariate

### haplotype group stratified by on_bp_med ###
summary(lm(dia_bp ~ hap + SEX + age, data = subset(df, df$on_bp_med==1))) # on_bp_med
summary(lm(dia_bp ~ hap + SEX + age, data = subset(df, df$on_bp_med==0))) # not on_bp_med

### pulse corrected ###
x <- df %>% group_by(group) %>% summarise(mean = mean(pul_cor_dia,na.rm = TRUE),n = n())

summary(lm(pul_cor_dia ~ hap, data = df))
summary(lm(pul_cor_dia ~ hap + SEX, data = df))
summary(lm(pul_cor_dia ~ hap + SEX + age, data = df))

############################################################################################################################
### Mean Blood Pressure 
############################################################################################################################

### rs699 alone ###
summary(lm(mbp ~ rs699_G, df))
summary(lm(mbp ~ rs699_G + SEX + age, df))

### covariates alone ###
summary(lm(mbp ~ SEX, df))
summary(lm(mbp ~ on_bp_med, df))
summary(lm(mbp ~ SEX+age, df))

### haplotype group ###
x <- df %>% group_by(hap) %>% summarise(mean = mean(mbp,na.rm = TRUE),n = n())

summary(lm(mbp ~ hap, data = df))
summary(lm(mbp ~ hap + SEX, data = df)) #p-value decreases
summary(lm(mbp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df)) 
summary(lm(mbp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df)) 
summary(lm(mbp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df)) 

##############################################
## plot for sys_bp
##############################################
hap_sys_bp <- ggplot(data = df, aes(x = as.factor(hap), y = sys_bp))+
  geom_jitter(color="black", size=0.4, alpha=0.5)+
  geom_boxplot(outlier.shape=NA)+theme_classic()

hap2_sys_bp <- ggplot(data = df, aes(x = as.factor(hap2), y = sys_bp))+
  geom_jitter(color="black", size=0.4, alpha=0.5)+
  geom_boxplot(outlier.shape=NA)+theme_classic()

rs5051_sys_bp <- ggplot(data = df, aes(x = as.factor(rs5051_T_round), y = sys_bp))+
  geom_jitter(color="black", size=0.4, alpha=0.5)+
  geom_boxplot(outlier.shape=NA)+theme_classic()

plot2 <- ggarrange(hap_sys_bp, hap2_sys_bp, rs5051_sys_bp,nrow = 3)

### haplotype group stratified by on_bp_med ###
summary(lm(mbp ~ hap + SEX + age, data = subset(df, df$on_bp_med==1))) # on_bp_med
summary(lm(mbp ~ hap + SEX + age, data = subset(df, df$on_bp_med==0))) # not on_bp_med

############################################################################################################################
### Z-score comparison
############################################################################################################################
z_rs5051 <- data.frame(z_score=c(rep(NA,4)), lowerCI=c(rep(NA,4)), upperCI=c(rep(NA,4)))
rownames(z_rs5051) <- c("systolic blood pressure","diastolic blood pressure","age of HTN diagnosis","mean blood pressure")

sd_sys_bp <- sd(df$sys_bp, na.rm = TRUE)
sd_dia_bp <- sd(df$dia_bp, na.rm = TRUE)
sd_mbp <- sd(df$mbp, na.rm = TRUE)
sd_age_htn <- sd(df$age_htn_dg, na.rm = TRUE)

z_rs5051$z_score[1] <-  (summary(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df))$coefficients[2,1])/sd_sys_bp
z_rs5051$lowerCI[1] <- confint(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "rs5051_T", level=0.95)[1,1]/sd_sys_bp
z_rs5051$upperCI[1] <- confint(lm(sys_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "rs5051_T", level=0.95)[1,2]/sd_sys_bp
z_rs5051$z_score[2] <-  (summary(lm(dia_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df))$coefficients[2,1])/sd_dia_bp
z_rs5051$lowerCI[2] <- confint(lm(dia_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "rs5051_T", level=0.95)[1,1]/sd_dia_bp
z_rs5051$upperCI[2] <- confint(lm(dia_bp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "rs5051_T", level=0.95)[1,2]/sd_dia_bp
z_rs5051$z_score[3] <-  (summary(lm(age_htn_dg ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df))$coefficients[2,1])/sd_age_htn
z_rs5051$lowerCI[3] <- confint(lm(age_htn_dg ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "rs5051_T", level=0.95)[1,1]/sd_age_htn
z_rs5051$upperCI[3] <- confint(lm(age_htn_dg ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "rs5051_T", level=0.95)[1,2]/sd_age_htn
z_rs5051$z_score[4] <-  (summary(lm(mbp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df))$coefficients[2,1])/sd_mbp
z_rs5051$lowerCI[4] <- confint(lm(mbp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "rs5051_T", level=0.95)[1,1]/sd_mbp
z_rs5051$upperCI[4] <- confint(lm(mbp ~ rs5051_T+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "rs5051_T", level=0.95)[1,2]/sd_mbp

rs5051_plot <- ggplot()+
  geom_point(data=z_rs5051,mapping=aes(x=z_score,y=rownames(z_rs5051)))+
  geom_errorbar(data=z_rs5051,mapping=aes(y=rownames(z_rs5051),x=z_score,xmin=lowerCI,xmax=upperCI))+
  theme_classic()+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    #  panel.background=element_blank(),
    # axis.line=element_line(colour="black"),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    axis.title=element_text(size=12),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=8),
    strip.text.y.right=element_text(angle=0,size=10))+
  ggtitle("rs5051_T")+
  geom_vline(xintercept = 0,linetype="dashed",color="firebrick",size=.6,alpha=0.85)+
  scale_y_discrete(expand = c(0,1))+
  scale_x_continuous(limits = c(-0.055,0.055))

z_hap <- data.frame(z_score=c(rep(NA,4)), lowerCI=c(rep(NA,4)), upperCI=c(rep(NA,4)))
rownames(z_hap) <- c("systolic blood pressure","diastolic blood pressure","age of HTN diagnosis","mean blood pressure")

z_hap$z_score[1] <-  (summary(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df))$coefficients[2,1])/sd_sys_bp
z_hap$lowerCI[1] <- confint(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap", level=0.95)[1,1]/sd_sys_bp
z_hap$upperCI[1] <- confint(lm(sys_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap", level=0.95)[1,2]/sd_sys_bp
z_hap$z_score[2] <-  (summary(lm(dia_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df))$coefficients[2,1])/sd_dia_bp
z_hap$lowerCI[2] <- confint(lm(dia_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap", level=0.95)[1,1]/sd_dia_bp
z_hap$upperCI[2] <- confint(lm(dia_bp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap", level=0.95)[1,2]/sd_dia_bp
z_hap$z_score[3] <-  (summary(lm(age_htn_dg ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df))$coefficients[2,1])/sd_age_htn
z_hap$lowerCI[3] <- confint(lm(age_htn_dg ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap", level=0.95)[1,1]/sd_age_htn
z_hap$upperCI[3] <- confint(lm(age_htn_dg ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap", level=0.95)[1,2]/sd_age_htn
z_hap$z_score[4] <-  (summary(lm(mbp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df))$coefficients[2,1])/sd_mbp
z_hap$lowerCI[4] <- confint(lm(mbp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap", level=0.95)[1,1]/sd_mbp
z_hap$upperCI[4] <- confint(lm(mbp ~ hap+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap", level=0.95)[1,2]/sd_mbp

hap_plot <- ggplot()+
  geom_point(data=z_hap,mapping=aes(x=z_score,y=rownames(z_hap)))+
  geom_errorbar(data=z_hap,mapping=aes(y=rownames(z_hap),x=z_score,xmin=lowerCI,xmax=upperCI))+
  theme_classic()+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    #  panel.background=element_blank(),
    # axis.line=element_line(colour="black"),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    axis.title=element_text(size=12),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=8),
    strip.text.y.right=element_text(angle=0,size=10))+
  ggtitle("hap")+
  geom_vline(xintercept = 0,linetype="dashed",color="firebrick",size=.6,alpha=0.85)+
  scale_y_discrete(expand = c(0,1))+
  scale_x_continuous(limits = c(-0.055,0.055))

z_hap2 <- data.frame(z_score=c(rep(NA,4)), lowerCI=c(rep(NA,4)), upperCI=c(rep(NA,4)))
rownames(z_hap2) <- c("systolic blood pressure","diastolic blood pressure","age of HTN diagnosis","mean blood pressure")

z_hap2$z_score[1] <-  (summary(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df))$coefficients[2,1])/sd_sys_bp
z_hap2$lowerCI[1] <- confint(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap2", level=0.95)[1,1]/sd_sys_bp
z_hap2$upperCI[1] <- confint(lm(sys_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap2", level=0.95)[1,2]/sd_sys_bp
z_hap2$z_score[2] <-  (summary(lm(dia_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df))$coefficients[2,1])/sd_dia_bp
z_hap2$lowerCI[2] <- confint(lm(dia_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap2", level=0.95)[1,1]/sd_dia_bp
z_hap2$upperCI[2] <- confint(lm(dia_bp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap2", level=0.95)[1,2]/sd_dia_bp
z_hap2$z_score[3] <-  (summary(lm(age_htn_dg ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df))$coefficients[2,1])/sd_age_htn
z_hap2$lowerCI[3] <- confint(lm(age_htn_dg ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap2", level=0.95)[1,1]/sd_age_htn
z_hap2$upperCI[3] <- confint(lm(age_htn_dg ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap2", level=0.95)[1,2]/sd_age_htn
z_hap2$z_score[4] <-  (summary(lm(mbp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df))$coefficients[2,1])/sd_mbp
z_hap2$lowerCI[4] <- confint(lm(mbp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap2", level=0.95)[1,1]/sd_mbp
z_hap2$upperCI[4] <- confint(lm(mbp ~ hap2+SEX+age+bmi+smk_stat+act_more+alc_less+salt_add+crp+apo_a+apo_b+edu, data = df), "hap2", level=0.95)[1,2]/sd_mbp

hap2_plot <- ggplot()+
  geom_point(data=z_hap2,mapping=aes(x=z_score,y=rownames(z_hap2)))+
  geom_errorbar(data=z_hap2,mapping=aes(y=rownames(z_hap2),x=z_score,xmin=lowerCI,xmax=upperCI))+
  theme_classic()+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    #  panel.background=element_blank(),
    # axis.line=element_line(colour="black"),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    axis.title=element_text(size=12),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=8),
    strip.text.y.right=element_text(angle=0,size=10))+
  ggtitle("hap2")+
  geom_vline(xintercept = 0,linetype="dashed",color="firebrick",size=.6,alpha=0.85)+
  scale_y_discrete(expand = c(0,1))+
  scale_x_continuous(limits = c(-0.055,0.055))

plot <- ggarrange(hap_plot, rs5051_plot, hap2_plot,nrow = 3)

ggsave(plot = plot,filename = "wing_plot_UKB.pdf",width = 90,height = 120, units = "mm")

ggsave(plot = plot2,filename = "box_plot_UKB.tiff",width = 90,height = 120, units = "mm")

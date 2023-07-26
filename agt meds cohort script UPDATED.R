library(lubridate)
library(stringr)
library(dplyr)
library(httr)
library(jsonlite)
library(xml2)

#starting with full dataset
AGT.Cohort.Med.Data.for.Nick <- read.csv("~/AGT Cohort Med Data for Nick.csv")
a <- AGT.Cohort.Med.Data.for.Nick
rm(AGT.Cohort.Med.Data.for.Nick)

#generic HTN drug list from lexicomp
htn.meds <- read.csv("~/htn meds.csv", header = F)
b <- htn.meds
rm(htn.meds)
b1 <- as.data.frame(str_split_fixed(b$V1," ",10))
b2 <- stack(b1)
b3 <- as.data.frame(gsub(",","",b2$values,ignore.case = T))
b4 <- as.data.frame(unique(b3$`gsub(",", "", b2$values, ignore.case = T)`))
b5 <- as.data.frame(b4[-c(75:79),])
b5[77,1] <- "furosemide"
colnames(b5) <- "V1"
b5 <- as.data.frame(b5[-c(43,48),]) # removes minoxidil and magnesium
colnames(b5) <- "V1"
write.csv(b5,"lexi_HTN.csv",row.names = F) 
rm(b1,b2,b3,b4,b)

#using openFDA to get brand names
c <- data.frame()
for(l in 1:nrow(b5)){
  print(l)
  y <- b5[l,1]
  d <- content(GET(paste0("https://api.fda.gov/drug/ndc.json?search=generic_name:\"",y,"\"&limit=1000")))
  n <- length(d[["results"]])
  b <- data.frame()
  
  tryCatch(for(x in 1:n){
    skip_to_next <- FALSE
    tryCatch(b[x,1] <- d[["results"]][[x]][["brand_name"]], error=function(e) { skip_to_next <<- TRUE })
    if(skip_to_next) { next }
  }, error=function(e) { skip_to_next <<- TRUE })
  if(skip_to_next) { next }
  aa <- na.omit(unique(b))
  aa[2] <- y
  c <- rbind(c,aa)
}

rm(d,l,x,aa,b,y,n,skip_to_next)
colnames(c) <- c("brand","generic")

#creating list to do initial datafilter # this will allow for manual curation of weird lines later
c1 <- stack(c)

c1$values <- gsub("-"," ",c1$values)
c1$values <- gsub("/"," ",c1$values)
c1$values <- gsub(",","",c1$values)
c1$values <- gsub("\\(","",c1$values)
c1$values <- gsub("\\)","",c1$values)

c2 <- as.data.frame(str_split_fixed(c1$values," ",10))
c3 <- stack(c2)
colnames(c3) <- "V1"
c4 <- as.data.frame(unique(tolower((c3$V1))))
colnames(c4) <- "V1"
write.csv(c4,"c4.csv",row.names = F) # cleanup weirdos in excel - removed victoza, rectiv, menopause, headache, 
#benzoyl, clindamycin, cosopt, combigan, timoptic, toxycology,liraglutide
#also removed hydrochloride, blanks, and other nonsense
c4 <- read.csv("~/c4.csv", sep="")

  #getting short generic names
  c6 <- as.data.frame(substr(b5$V1,1,5))
  colnames(c6) <- "V1"
  
  #adding manual search terms
  c66 <- as.data.frame(c("HCTZ","HCT"))
  colnames(c66) <- "V1"
  
c7 <- rbind(c6,c4)
c7 <- rbind(c7,c66)
colnames(c7) <- "V1"
c7 <- as.data.frame(unique(tolower(c7$V1)))
write.csv(c7,"HTN med list for filtering.csv")
rm(c1,c2,c3,c4,c5,c6,c66)


#doing initial datafilter
d1 <- as.data.frame(unique(tolower(a$ordered_med_name)))
colnames(d1) <- "V1"
d1$V1 <- gsub("/"," ",d1$V1)
d1$V1 <- gsub("-"," ",d1$V1)
d2 <- as.data.frame(str_split_fixed(d1$V1," ",10))
d3 <- stack(d2)
d4 <- as.data.frame(unique(tolower(d3$values))) # all of the possibilities to search for in the data
colnames(d4) <- "V1"
rm(d1,d2,d3)

e <- data.frame()
for(n in 1:nrow(c7)){
  print(n)
  p <- c7[n,1]
  h <- subset(d4,grepl(p,d4$V1,ignore.case = T))
  e <- rbind(e,h)
}

e1 <- as.data.frame(unique(tolower(e$V1)))
rm(e,n,p,h)
row.names(e1)=NULL
write.csv(e1,"e1.csv",row.names = F) # cleanup weirdos in excel 
#removed antivert, chlorhexidine,chlorzoxazone,prochlorperazine, hydroxychloroquine,chlorpheniramine,dichlor,diphenhydramine#hydroc
#hydrocerin#hydroco#hydrocod#hydrocodon#hydrocodone#hydrocort#hydrocortisone#hydromorphon#hydromorphone#hydroxyz#hydroxyzine#hydrobromide#dihydrochl#hydrocortiso#25mghydroxyzine
#325hydrocodone#325mghydrocodone#methylphenidate#flunisolide#nitrofurantoin#triam#triamcinolone!!!!kept "triamt" which is the short code for triamterene,#besylate#medoxomil

e1 <- read.csv("~/e1.csv")
e2 <- stack(e1)
e3 <- as.data.frame(unique(e2$values))
e4 <- as.data.frame(e3[-c(71,90),])
rm(e2,e3)


#AGT cohort filtering down to HTN meds
f <- data.frame()
d1 <- a
rm(n,d,f1)
for(n in 1:nrow(e1)){
  print(n)
  p <- e1[n,1]
  d <- subset(d1,grepl(p,d1$ordered_med_name,ignore.case = T))
  d1 <- subset(d1,!grepl(p,d1$ordered_med_name,ignore.case = T))
  f <- rbind(f,d)
}

d2 <- as.data.frame(unique(d1$ordered_med_name))
d3 <- as.data.frame(unique(d1$APG.ID))

#rm(n)
#d4 <- data.frame()
#for(n in 1:nrow(d3)){
#  if(!grepl(d3[n,1],f$APG.ID)) { d4[n,1] <- d3[n,1] }
#}
f1 <- as.data.frame(unique(f$APG.ID))
f1[2] <- "blah"
colnames(f1) <- 'APG_ID'
colnames(d4) <- 'APG_ID'
d4[2] <- "all"
f2 <- merge(d4,f1,by="APG_ID",all = T)
colnames(f2) <- c('APG_ID',"a","b")
f3 <- subset(f2,is.na(f2$b))
rm(f1,f2)
                    
library(stringr)
g <- f
g$ordered_med_name <- gsub("/"," ",g$ordered_med_name)
g$ordered_med_name <- gsub("-"," ",g$ordered_med_name)
g$ordered_med_name <- gsub("\\(","",g$ordered_med_name)
g$ordered_med_name <- gsub("\\)","",g$ordered_med_name)

g[5:14] <- str_split_fixed(g$ordered_med_name," ",10)

row.names(g)=NULL

h <- g
rm(s,n,p)
for(s in 5:14){
  print(s)
  for(n in 1:nrow(h)){
    if(!grepl(h[n,s],e4,ignore.case = T)) { h[n,s] <- "" }
  }
}

rm(s,n)
for(n in 1:nrow(h)){
  h[n,15] <- paste0(h[n,6],h[n,7],h[n,8],h[n,9],h[n,10],h[n,11],h[n,12],h[n,13],h[n,14])
}

write.csv(h,"h.csv",row.names = F) #fixed columns in excel - removed column V6 through V14, removed "ER","LA","M","TAR",and fixed hydroc to hycrochlorothiazide
h <- read.csv("~/h.csv")

#converting brand names to generic (this adds multiple ingredients if combo brand) in data frame
write.csv(e4,"e4.csv") #removed chlordiazepoxide which removed two rows.  
#then added generic names for the combo products with space delimiter
#then removed rows where no conversion is needed
e6 <- read.csv("~/e6.csv")
rm(e5)

i <- h
i$V5 <- tolower(i$V5)
i$V15 <- tolower(i$V15)
rm(n,n1,n2)
for(n in 1:nrow(e6)){
  print(n)
  n1 <- e6[n,1]
  n2 <- e6[n,2]
  for(s in 1:nrow(i)){
    if(i[s,5]==n1) { i[s,5] <- n2 }  
  }
}

rm(n,s,n1,n2)
for(n in 1:nrow(e6)){
  print(n)
  n1 <- e6[n,1]
  n2 <- e6[n,2]
  for(s in 1:nrow(i)){
    if(i[s,6]==n1) { i[s,6] <- n2 }  
  }
}

write.csv(i,"i.csv",row.names = F) 
# there was one instance of "losartan potassium hydroc" which was hydrochlorothiazide. this is a potential error since 
#hydroc was removed from the code since most hits would come back as hydrocodone. there are probably few occurances of hydroc coding for hydrochlorothiazide.
  hydroc_subset <- subset(a,grepl("hydroc",a$ordered_med_name,ignore.case = T)) #the losartan example is the only hydroc = HCTZ example.  all good!
  rm(hydroc_subset)
  i <- read.csv("~/i.csv") #changed hydroc to hydrochlorothiazide
  
j <- i
j[7:8] <- str_split_fixed(j$V5," ",2)
j <- as.data.frame(j[-c(5)])

#get combo drugs into same row
j1 <- j[c(1:5)]
j2 <- j[c(1:4,6)]
j3 <- j[c(1:4,7)]
colnames(j1)[5] <- "V5"
colnames(j2)[5] <- "V5"
colnames(j3)[5] <- "V5"

j4 <- rbind(j1,j2)
j4 <- rbind(j4,j3)

j5 <- j4[-which(j4$V5==""),] #removes blanks
rm(j,j1,j2,j3,j4)

j6 <- j5
j6[6] <- as_date(mdy(j5$OriginalOrder.Dts))
j6$V5 <- tolower(j6$V5)

#MAX concomitant
#making reference container
k <- as.data.frame(c(1996:2020))

rm(n,d,s,v,w)
l <- data.frame()
for(n in 1:25){
  d <- n*4
  l[(d),1] <- paste0("q1_",k[n,1])
  l[(d+1),1] <- paste0("q2_",k[n,1])
  l[(d+2),1] <- paste0("q3_",k[n,1])  
  l[(d+3),1] <- paste0("q4_",k[n,1])
}

l[2:13] <- NA
colnames(l)[1:13] <- c("quarter","drug1","drug2","drug3","drug4","drug5","drug6","drug7","drug8","drug9","drug10","drug11","drug12")
row.names(l)=NULL

#making dataframe to define quarters
l1 <- l
l <- as.data.frame(l[-c(1:3),])
l1 <- as.data.frame(l1[-c(5:11)])
colnames(l1)[1:4] <- c("quarter","start","end","year")
l1[3:4] <- str_split_fixed(l1$quarter,"_",2)
l1[3] <- NA

rm(d,n)
for(n in 1:25){
  d <- n*4
  l1[(d),2] <- as_date(mdy(paste0("1-","1-",l1[d,4])))
  l1[(d+1),2] <- as_date(mdy(paste0("4-","1-",l1[d,4])))
  l1[(d+2),2] <- as_date(mdy(paste0("7-","1-",l1[d,4])))  
  l1[(d+3),2] <- as_date(mdy(paste0("10-","1-",l1[d,4])))
}
l1$start <- as_date(l1$start)

rm(d,n)
for(n in 1:25){
  d <- n*4
  l1[(d),3] <- as_date(mdy(paste0("3-","31-",l1[d,4])))
  l1[(d+1),3] <- as_date(mdy(paste0("6-","30-",l1[d,4])))
  l1[(d+2),3] <- as_date(mdy(paste0("9-","30-",l1[d,4])))  
  l1[(d+3),3] <- as_date(mdy(paste0("12-","31-",l1[d,4])))
}
l1$end <- as_date(l1$end)

l1 <- as.data.frame(l1[-c(1:3),])

#assigning quarters to the cohort data
rm(n,d)
for(n in 1:nrow(j6)){
  m <- subset(l1,j6[n,6]>=l1$start&j6[n,6]<=l1$end)
  m1 <- m[1,1]
  j6[n,7] <- m1
}


#the loop
o <- as.data.frame(unique(j6$APG.ID))
o[2:20] <- NA
colnames(o) <- c("APG_ID","total_unique","total_fills","total_start","total_end","total_range","total_fills_per_year","total_max_concom","drug1","drug2","drug3","drug4","drug5","drug6","drug7","drug8","drug9","drug10","drug11","drug12")

rm(n,p,q,m,m1,x,y,r,s,t,d,z,qlist,n1,n2)
for(n in 1:nrow(o)){
  print(n)
  p <- subset(j6,j6$APG.ID==o[n,1])
  q <- as.data.frame(unique(tolower(p$V5)))
  o[n,2] <- nrow(q)
  o[n,3] <- nrow(p)
  o[n,4] <- min(p$V6)
  o[n,5] <- max(p$V6)
  o[n,6] <- o[n,5]-o[n,4]
  o[n,7] <- o[n,3]/o[n,6]*365
  
  ll <- l
  qlist <- list()
  for(x in 1:nrow(q)){
    qlist[[x]] <- as.data.frame(subset(p,q[x,1]==p$V5))
    o[n,(8+x)] <- q[x,1]
  }
  for(y in 1:length(qlist)){
    r <- as.data.frame(qlist[y])
    for(z in 1:nrow(l)){
      ll[z,y+1] <- any(grepl(ll[z,1],r$V7))
    }
  }
  for(t in 1:nrow(ll)){
    ll[t,14] <- length(ll[t,2:(nrow(q)+1)][ll[t,2:(nrow(q)+1)]==TRUE])
  }
  o[n,8] <- max(ll[14])
}

write.csv(o,"final_med_count_AGT_study.csv",row.names = F)
write.csv(f3,"subjects with zero HTN meds.csv",row.names = F)

#added subjects with zero HTN meds
final_med_count_AGT_study_w_zeros <- read.csv("~/final_med_count_AGT_study_w_zeros.csv")
o <- final_med_count_AGT_study_w_zeros

#merge with tylers haplotype and dx counts
u <- read.csv("~/u.csv")
u <- u[c(1:348),]
colnames(u)[1] <- "APG_ID"
v <- merge(u,o,by='APG_ID',all = T)
write.csv(v,"AGT_results_raw_data.csv")

v <- v[-c(200,317),]

rm(n,n1)
vlist1=list()
for(n in 1:9){
  vlist1[[n]] <- as.data.frame(subset(v,v$haplotype_group==n))
}

v3 <- data.frame(c("total","essential_HTN","essential_HTN_preC","Any_HTN","Any_HTN_preC"))
v3[2:11] <- NA
colnames(v3)[2:11] <- c("1","2","3","4","5","6","7","8","9","total")
rm(n,w,w1,w3)
for(n in 1:9){
  v4 <- vlist1[[n]]
  v3[1,n+1] <- nrow(v4)
  v3[2,n+1] <- length(v4[5][v4[5]==TRUE])
  v3[3,n+1] <- length(v4[6][v4[6]=="Yes"])
  v3[4,n+1] <- length(v4[7][v4[7]==TRUE])
  v3[5,n+1] <- length(v4[8][v4[8]=="Yes"])
}
v5 <- v3
v5[1,11] <- sum(v5[1,2:10])
v5[2,11] <- sum(v5[2,2:10])
v5[4,11] <- sum(v5[4,2:10])

v6 <- v5
v6[2:11] <- NA

rm(n)
for(n in 2:11){
  v6[1,n] <- v5[1,n]/v5[1,n]
  v6[2,n] <- v5[2,n]/v5[1,n]
  v6[3,n] <- v5[3,n]/v5[1,n]
  v6[4,n] <- v5[4,n]/v5[1,n]
  v6[5,n] <- v5[5,n]/v5[1,n]
}


#subset with HTN med data
v2 <- as.data.frame(subset(v,!is.na(v$total_max_concom)&!is.na(v$rs5051_copies)))
v2$total_fills_per_year <- gsub("Inf","0",v2$total_fills_per_year)
v2$total_fills_per_year <- as.numeric(v2$total_fills_per_year)

vlist=list()
rm(n,t,x,y,z)
for(n in 1:9){
  vlist[[n]] <- as.data.frame(subset(v2,v2$haplotype_group==n))
}


w2 <- data.frame(c("total","essential_HTN","essential_HTN_preC","Any_HTN","Any_HTN_preC","max_concomitant","total_fills_per_year"))
w2[2:10] <- NA
colnames(w2)[2:10] <- c("1","2","3","4","5","6","7","8","9")
rm(n,w,w1,w3)
for(n in 1:9){
  w3 <- vlist[[n]]
  w2[1,n+1] <- nrow(w3)
  w2[2,n+1] <- length(w3[5][w3[5]==TRUE])
  w2[3,n+1] <- length(w3[6][w3[6]=="Yes"])
  w2[4,n+1] <- length(w3[7][w3[7]==TRUE])
  w2[5,n+1] <- length(w3[8][w3[8]=="Yes"])
  w2[6,n+1] <- mean(w3$total_max_concom)
  w2[7,n+1] <- mean(na.omit(w3$total_fills_per_year))
}
w4 <- w2
w4[1,11] <- sum(w4[1,2:10])
w4[2,11] <- sum(w4[2,2:10])
w4[4,11] <- sum(w4[4,2:10])
w4[6,11] <- mean(v2$total_max_concom)
w4[7,11] <- mean(v2$total_fills_per_year,na.rm = T)

w5 <- w4
w5[2:11] <- NA
rm(n)
for(n in 2:11){
  w5[1,n] <- w4[1,n]/w4[1,n]
  w5[2,n] <- w4[2,n]/w4[1,n]
  w5[3,n] <- w4[3,n]/w4[1,n]
  w5[4,n] <- w4[4,n]/w4[1,n]
  w5[5,n] <- w4[5,n]/w4[1,n]
}

#looks like there are 19 people with med data but not in the list with haplotypes
v1 <- as.data.frame(subset(v,is.na(v$rs5051_copies)))

#ggplot
# max concom
w_means <- v2 %>% 
  group_by(haplotype_group) %>% 
  summarise(total_max_concom = mean(total_max_concom))

w_sd <- v2 %>% 
  group_by(haplotype_group) %>% 
  summarise(total_max_concom = sd(total_max_concom))
colnames(w_sd)[2] <- "sd"
w_means <- merge(w_means,w_sd,by="haplotype_group")

ggplot()+
  geom_point(data = v2,mapping = aes(x=haplotype_group,y=total_max_concom),col="black",size=1,position=position_jitter(h=.1,w=.1))+
  geom_point(data = w_means,mapping = aes(x=haplotype_group,y=total_max_concom),col="blue",size=3,position="identity")+
  geom_errorbar(data = w_means,mapping = aes(x=haplotype_group, ymin=total_max_concom-sd,ymax=total_max_concom+sd))

# total fills per year
fills_means <- v2 %>% 
  group_by(haplotype_group) %>% 
  summarise(total_fills_per_year = mean(total_fills_per_year,na.rm=T))

fills_sd <- v2 %>% 
  group_by(haplotype_group) %>% 
  summarise(total_fills_per_year = sd(total_fills_per_year,na.rm=T))
colnames(fills_sd)[2] <- "sd"
w_means <- merge(fills_means,fills_sd,by="haplotype_group")

ggplot()+
  geom_point(data = v2,mapping = aes(x=haplotype_group,y=total_fills_per_year),col="black",size=1,position=position_jitter(h=.1,w=.1))+
  geom_point(data = w_means,mapping = aes(x=haplotype_group,y=total_fills_per_year),col="blue",size=3,position="identity")+
  geom_errorbar(data = w_means,mapping = aes(x=haplotype_group, ymin=total_fills_per_year-sd,ymax=total_fills_per_year+sd))

  #t-test
  t1 <- as.data.frame(subset(v2,v2$haplotype_group==3))
  t1.5 <- as.data.frame(subset(v2,v2$haplotype_group==8))
  
  t.test(t1$total_fills_per_year,t1.5$total_fills_per_year)
  
#bar plots for rates
  #all
  v7 <- v6[2,3:9]
  v8 <- stack(v7)
  ggplot(data=v8,aes(x=ind,y=values))+
    geom_bar(stat="identity")
  
  v81 <- v8
  v81[3,1] <- .559  
  v81[4,1] <- .437  
  v81[5,1] <- .523
  
  ggplot(data=v81,aes(x=ind,y=values))+
    geom_bar(stat="identity")
  
  #with meds
  w7 <- w5[2,3:9]
  w8 <- stack(w7)
  ggplot(data=w8,aes(x=ind,y=values))+
    geom_bar(stat="identity")
  
  colnames(grouped) <- c("ind","values")
  ggplot(data=grouped,aes(x=ind,y=values))+
    geom_bar(stat="identity")

#plotting for manuscript
library(ggplot2)
library(ggpubr)
apg <- read.csv("~/AGT_results_raw_data.csv")

apg1 <- apg[c(5,8,15,16)]

a <- ggplot(data = subset(apg1, !is.na(apg1$haplotype_group)), aes(x = as.factor(haplotype_group), y = total_fills_per_year))+
  geom_jitter(color="black", size=0.5, alpha=0.5)+
  geom_boxplot(outlier.shape=NA, alpha=0.7)+theme_classic()

b <- ggplot(data = subset(apg1, !is.na(apg1$haplotype_group)), aes(x = as.factor(haplotype_group), y = total_max_concom))+
  geom_jitter(color="black", size=0.5, alpha=0.5)+
  geom_boxplot(outlier.shape=NA, alpha=0.7)+theme_classic()


ggarrange(a, b)

ggsave(plot = last_plot(),filename = "apg.pdf",units = "mm",width = 120, height = 90)
  
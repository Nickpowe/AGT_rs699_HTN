

setwd("/filepath/AGT")
AGT_hepg2_qPCR_results <- read.delim("/filepath/AGT/agt_hepg2_qPCR_results.txt", stringsAsFactors=FALSE)
AGT_HT29_qPCR_results <- read.delim("/filepath/AGT/AGT_colon_qPCR_results.txt")
AGT_HEK293_qPCR_results <- read.delim("/filepath/AGT/AGT_HEK293_qPCR_results.txt")

######################################################################
### HepG2
######################################################################
x <- subset(AGT_hepg2_qPCR_results,AGT_hepg2_qPCR_results$group=="WT")
AGT_hepg2_qPCR_results$Delta_Delta <- AGT_hepg2_qPCR_results$Delta-(mean(x$Delta))
AGT_hepg2_qPCR_results$fold_change <- 2^-AGT_hepg2_qPCR_results$Delta_Delta

t_set <- aggregate(AGT_hepg2_qPCR_results[c(3,4,8)],list(AGT_hepg2_qPCR_results$group2),mean)
t_set[5:6] <- str_split_fixed(t_set$Group.1,"-",2)
colnames(t_set)[5] <- "group"

var.test(AGT_CT ~ group,data = t_set,alternative="two.sided")
var.test(Delta ~ group,data = t_set,alternative="two.sided")

t1 <- t.test(AGT_CT ~ group,data = t_set)
t2 <- t.test(Delta ~ group,data = t_set)
t3 <- t.test(fold_change ~ group,data=t_set)

t11 <- t.test(AGT_CT ~ group,data = AGT_hepg2_qPCR_results)
t22 <- t.test(Delta ~ group,data = AGT_hepg2_qPCR_results)
t33 <- t.test(fold_change ~ group,data=AGT_hepg2_qPCR_results)

plot1 <- ggplot()+geom_boxplot(data=AGT_hepg2_qPCR_results,
                               aes(x=factor(group, levels = c("WT","Variant_T>C")),y=fold_change,
                                   fill=factor(group, levels = c("WT","Variant_T>C"))),outlier.shape = NA,alpha=0.5)+
  geom_quasirandom(data=AGT_hepg2_qPCR_results,
                   aes(x=factor(group, levels = c("WT","Variant_T>C")),
                       y=fold_change,
                       shape=group2,
                       color=group2),size=5,width = .4)+
  geom_quasirandom(data=AGT_hepg2_qPCR_results,
                   aes(x=factor(group, levels = c("WT","Variant_T>C")),
                       y=fold_change,
                       shape=group2,
                       color=group2),size=3,width = .4)+
  geom_quasirandom(data=AGT_hepg2_qPCR_results,
                   aes(x=factor(group, levels = c("WT","Variant_T>C")),
                       y=fold_change,
                       shape=group2,
                       color=group2),size=4,width = .4)+
  #scale_shape_manual(values = c(0,1,5,6,7,10,9,8))+
  scale_shape_manual(values = c(15,16,17,18,15,16,17,18))+
  #scale_color_manual(values=c("dodgerblue4","darkgreen","orangered3","grey23","goldenrod3","firebrick4","brown","purple"))+
  scale_color_manual(values=c("darkmagenta","salmon","dodgerblue","#009E73","goldenrod1","#0072B2","firebrick3","#CC79A7"))+
  scale_fill_manual(values=c("gray40","gray40"))+
  theme_classic()+
  theme(
    text = element_text(size=8),
    #axis.title.y = element_text(size=8),
    #axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    # panel.background = element_rect(fill="grey90"),
    # plot.background = element_rect(fill="grey90"),
    # legend.background = element_rect(fill="grey90"),
    axis.line.x = element_line(color="black"),
    axis.line.y = element_line(color="black"),
    legend.title = element_blank(),
    legend.position = "top")+
  guides(fill=F,shape=guide_legend(nrow = 1))+
  ylab("Fold Change (2^-DeltaDeltaCt)")+
  xlab(paste0("T-test: ",
              round(t3[["estimate"]][["mean in group WT"]],1),
              " vs. ",
              round(t3[["estimate"]][["mean in group VarC"]],1),
              "\nP-Value: ",formatC(t3[["p.value"]],format = "e",2),
              "\nT-test including pcr reps: ",
              round(t33[["estimate"]][["mean in group WT"]],1),
              " vs. ",
              round(t33[["estimate"]][["mean in group Variant_T>C"]],1),
              "\nP-Value: ",formatC(t33[["p.value"]],format = "e",2),
              "\nFold Change Mean: ",round(mean(subset(AGT_hepg2_qPCR_results$fold_change,
                                                       AGT_hepg2_qPCR_results$group=="Variant_T>C")),2)))+
  scale_y_continuous(limits = c(0,3))

plot1


######################################################################
### HT29
######################################################################
x <- subset(AGT_HT29_qPCR_results,AGT_HT29_qPCR_results$group=="WT")
AGT_HT29_qPCR_results$Delta_Delta <- AGT_HT29_qPCR_results$Delta-(mean(x$Delta))
AGT_HT29_qPCR_results$fold_change <- 2^-AGT_HT29_qPCR_results$Delta_Delta

t_set <- aggregate(AGT_HT29_qPCR_results[c(4,6,10)],list(AGT_HT29_qPCR_results$group2),mean)
t_set[5:6] <- str_split_fixed(t_set$Group.1,"-",2)
colnames(t_set)[5] <- "group"

var.test(AGT_CT ~ group,data = t_set,alternative="two.sided")
var.test(Delta ~ group,data = t_set,alternative="two.sided")

t1 <- t.test(AGT_CT ~ group,data = t_set)
t2 <- t.test(Delta ~ group,data = t_set)
t3 <- t.test(fold_change ~ group,data=t_set)

t11 <- t.test(AGT_CT ~ group,data = AGT_HT29_qPCR_results)
t22 <- t.test(Delta ~ group,data = AGT_HT29_qPCR_results)
t33 <- t.test(fold_change ~ group,data=AGT_HT29_qPCR_results)

plot2 <- ggplot()+geom_boxplot(data=AGT_HT29_qPCR_results,
                               aes(x=factor(group, levels = c("WT","Variant_T>C")),y=fold_change,
                                   fill=factor(group, levels = c("WT","Variant_T>C"))),outlier.shape = NA,alpha=0.5)+
  geom_quasirandom(data=AGT_HT29_qPCR_results,
                   aes(x=factor(group, levels = c("WT","Variant_T>C")),
                       y=fold_change,
                       shape=group2,
                       color=group2),size=5,width = .4)+
  geom_quasirandom(data=AGT_HT29_qPCR_results,
                   aes(x=factor(group, levels = c("WT","Variant_T>C")),
                       y=fold_change,
                       shape=group2,
                       color=group2),size=3,width = .4)+
  geom_quasirandom(data=AGT_HT29_qPCR_results,
                   aes(x=factor(group, levels = c("WT","Variant_T>C")),
                       y=fold_change,
                       shape=group2,
                       color=group2),size=4,width = .4)+
  #scale_shape_manual(values = c(0,1,5,6,7,10,9,8))+
  scale_shape_manual(values = c(15,16,17,18,15,16,17,18))+
  #scale_color_manual(values=c("dodgerblue4","darkgreen","orangered3","grey23","goldenrod3","firebrick4","brown","purple"))+
  scale_color_manual(values=c("darkmagenta","salmon","dodgerblue","#009E73","goldenrod1","#0072B2","firebrick3","#CC79A7"))+
  scale_fill_manual(values=c("gray40","gray40"))+
  theme_classic()+
  theme(
    text = element_text(size=8),
    #axis.title.y = element_text(size=8),
    #axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    # panel.background = element_rect(fill="grey90"),
    # plot.background = element_rect(fill="grey90"),
    # legend.background = element_rect(fill="grey90"),
    axis.line.x = element_line(color="black"),
    axis.line.y = element_line(color="black"),
    legend.title = element_blank(),
    legend.position = "top")+
  guides(fill=F,shape=guide_legend(nrow = 1))+
  ylab("Fold Change (2^-DeltaDeltaCt)")+
  xlab(paste0("T-test: ",
              round(t3[["estimate"]][["mean in group WT"]],1),
              " vs. ",
              round(t3[["estimate"]][["mean in group VarC"]],1),
              "\nP-Value: ",formatC(t3[["p.value"]],format = "e",2),
              "\nT-test including pcr reps: ",
              round(t33[["estimate"]][["mean in group WT"]],1),
              " vs. ",
              round(t33[["estimate"]][["mean in group Variant_T>C"]],1),
              "\nP-Value: ",formatC(t33[["p.value"]],format = "e",2),
              "\nFold Change Mean: ",round(mean(subset(AGT_HT29_qPCR_results$fold_change,
                                                       AGT_HT29_qPCR_results$group=="Variant_T>C")),2)))+
  scale_y_continuous(limits = c(0,3))

plot2

######################################################################
### Hek293
######################################################################
x <- subset(AGT_HEK293_qPCR_results,AGT_HEK293_qPCR_results$group=="WT")
AGT_HEK293_qPCR_results$Delta_Delta <- AGT_HEK293_qPCR_results$Delta.AGT-(mean(x$Delta.AGT))
AGT_HEK293_qPCR_results$fold_change <- 2^-AGT_HEK293_qPCR_results$Delta_Delta

t_set <- aggregate(AGT_HEK293_qPCR_results[c(4,6,10)],list(AGT_HEK293_qPCR_results$group2),mean)
t_set[5:6] <- str_split_fixed(t_set$Group.1,"-",2)
colnames(t_set)[5] <- "group"

# var.test(AGT_CT ~ group,data = t_set,alternative="two.sided")
# var.test(Delta.AGT ~ group,data = t_set,alternative="two.sided")
# 
# t1 <- t.test(AGT_CT ~ group,data = t_set)
# t2 <- t.test(Delta.AGT ~ group,data = t_set)

t11 <- t.test(AGT_CT ~ group,data = AGT_HEK293_qPCR_results)
t22 <- t.test(Delta.AGT ~ group,data = AGT_HEK293_qPCR_results)

plot3 <- ggplot()+geom_boxplot(data=AGT_HEK293_qPCR_results,
                               aes(x=factor(group, levels = c("WT","Variant_T>C")),y=fold_change,
                                   fill=factor(group, levels = c("WT","Variant_T>C"))),outlier.shape = NA,alpha=0.5)+
  geom_quasirandom(data=AGT_HEK293_qPCR_results,
                   aes(x=factor(group, levels = c("WT","Variant_T>C")),
                       y=fold_change,
                       shape=group2,
                       color=group2),size=5,width = .4)+
  geom_quasirandom(data=AGT_HEK293_qPCR_results,
                   aes(x=factor(group, levels = c("WT","Variant_T>C")),
                       y=fold_change,
                       shape=group2,
                       color=group2),size=3,width = .4)+
  geom_quasirandom(data=AGT_HEK293_qPCR_results,
                   aes(x=factor(group, levels = c("WT","Variant_T>C")),
                       y=fold_change,
                       shape=group2,
                       color=group2),size=4,width = .4)+
  #scale_shape_manual(values = c(0,1,5,6,7,10,9,8))+
  scale_shape_manual(values = c(15,16,17,18,15,16,17,18))+
  #scale_color_manual(values=c("dodgerblue4","darkgreen","orangered3","grey23","goldenrod3","firebrick4","brown","purple"))+
  scale_color_manual(values=c("darkmagenta",
                              #"salmon","dodgerblue","#009E73",
                              "goldenrod1"#,"#0072B2","firebrick3","#CC79A7"
                              ))+
  scale_fill_manual(values=c("gray40","gray40"))+
  theme_classic()+
  theme(
    text = element_text(size=8),
    #axis.title.y = element_text(size=8),
    #axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    # panel.background = element_rect(fill="grey90"),
    # plot.background = element_rect(fill="grey90"),
    # legend.background = element_rect(fill="grey90"),
    axis.line.x = element_line(color="black"),
    axis.line.y = element_line(color="black"),
    legend.title = element_blank(),
    legend.position = "top")+
  guides(fill=F,shape=guide_legend(nrow = 1))+
  ylab("Fold Change (2^-DeltaDeltaCt)")+
xlab(paste0("T-test: NA",
            "\nP-Value: NA",
            "\nT-test including pcr reps: NA",
            "\nP-Value: NA",
            "\nFold Change Mean: NA"))+
  scale_y_continuous(limits = c(0,3))

plot3

######################################################################
### combine
######################################################################
all_plots <- ggarrange(plot1,plot2,plot3,labels = "AUTO",common.legend = T,ncol = 3)
ggsave("/filepath/AGT/AGT_qPCR_combined.pdf",plot = all_plots,width = 300,height = 100,units = "mm")

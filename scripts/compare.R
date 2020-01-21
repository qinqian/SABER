args=commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)

####helper functions####
se <- function(x) sqrt(var(x, na.rm=TRUE)/length(x))

data_summary<-function(x){
  m<-median(x)
  ymin<-m-se(x)
  ymax<-m+se(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

data_summary1<-function(x){
  m<-median(x)
  ymin<-m-(IQR(x)/2)
  ymax<-m+(IQR(x)/2)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

data_summary2<-function(x){
  m<-mean(x)
  ymin<-m-se(x)
  ymax<-m+se(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

data_summary3<-function(x){
  m<-mean(x)
  ymin<-m-sd(x)
  ymax<-m+sd(x)
  return(c(y=m, ymin=ymin,ymax=ymax))
}

# diversity index functions
shannon_func<-function(x){
  prop<-x/sum(x)
  shannon_entropy<--1*sum(prop*log2(prop))
  return(shannon_entropy)
}

simpson_func<-function(x){
  prop<-x/sum(x)
  simpson_diversity<-1/sum(prop^2)
  return(simpson_diversity)
}

general_diversity_func<-function(x,q){
  prop<-x/sum(x)
  general_diversity<-(sum(prop^q))^(1/(1-q))
  return(general_diversity)
}

####main stat function####
output_stats<-function(condition1_data,condition1_label,condition2_data,condition2_label,output_folder,height,width){
  #load the data - data must be produced by analysis v4.0 or later script!!
  c1_df<-read.csv(paste0(condition1_data,"/SABER/data_out.csv"))
  c2_df<-read.csv(paste0(condition2_data,"/SABER/data_out.csv"))
  c1_df$class<-condition1_label
  c2_df$class<-condition2_label
  df_all<-rbind(c1_df,c2_df)
  
  df_vaf2<-df_all[which(df_all$stat=="vaf2"),]
  df_vaf1<-df_all[which(df_all$stat=="vaf1"),]
  df_simpson<-df_all[which(df_all$stat=="simpson"),]
  df_shannon<-df_all[which(df_all$stat=="shannon"),]
  
  dir.create(output_folder)
  #plot the data
  pvaf2<-ggplot(df_vaf2, aes(x = df_vaf2$class, y = df_vaf2$value))
  pvaf2<-pvaf2+
    geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill = df_vaf2$class)) + 
    scale_fill_manual(values=c("#DC0000","#3C5488")) + 
    scale_color_manual(values=c("black","black"))+
    stat_summary(fun.data=data_summary3, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
    theme_cowplot(font_size = 11)+
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank())+
    #scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
    scale_x_discrete(breaks = c(condition1_label,condition2_label), labels = c(condition1_label,condition2_label))+
    labs(y="Clones with VAF>0.02")
  pvaf2
  save_plot(plot = pvaf2, filename = paste0(output_folder,"/vaf2_plot.pdf"), base_height = height, base_width = width)
  
  pvaf1<-ggplot(df_vaf1, aes(x = df_vaf1$class, y = df_vaf1$value))
  pvaf1<-pvaf1+
    geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill = df_vaf1$class)) + 
    scale_fill_manual(values=c("#DC0000","#3C5488")) + 
    scale_color_manual(values=c("black","black"))+
    stat_summary(fun.data=data_summary3, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
    theme_cowplot(font_size = 11)+
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank())+
    #scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
    scale_x_discrete(breaks = c(condition1_label,condition2_label), labels = c(condition1_label,condition2_label))+
    labs(y="Clones with VAF>0.01")
  pvaf1
  save_plot(plot = pvaf1, filename = paste0(output_folder,"/vaf1_plot.pdf"), base_height = height, base_width = width)
  
  pshannon<-ggplot(df_shannon, aes(x = df_shannon$class, y = df_shannon$value))
  pshannon<-pshannon+
    geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill = df_shannon$class)) + 
    scale_fill_manual(values=c("#DC0000","#3C5488")) + 
    scale_color_manual(values=c("black","black"))+
    stat_summary(fun.data=data_summary3, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
    theme_cowplot(font_size = 11)+
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank())+
    #scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
    scale_x_discrete(breaks = c(condition1_label,condition2_label), labels = c(condition1_label,condition2_label))+
    labs(y="Shannon Entropy")
  pshannon
  save_plot(plot = pshannon, filename = paste0(output_folder,"/shannon_plot.pdf"), base_height = height, base_width = width)
  
  psimpson<-ggplot(df_simpson, aes(x = df_simpson$class, y = df_simpson$value))
  psimpson<-psimpson+
    geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill = df_simpson$class)) + 
    scale_fill_manual(values=c("#DC0000","#3C5488")) + 
    scale_color_manual(values=c("black","black"))+
    stat_summary(fun.data=data_summary3, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
    theme_cowplot(font_size = 11)+
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank())+
    #scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
    scale_x_discrete(breaks = c(condition1_label,condition2_label), labels = c(condition1_label,condition2_label))+
    labs(y="Simpson Diversity")
  psimpson
  save_plot(plot = psimpson, filename = paste0(output_folder,"/simpson_plot.pdf"), base_height = height, base_width = width)
  
  genstat<-function(x,y,stat,output_folder,condition1_label,condition2_label){
    xvec<-x[which(x$stat==stat),1]
    yvec<-y[which(y$stat==stat),1]
    ks<-ks.test(xvec,yvec)
    shapirox<-shapiro.test(xvec)
    shapiroy<-shapiro.test(yvec)
    t_equal<-t.test(xvec,yvec,alternative = "two.sided", var.equal = TRUE)
    t_unequal<-t.test(xvec,yvec,alternative = "two.sided", var.equal = FALSE)
    wilcox<-wilcox.test(xvec,yvec,alternative = "two.sided")
    returnlist<-list(paste0("X is ", condition1_label),xvec,paste0("SD of X is ",sd(xvec)),paste0("SEM of X is ",se(xvec)),
                     paste0("Y is ", condition2_label),yvec,paste0("SD of Y is ",sd(yvec)),paste0("SEM of Y is ",se(yvec)),
                     paste0("Statistic is ",stat),
                     shapirox,shapiroy,ks,t_equal,t_unequal,wilcox)
    sink(paste0(output_folder,"/stat_report_",stat,".txt"))
    print(returnlist)
  }
  genstat(x = c1_df, y = c2_df, stat = "vaf2", condition1_label = condition1_label, condition2_label = condition2_label, output_folder = output_folder);sink()
  genstat(x = c1_df, y = c2_df, stat = "vaf1", condition1_label = condition1_label, condition2_label = condition2_label,output_folder = output_folder);sink()
  genstat(x = c1_df, y = c2_df, stat = "shannon", condition1_label = condition1_label, condition2_label = condition2_label,output_folder = output_folder);sink()
  genstat(x = c1_df, y = c2_df, stat = "simpson", condition1_label = condition1_label, condition2_label = condition2_label,output_folder= output_folder);sink()
  
}
output_stats(condition1_data = args[1],
             condition1_label = args[2],
             condition2_data = args[3],
             condition2_label = args[4],
             output_folder = args[5],
             height = 2.3,
             width = 2.5)





####function calls####
# output_stats(condition1_data = "output_training_tvaf_10_tpaf_0.05",
#              condition1_label = "Training",
#              condition2_data = "output_validation_tvaf_10_tpaf_0.05",
#              condition2_label = "Validation",
#              output_folder = "comparison_training_validation_tvaf_10_tpaf_0.05_851ed",
#              height = 2.3,
#              width = 2.5)
# 
# output_stats(condition1_data = "output_training_tvaf_3_tpaf_0.05",
#              condition1_label = "Training",
#              condition2_data = "output_validation_tvaf_3_tpaf_0.05",
#              condition2_label = "Validation",
#              output_folder = "comparison_training_validation_tvaf_3_tpaf_0.05_4aa7f",
#              height = 2.3,
#              width = 2.5)
# 
# output_stats(condition1_data = "output_training_tvaf_1_tpaf_0.05",
#              condition1_label = "Training",
#              condition2_data = "output_validation_tvaf_1_tpaf_0.05",
#              condition2_label = "Validation",
#              output_folder = "comparison_training_validation_tvaf_1_tpaf_0.05_aa2eb",
#              height = 2.3,
#              width = 2.5)
# 
# output_stats(condition1_data = "output_training_tvaf_0.3_tpaf_0.05",
#              condition1_label = "Training",
#              condition2_data = "output_validation_tvaf_0.3_tpaf_0.05",
#              condition2_label = "Validation",
#              output_folder = "comparison_training_validation_tvaf_0.3_tpaf_0.05_1fe38",
#              height = 2.3,
#              width = 2.5)
# 
# output_stats(condition1_data = "output_training_tvaf_0.1_tpaf_0.05",
#              condition1_label = "Training",
#              condition2_data = "output_validation_tvaf_0.1_tpaf_0.05",
#              condition2_label = "Validation",
#              output_folder = "comparison_training_validation_tvaf_0.1_tpaf_0.05_cb098",
#              height = 2.3,
#              width = 2.5)
# 
# output_stats(condition1_data = "output_training_tvaf_0.03_tpaf_0.05",
#              condition1_label = "Training",
#              condition2_data = "output_validation_tvaf_0.03_tpaf_0.05",
#              condition2_label = "Validation",
#              output_folder = "comparison_training_validation_tvaf_0.03_tpaf_0.05_54af4",
#              height = 2.3,
#              width = 2.5)
# 
# output_stats(condition1_data = "output_gestaltAE_mcs_tvaf_0.3_tpaf_0.05",
#              condition1_label = "MCS",
#              condition2_data = "output_gestaltAE_prkcda_tvaf_0.3_tpaf_0.05",
#              condition2_label = "prkcda",
#              output_folder = "comparison_gestaltAE_6mo_tvaf_0.03_tpaf_0.05_bd85c",
#              height = 2.3,
#              width = 2.5)
# 
# output_stats(condition1_data = "output_gestaltAJ_mcs_tvaf_0.3_tpaf_0.1",
#              condition1_label = "MCS",
#              condition2_data = "output_gestaltAJ_prkcda_tvaf_0.3_tpaf_0.1",
#              condition2_label = "prkcda",
#              output_folder = "comparison_gestaltAJ_9mo_tvaf_0.03_tpaf_0.1_54a37",
#              height = 2.3,
#              width = 2.5)
# 
# output_stats(condition1_data = "output_gestalt_training_tvaf_0.3",
#              condition1_label = "Training",
#              condition2_data = "output_gestalt_validation_tvaf_0.3",
#              condition2_label = "Validation",
#              output_folder = "comparison_training_validation_20190829",
#              height = 2.3,
#              width = 2.5)
# 
# output_stats(condition1_data = "output_gestalt_training_tvaf_0.3",
#              condition1_label = "Training",
#              condition2_data = "output_gestalt_validation_tvaf_1",
#              condition2_label = "Validation",
#              output_folder = "comparison_training_validation_20190829_v2",
#              height = 2.3,
#              width = 2.5)
# 
# output_stats(condition1_data = "output_gestalt_training_tvaf_1",
#              condition1_label = "Training",
#              condition2_data = "output_gestalt_validation_tvaf_1",
#              condition2_label = "Validation",
#              output_folder = "comparison_training_validation_20190829_v3",
#              height = 2.3,
#              width = 2.5)

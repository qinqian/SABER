library(cowplot)

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
output_stats<-function(condition1_data,condition1_label,condition2_data,condition2_label, output_folder,height,width,font_size,title,sig){
  #load the data - data must be produced by analysis v4.0 or later script!!
  c1_vaf2_df<-read.csv(paste0(condition1_data,"/data_out_vaf2.csv"));c1_vaf2_df$class<-condition1_label
  c2_vaf2_df<-read.csv(paste0(condition2_data,"/data_out_vaf2.csv"));c2_vaf2_df$class<-condition2_label
  
  c1_shannon_df<-read.csv(paste0(condition1_data,"/data_out_shannon.csv"));c1_shannon_df$class<-condition1_label
  c2_shannon_df<-read.csv(paste0(condition2_data,"/data_out_shannon.csv"));c2_shannon_df$class<-condition2_label
  
  c1_simpson_df<-read.csv(paste0(condition1_data,"/data_out_simpson.csv"));c1_simpson_df$class<-condition1_label
  c2_simpson_df<-read.csv(paste0(condition2_data,"/data_out_simpson.csv"));c2_simpson_df$class<-condition2_label
  
  
  df_vaf2<-rbind(c1_vaf2_df,c2_vaf2_df)
  df_shannon<-rbind(c1_shannon_df,c2_shannon_df)
  df_simpson<-rbind(c1_simpson_df,c2_simpson_df)
  
  output_folder_old<-output_folder
  dt_stamp<-str_replace_all(string = Sys.time(), pattern = "([-\\s:ESTD])",replacement = "")
  output_folder<-paste0(output_folder_old,"_",dt_stamp)
  dir.create(output_folder)
  
  #plot the data
  pvaf2<-ggplot(df_vaf2, aes(x = df_vaf2$class, y = df_vaf2$value))
  pvaf2<-pvaf2+
    geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill = df_vaf2$class)) + 
    scale_fill_manual(values=c("#DC0000","#3C5488")) + 
    scale_color_manual(values=c("black","black"))+
    stat_summary(fun.data=data_summary3, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
    theme_cowplot(font_size = font_size)+
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank())+
    scale_y_continuous(limits = c(0,max(df_vaf2$value)*1.15))+
    annotate("text", x = 1.5, y = max(df_vaf2$value)*1.1, label = sig, size = 3)+
    scale_x_discrete(breaks = c(condition1_label,condition2_label), labels = c(condition1_label,condition2_label))+
    labs(y="Clones with VAF>0.02", title = title)+
    theme(plot.title = element_text(hjust = 0.5))
  pvaf2
  save_plot(plot = pvaf2, filename = paste0(output_folder,"/vaf2_plot.pdf"), base_height = height, base_width = width)
  
  pshannon<-ggplot(df_shannon, aes(x = df_shannon$class, y = df_shannon$value))
  pshannon<-pshannon+
    geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill = df_shannon$class)) + 
    scale_fill_manual(values=c("#DC0000","#3C5488")) + 
    scale_color_manual(values=c("black","black"))+
    stat_summary(fun.data=data_summary2, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
    theme_cowplot(font_size = font_size)+
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank())+
    scale_y_continuous(limits = c(0,max(df_shannon$value)*1.15))+
    #annotate("text", x = 1.5, y = max(df_shannon$value)*1.1, label = sig, size = 2)+
    scale_x_discrete(breaks = c(condition1_label,condition2_label), labels = c(condition1_label,condition2_label))+
    labs(y="Shannon Entropy", title = title)+
    theme(plot.title = element_text(hjust = 0.5))
  pshannon
  save_plot(plot = pshannon, filename = paste0(output_folder,"/shannon_plot.pdf"), base_height = height, base_width = width)
  
  psimpson<-ggplot(df_simpson, aes(x = df_simpson$class, y = df_simpson$value))
  psimpson<-psimpson+
    geom_dotplot(binaxis="y", stackdir="center", dotsize=1.5,position=position_jitter(w=0.03),aes(fill = df_simpson$class)) + 
    scale_fill_manual(values=c("#DC0000","#3C5488")) + 
    scale_color_manual(values=c("black","black"))+
    stat_summary(fun.data=data_summary2, color="black", size=0.5, width=0.3, fill=c("#DC0000","#3C5488"), alpha=0.3,geom="crossbar")+
    theme_cowplot(font_size = font_size)+
    theme(legend.position = "none")+
    theme(axis.title.x = element_blank())+
    scale_y_continuous(limits = c(0,max(df_simpson$value)*1.15))+
    #annotate("text", x = 1.5, y = max(df_simpson$value)*1.1, label = sig, size = 2)+
    scale_x_discrete(breaks = c(condition1_label,condition2_label), labels = c(condition1_label,condition2_label))+
    labs(y="Simpson Diversity", title = title)+
    theme(plot.title = element_text(hjust = 0.5))
  psimpson
  save_plot(plot = psimpson, filename = paste0(output_folder,"/simpson_plot.pdf"), base_height = height, base_width = width)
  
  genstat<-function(x,y,stat,output_folder,condition1_label,condition2_label){
    xvec<-x[,1]
    yvec<-y[,1]
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
  genstat(x = c1_vaf2_df, y = c2_vaf2_df, stat = "vaf2", condition1_label = condition1_label, condition2_label = condition2_label, output_folder = output_folder);sink()
  genstat(x = c1_shannon_df, y = c2_shannon_df, stat = "shannon", condition1_label = condition1_label, condition2_label = condition2_label,output_folder = output_folder);sink()
  genstat(x = c1_simpson_df, y = c2_simpson_df, stat = "simpson", condition1_label = condition1_label, condition2_label = condition2_label,output_folder= output_folder);sink()
  
}





####function calls####

output_stats(condition1_data = "/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/resubmisssion_analysis/prkcda_AU_12mo_mcs_tvaf_0.003_20200109083036",
             condition1_label = "mcs",
             condition2_data = "/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/resubmisssion_analysis/prkcda_AU_12mo_prkcda_tvaf_0.003_20200109083602",
             condition2_label = "prkcda",
             output_folder = "mcs_prkcda_12mpf",
             height = 2.25,
             width = 2.5,
             font_size = 11,
             title = "12 mpf",
             sig = "")

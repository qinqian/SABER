library(tidyverse)
library(cowplot)

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

long_analysis<-function(inlist,vector_of_mpf,vaf_threshold,title, h,w,fillcolor, fontsize,outfile){
  #initialize the data frame
  df<-read.csv(inlist[1])
  df$mpf<-as.numeric(rep(vector_of_mpf[1], times = nrow(df)))
  df<-df[c(-1,-2),c(-1,-4,-5,-6,-8)]
  colnames(df)<-c("read_count", "variant","specimen","mpf")
  df$vaf<-df$read_count/sum(df$read_count)
  df<-filter(df,df$vaf>=vaf_threshold)
  df$norm_vaf<-df$read_count/sum(df$read_count)
  #load in the data
  for (i in 2:length(inlist)) {
    df_new<-read.csv(inlist[i])
    df_new$mpf<-as.numeric(rep(vector_of_mpf[i], times = nrow(df_new)))
    df_new<-df_new[c(-1,-2),c(-1,-4,-5,-6,-8)]
    colnames(df_new)<-c("read_count", "variant","specimen","mpf")
    df_new$vaf<-df_new$read_count/sum(df_new$read_count)
    df_new<-filter(df_new,df_new$vaf>=vaf_threshold)
    df_new$norm_vaf<-df_new$read_count/sum(df_new$read_count)
    df<-bind_rows(df,df_new)
  }
  df$variant<-factor(df$variant)
  all_vars<-unique(df$variant)
  all_times<-unique(df$mpf)
  filler_df<-data.frame(variant = factor(rep(all_vars,times = length(all_times))),mpf = c(rep(all_times, each = length(all_vars))),vaf = rep(0,times = length(all_vars)*length(all_times)),norm_vaf = rep(0,times = length(all_vars)*length(all_times)))
  df<-bind_rows(df, filler_df)
  df<-df[order(df[,"variant"],-df[,"norm_vaf"]),]
  df$var_mpf<-paste0(df$variant,"_",df$mpf)
  df<-df[!duplicated(df$var_mpf),]
  df1<-filter(df,df$mpf==vector_of_mpf[1])
  df1$rank<-rank(df1$norm_vaf)
  df1$variant<-reorder(df1$variant,df1$rank)
  df$variant<-factor(df$variant, levels = levels(df1$variant))
  plot<-ggplot(df, aes(x = mpf, y = norm_vaf, alpha = variant))+
    geom_area(fill = fillcolor, color = "black")+
    theme_cowplot(font_size = fontsize)+
    theme(legend.position = "none")+
    labs(y = "Cumulative Frequency", x = "mpf", title = title)+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(breaks = vector_of_mpf)
  save_plot(plot = plot, filename = outfile, base_height = h, base_width = w)

}

#function call
long_analysis(inlist = c("/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/pipeline6/prkcda_bleeds/mcs_3mo_tvaf_0.3/thresh_0_tvaf_0.3/crispRvariants_AB048.csv",
                     "/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/pipeline6/prkcda_bleeds/mcs_6mo_tvaf_0.3/thresh_0_tvaf_0.3/crispRvariants_AE006.csv",
                     "/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/pipeline6/prkcda_bleeds/mcs_9mo_tvaf_0.3/thresh_0_tvaf_0.3/crispRvariants_AJ006.csv",
                     "/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/pipeline6/prkcda_bleeds/mcs_12mo_tvaf_0.3/thresh_0_tvaf_0.3/crispRvariants_AU016.csv"),
          vector_of_mpf = c(3,6,9,12),vaf_threshold = 0.02,h = 1.8, w = 1.25,fillcolor = "red",fontsize = 11,outfile = "longitudinal_analysis/testfunction.pdf",title = "")

long_analysis(inlist = c("/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/pipeline6/lgal_3mo_20191218_tvaf0.3/thresh_0_tvaf_0.3/crispRvariants_AB031.csv",
                     "/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/pipeline6/lgal_12mo_20191218_tvaf0.3/thresh_0_tvaf_0.3/crispRvariants_AR032.csv"),
          vector_of_mpf = c(3,12), vaf_threshold = 0.02, h = 2.25, w = 1.5, fillcolor = "red", fontsize = 11, outfile = "longitudinal_analysis_lgal/mcs2.pdf", title = "")

long_analysis(inlist = c("/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/pipeline6/lgal_3mo_20191218_tvaf0.3/thresh_0_tvaf_0.3/crispRvariants_AB040.csv",
                     "/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/pipeline6/lgal_12mo_20191218_tvaf0.3/thresh_0_tvaf_0.3/crispRvariants_AR039.csv"),
          vector_of_mpf = c(3,12),vaf_threshold = 0.02, h = 2.25, w = 1.55, fillcolor = "red", fontsize = 11, outfile = "longitudinal_analysis_lgal/mcs7.pdf", title = "")

long_analysis(inlist = c("/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/pipeline6/lgal_3mo_20191218_tvaf0.3/thresh_0_tvaf_0.3/crispRvariants_AB030.csv",
                     "/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/pipeline6/lgal_12mo_20191218_tvaf0.3/thresh_0_tvaf_0.3/crispRvariants_AR035.csv"),
          vector_of_mpf = c(3,12),vaf_threshold = 0.02, h = 2.25, w = 1.5, fillcolor = "red", fontsize = 11, outfile = "longitudinal_analysis_lgal/mcs4.pdf", title = "")

long_analysis(inlist = c("/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/pipeline6/lgal_3mo_20191218_tvaf0.3/thresh_0_tvaf_0.3/crispRvariants_AB036.csv",
                     "/home/OSUMC.EDU/blas02/pipelines/gestalt_barcodes/pipeline6/lgal_12mo_20191218_tvaf0.3/thresh_0_tvaf_0.3/crispRvariants_AR033.csv"),
          vector_of_mpf = c(3,12),vaf_threshold = 0.02, h = 2.25, w = 1.5, fillcolor = "red", fontsize = 11, outfile = "longitudinal_analysis_lgal/mcs3.pdf", title = "")


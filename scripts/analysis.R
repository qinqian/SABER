args=commandArgs(trailingOnly=TRUE)
####libraries####
library(CrispRVariants)
library("Rsamtools")
library(rtracklayer)
library(GenomicFeatures)
library(gridExtra)
library(GenomicRanges)
library(RColorBrewer)
library(rPython)
library(sqldf)
library(rmarkdown)
library(foreach)
library(doParallel)
library(stringr)
library(tidyverse)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(ggdendro)
library(scales)
library(viridis)
library(ggExtra)
library(ggrepel)
library(boot)
library(yaml)
config <- yaml.load_file(args[1])


####version notes####
# version 3.0 uses paired illumina reads from genewiz and needleall aligner\
# version 3.1 changes thresholding method for stacked bar plots and csvs to proportion of all reads
# version 3.2 parallelizes needleall
# no umis
# You need a genome file in a folder named "genome" and your fastqs in a folder named "fastq" for this version to work
# version 4.0 combines main pipeline with clone comparer
# version 4.1 spits out fraction informative data too
# version 4.2_normal_samples is designed to eval normal samples
# analysis_v1.0 derived from pipeline v4.2
# analysis_v2.0 cleans up analysis_v1.0
# analysis v6 comes from v5 but changes data_out.csv to include sample names and top clones

####input data##################################################################################
#Be sure this is correct before running!!!
#enter the experiment data

genome<-"SABER_pipeline4.fasta" # include .fasta.  Genome file has to be in genome folder.
genome_fp<-paste0(getwd(),"/references/",genome)
output_file<-"crispRvariants"#for the cripsrvariants plot and prefix on sample output csv and vaf plot
group_name<-"crispRvariants"
threshold<-200# number of reads below which they don't appear on the big crispr plot
dt_stamp<-str_replace_all(string = Sys.time(), pattern = "([-\\s:ESTD])",replacement = "")

#new_analysis<-T#if F, then load the R data file you want to work with and then source the script.  If T, be sure the output folder is set correctly

#output_folder<-"resub_vfinal_AG";dir.create(output_folder)#this will make a new output folder each time you run it

if ("tvaf" %in% names(config)) {
  tvaf_preset<-T
  final_tvaf<-config$tvaf
}else{
  tvaf_preset<-F
}
# make working subdirectories and create variables
output_folder<-args[3]
bam_folder<-args[2]
bam_files<-list.files(path = bam_folder,".*.bam$")
bam_fnames<-paste0(bam_folder,"/",bam_files)
mid_names<-tools::file_path_sans_ext(bam_files)


#graphical parameters
gp_crisprplot<-c("width" = 14, "height" = 3, "top.n" = 20)
gp_sparsity<-c("width" = 3.5, "height" = 2.25, "fontsize" = 11)
gp_heatmap<-c("width" = 3.6, "height" = 3.33, "read_cutoff" = 15000, "fontsize" = 11)
gp_sharing_curve_family<-c("width" = 3, "height" = 2.7, "fontsize" = 11)
gp_curve_plot<-c("width" = 3, "height" = 2.7, "fontsize" = 11)
gp_boot_corplot<-c("width" = 3, "height" = 2.7, "fontsize" = 11)
gp_boot_estimate<-c("width" = 3, "height" = 2.7, "fontsize" = 11)
gp_fidist<-c("width" = 3, "height" =2.7, "fontsize" = 11)
gp_histograms<-c("width" = 2.5, "height" = 2.3, "fontsize" = 11)
gp_informative_boxplot<-c("width" = 3, "height" = 2.7, "fontsize" = 11)
gp_sharing_curve_hifi<-c("width" = 3, "height" = 2.7, "fontsize" = 11)
gp_tvaf_product_plot<-c("width" = 3, "height" = 2.7, "fontsize" = 11)
##############################################################################################


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


#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


group_desig<-rep(group_name, times=length(bam_fnames))
md<-read.csv("references/blank_metadata.csv", header = TRUE)
newrow<-data.frame(bamfile=bam_fnames, directory=getwd(),Short.name=mid_names,Targeting.type="",sgRNA1="",sgRNA2="",Group=group_desig)
md<-rbind(md,newrow)


#create target region
gd <- rtracklayer::import("references/SABER2.bed")
gdl <- GenomicRanges::resize(gd, width(gd) + 0, fix = "center") #resize region for analysis
reference0<-read_file("references/SABER2_ref.txt")
reference1<-substr(reference0,1,310)
reference<-Biostrings::DNAString(reference1)
reference

# make the crispr set
crispr_set <- readsToTarget(bam_fnames,
                            target = gd,
                            reference = reference,
                            names = md$Short.name,
                            target.loc = 16,
                            collapse.pairs = FALSE,
                            split.snv=FALSE)#split.snv=FALSE adds SNVs into no variant count





#####generate informative, no variant and common variant tables with vaf and paf thresholds####
generate_v1<-function(x,thresh,vaf_thresh,paf_thresh,vc_prop,output_folder,output_file){
  df<-x#isolate df
  rnames<-row.names(df)#duplicate row names as new variable
  df$rnames<-rnames#add row name variabe to df
  df<-df[which(df$rnames!="Other"),]# removes other
  ##important:  CrispRVariants proportion generates PERCENT values, not FRACTION values.  So TVAF = 1 => 1% vaf, not 100% VAF.
  #common_vars<-names(which(rowSums(vc_prop>vaf_thresh)/ncol(vc_prop)>(paf_thresh)))#returns the names of rows with a propotion greater than vaf_thresh occuring in more than paf_thresh
  common_vars<-names(which(rowSums(vc_prop>vaf_thresh)>1))#returns the names of rows with a propotion greater than vaf_thresh occuring in more than paf_thresh
  ##end important
  remove<-c("no variant")# removes only no variant.  Leaves "Other" in.  Maybe need to pull it out.
  common_vars<-setdiff(common_vars,remove)
  for (i in 1:length(common_vars)) df$rnames<-gsub(paste0("^",common_vars[i],"$"),"common variant",df$rnames)#replaces all common variants as "common variant"
  #df<-df[which(df[1]/sum(df[1])>thresh | df$rnames=="no variant"),]# removes reads below threshold and keeps no variant##20190509 changed to use proportional threshold
  common_variant_sum<-sum(df[1][which(df$rnames=="common variant"),])#sums all common variant calls
  df$is_common_variant<-NA#adds a column that will hold the common variant marker
  df[nrow(df) + 1,] = list(common_variant_sum,"common variant sum","yes")#adds a new row to df with summed novariant calls.  Yes is a no variant marker
  df$is_no_variant<-NA#adds a column that will hold the no_variant marker
  df[which(df$rnames == "no variant"),4]<-"yes"
  df<-df[which(df$rnames!="common variant"),,drop=FALSE]#drop all constituents of no variant sum
  df<-df[which(df$rnames!="Other"),,drop=FALSE]#drop other reads
  df<-df[order(df$is_no_variant, df$is_common_variant,-df[1]),]# sort on no variant marker then on read counts
  df<-df[which(df[1]/sum(df[1])>thresh | df$is_no_variant=="yes" | df$is_common_variant=="yes"),]# removes reads below read count threshold, keeping no variant calls and common variant calls
  n_variants<-length(df$rnames)#total number of reads in each specimen
  total_reads<-sum(df[1])#sum of reads for each specimen
  allele_freq<-df[1]/total_reads# calculate allele freq by line
  names(allele_freq)[1] <- "allele_freq"# rename column for allele_freq
  df<-cbind(df,allele_freq)#add allele_freq column to the working df
  specimen_id<-names(df)[1]# pull the specimen name from the column title
  specimen<-rep_len(specimen_id, n_variants)# make a new column for the specimen id
  df<-cbind(df,specimen)# bind it to the working df
  #fix levels
  index<-as.factor(c(1:length(df$rnames)))#create index string as long as the number of variants, make it as a factor
  index<-factor(index, levels=rev(levels(index)))#reverse the levels
  df<-cbind(df,index)#bind it to the working df
  dir.create(paste0(output_folder,"/thresh_",thresh,"_tvaf_",vaf_thresh))
  write.csv(df,file = paste0(output_folder,"/thresh_",thresh,"_tvaf_",vaf_thresh,"/",output_file,"_",specimen_id,".csv"))
}

####generate variant tables####
vc_all <- as.data.frame(variantCounts(crispr_set), stringsAsFactors=FALSE)#big data frame of read counts
vc_list<-list()
for(i in 1:length(mid_names)){vc_list[[length(vc_list)+1]]<-vc_all[i]}#puts individual samples into a list of vectors
vc_prop<-variantCounts(crispr_set, result = "proportions")# generates a matrix of variant allele percents, similar to vc_all
vc_prop_decimal<-vc_prop/100

thresh<-0;tvaf<-1;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, vc_prop = vc_prop_decimal,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-0.3;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, vc_prop = vc_prop_decimal,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-0.1;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, vc_prop = vc_prop_decimal,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-0.03;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, vc_prop = vc_prop_decimal,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-0.01;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, vc_prop = vc_prop_decimal,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-0.003;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, vc_prop = vc_prop_decimal,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-0.001;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, vc_prop = vc_prop_decimal,output_folder = output_folder, output_file = output_file)
thresh<-0;tvaf<-0.0003;parLapply(cl,vc_list,generate_v1,thresh = thresh, vaf_thresh = tvaf, vc_prop = vc_prop_decimal,output_folder = output_folder, output_file = output_file)

#stop cluster
stopCluster(cl)


####sparsity####
#calculate sparsity of the original dataframe
#load the data
sparsity<-as.data.frame(colSums(vc_all==0)/(colSums(vc_all==0) + colSums(vc_all!=0)))#95-99% sparse
colnames(sparsity)<-"Sparsity"
sparsity<-sparsity %>%
  rownames_to_column("Sample")
sparsity$nonsparse<-1-sparsity$Sparsity
sparsity$index<-factor(rank(sparsity$nonsparse, ties.method = "first"))

# plot sparsity
sparse_plot<-ggplot(sparsity, aes(x = sparsity$index, y = sparsity$Sparsity))
sparse_plot<-sparse_plot+
  geom_bar(stat = "identity", color = "black", fill = "white")+
  scale_x_discrete(breaks = c(5,10,15,20,25,30,35))+
  theme_cowplot(font_size = unname(gp_sparsity["fontsize"]))+
  ylab(label = "Sparsity")+
  xlab(label = "Sample Index")
sparse_plot
save_plot(plot = sparse_plot, filename = paste0(output_folder,"/sparsity_plot.pdf"), base_height = unname(gp_sparsity["height"]), base_width = unname(gp_sparsity["width"]))


####big heatmap####
# plot heatmap of all variants present over 20000 reads in any sample.  i.e. modify the crispr plots heatmap
heat1_df<-vc_all
heat1_df<-heat1_df %>%
  rownames_to_column("Variants") %>%
  filter_if(is.numeric, any_vars(. >unname(gp_heatmap["read_cutoff"])))%>%
  column_to_rownames("Variants")

#cluster
heat1_matrix<-as.matrix(heat1_df)
heat1_dendro<-as.dendrogram(hclust(d=dist(x=heat1_matrix)))

x <- heat1_matrix
dd.col <- as.dendrogram(hclust(dist(x)))
dd.row <- as.dendrogram(hclust(dist(t(x))))
dx <- dendro_data(dd.row)
dy <- dendro_data(dd.col)

# helper function for creating dendograms
ggdend <- function(df) {
  ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}

# x/y dendograms
px <- ggdend(dx$segments)
py <- ggdend(dy$segments) + coord_flip()


heat1_df<-heat1_df %>%
  rownames_to_column("Variants")
heat1_df_long<-melt(heat1_df)

# #reorder heatmap
heat1_orderx<-order.dendrogram(dd.col)
heat1_ordery<-order.dendrogram(dd.row)
heat1_df_long$Variants<-factor(x=heat1_df_long$Variants, levels = heat1_df$Variants[rev(heat1_orderx)], ordered = TRUE)
heat1_df_long$variable<-factor(x=heat1_df_long$variable, levels = colnames(vc_all)[heat1_ordery], ordered = TRUE)
heat1_df_long$rank<-rank(heat1_df_long$variable, ties.method = "min")
heat2_df_long<-unique(heat1_df_long[,c(2,4)])
heat2_df_long<-arrange(heat2_df_long,rank)
heat2_df_long$index<-factor(seq.int(nrow(heat2_df_long)))
heat3_df_long<-left_join(x = heat1_df_long, y = heat2_df_long, by = "variable")

heat3_df_long$scaled_value<-scale(heat3_df_long$value)

heat1<-ggplot(data = heat3_df_long, aes(x = variable, y=Variants, fill = value))
heat1<-heat1+
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu")+
  theme_cowplot(font_size = unname(gp_heatmap["fontsize"]))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 8))+
  theme(axis.text.y = element_blank())+
  theme(axis.ticks.y = element_blank())+
  theme(axis.line.y = element_blank())+
  guides(fill=guide_colourbar(title="Reads", reverse = F))+
  xlab("Sample")+
  ylab("Variant")+
  coord_fixed()
heat1
save_plot(plot = heat1, filename = paste0(output_folder,"/full_heatmap.pdf"), base_width = unname(gp_heatmap["width"]), base_height = unname(gp_heatmap["height"]))##check to be sure what the cutoff is!!


#####analysis of tvaf, cvsa settings####
cor_func3<-function(include_uninformative, thresh, tvaf, method, output_folder){
  input_directory<-paste0(output_folder,"/thresh_",thresh,"_tvaf_",tvaf)
  inlist<-list.files(input_directory,full.names = TRUE, pattern = "\\d.csv$")
  df_cor<-data.frame(Variant = character(), Count = integer(), Specimen = character())
  for (i in 1:length(inlist)){
    df_new<-read.csv(file = inlist[i], header = TRUE)
    df_new<-df_new[,c(3,2,7)]
    colnames(df_new)<-c("Variant", "Count", "Specimen")
    df_cor<-rbind(df_cor,df_new)
  }
  df_cor_wide0<-reshape(df_cor, idvar = "Variant", timevar = "Specimen", direction = "wide")
  df_cor_wide<-df_cor_wide0
  df_cor_wide[is.na(df_cor_wide)]<-0
  colnames(df_cor_wide)<-substr(colnames(df_cor_wide),7,11)
  rownames(df_cor_wide)<-df_cor_wide[,1]
  df_cor_wide<-df_cor_wide[,-1]
  if (include_uninformative==TRUE) {df_cor_wide1<-df_cor_wide} else {df_cor_wide1<-df_cor_wide[-c(1,2),]}


  stat_matrix<-matrix(vector(),nrow = ncol(df_cor_wide1), ncol = ncol(df_cor_wide1))
  colnames(stat_matrix)<-colnames(df_cor_wide1)
  rownames(stat_matrix)<-colnames(df_cor_wide1)

  for (i in 1:ncol(df_cor_wide1)) {
    for (j in 1:ncol(df_cor_wide1)) {
      a<-pmin(df_cor_wide1[,i],df_cor_wide1[,j])
      b<-pmax(df_cor_wide1[,i],df_cor_wide1[,j])
      c<-(sum(a)*2)/(sum(a)+sum(b))
      stat_matrix[i,j]<-c

    }

  }

  reorder_cormat <- function(x){
    # Use correlation between variables as distance
    dd <- as.dist((1-x)/2)
    hc <- hclust(dd)
    x <-x[hc$order, hc$order]
  }

  get_lower_tri <- function(x){
    x[upper.tri(x)]<- NA
    return(x)
  }

  sm_reordered<-reorder_cormat(stat_matrix)
  sm_lower<-get_lower_tri(sm_reordered)

  msm<-melt(sm_lower)

  msm_plot<-ggplot(data = msm, aes(x = msm$Var1, y = msm$Var2, fill = msm$value))
  msm_plot<-msm_plot+
    geom_tile()+
    scale_fill_distiller(palette = "RdYlBu", na.value = "transparent",guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE),limits = c(0,1))+
    theme_cowplot(font_size = 11)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))+
    theme(axis.text.y = element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(axis.line.y = element_blank())+
    theme(axis.ticks.y = element_blank())+
    xlab("Sample")+
    theme(axis.title.y = element_blank())+
    theme(plot.title = element_blank())+
    theme(plot.subtitle = element_blank())+
    guides(fill=guide_colorbar(title="Sharing\nFactor", reverse = F))+
    theme(legend.position = c(0,0.7))+
    coord_fixed()
  msm_plot

  msm1<-msm[complete.cases(msm),]
  msm1$rank<-rank(msm1$value,ties.method = "random")
  msm1$tvaf<-tvaf
  msm1$include_uninformative<-include_uninformative

  #calculate Fraction informative
  FI<-as.data.frame(colSums(df_cor_wide[-c(1,2),])/colSums(df_cor_wide))#ratio of reads excluding no variant and common variant to all reads in each sample.
  colnames(FI)<-"Fraction_Informative"
  FI$sample_x<-rownames(FI)

  avg_share<-round(mean(msm[which(msm$Var1!=msm$Var2),3],na.rm = TRUE),5)#averages sharing factor excluding identical comparisons

# avg_share<-round(mean(msm[which(msm$Var1!=msm$Var2),3],na.rm = TRUE),4)#averages pearson correlation excluding identical comparisons

  mean_inf<-round(mean(FI$Fraction_Informative),3)

  save_plot(plot = msm_plot, filename = paste0(input_directory,"/msm_plot_tvaf_",tvaf,"_uninform_",include_uninformative,".pdf"), base_width = 2.75)
  nubbin<-data.frame(threshold = thresh, tvaf = tvaf, mean_inf = mean_inf, avg_share = avg_share)
  write.csv(x = nubbin, file = paste0(input_directory,"/nubbin_",include_uninformative,"_",method,".csv"))
  return(msm1)
}

msm_1_true<-cor_func3(include_uninformative = TRUE, thresh = 0, tvaf = 1, method = "sf", output_folder = output_folder)

msm_1_false<-cor_func3(include_uninformative = FALSE, thresh = 0, tvaf = 1, method = "sf", output_folder = output_folder)
msm_0.3_false<-cor_func3(include_uninformative = FALSE, thresh = 0, tvaf = 0.3, method = "sf", output_folder = output_folder)
msm_0.1_false<-cor_func3(include_uninformative = FALSE, thresh = 0, tvaf = 0.1, method = "sf", output_folder = output_folder)
msm_0.03_false<-cor_func3(include_uninformative = FALSE, thresh = 0, tvaf = 0.03, method = "sf", output_folder = output_folder)
msm_0.01_false<-cor_func3(include_uninformative = FALSE, thresh = 0, tvaf = 0.01, method = "sf", output_folder = output_folder)
msm_0.003_false<-cor_func3(include_uninformative = FALSE, thresh = 0, tvaf = 0.003, method = "sf", output_folder = output_folder)
msm_0.001_false<-cor_func3(include_uninformative = FALSE, thresh = 0, tvaf = 0.001, method = "sf", output_folder = output_folder)
msm_0.0003_false<-cor_func3(include_uninformative = FALSE, thresh = 0, tvaf = 0.0003, method = "sf", output_folder = output_folder)

####curve plot with sharing statistic####
colors<-c("1_TRUE" = "#1F77B4","1_FALSE" = "#FF7F0E","0.3_FALSE" = "#2CA02C","0.1_FALSE" = "#D62728","0.03_FALSE" = "#9467BD","0.01_FALSE" = "#8C564B","0.003_FALSE" = "#E377C2","0.001_FALSE" = "#7F7F7F","3e-04_FALSE" = "#BCBD22")
labels<-c("1_TRUE" = "All Vars","1_FALSE" = "All Edited Vars","0.3_FALSE" = expression(paste(theta[V],"= 0.3")),"0.1_FALSE" = expression(paste(theta[V],"= 0.1")),"0.03_FALSE" = expression(paste(theta[V],"= 0.03")),"0.01_FALSE" = expression(paste(theta[V],"= 0.01")),"0.003_FALSE" = expression(paste(theta[V],"= 0.003")),"0.001_FALSE" = expression(paste(theta[V],"= 0.001")),"3e-04_FALSE" = expression(paste(theta[V],"= 0.0003")))

sharing_inlist<-list.files(output_folder, recursive = TRUE, pattern = "nubbin_FALSE_sf.csv",full.names = TRUE)
sharing_inlist1<-as.character(sort(factor(sharing_inlist, levels = factor(sharing_inlist[c(6,8,5,7,4,3,2,1)]))))#resorts

curve_df0<-read.csv(sharing_inlist1[1])
for (i in 2:length(sharing_inlist1)){
  curve_df0<-rbind(curve_df0,read.csv(sharing_inlist1[i]))
}

curve_df1<-curve_df0

curve_df1.1<-curve_df1 %>% arrange(desc(tvaf)) %>% rownames_to_column("condition")
curve_df1.2<-curve_df1.1[,-2]
curve_df1.2[,1] <- as.factor(curve_df1.2[,1])
curve_df2<-curve_df1.2

curve_df4<-curve_df2

curve_df4$condition<-dplyr::recode(curve_df4$condition,
                                   "2" = "0.3_FALSE",
                                   "8" = "3e-04_FALSE",
                                   "3" = "0.1_FALSE",
                                   "1" = "1_FALSE",
                                   "4" = "0.03_FALSE",
                                   "5" = "0.01_FALSE",
                                   "6" = "0.003_FALSE",
                                   "7" = "0.001_FALSE")
curve_df5<-curve_df4[!duplicated(curve_df4$avg_share),]
curve_df6<-curve_df5 %>% mutate(mean_inf_scaled = scale(mean_inf),avg_unshared_scaled = scale(1-avg_share),product = mean_inf_scaled*avg_unshared_scaled)

curve_plot5<-ggplot(curve_df5, aes(x = avg_share, y = mean_inf, fill = condition))
curve_plot5<-curve_plot5+
  geom_point(shape = 21,size = 3, color = "black", alpha = 0.7)+
  theme_cowplot(font_size = unname(gp_curve_plot["fontsize"]))+
  scale_fill_manual(values=  colors, labels = labels)+
  theme(legend.text = element_text(size = 8))+
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(keyheight = 0.1, keywidth = 0.025, default.unit = "inch",ncol = 3, override.aes = list(size = 2)))+
  labs(x = "Mean Sharing Factor", y = expression(paste("Mean ",Phi)), fill ="")
curve_plot5
save_plot(plot = curve_plot5, filename = paste0(output_folder,"/curve_plot.pdf"), base_height = unname(gp_curve_plot["height"]), base_width = unname(gp_curve_plot["width"]))

curve_plot6<-ggplot(curve_df6, aes(x = tvaf, y = product, fill = condition))
curve_plot6<-curve_plot6+
  geom_point(shape = 22, size = 3, color = "black", alpha = 0.7)+
  theme_cowplot(font_size = unname(gp_tvaf_product_plot["fontsize"]))+
  scale_fill_manual(values = colors, labels = labels)+
  theme(legend.text = element_text(size = 8))+
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(keyheight = 0.1, keywidth = 0.025, default.unit = "inch",ncol = 3, override.aes = list(size = 2)))+
  scale_x_log10()+
  labs(x = expression(paste(theta[V])),y = "Z-score Product", fill ="")
curve_plot6
save_plot(plot = curve_plot6, filename = paste0(output_folder,"/tvaf_product_plot.pdf"), base_height = unname(gp_tvaf_product_plot["height"]), base_width = unname(gp_tvaf_product_plot["width"]))


####make family of curves for sharing factor at each tvaf###
msm_list<-list(msm_1_true,msm_1_false,msm_0.3_false,msm_0.1_false,msm_0.03_false,msm_0.01_false,msm_0.003_false,msm_0.001_false,msm_0.0003_false)

msm_full<-dplyr::bind_rows(msm_list)
msm_full$condition<-paste0(msm_full$tvaf,"_",msm_full$include_uninformative)
unique(msm_full$condition)
msm_full$condition<-factor(msm_full$condition, levels = unique(msm_full$condition))

sf_curve_plot<-msm_full %>% filter(msm_full$tvaf %in% curve_df5$tvaf)
paired_share_factors<-ggplot(sf_curve_plot, aes(x = rank, y = value, fill = condition, group = condition))+#reorder(rank,-rank)
  geom_point(shape=21, color = "grey20", size = 1, alpha = 0.4, stroke = 0.25)+
  scale_x_reverse()+
  theme_cowplot(font_size = unname(gp_sharing_curve_family["fontsize"]))+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  scale_y_log10()+
  scale_fill_manual(values = colors, labels = labels)+
  theme(legend.position = c(0.025,0.25))+
  theme(legend.text = element_text(size = 8))+
  guides(fill = guide_legend(keyheight = 0.1, default.unit = "inch",ncol = 1, override.aes = list(alpha = 0.7,size = 2)))+
  labs(x = "Ranked Sample Pair", y = "Sharing Factor", fill = "")
  paired_share_factors
save_plot(paired_share_factors, filename = paste0(output_folder,"/sharing_factor_curve_family.pdf"), base_width = unname(gp_sharing_curve_family["width"]), base_height = unname(gp_sharing_curve_family["height"]))


####control point to bring in preset threshold if you want to do that####
if (tvaf_preset==F) {
  final_tvaf<-curve_df6[order(-curve_df6$product),3][1]#tvaf at max z score product
  output_folder_old<-output_folder
  output_folder<-paste0(output_folder_old,"_tvaf_",final_tvaf,"_",dt_stamp)
  file.rename(from = output_folder_old,to = output_folder)
}

if (tvaf_preset==T) {
  output_folder_old<-output_folder
  output_folder<-paste0(output_folder_old,"_tvaf_",final_tvaf,"_",dt_stamp)
  file.rename(from = output_folder_old,to = output_folder)
}

####bootstrap analysis on whole dataset ####
#helper function for bootstrap analysis
bootf<-function(data, indices){
  dt<-data[indices,]
  c(mean(dt[,2]),sd(dt[,3]))
}

#generate dataframe of fi, clones over 2% and sample
#load the data
infiles<-list.files(paste0(output_folder,"/thresh_0_tvaf_",final_tvaf), pattern = "\\d.csv", full.names = TRUE)

bootdata_0<-read.csv(infiles[1])
fi<-1-sum(c(bootdata_0[1,2],bootdata_0[2,2]))/sum(bootdata_0[2])
vaf<-c(NA,NA,bootdata_0[-c(1:2),2]/sum(bootdata_0[-c(1:2),2]))
bootdata_1<-bootdata_0[,c(7,2,3)]
colnames(bootdata_1)<-c("specimen", "counts", "variant")
bootdata_1$FI<-rep(x = fi,times = nrow(bootdata_1))
bootdata_1$vaf<-vaf
bootdata_1$vaf2<-rep(sum(vaf>0.02, na.rm = TRUE), times = nrow(bootdata_1))
for (i in 2:length(infiles)){
  newdf<-read.csv(infiles[i])
  fi<-1-sum(c(newdf[1,2],newdf[2,2]))/sum(newdf[2])
  vaf<-c(NA,NA,newdf[-c(1:2),2]/sum(newdf[-c(1:2),2]))
  newdf<-newdf[,c(7,2,3)]
  colnames(newdf)<-c("specimen", "counts", "variant")
  newdf$FI<-rep(x = fi,times = nrow(newdf))
  newdf$vaf<-vaf
  newdf$vaf2<-rep(sum(vaf>0.02, na.rm = TRUE), times = nrow(newdf))
  bootdata_1<-rbind(bootdata_1,newdf)
}
bootdata_2<-bootdata_1[,c(1,4,6)]#gets rid of the variant-level data for each sample since you don't need it here
bootdata_3<-unique(bootdata_2)#gets rid of duplicate rows

set.seed(12345)
bootstrap<-boot(bootdata_3, bootf, R=1000)

boot_df<-as.data.frame(bootstrap$t)
linear_mod<-lm(boot_df, formula = (boot_df$V2 ~ boot_df$V1))
cor(boot_df$V1, boot_df$V2, method = "pearson")
cor_return<-cor.test(boot_df$V1, boot_df$V2, method = "pearson")
summary(linear_mod)


boot_cor_plot<-ggplot(boot_df, aes(x = boot_df$V1, y = boot_df$V2))
boot_cor_plot<-boot_cor_plot+
  geom_density2d(color = "black")+
  theme_cowplot(font_size = unname(gp_boot_corplot["fontsize"]))+
  labs(x = expression(paste("Mean ",Phi)), y = expression(paste("S.D. ",B["0.02"],sep = "")))
boot_cor_plot
save_plot(plot = boot_cor_plot, filename = paste0(output_folder,"/bootstrap_correlation_sf.pdf"), base_height = unname(gp_boot_corplot["height"]), base_width = unname(gp_boot_corplot["width"]))


boot_ci<-boot.ci(bootstrap, index = 1)#makes a list of the bootstrap confidence intervals
boot_df1<-boot_df
boot_df1$v1bca<-rep(boot_ci$bca[1,4], times = nrow(boot_df1))#appends the lower bound of the 95%ci for the V1 value (mean inf) to the df
boot_df1$lower<-boot_df1$V1<boot_df1$v1bca

boot_plot<-ggplot(boot_df, aes(x = boot_df$V1))
boot_plot<-boot_plot+
  geom_histogram(binwidth = 0.01, fill = "white", color = "black")+
  geom_vline(aes(xintercept = mean(boot_df$V1), color = "Estimate"), linetype = "dashed")+#plots bootstrap estimate of mean informative fraction
  geom_vline(aes(xintercept = boot_ci$bca[1,4], color = "Lower\nBound"), linetype = "dashed")+#plots lower bound of 95% ci
  scale_color_manual(name = "", values = c('Lower\nBound' = "#DC0000",Estimate = "#3c5488"))+
  theme_cowplot(font_size = unname(gp_boot_estimate["fontsize"]))+
  theme(legend.position = c(0.65,0.95))+
  theme(legend.text = element_text(size = 8))+
  guides(color = guide_legend(label.position = "left", label.hjust = 1))+
  labs(y = "Replicates", x = expression(paste("Mean ",Phi)))
boot_plot
save_plot(plot = boot_plot, filename = paste0(output_folder,"/bootstrap_estimate_mif_sf.pdf"), base_height = unname(gp_boot_estimate["height"]), base_width = unname(gp_boot_estimate["width"]))

returnlist<-list("bootstrap report:  " = bootstrap,
                 "correlation of FI and SD over replicates:  " = cor_return,
                 "95% ci of bootstrap values:  " = boot_ci,
                 "estimate:  " = mean(boot_df$V1))
capture.output(print(returnlist), file = paste0(output_folder,"/bootstrap_report.txt"))

####re-run training samples using final tvaf condition and identify informative vs uninformative samples####
#load the data
infiles<-list.files(paste0(output_folder,"/thresh_0_tvaf_",final_tvaf), pattern = "\\d.csv", full.names = TRUE)

fidf_0<-read.csv(infiles[1])
fi<-1-sum(c(fidf_0[1,2],fidf_0[2,2]))/sum(fidf_0[2])
fidf_1<-fidf_0[,c(7,2,3)]
colnames(fidf_1)<-c("specimen", "counts", "variant")
fidf_1$FI<-rep(x = fi,times = nrow(fidf_0))
for (i in 2:length(infiles)){
  newdf<-read.csv(infiles[i])
  fi<-1-sum(c(newdf[1,2],newdf[2,2]))/sum(newdf[2])
  newdf<-newdf[,c(7,2,3)]
  colnames(newdf)<-c("specimen", "counts", "variant")
  newdf$FI<-rep(x = fi,times = nrow(newdf))
  fidf_1<-rbind(fidf_1,newdf)
}
fidf_2<-fidf_1[,c(1,4)]#gets rid of the variant-level data for each sample since you don't need it here
fidf_3<-unique(fidf_2)#gets rid of duplicate rows
fidf_3$rank<-rank(-(fidf_3$FI), ties.method = "first")#ranks samples based on fraction informative with 1 = most informative

#use bootstrap to define cutoff
cutoff<-boot_ci$bca[1,4]
top_boot_clust<-fidf_3[which(fidf_3$FI>=cutoff),]
top_boot_clust$boot_clust<-rep("High FI", times = nrow(top_boot_clust))
bottom_boot_clust<-fidf_3[which(fidf_3$FI<cutoff),]
bottom_boot_clust$boot_clust<-rep("Low FI", times = nrow(bottom_boot_clust))
fidf_3.1<-rbind(top_boot_clust,bottom_boot_clust)

#plot the fi distribution
fi_dist<-ggplot(fidf_3.1, aes(x=reorder(fidf_3.1$specimen,fidf_3.1$rank), y=fidf_3.1$FI, fill=fidf_3.1$boot_clust))
fi_dist<-fi_dist+
  geom_bar(stat = "identity", color = "black")+
  scale_fill_manual(values = c("#72ce55","#fde725"),name = "",labels = c(expression(paste("High ",Phi),paste("Low ",Phi))))+
  geom_hline(yintercept = cutoff, color = "#DC0000", linetype = "dashed")+
  theme_cowplot(font_size = unname(gp_fidist["fontsize"]))+
  theme(legend.position = c(0.7,0.95))+
  theme(legend.text = element_text(size = 8))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 8))+
  labs(x="Sample", y=expression(paste(Phi)))
fi_dist
save_plot(plot = fi_dist, filename = paste0(output_folder,"/fidist_sf.pdf"), base_width = unname(gp_fidist["width"]), base_height = unname(gp_fidist["height"]))

# save_plot(plot = fi_dist, filename = paste0(output_folder,"/fidist_sharing.pdf"), base_width = 3.5, base_height = 2.75)

informative_df<-top_boot_clust
informative_spec_vec<-informative_df$specimen

all_samples<-fidf_3.1$specimen
informative_samples<-filter(fidf_3.1, fidf_3.1$boot_clust=="High FI")

####quantifying alleles from thresholded data in informative samples####
#reassemble the long form data frame of read counts using only informative samples
inlist<-as.data.frame(cbind(as.character(sort(all_samples)), infiles))
colnames(inlist)<-c("specimen", "path")
inlist_informative<-left_join(informative_samples,inlist, by = "specimen")
informative_path<-as.character(inlist_informative$path)
qdf0<-read.csv(informative_path[1])
colnames(qdf0)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
for (i in 2:length(informative_path)) {
  qdf_new<-read.csv(informative_path[i])
  colnames(qdf_new)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
  qdf0<-rbind(qdf0,qdf_new)
}
qdf1<-qdf0[,-1]

# print out an example for the paper
# qdf_AB031<-filter(qdf1,qdf1$specimen=="AB031")
# qdf_AB031$log10_read_count<-log10(qdf_AB031$read_count)
# qdf_AB031_inf<-qdf_AB031[-c(1,2),];qdf_AB031_inf$index<-qdf_AB031_inf$index-2
# qdf_AB031_inf_filt<-filter(qdf_AB031_inf, qdf_AB031_inf$index<26)
# qdf_AB031_inf_filt$index<-as.factor(qdf_AB031_inf_filt$index)
# AB031_plot<-ggplot(data = qdf_AB031_inf_filt, aes(x = qdf_AB031_inf_filt$index, y = qdf_AB031_inf_filt$log10_read_count))
# AB031_plot<-AB031_plot+
#   theme_cowplot(font_size = 11)+
#   geom_bar(stat = "identity", color = "black", fill = "white")+
#   scale_x_discrete(breaks = c(5,10,15,20,25))+
#   labs(x="Variant Index", y = expression(paste(Log[10]," Read Count",sep="")))#does not recalculate VAF with exclusion of no variant and common variant
# AB031_plot
# save_plot(filename = paste0(output_folder,"/AB031_example_sf_2e341.pdf"), plot = AB031_plot, base_width = 2.5, base_height = 2.3)


#generate a vector of shannon entropy values for the informative training set
spec_list<-unique(qdf1$specimen)

shannon_df <- data.frame(value=numeric(),value_type = character(),sample_name=character(),clone1=character(),clone2=character(),clone3=character(),clone4=character(),clone5=character(),stringsAsFactors=FALSE)
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  shannon_new<-data.frame(value=numeric(),value_type = character(),sample_name=character(),clone1=character(),clone2=character(),clone3=character(),clone4=character(),clone5=character(),stringsAsFactors=FALSE)
  shannon_new[1,1]<-shannon_func(qdf_filt_inf[1])
  shannon_new$value_type<-"shannon"
  shannon_new$sample_name<-unique(qdf_filt_inf$specimen)
  shannon_new$clone1<-qdf_filt_inf[1,2]
  shannon_new$clone2<-qdf_filt_inf[2,2]
  shannon_new$clone3<-qdf_filt_inf[3,2]
  shannon_new$clone4<-qdf_filt_inf[4,2]
  shannon_new$clone5<-qdf_filt_inf[5,2]
  shannon_df<-rbind(shannon_df,shannon_new)
}
shannon_df$value
shapiro.test(shannon_df$value)

#generate a vector of simpson diversity indices for the informative training set
simpson_df <- data.frame(value=numeric(),value_type = character(),sample_name=character(),clone1=character(),clone2=character(),clone3=character(),clone4=character(),clone5=character(),stringsAsFactors=FALSE)
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  simpson_new<-data.frame(value=numeric(),value_type = character(),sample_name=character(),clone1=character(),clone2=character(),clone3=character(),clone4=character(),clone5=character(),stringsAsFactors=FALSE)
  simpson_new[1,1]<-simpson_func(qdf_filt_inf[1])
  simpson_new$value_type<-"simpson"
  simpson_new$sample_name<-unique(qdf_filt_inf$specimen)
  simpson_new$clone1<-qdf_filt_inf[1,2]
  simpson_new$clone2<-qdf_filt_inf[2,2]
  simpson_new$clone3<-qdf_filt_inf[3,2]
  simpson_new$clone4<-qdf_filt_inf[4,2]
  simpson_new$clone5<-qdf_filt_inf[5,2]
  simpson_df<-rbind(simpson_df,simpson_new)
}
simpson_df$value
shapiro.test(simpson_df$value)

#counts barcodes over 2%
vaf2_df <- data.frame(value=numeric(),value_type = character(),sample_name=character(),clone1=character(),clone2=character(),clone3=character(),clone4=character(),clone5=character(),stringsAsFactors=FALSE)
for (i in 1:length(spec_list)){
  qdf_filt<-filter(qdf1, qdf1[,6]==spec_list[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  qdf_filt_inf$adj_allele_freq<-qdf_filt_inf$read_count/sum(qdf_filt_inf$read_count)#recalculates allele frequencies from read counts without no variant and common variant
  qdf_2percent<-filter(qdf_filt_inf, qdf_filt_inf$adj_allele_freq>=0.02)
  value<-nrow(qdf_2percent)
  vaf2_new<-data.frame(value=numeric(),value_type = character(),sample_name=character(),clone1=character(),clone2=character(),clone3=character(),clone4=character(),clone5=character(),stringsAsFactors=FALSE)
  vaf2_new[1,1]<-value
  vaf2_new$value_type<-"vaf2"
  vaf2_new$sample_name<-unique(qdf_filt_inf$specimen)
  vaf2_new$clone1<-qdf_filt_inf[1,2]
  vaf2_new$clone2<-qdf_filt_inf[2,2]
  vaf2_new$clone3<-qdf_filt_inf[3,2]
  vaf2_new$clone4<-qdf_filt_inf[4,2]
  vaf2_new$clone5<-qdf_filt_inf[5,2]
  vaf2_df<-rbind(vaf2_df,vaf2_new)
}
vaf2_df$value
shapiro.test(vaf2_df$value)


#compile and summarize data

p2<-ggplot(vaf2_df, aes(x = vaf2_df$value))
p2<-p2+
  geom_histogram(aes(y = ..density..), color = "black", fill = "#DC0000", binwidth = 1) +
  geom_density(alpha=.2, fill="#DC0000")+
  theme_cowplot(font_size = unname(gp_histograms["fontsize"]))+
  theme(legend.position = "none")+
  scale_x_continuous(limits = c(0,max(vaf2_df$value)*1.2))+
  labs(x="Clones with VAF>0.02", y="Density")
p2
save_plot(plot = p2, filename = paste0(output_folder,"/vaf_0.02.pdf"), base_width = unname(gp_histograms["width"]), base_height = unname(gp_histograms["height"]))

p3<-ggplot(shannon_df, aes(x = shannon_df$value))
p3<-p3+
  geom_histogram(aes(y = ..density..), color = "black", fill = "#DC0000", binwidth = 0.5) +
  geom_density(alpha=.2, fill="#DC0000")+
  theme_cowplot(font_size = unname(gp_histograms["fontsize"]))+
  theme(legend.position = "none")+
  scale_x_continuous(limits = c(0,max(shannon_df$value)*1.2))+
  labs(x="Shannon Entropy", y="Density")
p3
save_plot(plot = p3, filename = paste0(output_folder,"/shannon.pdf"), base_width = unname(gp_histograms["width"]), base_height = unname(gp_histograms["height"]))

p4<-ggplot(simpson_df, aes(x = simpson_df$value))
p4<-p4+
  geom_histogram(aes(y = ..density..), color = "black", fill = "#DC0000", binwidth = 1) +
  geom_density(alpha=.2, fill="#DC0000")+
  theme_cowplot(font_size = unname(gp_histograms["fontsize"]))+
  theme(legend.position = "none")+
  scale_x_continuous(limits = c(0,max(simpson_df$value)*1.2))+
  labs(x="Simpson Diversity", y="Density")
p4
save_plot(plot = p4, filename = paste0(output_folder,"/simpson.pdf"), base_width = unname(gp_histograms["width"]), base_height = unname(gp_histograms["height"]))

####go back to whole training set (informative + uninformative) and plot fraction informative vs clones over 2%
inlist_informative2<-left_join(fidf_3.1,inlist, by = "specimen")#all samples included
informative_path2<-as.character(inlist_informative2$path)
qdf2<-read.csv(informative_path2[1])
colnames(qdf2)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
for (i in 2:length(informative_path2)) {
  qdf_new<-read.csv(informative_path2[i])
  colnames(qdf_new)<-c("X","read_count","rnames","is_common_variant","is_no_variant","allele_freq","specimen","index")
  qdf2<-rbind(qdf2,qdf_new)
}
qdf3<-qdf2[,-1]

#generate a vector of the number of variants over 2%
vaf2_all<-numeric()
for (i in 1:length(all_samples)){
  qdf_filt<-filter(qdf3, qdf3[,6]==all_samples[i])#loads in dfs one sample at a time
  qdf_filt_inf<-qdf_filt[-c(1,2),]# removes no variant and common variant
  qdf_filt_inf$index<-qdf_filt_inf$index-2#resets index count
  qdf_filt_inf$adj_allele_freq<-qdf_filt_inf$read_count/sum(qdf_filt_inf$read_count)#recalculates allele frequencies from read counts without no variant and common variant
  qdf_2percent<-filter(qdf_filt_inf, qdf_filt_inf$adj_allele_freq>=0.02)
  vaf2_all<-c(vaf2_all,nrow(qdf_2percent))
}
vaf2_all
fidf_5<-fidf_3.1
fidf_5$vars_over_0.02<-vaf2_all

colorder<-c("Low FI","High FI")
p5<-ggplot(fidf_5, aes(x = fidf_5$boot_clust, y = fidf_5$vars_over_0.02, fill = fidf_5$boot_clust))
p5<-p5+
  geom_boxplot(notch = FALSE)+
  #geom_violin()+
  theme_cowplot(font_size = unname(gp_informative_boxplot["fontsize"]))+
  scale_x_discrete(limits = colorder, labels = c(expression(paste("Low ",Phi),paste("High ",Phi))))+
  scale_fill_manual(values = c("#72ce55","#fde725"))+
  theme(legend.position = "none")+
  labs(x = "", y = expression(paste(B["0.02"],sep = "")))
p5
save_plot(plot = p5, filename = paste0(output_folder,"/informative_boxplot.pdf"), base_width = unname(gp_informative_boxplot["width"]), base_height = unname(gp_informative_boxplot["height"]))

fidf_5_highfi<-filter(fidf_5, fidf_5$boot_clust=="High FI")
fidf_5_lowfi<-filter(fidf_5, fidf_5$boot_clust=="Low FI")

sink(paste0(output_folder,"/fi_report.txt"))
# fidf_5_highfi<-filter(fidf_5, fidf_5$boot_clust=="High FI")
fidf_5_highfi
# fidf_5_lowfi<-filter(fidf_5, fidf_5$boot_clust=="Low FI")
fidf_5_lowfi
mean(fidf_5_highfi$FI)
se(fidf_5_highfi$FI)
sd(fidf_5_highfi$FI)
mean(fidf_5_lowfi$FI)
se(fidf_5_lowfi$FI)
sd(fidf_5_lowfi$FI)
ks.test(fidf_5_highfi$FI,fidf_5_lowfi$FI)
sink()

fileConn<-(paste0(output_folder,"/data_out.txt"))
l1<-paste("The list of informative samples after thresholding is:",collapse = "")
l2<-paste(vaf2_df$sample_name, collapse = ",")
l3<-paste("\nThe number of clones with Vaf>0.02 in each sample is:",collapse = "")
l4<-paste(vaf2_df$value, collapse = ",")
l5<-paste("\n",shapiro.test(vaf2_df$value)[3][[1]],"for VAF>0.02:  W = ",shapiro.test(vaf2_df$value)[1][[1]],", p = ",shapiro.test(vaf2_df$value)[2][[1]], collapse = "")
l6<-paste("\nThe Shannon entropy for each sample is:",collapse = "")
l7<-paste(round(shannon_df$value,4), collapse = ",")
l8<-paste("\n",shapiro.test(shannon_df$value)[3][[1]],"for Shannon entropy:  W = ",shapiro.test(shannon_df$value)[1][[1]],", p = ",shapiro.test(shannon_df$value)[2][[1]], collapse = "")
l9<-paste("\nThe Simpson diversity index for each sample is:",collapse = "")
l10<-paste(round(simpson_df$value,4), collapse = ",")
l11<-paste("\n",shapiro.test(simpson_df$value)[3][[1]],"for Simpson diversity index:  W = ",shapiro.test(simpson_df$value)[1][[1]],", p = ",shapiro.test(simpson_df$value)[2][[1]], collapse = "")
writeLines(c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11), fileConn)

#output data in csv
write_csv(vaf2_df, path = paste0(output_folder,"/data_out_vaf2.csv"))
write_csv(shannon_df, path = paste0(output_folder,"/data_out_shannon.csv"))
write_csv(simpson_df, path = paste0(output_folder,"/data_out_simpson.csv"))

#####generate sharing factor plot for informative only####
cor_func4<-function(fidf_5_highfi, thresh, tvaf, method, output_folder){
  input_directory<-paste0(output_folder,"/thresh_",thresh,"_tvaf_",tvaf)
  #inlist<-list.files(input_directory,full.names = TRUE, pattern = "\\d.csv$")
  inlist<-list.files(input_directory,full.names = TRUE)
  insample<-str_sub(inlist,-9,-5)
  inframe<-as.data.frame(cbind(inlist,insample))
  inlist1<-left_join(fidf_5_highfi,inframe,by = c("specimen" = "insample"))
  inlist1$inlist<-as.character(inlist1$inlist)
  df_cor<-data.frame(Variant = character(), Count = integer(), Specimen = character())
  for (i in 1:nrow(inlist1)){
    df_new<-read.csv(file = inlist1$inlist[i], header = TRUE)
    df_new<-df_new[,c(3,2,7)]
    colnames(df_new)<-c("Variant", "Count", "Specimen")
    df_cor<-rbind(df_cor,df_new)
  }
  df_cor_wide0<-reshape(df_cor, idvar = "Variant", timevar = "Specimen", direction = "wide")
  df_cor_wide<-df_cor_wide0
  df_cor_wide[is.na(df_cor_wide)]<-0
  colnames(df_cor_wide)<-substr(colnames(df_cor_wide),7,11)
  rownames(df_cor_wide)<-df_cor_wide[,1]
  df_cor_wide<-df_cor_wide[,-1]
  df_cor_wide1<-df_cor_wide[-c(1,2),]

  stat_matrix<-matrix(vector(),nrow = ncol(df_cor_wide1), ncol = ncol(df_cor_wide1))
  colnames(stat_matrix)<-colnames(df_cor_wide1)
  rownames(stat_matrix)<-colnames(df_cor_wide1)

  for (i in 1:ncol(df_cor_wide1)) {
    for (j in 1:ncol(df_cor_wide1)) {
      a<-pmin(df_cor_wide1[,i],df_cor_wide1[,j])
      b<-pmax(df_cor_wide1[,i],df_cor_wide1[,j])
      c<-(sum(a)*2)/(sum(a)+sum(b))
      stat_matrix[i,j]<-c

    }

  }

  reorder_cormat <- function(x){
    # Use correlation between variables as distance
    dd <- as.dist((1-x)/2)
    hc <- hclust(dd)
    x <-x[hc$order, hc$order]
  }

  get_lower_tri <- function(x){
    x[upper.tri(x)]<- NA
    return(x)
  }

  sm_reordered<-reorder_cormat(stat_matrix)
  sm_lower<-get_lower_tri(sm_reordered)

  msm<-melt(sm_lower)

  msm_plot<-ggplot(data = msm, aes(x = msm$Var1, y = msm$Var2, fill = msm$value))
  msm_plot<-msm_plot+
    geom_tile()+
    #scale_fill_viridis_c(begin = 0, end = 1, guide = "colourbar", aesthetics = "fill", na.value = "white")+
    scale_fill_distiller(palette = "RdYlBu", na.value = "transparent",guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE),limits = c(0,1))+
    #scale_x_discrete(breaks = c(5,10,15,20,25,30,35))+
    theme_cowplot(font_size = 11)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))+
    theme(axis.text.y = element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(axis.line.y = element_blank())+
    theme(axis.ticks.y = element_blank())+
    xlab("Sample")+
    theme(axis.title.y = element_blank())+
    theme(plot.title = element_blank())+
    theme(plot.subtitle = element_blank())+
    guides(fill=guide_colorbar(title="Sharing\nFactor", reverse = F))+
    theme(legend.position = c(0,0.7))+
    coord_fixed()
  msm_plot

  #calculate Fraction informative
  FI<-as.data.frame(colSums(df_cor_wide[-c(1,2),])/colSums(df_cor_wide))#ratio of reads excluding no variant and common variant to all reads in each sample.
  colnames(FI)<-"Fraction_Informative"
  FI$sample_x<-rownames(FI)

  avg_share<-round(mean(msm[which(msm$Var1!=msm$Var2),3],na.rm = TRUE),4)#averages pearson correlation excluding identical comparisons

  mean_inf<-round(mean(FI$Fraction_Informative),3)

  save_plot(plot = msm_plot, filename = paste0(output_folder,"/share_factor_highfi.pdf"), base_width = 2.5)
  nubbin<-data.frame(threshold = thresh, tvaf = tvaf, mean_inf = mean_inf, avg_share = avg_share)
  write.csv(x = nubbin, file = paste0(input_directory,"/nubbin_informative_only_",method,".csv"))
  return(nubbin)
}

informative_only_nubbin<-cor_func4(fidf_5_highfi = fidf_5_highfi, thresh = 0, tvaf = final_tvaf, method = "sf", output_folder = output_folder)


select_hifi<-function(data_in,keeplist) {
  data<-data_in[which(data_in[,1] %in% keeplist & data_in[,2] %in% keeplist),]
  data$rerank<-rank(data[,3], ties.method = "random")
  return(data)
}

msm_1_true_hifi<-select_hifi(data_in = msm_1_true, keeplist = fidf_5_highfi$specimen)
msm_1_false_hifi<-select_hifi(data_in = msm_1_false, keeplist = fidf_5_highfi$specimen)
msm_0.3_false_hifi<-select_hifi(data_in = msm_0.3_false, keeplist = fidf_5_highfi$specimen)
msm_0.1_false_hifi<-select_hifi(data_in = msm_0.1_false, keeplist = fidf_5_highfi$specimen)
msm_0.03_false_hifi<-select_hifi(data_in = msm_0.03_false, keeplist = fidf_5_highfi$specimen)
msm_0.01_false_hifi<-select_hifi(data_in = msm_0.01_false, keeplist = fidf_5_highfi$specimen)
msm_0.003_false_hifi<-select_hifi(data_in = msm_0.003_false, keeplist = fidf_5_highfi$specimen)
msm_0.001_false_hifi<-select_hifi(data_in = msm_0.001_false, keeplist = fidf_5_highfi$specimen)
msm_0.0003_false_hifi<-select_hifi(data_in = msm_0.0003_false, keeplist = fidf_5_highfi$specimen)




msm_hifi_list<-list(msm_1_true_hifi,msm_1_false_hifi,msm_0.3_false_hifi,msm_0.1_false_hifi,msm_0.03_false_hifi,msm_0.01_false_hifi,msm_0.003_false_hifi,msm_0.001_false_hifi,msm_0.0003_false_hifi)
msm_hifi<-dplyr::bind_rows(msm_hifi_list)
msm_hifi$condition<-paste0(msm_hifi$tvaf,"_",msm_hifi$include_uninformative)
msm_hifi$condition<-factor(msm_hifi$condition, levels = unique(msm_hifi$condition))

msm_hifi<-msm_hifi[which(msm_hifi$tvaf==final_tvaf | msm_hifi$tvaf==1),]

paired_share_factors2<-ggplot(msm_hifi, aes(x = rerank, y = value, fill = condition, group = condition))+
  geom_point(shape=21, color = "grey20", size = 1, alpha = 0.6, stroke = 0.25)+
  scale_x_reverse()+
  theme_cowplot(font_size = unname(gp_sharing_curve_hifi["fontsize"]))+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())+
  scale_y_log10()+
  scale_fill_manual(values = colors, labels = labels)+
  theme(legend.position = c(0.1,0.2))+
  theme(legend.text = element_text(size = 8))+
  guides(fill = guide_legend(keyheight = 0.1, default.unit = "inch",ncol = 1, override.aes = list(shape = 21,size = 2, alpha = 1)))+
  labs(x = "Ranked Sample Pair", y = "Sharing Factor", fill = "")
paired_share_factors2
save_plot(paired_share_factors2, filename = paste0(output_folder,"/sharing_factor_curves_hifi.pdf"), base_width = unname(gp_sharing_curve_hifi["width"]), base_height = unname(gp_sharing_curve_hifi["height"]))

#print tvaf qc report

fileConn<-file(paste0(output_folder,"/tvaf_qc_report.txt"))
l1<-"Theta (tvaf) QC Report\n\n"
if (tvaf_preset == F){
  l2<-(paste0("The value automatically selected for theta is ",final_tvaf,"\n\n"))
} else {
  l2<-(paste0("The value you have selected for theta is ",final_tvaf,"\n\n"))
}
if(final_tvaf>0.03){lwarning1<-"**Warning - the selected theta (tvaf) is over 3 percent.  Consider manually choosing a lower value.**\n\n"} else {lwarning1<-""}
l3<-(paste0("This means that variants present at greater than ",final_tvaf*100," percent in more than 1 sample are considered shared variants.\n\n"))
l4<-(paste0("At this value of theta, the mean sharing factor in all submitted samples is ",curve_df6[order(-curve_df6$product),5][1],"\n\n"))
l5<-(paste0("This means that prior to removal of uninformative samples, the average degree of barcode sharing is ",curve_df6[order(-curve_df6$product),5][1]*100," percent.\n\n"))
l6<-(paste0("And the average percent of informative reads per sample is ",curve_df6[order(-curve_df6$product),4][1]*100," percent.\n\n"))
l7<-(paste0("After removal of uninformative samples, the average degree of barcode sharing is ",informative_only_nubbin[1,4]*100," percent.\n\n"))
if(informative_only_nubbin[1,4]>0.01){lwarning2<-"**Warning - the average sharing factor is over 1%.  This dataset may be unreliable.**\n\n"} else {lwarning2<-""}
l8<-(paste0("And the average percent of informative reads is ",informative_only_nubbin[1,3]*100," percent.\n\n"))
if(informative_only_nubbin[1,3]<0.6){lwarning3<-"**Warning - the percent of informative reads in the final sample set is lower than 50%.  This dataset may be unreliable.**\n\n"} else {lwarning3<-""}
writeLines(c(l1,l2,lwarning1,l3,l4,l5,l6,l7,lwarning2,l8,lwarning3), fileConn)
close(fileConn)

#calculate total number of variants in informative samples, excluding no variant and common variant:
final_var_counts<-tibble(qdf3)
summary_table<-qdf3 %>% filter(specimen %in% fidf_5_highfi$specimen) %>% filter(rnames!="no variant" & rnames!="common variant sum") %>% group_by(specimen) %>% summarise(n = n())
var_summary<-summary_table %>% mutate(mean = mean(n), sd = sd(n), min = min(n), max = max(n)) %>% select(-specimen,-n) %>% distinct()

#make the crisprvariants plot using colors to identify High FI and Low FI samples
grp<-fidf_5 %>% arrange(specimen) %>% pull(boot_clust)# %>% recode("High FI" = 1, "Low FI" = 2)
grp<-factor(grp, levels = c("High FI", "Low FI"))
#grp_colors<-c()
# plot the variants
ps<-37
pam_seq<-seq(ps,280,27)

#while (!is.null(dev.list())) dev.off()
pdf(file=paste0(output_folder,"/",output_file,".pdf"), width = unname(gp_crisprplot["width"]), height = unname(gp_crisprplot["height"]))
plotVariants(crispr_set,
             col.wdth.ratio = c(1,1),
             plotAlignments.args = list(pam.start = pam_seq, #c(37,64), #draws a 3-nt box starting including the position noted
                                        target.loc = pam_seq-3, #draws a vertical line after the position noted
                                        guide.loc = IRanges::IRanges(pam_seq-20,pam_seq+2), #first parameter - beginning of target sequence, second - end of target sequence
                                        min.count = threshold,
                                        tile.height = 0.9,
                                        xtick.labs = c("0","100","200","300"),
                                        xtick.breaks = c(0,100,200,300),
                                        plot.text.size = 0,
                                        top.n = unname(gp_crisprplot["top.n"])),
             plotFreqHeatmap.args = list(min.count = threshold,
                                         plot.text.size = 2,
                                         x.size = 8,
                                         group = grp,
                                         legend.text.size = 8,
                                         top.n = unname(gp_crisprplot["top.n"]),
                                         header = "counts",
                                         type = "proportions",
                                         #legend.position = "none",
                                         legend.key.height = grid::unit(0.5, "lines")))
#dev.copy2pdf(file=paste0(output_folder,"/",output_file,".pdf"), width = unname(gp_crisprplot["width"]), height = unname(gp_crisprplot["height"]))  #for 10 samples use 24 x 24, for 30-40 samples use 48 x 48
dev.off()

# clean up directories
# unlink("temp", recursive = TRUE)
# unlink("fastq", recursive = TRUE)
# dir.create("fastq")
# system(paste0("cp analysis_resub_vfinal.R ",output_folder,"/analysis_resub_vfinal.R"))
save.image(paste0(output_folder,"/rdata.RData"))

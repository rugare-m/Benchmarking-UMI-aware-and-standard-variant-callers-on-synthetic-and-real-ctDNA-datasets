library(ggrepel)

vc_colors <- c("bcftools" = "#FF5733",
               "FreeBayes" = "#00BFFF",
               "LoFreq" = "#FFC300",
               "Mutect2" = "#00FF7F",
               "UMI-VarCal" = "#8A2BE2",
               "UMIErrorCorrect" = "#FF1493")

#Load datasets 
replicates <- read.csv("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/replicate_comps.csv")
comet <- read.csv("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/comet2.csv")
real_depths <- read.csv("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/real_depths.csv")
vafs <- read.csv("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/vafs_complete.csv")
af_vs_depth <- read.csv("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/depth_v_af.csv") 




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+REAL DATA ALLELE FREQUENCIES
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
allele_depth_order <-  af_vs_depth
allele_depth_order$vc <- factor(allele_depth_order$vc, levels = c("bcftools","freebayes", "lofreq", "mutect2","mageri", "umivarcal", "umierrorcorrect"))

allele_depth.mod <- allele_depth_order%>% 
  mutate(vc = recode(vc, "freebayes" = "FreeBayes", "lofreq" = "LoFreq","bcftools" = "bcftools", "mutect2" = "Mutect2", "mageri"= "MAGERI", "umivarcal" = "UMI-VarCal", "umierrorcorrect" = "UMIErrorCorrect"))


point_plots <- allele_depth.mod %>% filter(status == "umi") %>%
  ggplot(aes(af, depth))+
  geom_point(size = 0.1)+
  facet_grid(cols = vars(status), rows = vars(vc))+
  xlim (0,1)
point_plots


umi_dist_plot <- allele_depth.mod %>% filter(status == "umi") %>%
  filter(af <= 0.9) %>%
  filter(vc != "MAGERI") %>%
  ggplot(aes(x = af))+
  facet_grid(rows = vars(vc))+
  geom_histogram(bins = 75)+
  ylim(0,75)+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))+
  xlab("Allele Frequency")+
  theme_classic()+
  theme(strip.text.x = element_text(size = 14, face="bold"), 
      strip.text.y = element_text(size = 14, face="bold"),
      axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0.5, face="bold", size = 20),
      axis.text.y = element_text(angle = 0, hjust=0.5, vjust=0.5, face="bold", size = 15),
      axis.text=element_text(size=20, face="bold"), 
      axis.title=element_text(size=20, face="bold"),
      axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")))
umi_dist_plot

ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/publications/Chapter One/figures/real/umi/distribution_af.png",
  plot = umi_dist_plot,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#SYNTHETIC DATA 
synth_umi_only <- replicates %>% filter(status=='umi'  & vc != "octopus") %>% mutate(appx_af = as.character(appx_af)) %>%  mutate(actual_mean_af = as.character(actual_mean_af))
#code to reorder the facet grid - umi_only dataset
synth_umi_only_order <-  synth_umi_only
synth_umi_only_order$vc <- factor(synth_umi_only_order$vc, levels = c("bcftools","freebayes", "lofreq", "mutect2","mageri", "umivarcal", "umierrorcorrect"))
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
synth_umi_only_order.mod <- synth_umi_only_order%>% 
  mutate(vc = recode(vc, "freebayes" = "FreeBayes", 
                     "lofreq" = "LoFreq",
                     "bcftools" = "bcftools", 
                     "mutect2" = "Mutect2", 
                     "mageri"= "MAGERI", 
                     "umivarcal" = "UMI-VarCal", 
                     "umierrorcorrect" = "UMIErrorCorrect")) %>% mutate(appx_depth = recode(appx_depth, "200" = "200", 
                                                                                            "440" = "450", 
                                                                                            "830" = "850")) %>%filter (vc != "MAGERI")


#http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data: a dataframe
# varname: the name of a column containing the variable to be summarised
# groupnames: vector of column names to be used as grouping variables

library(plyr)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

########### Total Calls  ####################### 
##++++++++++++++++++++++++++++++++++++++
##+SYNTHETIC DATASET W/UMIs
##++++++++++++++++++++++++++++++++++++++

# make the bar plot w/umis
df.barplot <- data_summary(synth_umi_only_order.mod, varname="total", 
                    groupnames=c("appx_af", "appx_depth", "vc"))

# Convert appx_af to a factor variable
df.barplot$appx_af=as.factor(df.barplot$appx_af)
head(df.barplot)

df.barplot.mod <- df.barplot %>% filter (vc != "octopus") %>% filter (vc != "mageri") %>% mutate(appx_depth = recode(appx_depth, "200" = "200", "440" = "450", "830" = "850")) %>% 
  mutate(vc = recode(vc, "freebayes" = "FreeBayes", "lofreq" = "LoFreq","bcftools" = "bcftools", "mutect2" = "Mutect2", "umivarcal" = "UMI-VarCal", "umierrorcorrect" = "UMIErrorCorrect")) 

barplot.synth.umis <- ggplot(df.barplot.mod, aes(x=appx_af, y=total)) + 
  geom_bar(stat="identity", position=position_dodge(), fill="white", colour = "black", linewidth=0.7) +
  geom_errorbar(aes(ymin=total-sd, ymax=total+sd), width=.2, position=position_dodge(.9),  colour = "black", linewidth=1.2)+
  facet_grid(rows = vars(appx_depth), cols = vars(vc), labeller=labeller(depth = labels))+
  theme_classic()+
  theme(strip.text.x = element_text(size = 12, face="bold"), 
        strip.text.y = element_text(size = 12, face="bold"),
        axis.text.x = element_text(size = 15,angle = 45, hjust=1, face="bold"),
        axis.text.y = element_text(size = 15,face="bold"),
        axis.text=element_text(size=11, face="bold"), 
        axis.title=element_text(size=14, face="bold"))+
  xlab("Binned Allele Frequency")+
  ylab("Mean Total Calls")
barplot.synth.umis

ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 2/Figures/total_synthetic_bar.png",
  plot = barplot.synth.umis,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")


#Total Number of variants detected - box
total_box <- ggplot (synth_umi_only_order.mod, aes(x=appx_af, y=total))+
  geom_boxplot(aes())+
  facet_grid(rows = vars(appx_depth), cols = vars(vc))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Mean Total Variants")+
  xlab("Approximate Allele Frequency")+
  theme_linedraw()+
  theme(strip.text.x = element_text(size = 14, face="bold"), 
        strip.text.y = element_text(size = 14, face="bold"),
        axis.text.x = element_text(angle = 45, hjust=0.5, vjust=0.5, face="bold", size = 15),
        axis.text.y = element_text(angle = 0, hjust=0.5, vjust=0.5, face="bold", size = 15),
        axis.text=element_text(size=20, face="bold"), 
        axis.title=element_text(size=20, face="bold"),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
total_box

ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 2/Figures/total_synthetic_box.png",
  plot = total_box,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")


###########################True Positive calls detected###############################

##++++++++++++++++++++++++++++++++++++++
##+Synthetic Dataset with UMIs
##++++++++++++++++++++++++++++++++++++++
df.umi.tp <- synth_umi_only_order %>% filter (vc != "octopus") %>%  filter (vc != "mageri") %>% mutate(appx_depth = recode(appx_depth, "200" = "200x", "440" = "450x", "830" = "850x")) %>% 
  mutate(vc = recode(vc,"bcftools" = "bcftools", "freebayes" = "FreeBayes","lofreq" = "LoFreq", "mutect2" = "Mutect2", "umivarcal" = "UMI-VarCal", "umierrorcorrect" = "UMIErrorCorrect")) %>%
  mutate(actual_mean_af = recode (actual_mean_af, "0.0059" = "0.006", "0.016" = "0.015", "0.033" = "0.03", "0.056" = "0.055", "0.058" = "0.06", "0.078" = "0.08", "0.079" = "0.08", "0.083" = "0.08"))





boxplot.umi.tp <- ggplot(df.umi.tp, aes(x=actual_mean_af, y=tp, fill = vc))+
  stat_summary(fun = mean, geom = "line", linewidth = 1,aes(group = vc,colour = vc))+
  #geom_boxplot(lwd=0.7)+
  #stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(0.9)) +
  facet_grid(cols = vars(appx_depth))+ 
  theme_linedraw()+
  theme(strip.text.x = element_text(size = 12,face="bold"),
        strip.text.y = element_text(size = 40),
        axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14, face="bold"),
        axis.text.x = element_text(size = 14,angle = 45, hjust=1, vjust = 1, face="bold"),
        axis.text.y = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = rel(1.2), angle = 90),
        axis.title.x = element_text(size = rel(1.5), angle = 0))+
  theme(legend.position="top")+
  ylab ("True Positive Calls")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))+
  xlab ("Mean Allele Frequency")+
  labs(colour = "Variant Caller")+
  guides(fill=guide_legend(title="Variant Caller", shape=FALSE),
         guide_colorbar(barwidth = 100))+
  scale_color_manual(values = vc_colors) +
  theme(legend.title = element_text(face = "bold", size = 14))+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(face = "bold", size = 12))
boxplot.umi.tp

ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 2/Figures/tp_synthetic.png",
  plot = boxplot.umi.tp,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")


###########################True Positive calls detected###############################

########## False Positive ##############################

##++++++++++++++++++++++++++++++++++++++
##+Synthetic Dataset with  UMIs
##++++++++++++++++++++++++++++++++++++++

df.umi.fp <- synth_umi_only_order %>% filter (vc != "octopus") %>% filter (vc != "mageri") %>% mutate(appx_depth = recode(appx_depth, "200" = "200x", "440" = "450x", "830"="850x")) %>% 
  mutate(vc = recode(vc,"bcftools" = "bcftools", "freebayes" = "FreeBayes","lofreq" = "LoFreq", "mutect2" = "Mutect2", "umivarcal" = "UMI-VarCal", "umierrorcorrect" = "UMIErrorCorrect")) %>%
  mutate(actual_mean_af = recode (actual_mean_af, "0.0059" = "0.006","0.0057" = "0.006","0.009" = "0.01","0.018" = "0.02","0.019" = "0.02","0.016" = "0.015", "0.033" = "0.03", "0.056" = "0.055", "0.058" = "0.06", "0.078" = "0.08", "0.079" = "0.08",
                                                                                                                                                                                                                           "0.083" = "0.08", "0.036" = "0.035","0.039" = "0.04", "0.051" = "0.05", "0.059" = "0.06", "0.076" = "0.075"))


boxplot.umi.fp <- ggplot(df.umi.fp, aes(x=actual_mean_af , y=total_minus_tp, fill = vc))+
  stat_summary(fun = mean, geom = "line", linewidth = 1, aes(group = vc, colour = vc))+
  #geom_boxplot(lwd=0.7)+
  facet_grid(cols = vars(appx_depth))+ 
  theme_linedraw()+
  theme(strip.text.x = element_text(size = 12,face="bold"),
        strip.text.y = element_text(size = 40),
        axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14, face="bold"),
        axis.text.x = element_text(size = 14,angle = 45, hjust=1, vjust = 1, face="bold"),
        axis.text.y = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = rel(1.5), angle = 90),
        axis.title.x = element_text(size = rel(1.5), angle = 0))+
  ylab ("Putative False Positive")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))+
  xlab ("Mean Allele Frequency")+
  labs(colour = "Variant Caller")+
  guides(fill=guide_legend(title="Depth"))+
  scale_color_manual(values = vc_colors) +
  theme(legend.title = element_text(face = "bold", size = 14))+
  theme(legend.text = element_text(face = "bold", size = 12))+
  theme(legend.position = "top")
boxplot.umi.fp

ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 2/Figures/fp_synthetic.png",
  plot = boxplot.umi.fp,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")



########## False Positive ##############################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Real Datasets  - #subset data
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
real_umi_only <- comet %>% filter(status=='umi')
#code to reorder the facet grid - umi_only dataset
real_umi_only_order <-  real_umi_only
real_umi_only_order$vc <- factor(real_umi_only_order$vc, levels = c("bcftools","freebayes", "lofreq", "mutect2","mageri", "umi-varcal", "umierrorcorrect"))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Real Dataset - umis
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sra.umi.mod <- real_umi_only_order %>% filter(status == "umi")%>% filter(vc != "mageri")%>%
  filter(accession != "SRR15081471", accession != "SRR15081464")%>%
  mutate(vc = recode(vc, "umi-varcal" = "UMI-VarCal", "lofreq" = "LoFreq", 
                     "mutect2" = "Mutect2", "freebayes"= "FreeBayes",
                     "umierrorcorrect"="UMIErrorCorrect")) 

#Total number of variants detected by callers across all 10 patients 
barplot.sra.umi <- ggplot(data = sra.umi.mod, aes(x = accession, y=total))+
  geom_bar(stat="identity", fill="white", colour = "black",linewidth=0.8)+
  facet_grid(cols = vars(vc))+
  theme_linedraw()+
  theme(strip.text.x = element_text(size = 12,face="bold"),
        strip.text.y = element_text(size = 40),
        axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14, face="bold"),
        axis.text.x = element_text(size = 10,angle = 45, hjust=1, vjust = 1, face="bold"),
        axis.text.y = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = rel(1.5), angle = 90),
        axis.title.x = element_text(size = rel(1.5), angle = 0))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))+
  xlab (NULL)+
  ylab ("Total Calls")
barplot.sra.umi
ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 2/Figures/total_calls_real.png",
  plot = barplot.sra.umi,
  dpi = 600,
  width=10, height=10,
  limitsize = FALSE,
  bg = "transparent")

#Percentage Pathogenic variants  - UMIs
pointplot.sra.umi <-sra.umi.mod  %>%  filter(vc != "mageri", accession != "SRR15081471", accession != "SRR15081464") %>%
  ggplot(aes(x = vc, y=percentage_pathogenic))+
  geom_point(aes(colour = accession), size = 5)+
  geom_text_repel(box.padding = 0.5,aes(label = total, colour = accession, fontface = "bold"), 
                  vjust = 0.5, hjust = -1, size = 4, force = 4) +
  theme_linedraw()+
  theme(strip.text.x = element_text(size = 12, face="bold"),
        strip.text.y = element_text(size = 16, face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=12, face="bold"),
        axis.text=element_text(size=14, face="bold"),
        axis.title=element_text(size=16, face="bold"),
        axis.text.x = element_text(angle = 0, size = 16),
        axis.text.y = element_text(angle = 0, size = 14),
        legend.title = element_text(size=15, face="bold"),
        legend.position="top")+
  ylim(0,110)+
  labs(fill = "Dose (mg)")+
  xlab (NULL)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))+
  ylab ("Percentage of COSMIC Variants Detected")+
  labs(color=NULL)
pointplot.sra.umi
ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 2/Figures/percentage_cosmic_real.png",
  plot = pointplot.sra.umi,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")


#Concordance of variants <- UpSetR <- umi only variants  - TBD
all.calls<- list(
  UMIVarCal = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/umis/umivarcal.txt"),
  UMIErrorCorrect = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/umis/umierrorcorrect.txt"),
  bcftools = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/umis/bcftools.txt"),
  FreeBayes = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/umis/fbayes.txt"),
  LoFreq = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/umis/lofreq.txt"),
  Mutect2 = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/umis/mut2.txt")
  )

png("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 2/Figures/concordance_real.png", width = 465, height = 225, units='mm', res = 500)
upset(fromList(all.calls), nsets = 7, number.angles = 0, point.size = 6, line.size = 3, 
      mainbar.y.label = "Intersects", sets.x.label = "Total number of variants called", 
      text.scale = c(2, 2, 2, 2, 2.5, 2),order.by = "degree", nintersects = 60, 
      set_size.show = TRUE, set_size.numbers_size = 8, set_size.scale_max = 7500,
      show.numbers = "yes")
dev.off()



#compare true positive, low frequency variants called between UMI and non-UMI datasets, across variant callers
comps <- replicates %>% filter (vc != "octopus") %>% filter (vc != "mageri") %>% filter (vc != "umivarcal") %>% filter (vc != "umierrorcorrect")

x <- ggplot(comps, aes(x=vc, y=tp, fill=status)) + 
  geom_boxplot()+
  facet_wrap(~appx_af, scale="free")



my_comparisons <- list( c("umi", "no_umi"))

x+ stat_compare_means(aes(label = ..p.signif..), 
                      label.x = 1.5, label.y = 40)









#False positives?
ggplot(comps, aes(x=vc, y=total_minus_tp, fill=status)) + 
  geom_boxplot()+
  facet_wrap(~appx_af, scale="free")

library(ggpubr)

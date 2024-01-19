#Datasets1 > Dietz et al Real dataset 69x
#Datasets2 > Synthetic Dataset SRR10296599 with spiked in variants, 150x 350x 450x
#Dataset3 > Foll et al

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
dietz <- read_excel("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/dietz.xlsx")
real_depths <- read.csv("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/real_depths.csv")
vafs <- read.csv("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/vafs_complete.csv")
af_vs_depth <- read.csv("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/depth_v_af.csv") 
daf <- read_excel("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/dietz_af.xlsx")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+REAL COMET DATA ALLELE FREQUENCIES
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
allele_depth_order <-  af_vs_depth
allele_depth_order$vc <- factor(allele_depth_order$vc, levels = c("bcftools","freebayes", "lofreq", "mutect2","mageri", "umivarcal", "umierrorcorrect"))

allele_depth.mod <- allele_depth_order%>% 
  mutate(vc = recode(vc, "freebayes" = "FreeBayes", "lofreq" = "LoFreq","bcftools" = "bcftools", "mutect2" = "Mutect2", "mageri"= "MAGERI", "umivarcal" = "UMI-VarCal", "umierrorcorrect" = "UMIErrorCorrect"))

no_umi_dist <- allele_depth.mod %>% filter(status == "no_umi") %>%
  filter(af < 1.0) %>%
  ggplot(aes(x = af))+
  facet_grid(rows = vars(vc))+
  geom_histogram(bins = 75, color="black", fill = "grey")+
  ylim (0,295)+
  xlab("\nAllele Frequency")+
  theme_linedraw()+
  theme(strip.text.x = element_text(size = 14, face="bold"), 
        strip.text.y = element_text(size = 14, face="bold"),
        axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0.5, face="bold", size = 15),
        axis.text.y = element_text(angle = 0, hjust=0.5, vjust=0.5, face="bold", size = 15),
        axis.text=element_text(size=20, face="bold"), 
        axis.title=element_text(size=20, face="bold"),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  xlim(0, 1)+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
no_umi_dist

ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/real_mbc_distr.png",
  plot = no_umi_dist,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+REAL WES DATA ALLELE FREQUENCIES
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


allele_depth <- daf%>% filter(AF < 1.0) %>%
  mutate(vc = recode(vc, "freebayes" = "FreeBayes", "lofreq" = "LoFreq","bcftools" = "bcftools", "mutect2" = "Mutect2", "mageri"= "MAGERI", "umivarcal" = "UMI-VarCal", "umierrorcorrect" = "UMIErrorCorrect"))

no_umi_dist <- allele_depth %>%
  ggplot(aes(x = AF))+
  facet_grid(rows = vars(vc))+
  geom_histogram(bins = 75, color="black", fill = "grey")+
  ylim (0,100)+
  xlab("\nAllele Frequency")+
  theme_linedraw()+
  theme(strip.text.x = element_text(size = 14, face="bold"), 
        strip.text.y = element_text(size = 14, face="bold"),
        axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0.5, face="bold", size = 15),
        axis.text.y = element_text(angle = 0, hjust=0.5, vjust=0.5, face="bold", size = 15),
        axis.text=element_text(size=20, face="bold"), 
        axis.title=element_text(size=20, face="bold"),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1))
no_umi_dist

ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/real_wes_distr.png",
  plot = no_umi_dist,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#SYNTHETIC DATA TOTAL VARIANTS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
synth_no_umi_only <- replicates %>% filter(status=='no_umi' & vc != "octopus") %>% mutate(appx_af = as.character(appx_af)) %>%  mutate(actual_mean_af = as.character(actual_mean_af))

#code to reorder the facet grid - no umi_only dataset
synth_no_umi_only_order <-  synth_no_umi_only
synth_no_umi_only_order$vc <- factor(synth_no_umi_only_order$vc, levels = c("bcftools","freebayes", "lofreq", "mutect2"))

synth_no_umi_only_order.mod <- synth_no_umi_only_order%>% 
  mutate(vc = recode(vc, "freebayes" = "FreeBayes", 
                         "lofreq" = "LoFreq",
                         "bcftools" = "bcftools", 
                         "mutect2" = "Mutect2", 
                         "mageri"= "MAGERI", 
                         "umivarcal" = "UMI-VarCal", 
                         "umierrorcorrect" = "UMIErrorCorrect")) %>% mutate(appx_depth = recode(appx_depth, "200" = "200", 
                                                                                                            "440" = "450", 
                                                                                                            "830" = "850"))
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Total Number of variants detected 
total_box <- ggplot (synth_no_umi_only_order.mod, aes(x=appx_af, y=total))+
  geom_boxplot(aes())+
  facet_grid(rows = vars(appx_depth), cols = vars(vc))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Mean Total Variants")+
  xlab("Approximate Allele Frequency")+
  theme_linedraw()+
  theme(strip.text.x = element_text(size = 14, face="bold"), 
        strip.text.y = element_text(size = 14, face="bold"),
        axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0.5, face="bold", size = 15),
        axis.text.y = element_text(angle = 0, hjust=0.5, vjust=0.5, face="bold", size = 15),
        axis.text=element_text(size=20, face="bold"), 
        axis.title=element_text(size=20, face="bold"),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
total_box

ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/total_synthetic.png",
  plot = total_box,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")

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


########### Synthetic Total Calls  ####################### 


##++++++++++++++++++++++++++++++++++++++
##+SYNTHETIC DATA TOTAL VARIANTS - BARPLOT
##++++++++++++++++++++++++++++++++++++++

# make the bar plot w/no_umis
df2.barplot <- data_summary(synth_no_umi_only_order.mod, varname="total", 
                    groupnames=c("appx_af", "appx_depth", "vc"))
# Convert appx_af to a factor variable - no/umi
df2.barplot$appx_af=as.factor(df2.barplot$appx_af)
head(df2.barplot)

barplot.synth.no.umis<- ggplot(df2.barplot, aes(x=appx_af, y=total)) + 
  geom_bar(stat="identity", position=position_dodge(), fill="white", colour = "black", linewidth=0.7) +
  geom_errorbar(aes(ymin=total-sd, ymax=total+sd), width=.2, position=position_dodge(.9),  colour = "black", linewidth=1.2)+
  facet_grid(rows = vars(appx_depth), cols = vars(vc), labeller=labeller(depth = labels))+
  theme_linedraw()+
  theme(strip.text.x = element_text(size = 14, face="bold"), 
        strip.text.y = element_text(size = 14, face="bold"),
        axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0.5, face="bold", size = 15),
        axis.text.y = element_text(angle = 0, hjust=0.5, vjust=0.5, face="bold", size = 15),
        axis.text=element_text(size=20, face="bold"), 
        axis.title=element_text(size=20, face="bold"),
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  xlab("Approximate Allele Frequency")+
  ylab("Mean Total Variants")
barplot.synth.no.umis


ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/total_synthetic_bar.png",
  plot = barplot.synth.no.umis,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")


###########################True Positive calls - comet###############################
##++++++++++++++++++++++++++++++++++++++
##+Synthetic Dataset with no UMIs
##++++++++++++++++++++++++++++++++++++++
df.no.umi.tp <- synth_no_umi_only_order %>% filter (vc != "octopus") %>% mutate(appx_depth = recode(appx_depth, "200" = "200x", "440" = "450x", "830"="850x")) %>% 
  mutate(vc = recode(vc,"bcftools" = "bcftools", "freebayes" = "FreeBayes","lofreq" = "LoFreq", "mutect2" = "Mutect2")) %>%
  mutate(actual_mean_af = recode (actual_mean_af, "0.0059" = "0.006","0.0057" = "0.006","0.009" = "0.01","0.018" = "0.02","0.019" = "0.02",
                                  "0.016" = "0.015", "0.033" = "0.03", "0.056" = "0.055", "0.058" = "0.06", "0.078" = "0.08", "0.079" = "0.08",
                                  "0.083" = "0.08", "0.036" = "0.035","0.039" = "0.04", "0.051" = "0.05", "0.059" = "0.06", "0.076" = "0.075"))


boxplot.no.umi.tp <- ggplot(df.no.umi.tp, aes(x=actual_mean_af, y=tp, fill = vc))+
  stat_summary(fun = mean, geom = "line",linewidth = 1, aes(group = vc, colour = vc), position = position_dodge(0.75))+
  geom_boxplot(aes(x=actual_mean_af, y=tp, colour = vc, fill=NULL ), lwd=0.5, position = position_dodge(0.75), show.legend = FALSE)+
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
  ylab ("True Positive Calls")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))+
  xlab ("Mean Allele Frequency")+
  labs(colour = "Variant Caller")+
  labs(colour = NULL)+
  guides(fill = guide_legend(title = "Depth", override.aes = list(fill = NA))) +
  scale_color_manual(values = vc_colors) +
  theme(legend.title = element_text(face = "bold", size = 14))+
  theme(legend.text = element_text(face = "bold", size = 15))+
  theme(legend.position = "top")
boxplot.no.umi.tp


ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/tp_synthetic.png",
  plot = boxplot.no.umi.tp,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")

df.no.umi.tp

########################### True Positive calls - dietz ###############################

cosmics.pointplot.dietz <- dietz %>%
  filter (accession != "SRR3401418") %>%
  ggplot(aes(x = vc, y = percentage_cosmic)) +
  geom_point(aes(colour = accession), size = 8, ) +
  #geom_text_repel(aes(label = percentage_cosmic, colour = accession, fontface = "bold"), vjust = 0.5, hjust = -1, 
  #                size = 4, force = 3) +  # Adding labels
  theme_linedraw() +
  theme(strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 0, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.position = "top") +
  labs(fill = "Dose (mg)") +
  ylab("COSMIC Variants") +
  xlab(NULL)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20), limits = c(2, 44)) +
  #ylab("Percentage of COSMIC Variants Detected") +
  labs(color = NULL)

cosmics.pointplot.dietz

ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/percentage_cosmic_dietz.png",
  plot = cosmics.pointplot.dietz,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")



#### test ####
########## Synthetic False Positive ##############################

df.no.umi.fp <- synth_no_umi_only_order %>% filter (vc != "octopus") %>% mutate(appx_depth = recode(appx_depth, "200" = "200x", "440" = "450x", "830"="850x")) %>% 
  mutate(vc = recode(vc,"bcftools" = "bcftools", "freebayes" = "FreeBayes","lofreq" = "LoFreq", "mutect2" = "Mutect2")) %>%
  mutate(actual_mean_af = recode (actual_mean_af, "0.0059" = "0.006","0.0057" = "0.006","0.009" = "0.01","0.018" = "0.02","0.019" = "0.02",
                                  "0.016" = "0.015", "0.033" = "0.03", "0.056" = "0.055", "0.058" = "0.06", "0.078" = "0.08", "0.079" = "0.08",
                                  "0.083" = "0.08", "0.036" = "0.035","0.039" = "0.04", "0.051" = "0.05", "0.059" = "0.06", "0.076" = "0.075"))

boxplot.no.umi.fp <- ggplot(df.no.umi.fp, aes(x=actual_mean_af , y=total_minus_tp, fill = vc))+
  stat_summary(fun = mean, geom = "line", linewidth = 1,aes(group = vc, colour = vc), position = position_dodge(0.75))+
  geom_boxplot(aes(x=actual_mean_af, y=total_minus_tp, colour = vc, fill=NULL ), lwd=0.45, position = position_dodge(0.75), show.legend = FALSE)+
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
  xlab ("Mean Allele Frequency")+
  labs(colour = "Variant Caller")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))+
  ylab ("Putative False Positive")+
  theme(legend.title = element_text(face = "bold", size = 14))+
  scale_color_manual(values = vc_colors) +
  theme(legend.text = element_text(face = "bold", size = 12))+
  guides(fill=guide_legend(title="Depth"))
boxplot.no.umi.fp

ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/fp_synthetic.png",
  plot = boxplot.no.umi.fp,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############### subset data real dataset (comet)  ###################
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
real_no_umi <- comet %>% filter(status=='no_umi')

#code to reorder the facet grid - no umi_only dataset
real_no_umi_only_order <-  real_no_umi
real_no_umi_only_order$vc <- factor(real_no_umi_only_order$vc, levels = c("bcftools","freebayes", "lofreq", "mutect2"))



##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############### Total Calls - comet ##########################
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Give callers proper names
sra.no.umi.mod <- real_no_umi_only_order %>%  filter(vc != "octopus",  accession != "SRR15081471", accession != "SRR15081464") %>%  filter(vc != "platypus")%>% filter(status == "no_umi")%>% 
  mutate(vc = recode(vc, "octopus" = "Octopus", "lofreq" = "LoFreq", "mutect2" = "Mutect2", "freebayes"= "FreeBayes")) 

#Total number of variants detected by callers across all 10 patients 
barplot.sra.no.umi <- ggplot(data = sra.no.umi.mod, aes(x = accession, y=total))+
  geom_bar(stat="identity", fill="white", colour = "black",linewidth=0.8)+
  facet_grid(cols = vars(vc))+
  theme_linedraw()+
  theme(strip.text.x = element_text(size = 12,face="bold"),
        strip.text.y = element_text(size = 40),
        axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14, face="bold"),
        axis.text.x = element_text(size = 14,angle = 45, hjust=1, vjust = 1, face="bold"),
        axis.text.y = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = rel(1.5), angle = 90),
        axis.title.x = element_text(size = rel(1.5), angle = 0))+
  ylab ("Total Calls")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))+
  xlab ("Accession")
barplot.sra.no.umi

ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/real_comet_total.png",
  plot = barplot.sra.no.umi,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")



################## Total calls - dietz ##########################
dietz <- read.csv("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/dietz.csv")

#Total number of variants detected by callers across all 10 patients 
barplot.dietz <- ggplot(data = dietz, aes(x = accession, y=total))+
  geom_bar(stat="identity", fill="white", colour = "black",linewidth=0.8)+
  facet_grid(cols = vars(vc))+
  theme_linedraw()+
  theme(strip.text.x = element_text(size = 12,face="bold"),
        strip.text.y = element_text(size = 40),
        axis.text=element_text(size=11, face="bold"),
        axis.title=element_text(size=14, face="bold"),
        axis.text.x = element_text(size = 14,angle = 45, hjust=1, vjust = 1, face="bold"),
        axis.text.y = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = rel(1.5), angle = 90),
        axis.title.x = element_text(size = rel(1.5), angle = 0))+
  ylab ("Total Calls")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))+
  xlab ("Accession")
barplot.dietz

ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/real_comet_total.png",
  plot = barplot.sra.no.umi,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")

############################## Percentage pathogenic variants - comet ###################################################
cosmics.pointplot.no.umi <- sra.no.umi.mod %>%
  filter(vc != "mageri", accession != "SRR15081471", accession != "SRR15081464") %>%
  ggplot(aes(x = vc, y = percentage_pathogenic)) +
  geom_point(aes(colour = accession), size = 8) +
  geom_text_repel(aes(label = total, colour = accession, fontface = "bold"), vjust = 0.5, hjust = -1, 
                  size = 4,nudge_x = 0.1, direction = "y", force = 1, max.iter = 1000) +  # Adding labels
  theme_linedraw() +
  theme(strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 0, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.position = "top") +
  ylim(0, 175) +
  labs(fill = "Dose (mg)") +
  xlab(NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20), limits = c(1, 46)) +
  ylab("Percentage of COSMIC Variants Detected") +
  labs(color = NULL)

cosmics.pointplot.no.umi


ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/percentage_cosmic_calls.png",
  plot = cosmics.pointplot.no.umi,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")


# Group by 'vc' and calculate mean cosmic percentage
result <- sra.no.umi.mod %>%
  filter(vc != "mageri", accession != "SRR15081471", accession != "SRR15081464") %>%
  group_by(vc) %>%
  summarise(mean_cosmic_percentage = mean(percentage_, na.rm = TRUE))

# Print the result
print(result)


############################## Percentage pathogenic variants - WES ###################################################
dietz_pathogenic <- dietz %>%
  filter (accession != "SRR3401418") %>%
  ggplot(aes(x = vc, y = percentage_cosmic)) +
  geom_point(aes(colour = accession), size = 8) +
  geom_text_repel(aes(label = total, colour = accession, fontface = "bold"), vjust = 0.5, hjust = -1, 
                  size = 4, force = 3) +  # Adding labels
  theme_linedraw() +
  theme(strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 0, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.position = "top") +
  ylim(0, 175) +
  xlab(NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20), limits = c(1, 44)) +
  ylab("Percentage of COSMIC Variants Detected") +
  labs(color = NULL)

dietz_pathogenic

ggsave(
  "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/percentage_cosmic_calls.png",
  plot = dietz_pathogenic,
  dpi = 500,
  limitsize = FALSE,
  bg = "transparent")

# Assuming your dataframe is called df
library(dplyr)




# Group by 'vc' and calculate mean cosmic percentage
result <- dietz %>%
  filter (accession != "SRR3401418")%>%
  group_by(vc) %>%
  summarise(mean_cosmic_percentage = mean(percentage_cosmic, na.rm = TRUE))

# Print the result
print(result)


###### Concordance in comet dataset #########
#Concordance of variants <- UpSetR - no_umis

#present just SRR15081472 intersects
myVariantSets <- list(
  bcftools = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/no_umis/COSMIC_SRR15081472.marked_duplicates.bam_bcf.txt"),
  FreeBayes = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/no_umis/COSMIC_SRR15081472.marked_duplicates.bam_fbayes.txt"),
  LoFreq = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/no_umis/COSMIC_SRR15081472.marked_duplicates.bam_lofreq.txt"),
  Mutect2 = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/no_umis/COSMIC_SRR15081472.marked_duplicates.bam_mut2.txt"))

png("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/comet_concordance.png", width = 500, height = 225, units='mm', res = 500)
upset(fromList(myVariantSets), nsets = 6, number.angles = 0, point.size = 6, line.size = 3, 
             mainbar.y.label = "Intersects", sets.x.label = "Total number of variants called", 
             text.scale = c(2, 2, 2, 2, 2, 2.5),order.by = "degree", nintersects = 60, 
             set_size.show = TRUE, set_size.numbers_size = 8, set_size.scale_max = 1500,
             show.numbers = "yes")
dev.off()


###### Concordance in dietz dataset #########

#present just SRR3401415 intersects
myVariantSetsDietz <- list(
  bcftools = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/no_umis/COSMIC_bcftools.SRR3401415.txt"),
  FreeBayes = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/no_umis/COSMIC_gatkfbayes.SRR3401415.txt"),
  LoFreq = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/no_umis/COSMIC_lofreq.SRR3401415.txt"),
  Mutect2 = readLines("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/upsetr/no_umis/COSMIC_mutect2.SRR3401415.txt"))



png("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/dietz_concordance.png", width = 500, height = 225, units='mm', res = 500)
upset(fromList(myVariantSetsDietz), nsets = 6, number.angles = 0, point.size = 6, line.size = 3, 
      mainbar.y.label = "Intersects", sets.x.label = "Total number of variants called", 
      text.scale = c(2, 2, 2, 2, 2, 2.5),order.by = "degree", nintersects = 60, 
      set_size.show = TRUE, set_size.numbers_size = 8, set_size.scale_max = 2750,
      show.numbers = "yes")
dev.off()

mutate_all(~ ifelse(. == "Never", 0, .))


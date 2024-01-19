library(gt)
library(tidyverse)
library(glue)
library(kableExtra)
library(gtExtras)
library("scales")
library(readxl)


data_summaries <- read.csv("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/data_summaries.csv")
variant_caller_summaries<- read.csv("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/variant_caller_summaries.csv")
dietz_dat_sum <- read_excel("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/dietz_data_sum.xlsx")
bam_qc <- read_excel("/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/quality control/bamqc.xlsx")

gt(data_summaries)

#    


data_sum <- data_summaries %>%
  select(-Post.UMI.Correction.Depth)%>%
  filter(Accession != "SRR3401418")%>%
  gt()%>%
  cols_hide (columns = c(Status2)) %>%
  tab_spanner(label = "UMI", columns = matches("Depths"))%>%
  cols_label(Original.Depth = "Raw Depth",
             #Post.UMI.Correction.Depth = "After UMI Read Collapse",
             Post.MarkDuplicates.Depth = "After gatk MarkDuplicates") %>%
  tab_header(title = md("Description of benchmarking datasets")) %>%
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      #Give a thick border below
      cell_borders(sides = "bottom", weight = px(3)),
      #Make text bold
      cell_text(weight = "bold")
    )) %>% 
  cols_align(
    align = c("center"),
    columns = everything()
  ) %>%
  tab_style(
    cells_column_spanners(spanners = everything()),
    style     = list(
      #Give a thick border below
      cell_borders(sides = "bottom", weight = px(3)),
      #Make text bold
      cell_text(weight = "bold")
    )) %>%
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  )

data_sum
gtsave(data_sum, expand = 10, zoom = 10, filename = "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/data_description.png")


dietz_dat_sum
ddata_sum <- dietz_dat_sum %>%
  gt()%>%
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      #Give a thick border below
      cell_borders(sides = "bottom", weight = px(3)),
      #Make text bold
      cell_text(weight = "bold")
    )) %>% 
  cols_align(
    align = c("center"),
    columns = everything()
  ) %>%
  tab_style(
    cells_column_spanners(spanners = everything()),
    style     = list(
      #Give a thick border below
      cell_borders(sides = "bottom", weight = px(3)),
      #Make text bold
      cell_text(weight = "bold")
    )) %>%
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  )
ddata_sum
gtsave(ddata_sum,expand = 10, zoom = 10, filename = "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/ddata_discription.png")



bam_qc
bamqc <- bam_qc %>% mutate(Dataset = recode(Dataset, "Dietz" = "real WES", "COMET" = "real mBC")) %>% 
  filter("Sample Name" != "SRR3401418")  %>%
  select(-"1X", -"5X", -"10X", -"30X", -"≥ 50X",
         -"Median cov", -"Mean cov")%>%
  gt()%>%
  tab_header(title = md("BAM file quality control metrics")) %>%
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      #Give a thick border below
      cell_borders(sides = "bottom", weight = px(3)),
      #Make text bold
      cell_text(weight = "bold")
    )) %>% 
  cols_align(
    align = c("center"),
    columns = everything()
  ) %>%
  tab_style(
    cells_column_spanners(spanners = everything()),
    style     = list(
      #Give a thick border below
      cell_borders(sides = "bottom", weight = px(3)),
      #Make text bold
      cell_text(weight = "bold")
    )) %>%
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  )
bamqc

gtsave(bamqc,expand = 10, zoom = 10, filename = "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Thesis/Chapter 1/Figures/ddata_discription.png")



var_cal <- variant_caller_summaries %>%dplyr::filter(Variant.Caller != "MAGERI") %>%
  #dplyr::filter(Variant.Caller != "UMIErrorCorrect") %>%dplyr::filter(Variant.Caller != "UMI-VarCal") %>%
  #dplyr::filter(Variant.Caller != "UMI-VarCal") %>%
  gt()%>%
  cols_label(Variant.Caller = "Variant Caller", UMI.Awareness = "Native UMI Support", 
             Type.of.Variant = "Supported Variants", Input.Files = "Input Files",
             Version = "Version", Version.Release = "Version Release") %>%
  tab_header(title = md("Variant Caller Overview")) %>%
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      #Give a thick border below
      cell_borders(sides = "bottom", weight = px(3)),
      #Make text bold
      cell_text(weight = "bold")
    )) %>% 
  cols_align(
    align = c("center"),
    columns = everything()
  ) %>%
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  )

var_cal

gtsave(var_cal,expand = 10, zoom = 10, filename = "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/PhD/Publications/Chapter One/bmc_genomics/Figures/variant_caller_summaries.png")




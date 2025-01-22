# File Name: Haplotype_analysis.R
# Author: Kenan Tan
# Date: 2024-01-08
# Description: Description of the script


# Library Imports -------------------------
library(tidyverse)
library(geneHapR)
library(ggsignif)
library(RColorBrewer)
library(ggbeeswarm)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(maps)
library(mapdata)
library(scatterpie)

# Data Import and Setup --------------------
setwd("U:/project/HvCLV(assumption)/vcf")
pheno <- import_AccINFO("//filer-5/agruppen/PBP/tan/gwas/pheno.fam")
vcf <- import_vcf("../../HvCMF6ab/vcf/snp_id_list_res.recode.vcf")

# Haplotypes Processing --------------------
hapResult <- vcf2hap(vcf, hapPrefix = "H", na_drop = TRUE, hetero_remove = TRUE)
hapSummary <- hap_summary(hapResult, hapPrefix = "H")

# Plot Haplotypes Summary --------------------
plotHapTable(hapSummary, hapPrefix = "H")

# Haplotypes Network ------------------------
hapnet <- get_hapNet(hapSummary, AccINFO = pheno, groupName = "continent", na.label = "Unknown")
plotHapNet(hapnet, scale = "log2", show.mutation = 2, col.link = 2, link.width = 2, 
           pie.lim = c(0.5, 2), legend_version = 1, labels = TRUE, legend = c(12, 0), cex.legend = 0.6)

# Haplotypes vs Phenotype --------------------
hapVsPheno(hap = hapResult, phenoName = "FSN_all", pheno = pheno, hapPrefix = "H", 
           minAcc = 5, symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("***", "**", "*", "ns")), mergeFigs = TRUE)

# Plot Phenotype Comparison -----------------
pheno <- read.csv("pheno.fam",sep = "\t")
phenotype_Plot <- function(data, x, y) {
  test_result <- t.test(as.formula(paste(y, "~", x)), data = data[!is.na(data[[x]]), ])
  p_value <- test_result$p.value
  
  ggplot(data[!is.na(data[[x]]),], aes_string(x = x, fill = x, y = y)) +
    geom_boxplot(size = 1, alpha = 0.7) +
    geom_beeswarm(aes_string(color = x), size = 3, alpha = 0.3) +
    scale_color_manual(values = c("Hap1" = "#1f78b4", "Hap2" = "#e31a1c")) +
    scale_fill_manual(values = c("Hap1" = "#1f78b4", "Hap2" = "#e31a1c")) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none") +
    annotate("text", x = 1.5, y = max(data[[y]], na.rm = TRUE) + 2, 
             label = paste("p =", format(p_value, digits = 2)), size = 4)
}

# Call the phenotype plot function
phenotype_Plot(data = pheno, x = "CLV", y = "PTD_field")

# Map and Pie Plot -------------------------
# Data preparation for scatter pie plot
data_summary <- pheno %>%
  group_by(Longitude, Latitude) %>%
  summarise(Hap1 = sum(CMF6ab == "Hap1"), Hap2 = sum(CMF6ab == "Hap2")) %>%
  mutate(Total_count = Hap1 + Hap2) %>%
  filter(!is.na(Longitude) & Total_count > 0)

# World map with scatter pie plot
world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf(fill = "lightgrey", color = "white", size = 0.2) +
  geom_scatterpie(data = data_summary, aes(x = Longitude, y = Latitude, r = Total_count / max(Total_count) * 10), 
                  cols = c("Hap1", "Hap2"), color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("Hap2" = "#1f78b4", "Hap1" = "#e31a1c")) +
  coord_sf(expand = FALSE) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

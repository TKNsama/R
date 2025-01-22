# File Name: GWAS_plot.R
# Author: Kenan Tan
# Date: 2024-05-20
# Description: Description of the script

# 加载必要的库
library(tidyverse)
library(CMplot)
library(patchwork)

setwd("//filer-5/agruppen/PBP/tan/gwas")

# 定义绘制CMplot的函数 --------------------
plot_gwas_results <- function(file_prefix, threshold_logp = 2, af_threshold = 0.1) {
  df <- read_tsv(paste0("output/", file_prefix, "_res.assoc.txt"))
  df$log_p_value <- -log10(df$p_wald)
  
  # 绘制 QQ plot
  df_qq <- df[df$af > af_threshold, c(2, 1, 3, 12)]
  CMplot(df_qq, plot.type = "q", file.name = paste0(file_prefix, "_QQ_plot"))
  
  # 绘制 Manhattan plot
  df_manhattan <- df[df$log_p_value > threshold_logp, c(2, 1, 3, 12)]
  CMplot(df_manhattan, plot.type = "m", 
         threshold = c(0.01, 0.05) / nrow(df_manhattan),
         threshold.col = c('grey', 'black'),
         threshold.lty = c(1, 2), 
         threshold.lwd = c(1, 1), 
         amplify = TRUE, 
         signal.cex = c(1, 1), 
         signal.pch = c(20, 20), 
         signal.col = c("red", "orange"),
         file.name = paste0(file_prefix, "_plot"))
}

# 批量绘制CMplot --------------------
for (i in 1:7) {
  plot_gwas_results(i)
}

# 过滤GWAS结果 --------------------
filter_gwas_results <- function(file_prefix, logp_threshold = 2, af_threshold = 0.1) {
  df <- read_tsv(paste0("output/roop_gaws/", file_prefix, "_res.assoc.txt"))
  df$log_p_value <- -log10(df$p_wald)
  df_filtered <- df[df$log_p_value > logp_threshold & df$af > af_threshold, ]
  write_tsv(df_filtered, paste0(file_prefix, "_filter.tsv"))
}

# 过滤第5组结果
filter_gwas_results(5)

# ggplot2 自定义绘图 --------------------
gwas_ggplot <- function(df, title = "", xlim_range = NULL, y_hline = 5, color_low = "grey", color_high = "red") {
  ggplot(data = df, aes(x = BP, y = -log10(P), color = -log10(P))) +
    geom_point() +
    geom_hline(yintercept = y_hline, size = 1, color = "grey", linetype = "dashed") +
    scale_color_gradient(low = color_low, high = color_high) + 
    labs(x = "", y = "-log10(p)", title = title) +
    facet_wrap(~CHR, ncol = 7, scales = "free_x") +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
          axis.line = element_line(color = "black", size = 0.8),
          legend.position = "none",
          panel.spacing = unit(0.2, "lines"))
}

# 使用函数绘制GWAS结果
df <- read.csv("3_res_3.txt", sep = "\t")
gwas_ggplot(df)

# 局部绘图函数 --------------------
localized_plot <- function(df, xlim_range, y_hline = 5, vlines = NULL, color_low = "white", color_high = "blue") {
  p <- ggplot(data = df, aes(x = BP, y = -log10(P), color = R2_s1_521753849)) +
    geom_point() +
    xlim(xlim_range) +
    geom_hline(yintercept = y_hline, size = 1, color = "grey", linetype = "dashed") +
    scale_color_gradient(low = color_low, high = color_high) + 
    labs(x = "", y = "-log10(p)") +
    theme_bw() +
    theme(panel.border = element_rect(color = "black", size = 1),
          axis.text.x = element_blank(),
          legend.position = "none")
  
  if (!is.null(vlines)) {
    for (vline in vlines) {
      p <- p + geom_vline(xintercept = vline, size = 1, color = "black", linetype = "dashed")
    }
  }
  return(p)
}

# 绘制局部GWAS图
a <- localized_plot(df, xlim_range = c(521600000, 521800000))

# 绘制PPR区域的矩形图 --------------------
plot_region <- function(df, xlim_range) {
  ggplot(df, aes(fill = type)) +
    xlim(xlim_range) +
    geom_rect(aes(xmin = Start.position, xmax = End.position, ymin = -1, ymax = 1)) +
    theme_bw() +
    theme(panel.border = element_rect(color = "black", size = 1),
          legend.position = "none",
          axis.text.y = element_blank()) +
    scale_fill_manual(values = c(cmf = "red", Other = "black", clv = "blue"))
}

# 读取PPR数据并绘制
df <- read.csv("cmf.csv", sep = ",")
b <- plot_region(df, xlim_range = c(521600000, 521800000))

# 组合局部GWAS图和PPR区域图 --------------------
combined_plot <- a / b + plot_layout(heights = c(20, 1))
combined_plot

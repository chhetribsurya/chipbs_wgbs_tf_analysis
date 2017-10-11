library(data.table)
library(ggplot2)
library(dplyr)
library(ggExtra)
library(grid)
#library(MASS)


#####################
#####################
#####################

#output_dir = "/Users/suryachhetri/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_wgbs_ratio_plots"
output_dir = "~/Dropbox/encode_3/chip_bs_wgbs_analysis/plots/motifs_binned_meth_plots"
file_dir <- "~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect" 
file_list <- Sys.glob(file.path(file_dir, "*cpg*.txt"))

### JunD factor:
### Read the file:
meth_table <- fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/JunD_cpg_grouped_chipbs_wgbs_ratio.txt", sep="\t", header= TRUE)
meth_table[meth_table$chip_coverage >100]

meth_table$coverage_anno <- ""
meth_table[meth_table$chip_coverage < 5]$coverage_anno <-  "<5"
meth_table[meth_table$chip_coverage >= 5 & meth_table$chip_coverage < 10]$coverage_anno <- "5-10"
meth_table[meth_table$chip_coverage >= 10 & meth_table$chip_coverage < 20]$coverage_anno <- "10-20"
meth_table[meth_table$chip_coverage >= 20 & meth_table$chip_coverage < 50]$coverage_anno <- "20-50"
meth_table[meth_table$chip_coverage >= 50 & meth_table$chip_coverage < 100]$coverage_anno <- "50-100"
meth_table[meth_table$chip_coverage >= 100 ]$coverage_anno <- ">100"

#color_manual <-  c("red", "yellow", "black", "green", "blue", "darkblue")
pdf(file.path(output_dir,"JUND_cpg_methylation_profile_with_meth_fraction_plot1.pdf"))

meth_table$coverage_anno <- factor(meth_table$coverage_anno, levels=c("<5","5-10","10-20","20-50","50-100",">100"))
#ggplot(meth_table, aes(chip_coverage,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
ggplot(meth_table, aes(wgbs_meth_perc,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
  geom_point(alpha = 0.5) +
  geom_density2d(colour="black", h=0.3) + 
  labs(x="WGBS methylation fraction", y="log10 (ChIPBS:WGBS ratio)", color="Coverage") + 
  theme_bw() + 
  ggtitle("JUND binding sites CpG methylation distribution") + 
  theme(
    axis.text.y = element_text(size=8),
    plot.title=element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
    legend.title=element_text(face="bold")
  ) + scale_y_continuous(limits=c(-5,5), expand=c(0,0))
#guides(fill=guide_legend(title="TF Category"))
#ggMarginal(this_plot + theme_gray(), col = "darkblue")  #scale_color_manual(values=color_manual)

dev.off()


meth_table <- fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/JunD_cpg_grouped_chipbs_wgbs_ratio.txt", sep="\t", header= TRUE)
meth_table[meth_table$chip_coverage >100]

meth_table$coverage_anno <- ""
meth_table[meth_table$chip_coverage < 5]$coverage_anno <-  "<5"
meth_table[meth_table$chip_coverage >= 5 & meth_table$chip_coverage < 10]$coverage_anno <- "5-10"
meth_table[meth_table$chip_coverage >= 10 & meth_table$chip_coverage < 20]$coverage_anno <- "10-20"
meth_table[meth_table$chip_coverage >= 20 & meth_table$chip_coverage < 50]$coverage_anno <- "20-50"
meth_table[meth_table$chip_coverage >= 50 & meth_table$chip_coverage < 100]$coverage_anno <- "50-100"
meth_table[meth_table$chip_coverage >= 100 ]$coverage_anno <- ">100"

#color_manual <-  c("red", "yellow", "black", "green", "blue", "darkblue")
pdf(file.path(output_dir,"JUND_cpg_methylation_profile_with_chipbs_coverage_plot2.pdf"))

meth_table$coverage_anno <- factor(meth_table$coverage_anno, levels=c("<5","5-10","10-20","20-50","50-100",">100"))
ggplot(meth_table, aes(chip_coverage,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
#ggplot(meth_table, aes(wgbs_meth_perc,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
  geom_point(alpha = 0.5) +
  geom_density2d(colour="black", h=0.3) + 
  labs(x="ChIPBS Coverage", y="log10 (ChIPBS:WGBS ratio)", color="Coverage") + 
  theme_bw() + 
  ggtitle("JUND binding sites CpG methylation distribution") + 
  theme(
    axis.text.y = element_text(size=8),
    plot.title=element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
    legend.title=element_text(face="bold")
  ) + scale_y_continuous(limits=c(-5,5), expand=c(0,0))
#guides(fill=guide_legend(title="TF Category"))
#ggMarginal(this_plot + theme_gray(), col = "darkblue")  #scale_color_manual(values=color_manual)

dev.off()






### Pol2 factor:
### Read the file:
meth_table <- fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/Pol2_cpg_grouped_chipbs_wgbs_ratio.txt", sep="\t", header= TRUE)
meth_table[meth_table$chip_coverage >100]

meth_table$coverage_anno <- ""
meth_table[meth_table$chip_coverage < 5]$coverage_anno <-  "<5"
meth_table[meth_table$chip_coverage >= 5 & meth_table$chip_coverage < 10]$coverage_anno <- "5-10"
meth_table[meth_table$chip_coverage >= 10 & meth_table$chip_coverage < 20]$coverage_anno <- "10-20"
meth_table[meth_table$chip_coverage >= 20 & meth_table$chip_coverage < 50]$coverage_anno <- "20-50"
meth_table[meth_table$chip_coverage >= 50 & meth_table$chip_coverage < 100]$coverage_anno <- "50-100"
meth_table[meth_table$chip_coverage >= 100 ]$coverage_anno <- ">100"

pdf(file.path(output_dir,"Pol2_cpg_methylation_profile_with_meth_fraction_plot1.pdf"))

meth_table$coverage_anno <- factor(meth_table$coverage_anno, levels=c("<5","5-10","10-20","20-50","50-100",">100"))
ggplot(meth_table, aes(wgbs_meth_perc,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
  geom_point(alpha = 0.5) +
  geom_density2d(colour="black", h=0.3) + 
  labs(x="WGBS methylation fraction", y="log10 (ChIPBS:WGBS ratio)", color="Coverage") + 
  theme_bw() + 
  ggtitle("Pol2 binding sites CpG methylation distribution") + 
  theme(
    axis.text.y = element_text(size=8),
    plot.title=element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
    legend.title=element_text(face="bold")
  ) + scale_y_continuous(limits=c(-5,5), expand=c(0,0))
dev.off()


### Read the file:
meth_table <- fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/Pol2_cpg_grouped_chipbs_wgbs_ratio.txt", sep="\t", header= TRUE)
meth_table[meth_table$chip_coverage >100]

meth_table$coverage_anno <- ""
meth_table[meth_table$chip_coverage < 5]$coverage_anno <-  "<5"
meth_table[meth_table$chip_coverage >= 5 & meth_table$chip_coverage < 10]$coverage_anno <- "5-10"
meth_table[meth_table$chip_coverage >= 10 & meth_table$chip_coverage < 20]$coverage_anno <- "10-20"
meth_table[meth_table$chip_coverage >= 20 & meth_table$chip_coverage < 50]$coverage_anno <- "20-50"
meth_table[meth_table$chip_coverage >= 50 & meth_table$chip_coverage < 100]$coverage_anno <- "50-100"
meth_table[meth_table$chip_coverage >= 100 ]$coverage_anno <- ">100"

pdf(file.path(output_dir,"Pol2_cpg_methylation_profile_with_chipbs_coverage_plot2.pdf"))

meth_table$coverage_anno <- factor(meth_table$coverage_anno, levels=c("<5","5-10","10-20","20-50","50-100",">100"))
ggplot(meth_table, aes(chip_coverage,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
  geom_point(alpha = 0.5) +
  geom_density2d(colour="black", h=0.3) + 
  labs(x="ChIPBS Coverage", y="log10 (ChIPBS:WGBS normalized ratio)", color="Coverage") + 
  theme_bw() + 
  ggtitle("Pol2 binding sites CpG methylation distribution") + 
  theme(
    axis.text.y = element_text(size=8),
    plot.title=element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
    legend.title=element_text(face="bold")
  ) + scale_y_continuous(limits=c(-5,5), expand=c(0,0))
dev.off()


### TAF1 factor:
### Read the file:
meth_table <- fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/TAF1_cpg_grouped_chipbs_wgbs_ratio.txt", sep="\t", header= TRUE)
meth_table[meth_table$chip_coverage >100]

meth_table$coverage_anno <- ""
meth_table[meth_table$chip_coverage < 5]$coverage_anno <-  "<5"
meth_table[meth_table$chip_coverage >= 5 & meth_table$chip_coverage < 10]$coverage_anno <- "5-10"
meth_table[meth_table$chip_coverage >= 10 & meth_table$chip_coverage < 20]$coverage_anno <- "10-20"
meth_table[meth_table$chip_coverage >= 20 & meth_table$chip_coverage < 50]$coverage_anno <- "20-50"
meth_table[meth_table$chip_coverage >= 50 & meth_table$chip_coverage < 100]$coverage_anno <- "50-100"
meth_table[meth_table$chip_coverage >= 100 ]$coverage_anno <- ">100"

pdf(file.path(output_dir,"TAF1_cpg_methylation_profile_with_meth_fraction_plot1.pdf"))

meth_table$coverage_anno <- factor(meth_table$coverage_anno, levels=c("<5","5-10","10-20","20-50","50-100",">100"))
ggplot(meth_table, aes(wgbs_meth_perc,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
  geom_point(alpha = 0.5) +
  geom_density2d(colour="black", h=0.3) + 
  labs(x="WGBS methylation fraction", y="log10 (ChIPBS:WGBS ratio)", color="Coverage") + 
  theme_bw() + 
  ggtitle("TAF1 binding sites CpG methylation distribution") + 
  theme(
    axis.text.y = element_text(size=8),
    plot.title=element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
    legend.title=element_text(face="bold")
  ) + scale_y_continuous(limits=c(-5,5), expand=c(0,0))
dev.off()


### Read the file:
meth_table <- fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/TAF1_cpg_grouped_chipbs_wgbs_ratio.txt", sep="\t", header= TRUE)
meth_table[meth_table$chip_coverage >100]

meth_table$coverage_anno <- ""
meth_table[meth_table$chip_coverage < 5]$coverage_anno <-  "<5"
meth_table[meth_table$chip_coverage >= 5 & meth_table$chip_coverage < 10]$coverage_anno <- "5-10"
meth_table[meth_table$chip_coverage >= 10 & meth_table$chip_coverage < 20]$coverage_anno <- "10-20"
meth_table[meth_table$chip_coverage >= 20 & meth_table$chip_coverage < 50]$coverage_anno <- "20-50"
meth_table[meth_table$chip_coverage >= 50 & meth_table$chip_coverage < 100]$coverage_anno <- "50-100"
meth_table[meth_table$chip_coverage >= 100 ]$coverage_anno <- ">100"

pdf(file.path(output_dir,"TAF1_cpg_methylation_profile_with_chipbs_coverage_plot2.pdf"))

meth_table$coverage_anno <- factor(meth_table$coverage_anno, levels=c("<5","5-10","10-20","20-50","50-100",">100"))
ggplot(meth_table, aes(chip_coverage,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
  geom_point(alpha = 0.5) +
  geom_density2d(colour="black", h=0.3) + 
  labs(x="ChIPBS Coverage", y="log10 (ChIPBS:WGBS normalized ratio)", color="Coverage") + 
  theme_bw() + 
  ggtitle("TAF1 binding sites CpG methylation distribution") + 
  theme(
    axis.text.y = element_text(size=8),
    plot.title=element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
    legend.title=element_text(face="bold")
  ) + scale_y_continuous(limits=c(-5,5), expand=c(0,0))
dev.off()



### YY1 factor:
### Read the file:
meth_table <- fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/YY1_cpg_grouped_chipbs_wgbs_ratio.txt", sep="\t", header= TRUE)
meth_table[meth_table$chip_coverage >100]

meth_table$coverage_anno <- ""
meth_table[meth_table$chip_coverage < 5]$coverage_anno <-  "<5"
meth_table[meth_table$chip_coverage >= 5 & meth_table$chip_coverage < 10]$coverage_anno <- "5-10"
meth_table[meth_table$chip_coverage >= 10 & meth_table$chip_coverage < 20]$coverage_anno <- "10-20"
meth_table[meth_table$chip_coverage >= 20 & meth_table$chip_coverage < 50]$coverage_anno <- "20-50"
meth_table[meth_table$chip_coverage >= 50 & meth_table$chip_coverage < 100]$coverage_anno <- "50-100"
meth_table[meth_table$chip_coverage >= 100 ]$coverage_anno <- ">100"

pdf(file.path(output_dir,"YY1_cpg_methylation_profile_with_meth_fraction_plot1.pdf"))

meth_table$coverage_anno <- factor(meth_table$coverage_anno, levels=c("<5","5-10","10-20","20-50","50-100",">100"))
ggplot(meth_table, aes(wgbs_meth_perc,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
  geom_point(alpha = 0.5) +
  geom_density2d(colour="black", h=0.3) + 
  labs(x="WGBS methylation fraction", y="log10 (ChIPBS:WGBS ratio)", color="Coverage") + 
  theme_bw() + 
  ggtitle("YY1 binding sites CpG methylation distribution") + 
  theme(
    axis.text.y = element_text(size=8),
    plot.title=element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
    legend.title=element_text(face="bold")
  ) + scale_y_continuous(limits=c(-5,5), expand=c(0,0))
dev.off()


### Read the file:
meth_table <- fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/YY1_cpg_grouped_chipbs_wgbs_ratio.txt", sep="\t", header= TRUE)
meth_table[meth_table$chip_coverage >100]

meth_table$coverage_anno <- ""
meth_table[meth_table$chip_coverage < 5]$coverage_anno <-  "<5"
meth_table[meth_table$chip_coverage >= 5 & meth_table$chip_coverage < 10]$coverage_anno <- "5-10"
meth_table[meth_table$chip_coverage >= 10 & meth_table$chip_coverage < 20]$coverage_anno <- "10-20"
meth_table[meth_table$chip_coverage >= 20 & meth_table$chip_coverage < 50]$coverage_anno <- "20-50"
meth_table[meth_table$chip_coverage >= 50 & meth_table$chip_coverage < 100]$coverage_anno <- "50-100"
meth_table[meth_table$chip_coverage >= 100 ]$coverage_anno <- ">100"

pdf(file.path(output_dir,"YY1_cpg_methylation_profile_with_chipbs_coverage_plot2.pdf"))

meth_table$coverage_anno <- factor(meth_table$coverage_anno, levels=c("<5","5-10","10-20","20-50","50-100",">100"))
ggplot(meth_table, aes(chip_coverage,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
  geom_point(alpha = 0.5) +
  geom_density2d(colour="black", h=0.3) + 
  labs(x="ChIPBS Coverage", y="log10 (ChIPBS:WGBS normalized ratio)", color="Coverage") + 
  theme_bw() + 
  ggtitle("YY1 binding sites CpG methylation distribution") + 
  theme(
    axis.text.y = element_text(size=8),
    plot.title=element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
    legend.title=element_text(face="bold")
  ) + scale_y_continuous(limits=c(-5,5), expand=c(0,0))
dev.off()


### ZBTB33 factor:
### Read the file:
meth_table <- fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/ZBTB33_cpg_grouped_chipbs_wgbs_ratio.txt", sep="\t", header= TRUE)
meth_table[meth_table$chip_coverage >100]

meth_table$coverage_anno <- ""
meth_table[meth_table$chip_coverage < 5]$coverage_anno <-  "<5"
meth_table[meth_table$chip_coverage >= 5 & meth_table$chip_coverage < 10]$coverage_anno <- "5-10"
meth_table[meth_table$chip_coverage >= 10 & meth_table$chip_coverage < 20]$coverage_anno <- "10-20"
meth_table[meth_table$chip_coverage >= 20 & meth_table$chip_coverage < 50]$coverage_anno <- "20-50"
meth_table[meth_table$chip_coverage >= 50 & meth_table$chip_coverage < 100]$coverage_anno <- "50-100"
meth_table[meth_table$chip_coverage >= 100 ]$coverage_anno <- ">100"

pdf(file.path(output_dir,"ZBTB33_cpg_methylation_profile_with_meth_fraction_plot1.pdf"))

meth_table$coverage_anno <- factor(meth_table$coverage_anno, levels=c("<5","5-10","10-20","20-50","50-100",">100"))
ggplot(meth_table, aes(wgbs_meth_perc,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
  geom_point(alpha = 0.5) +
  geom_density2d(colour="black", h=0.3) + 
  labs(x="WGBS methylation fraction", y="log10 (ChIPBS:WGBS ratio)", color="Coverage") + 
  theme_bw() + 
  ggtitle("ZBTB33 binding sites CpG methylation distribution") + 
  theme(
    axis.text.y = element_text(size=8),
    plot.title=element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
    legend.title=element_text(face="bold")
  ) + scale_y_continuous(limits=c(-5,5), expand=c(0,0))
dev.off()


### Read the file:
meth_table <- fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/ZBTB33_cpg_grouped_chipbs_wgbs_ratio.txt", sep="\t", header= TRUE)
meth_table[meth_table$chip_coverage >100]

meth_table$coverage_anno <- ""
meth_table[meth_table$chip_coverage < 5]$coverage_anno <-  "<5"
meth_table[meth_table$chip_coverage >= 5 & meth_table$chip_coverage < 10]$coverage_anno <- "5-10"
meth_table[meth_table$chip_coverage >= 10 & meth_table$chip_coverage < 20]$coverage_anno <- "10-20"
meth_table[meth_table$chip_coverage >= 20 & meth_table$chip_coverage < 50]$coverage_anno <- "20-50"
meth_table[meth_table$chip_coverage >= 50 & meth_table$chip_coverage < 100]$coverage_anno <- "50-100"
meth_table[meth_table$chip_coverage >= 100 ]$coverage_anno <- ">100"

pdf(file.path(output_dir,"ZBTB33_cpg_methylation_profile_with_chipbs_coverage_plot2.pdf"))

meth_table$coverage_anno <- factor(meth_table$coverage_anno, levels=c("<5","5-10","10-20","20-50","50-100",">100"))
ggplot(meth_table, aes(chip_coverage,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
  geom_point(alpha = 0.5) +
  geom_density2d(colour="black", h=0.3) + 
  labs(x="ChIPBS Coverage", y="log10 (ChIPBS:WGBS normalized ratio)", color="Coverage") + 
  theme_bw() + 
  ggtitle("ZBTB33 binding sites CpG methylation distribution") + 
  theme(
    axis.text.y = element_text(size=8),
    plot.title=element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
    legend.title=element_text(face="bold")
  ) + scale_y_continuous(limits=c(-5,5), expand=c(0,0))
dev.off()


# smoothScatter(wgbs_meth_perc~chip_wgbs_ratio_norm_log)
# with(meth_table, smoothScatter(chip_meth_perc~wgbs_meth_perc))


###########################################
### Methylation profile in each bins 20bp up and downstream of motifs:
###########################################
### JunD:

### Type 1 plot: Barplot using facet_grid() for meth distribution and CpG counts in bin:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/JunD_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]

pdf(file.path(output_dir,"JunD_cpg_methylation_plot_1.pdf"))
df_1 <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))
names(df_1) <- c("bin_mid", "cpg_count_in_bin", "assay_group", "meth_perc")

# df_1 <- subset(df_1, cpg_count_in_bin >=20)
df_2 <- melt(subset(df_1, select=c(bin_mid, meth_perc, cpg_count_in_bin, assay_group)), id=c("bin_mid", "assay_group"))

meth_table_names <- c("meth_perc" = "Methylation %", "cpg_count_in_bin" = "CpG counts")
ggplot(df_2, aes(x=bin_mid, y=value, fill=assay_group))+geom_bar(stat='identity', position='dodge') + 
          theme_gray() +
          facet_grid(variable~., scales = "free_y", labeller = as_labeller(meth_table_names) ) +
          labs(x="Distance from Center (bp)", y="values", fill="Assay Type")+
          theme(strip.text = element_text(face="bold", size=9),
          strip.background = element_rect(fill="lightblue", colour="black",size=1)) 
         # scale_y_continuous("Value") + 
         # scale_x_continuous("Distance from Center (bp)", limits=c(-20, 20))

dev.off()


### Type 2 plot : Combining barplot and boxplot across 20bp bins:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/JunD_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))

pdf(file.path(output_dir,"JunD_cpg_methylation_plot_2.pdf"))

### Barplot for meth distribution:
plot1 <- ggplot(mm, aes(x=bin_mid, y=value, fill=variable))+geom_bar(stat='identity', position='dodge') + theme_gray() +
          ##geom_text(aes(label=cpg_count_in_bin), size=1.1, hjust=0.5) + 
          ggtitle("JunD binding sites meth profile") +
          scale_x_continuous("Distance from Center (bp)") +
          scale_y_continuous("Methylation percent",limits=c(0,1)) +
          labs(fill="Assay Type") +
          theme(
            plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
            #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
            #axis.title.y = element_text(angle=0),
            legend.title=element_text(face="bold")
          )

### Barplot for CpG count in bins:
# plot2 <- ggplot(mm, aes(x=bin_mid, y=cpg_count_in_bin,fill=variable))+
#           geom_bar(stat='identity') + theme_gray() +
#           scale_x_continuous("Distance from Center (bp)") +
#           scale_y_continuous("CpG counts in bin") +
#           labs(fill="Assay") +
#           theme(
#             plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
#             #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
#             #axis.title.y = element_text(angle=0),
#             legend.title=element_text(face="bold")
#           )

### Box plot meth distribution for each CpG:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/JunD_chip_wgbs_motifs_intersect_20bp_each_cpgs.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc)), id=c("bin_mid"))

plot2 <- ggplot(mm, aes(x=as.factor(bin_mid), y=value, fill=variable)) + 
      geom_boxplot(outlier.shape = NA) + theme_gray()+
      theme(
        legend.title=element_text(face="bold"), 
        axis.title.x = element_text(size=10, hjust = 0.6, colour = "black"),
        axis.ticks.x =element_line(size=0.2, color="black")) +
        scale_x_discrete(name="Distance from Center (bp)", breaks=c(-19.5,-10.5,0.5,10.5,19.5), labels=c("-20", "-10", "0", "10", "20")) +
        scale_y_continuous(limits=c(0,1)) +
        labs(y="Methylation percent", fill="Assay Type")


grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
dev.off()



### Type 3 plot : Combining barplot and violinplot across 20bp bins:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/JunD_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))

pdf(file.path(output_dir,"JunD_cpg_methylation_plot_3.pdf"))

### Barplot for meth distribution:
plot1 <- ggplot(mm, aes(x=bin_mid, y=value, fill=variable))+geom_bar(stat='identity', position='dodge') + theme_gray() +
          #geom_text(aes(label=cpg_count_in_bin), size=1.1, hjust=0.5) + 
          ggtitle("JunD binding sites meth profile") +
          scale_x_continuous("Distance from Center (bp)") +
          scale_y_continuous("Methylation percent",limits=c(0,1)) +
          labs(fill="Assay Type") +
          theme(
            plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
            #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
            #axis.title.y = element_text(angle=0),
            legend.title=element_text(face="bold")
          )

### Barplot for CpG count in bins:
# plot2 <- ggplot(mm, aes(x=bin_mid, y=cpg_count_in_bin,fill=variable))+
#           geom_bar(stat='identity') + theme_gray() +
#           scale_x_continuous("Distance from Center (bp)") +
#           scale_y_continuous("CpG counts in bin") +
#           labs(fill="Assay") +
#           theme(
#             plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
#             #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
#             #axis.title.y = element_text(angle=0),
#             legend.title=element_text(face="bold")
#           )

### violin plot meth distribution for each CpG:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/JunD_chip_wgbs_motifs_intersect_20bp_each_cpgs.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc)), id=c("bin_mid"))

plot2 <- ggplot(mm, aes(x=as.factor(bin_mid), y=value, fill=variable)) + 
      geom_violin(trim=FALSE) + theme_gray()+
      theme(
        legend.title=element_text(face="bold"), 
        axis.title.x = element_text(size=10, hjust = 0.6, colour = "black"),
        axis.ticks.x =element_line(size=0.2, color="black")) +
        scale_x_discrete(name="Distance from Center (bp)", breaks=c(-19.5,-10.5,0.5,10.5,19.5), labels=c("-20", "-10", "0", "10", "20")) +
        scale_y_continuous(limits=c(0,1)) +
        labs(y="Methylation percent", fill="Assay Type")

grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
dev.off()




### Pol2:
###

### Type 1 plot: Barplot using facet_grid() for meth distribution and CpG counts in bin:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/Pol2_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]

pdf(file.path(output_dir,"Pol2_cpg_methylation_plot_1.pdf"))
df_1 <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))
names(df_1) <- c("bin_mid", "cpg_count_in_bin", "assay_group", "meth_perc")

# df_1 <- subset(df_1, cpg_count_in_bin >=20)
df_2 <- melt(subset(df_1, select=c(bin_mid, meth_perc, cpg_count_in_bin, assay_group)), id=c("bin_mid", "assay_group"))

meth_table_names <- c("meth_perc" = "Methylation %", "cpg_count_in_bin" = "CpG counts")
ggplot(df_2, aes(x=bin_mid, y=value, fill=assay_group))+geom_bar(stat='identity', position='dodge') + 
          theme_gray() +
          facet_grid(variable~., scales = "free_y", labeller = as_labeller(meth_table_names) ) +
          labs(x="Distance from Center (bp)", y="values", fill="Assay Type")+
          theme(strip.text = element_text(face="bold", size=9),
          strip.background = element_rect(fill="lightblue", colour="black",size=1)) 
         # scale_y_continuous("Value") + 
         # scale_x_continuous("Distance from Center (bp)", limits=c(-20, 20))

dev.off()


### Type 2 plot : Combining barplot and boxplot across 20bp bins:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/Pol2_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))

pdf(file.path(output_dir,"Pol2_cpg_methylation_plot_2.pdf"))

### Barplot for meth distribution:
plot1 <- ggplot(mm, aes(x=bin_mid, y=value, fill=variable))+geom_bar(stat='identity', position='dodge') + theme_gray() +
          #geom_text(aes(label=cpg_count_in_bin), size=1.1, hjust=0.5) + 
          ggtitle("Pol2 binding sites meth profile") +
          scale_x_continuous("Distance from Center (bp)") +
          scale_y_continuous("Methylation percent",limits=c(0,1)) +
          labs(fill="Assay Type") +
          theme(
            plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
            #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
            #axis.title.y = element_text(angle=0),
            legend.title=element_text(face="bold")
          )
### Barplot for CpG count in bins:
# plot2 <- ggplot(mm, aes(x=bin_mid, y=cpg_count_in_bin,fill=variable))+
#           geom_bar(stat='identity') + theme_gray() +
#           scale_x_continuous("Distance from Center (bp)") +
#           scale_y_continuous("CpG counts in bin") +
#           labs(fill="Assay") +
#           theme(
#             plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
#             #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
#             #axis.title.y = element_text(angle=0),
#             legend.title=element_text(face="bold")
#           )

### Box plot meth distribution for each CpG:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/Pol2_chip_wgbs_motifs_intersect_20bp_each_cpgs.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc)), id=c("bin_mid"))

plot2 <- ggplot(mm, aes(x=as.factor(bin_mid), y=value, fill=variable)) + 
      geom_boxplot(outlier.shape = NA) + theme_gray()+
      theme(
        legend.title=element_text(face="bold"), 
        axis.title.x = element_text(size=10, hjust = 0.6, colour = "black"),
        axis.ticks.x =element_line(size=0.2, color="black")) +
        scale_x_discrete(name="Distance from Center (bp)", breaks=c(-19.5,-10.5,0.5,10.5,19.5), labels=c("-20", "-10", "0", "10", "20")) +
        scale_y_continuous(limits=c(0,1)) +
        labs(y="Methylation percent", fill="Assay Type")

grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
dev.off()



### Type 3 plot : Combining barplot and violin plot across 20bp bins:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/Pol2_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))

pdf(file.path(output_dir,"Pol2_cpg_methylation_plot_3.pdf"))

### Barplot for meth distribution:
plot1 <- ggplot(mm, aes(x=bin_mid, y=value, fill=variable))+geom_bar(stat='identity', position='dodge') + theme_gray() +
          #geom_text(aes(label=cpg_count_in_bin), size=1.1, hjust=0.5) + 
          ggtitle("Pol2 binding sites meth profile") +
          scale_x_continuous("Distance from Center (bp)") +
          scale_y_continuous("Methylation percent",limits=c(0,1)) +
          labs(fill="Assay Type") +
          theme(
            plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
            #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
            #axis.title.y = element_text(angle=0),
            legend.title=element_text(face="bold")
          )
### Barplot for CpG count in bins:
# plot2 <- ggplot(mm, aes(x=bin_mid, y=cpg_count_in_bin,fill=variable))+
#           geom_bar(stat='identity') + theme_gray() +
#           scale_x_continuous("Distance from Center (bp)") +
#           scale_y_continuous("CpG counts in bin") +
#           labs(fill="Assay") +
#           theme(
#             plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
#             #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
#             #axis.title.y = element_text(angle=0),
#             legend.title=element_text(face="bold")
#           )

### violin plot meth distribution for each CpG:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/Pol2_chip_wgbs_motifs_intersect_20bp_each_cpgs.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc)), id=c("bin_mid"))

plot2 <- ggplot(mm, aes(x=as.factor(bin_mid), y=value, fill=variable)) + 
      geom_violin(trim=FALSE) + theme_gray()+
      theme(
        legend.title=element_text(face="bold"), 
        axis.title.x = element_text(size=10, hjust = 0.6, colour = "black"),
        axis.ticks.x =element_line(size=0.2, color="black")) +
        scale_x_discrete(name="Distance from Center (bp)", breaks=c(-19.5,-10.5,0.5,10.5,19.5), labels=c("-20", "-10", "0", "10", "20")) +
        scale_y_continuous(limits=c(0,1)) +
        labs(y="Methylation percent", fill="Assay Type")


grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
dev.off()





### TAF1:
###

### Type 1 plot: Barplot using facet_grid() for meth distribution and CpG counts in bin:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/TAF1_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]

pdf(file.path(output_dir,"TAF1_cpg_methylation_plot_1.pdf"))
df_1 <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))
names(df_1) <- c("bin_mid", "cpg_count_in_bin", "assay_group", "meth_perc")

# df_1 <- subset(df_1, cpg_count_in_bin >=20)
df_2 <- melt(subset(df_1, select=c(bin_mid, meth_perc, cpg_count_in_bin, assay_group)), id=c("bin_mid", "assay_group"))

meth_table_names <- c("meth_perc" = "Methylation %", "cpg_count_in_bin" = "CpG counts")
ggplot(df_2, aes(x=bin_mid, y=value, fill=assay_group))+geom_bar(stat='identity', position='dodge') + 
          theme_gray() +
          facet_grid(variable~., scales = "free_y", labeller = as_labeller(meth_table_names) ) +
          labs(x="Distance from Center (bp)", y="values", fill="Assay Type")+
          theme(strip.text = element_text(face="bold", size=9),
          strip.background = element_rect(fill="lightblue", colour="black",size=1)) 
         # scale_y_continuous("Value") + 
         # scale_x_continuous("Distance from Center (bp)", limits=c(-20, 20))

dev.off()


### Type 2 plot : Combining barplot and boxplot across 20bp bins:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/TAF1_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))

pdf(file.path(output_dir,"TAF1_cpg_methylation_plot_2.pdf"))

### Barplot for meth distribution:
plot1 <- ggplot(mm, aes(x=bin_mid, y=value, fill=variable))+geom_bar(stat='identity', position='dodge') + theme_gray() +
          #geom_text(aes(label=cpg_count_in_bin), size=1.1, hjust=0.5) + 
          ggtitle("TAF1 binding sites meth profile") +
          scale_x_continuous("Distance from Center (bp)") +
          scale_y_continuous("Methylation percent",limits=c(0,1)) +
          labs(fill="Assay Type") +
          theme(
            plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
            #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
            #axis.title.y = element_text(angle=0),
            legend.title=element_text(face="bold")
          )
### Barplot for CpG count in bins:
# plot2 <- ggplot(mm, aes(x=bin_mid, y=cpg_count_in_bin,fill=variable))+
#           geom_bar(stat='identity') + theme_gray() +
#           scale_x_continuous("Distance from Center (bp)") +
#           scale_y_continuous("CpG counts in bin") +
#           labs(fill="Assay") +
#           theme(
#             plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
#             #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
#             #axis.title.y = element_text(angle=0),
#             legend.title=element_text(face="bold")
#           )

### Box plot meth distribution for each CpG:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/TAF1_chip_wgbs_motifs_intersect_20bp_each_cpgs.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc)), id=c("bin_mid"))

plot2 <- ggplot(mm, aes(x=as.factor(bin_mid), y=value, fill=variable)) + 
      geom_boxplot(outlier.shape = NA) + theme_gray()+
      theme(
        legend.title=element_text(face="bold"), 
        axis.title.x = element_text(size=10, hjust = 0.6, colour = "black"),
        axis.ticks.x =element_line(size=0.2, color="black")) +
        scale_x_discrete(name="Distance from Center (bp)", breaks=c(-19.5,-10.5,0.5,10.5,19.5), labels=c("-20", "-10", "0", "10", "20")) +
        scale_y_continuous(limits=c(0,1)) +
        labs(y="Methylation percent", fill="Assay Type")

grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
dev.off()



### Type 3 plot : Combining barplot and violin plot across 20bp bins:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/TAF1_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))

pdf(file.path(output_dir,"TAF1_cpg_methylation_plot_3.pdf"))

### Barplot for meth distribution:
plot1 <- ggplot(mm, aes(x=bin_mid, y=value, fill=variable))+geom_bar(stat='identity', position='dodge') + theme_gray() +
          #geom_text(aes(label=cpg_count_in_bin), size=1.1, hjust=0.5) + 
          ggtitle("TAF1 binding sites meth profile") +
          scale_x_continuous("Distance from Center (bp)") +
          scale_y_continuous("Methylation percent",limits=c(0,1)) +
          labs(fill="Assay Type") +
          theme(
            plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
            #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
            #axis.title.y = element_text(angle=0),
            legend.title=element_text(face="bold")
          )
### Barplot for CpG count in bins:
# plot2 <- ggplot(mm, aes(x=bin_mid, y=cpg_count_in_bin,fill=variable))+
#           geom_bar(stat='identity') + theme_gray() +
#           scale_x_continuous("Distance from Center (bp)") +
#           scale_y_continuous("CpG counts in bin") +
#           labs(fill="Assay") +
#           theme(
#             plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
#             #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
#             #axis.title.y = element_text(angle=0),
#             legend.title=element_text(face="bold")
#           )

### violin plot meth distribution for each CpG:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/TAF1_chip_wgbs_motifs_intersect_20bp_each_cpgs.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc)), id=c("bin_mid"))

plot2 <- ggplot(mm, aes(x=as.factor(bin_mid), y=value, fill=variable)) + 
      geom_violin(trim=FALSE) + theme_gray()+
      theme(
        legend.title=element_text(face="bold"), 
        axis.title.x = element_text(size=10, hjust = 0.6, colour = "black"),
        axis.ticks.x =element_line(size=0.2, color="black")) +
        scale_x_discrete(name="Distance from Center (bp)", breaks=c(-19.5,-10.5,0.5,10.5,19.5), labels=c("-20", "-10", "0", "10", "20")) +
        scale_y_continuous(limits=c(0,1)) +
        labs(y="Methylation percent", fill="Assay Type")


grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
dev.off()






### YY1:
###

### Type 1 plot: Barplot using facet_grid() for meth distribution and CpG counts in bin:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/YY1_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]

pdf(file.path(output_dir,"YY1_cpg_methylation_plot_1.pdf"))
df_1 <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))
names(df_1) <- c("bin_mid", "cpg_count_in_bin", "assay_group", "meth_perc")

# df_1 <- subset(df_1, cpg_count_in_bin >=20)
df_2 <- melt(subset(df_1, select=c(bin_mid, meth_perc, cpg_count_in_bin, assay_group)), id=c("bin_mid", "assay_group"))

meth_table_names <- c("meth_perc" = "Methylation %", "cpg_count_in_bin" = "CpG counts")
ggplot(df_2, aes(x=bin_mid, y=value, fill=assay_group))+geom_bar(stat='identity', position='dodge') + 
          theme_gray() +
          facet_grid(variable~., scales = "free_y", labeller = as_labeller(meth_table_names) ) +
          labs(x="Distance from Center (bp)", y="values", fill="Assay Type")+
          theme(strip.text = element_text(face="bold", size=9),
          strip.background = element_rect(fill="lightblue", colour="black",size=1)) 
         # scale_y_continuous("Value") + 
         # scale_x_continuous("Distance from Center (bp)", limits=c(-20, 20))

dev.off()


### Type 2 plot : Combining barplot and boxplot across 20bp bins:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/YY1_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))

pdf(file.path(output_dir,"YY1_cpg_methylation_plot_2.pdf"))

### Barplot for meth distribution:
plot1 <- ggplot(mm, aes(x=bin_mid, y=value, fill=variable))+geom_bar(stat='identity', position='dodge') + theme_gray() +
          ##geom_text(aes(label=cpg_count_in_bin), size=1.1, hjust=0.5) + 
          ggtitle("YY1 binding sites meth profile") +
          scale_x_continuous("Distance from Center (bp)") +
          scale_y_continuous("Methylation percent",limits=c(0,1)) +
          labs(fill="Assay Type") +
          theme(
            plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
            #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
            #axis.title.y = element_text(angle=0),
            legend.title=element_text(face="bold")
          )
### Barplot for CpG count in bins:
# plot2 <- ggplot(mm, aes(x=bin_mid, y=cpg_count_in_bin,fill=variable))+
#           geom_bar(stat='identity') + theme_gray() +
#           scale_x_continuous("Distance from Center (bp)") +
#           scale_y_continuous("CpG counts in bin") +
#           labs(fill="Assay") +
#           theme(
#             plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
#             #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
#             #axis.title.y = element_text(angle=0),
#             legend.title=element_text(face="bold")
#           )

### Box plot meth distribution for each CpG:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/YY1_chip_wgbs_motifs_intersect_20bp_each_cpgs.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc)), id=c("bin_mid"))

plot2 <- ggplot(mm, aes(x=as.factor(bin_mid), y=value, fill=variable)) + 
      geom_boxplot(outlier.shape = NA) + theme_gray()+
      theme(
        legend.title=element_text(face="bold"), 
        axis.title.x = element_text(size=10, hjust = 0.6, colour = "black"),
        axis.ticks.x =element_line(size=0.2, color="black")) +
        scale_x_discrete(name="Distance from Center (bp)", breaks=c(-19.5,-10.5,0.5,10.5,19.5), labels=c("-20", "-10", "0", "10", "20")) +
        scale_y_continuous(limits=c(0,1)) +
        labs(y="Methylation percent", fill="Assay Type")

grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
dev.off()



### Type 3 plot : Combining barplot and violin plot across 20bp bins:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/YY1_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))

pdf(file.path(output_dir,"YY1_cpg_methylation_plot_3.pdf"))

### Barplot for meth distribution:
plot1 <- ggplot(mm, aes(x=bin_mid, y=value, fill=variable))+geom_bar(stat='identity', position='dodge') + theme_gray() +
          ##geom_text(aes(label=cpg_count_in_bin), size=1.1, hjust=0.5) + 
          ggtitle("YY1 binding sites meth profile") +
          scale_x_continuous("Distance from Center (bp)") +
          scale_y_continuous("Methylation percent",limits=c(0,1)) +
          labs(fill="Assay Type") +
          theme(
            plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
            #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
            #axis.title.y = element_text(angle=0),
            legend.title=element_text(face="bold")
          )
### Barplot for CpG count in bins:
# plot2 <- ggplot(mm, aes(x=bin_mid, y=cpg_count_in_bin,fill=variable))+
#           geom_bar(stat='identity') + theme_gray() +
#           scale_x_continuous("Distance from Center (bp)") +
#           scale_y_continuous("CpG counts in bin") +
#           labs(fill="Assay") +
#           theme(
#             plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
#             #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
#             #axis.title.y = element_text(angle=0),
#             legend.title=element_text(face="bold")
#           )

### Box plot meth distribution for each CpG:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/YY1_chip_wgbs_motifs_intersect_20bp_each_cpgs.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc)), id=c("bin_mid"))

plot2 <- ggplot(mm, aes(x=as.factor(bin_mid), y=value, fill=variable)) + 
      geom_violin(trim=FALSE) + theme_gray()+
      theme(
        legend.title=element_text(face="bold"), 
        axis.title.x = element_text(size=10, hjust = 0.6, colour = "black"),
        axis.ticks.x =element_line(size=0.2, color="black")) +
        scale_x_discrete(name="Distance from Center (bp)", breaks=c(-19.5,-10.5,0.5,10.5,19.5), labels=c("-20", "-10", "0", "10", "20")) +
        scale_y_continuous(limits=c(0,1)) +
        labs(y="Methylation percent", fill="Assay Type")


grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
dev.off()




### ZBTB33:
###

### Type 1 plot: Barplot using facet_grid() for meth distribution and CpG counts in bin:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/ZBTB33_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]

pdf(file.path(output_dir,"ZBTB33_cpg_methylation_plot_1.pdf"))
df_1 <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))
names(df_1) <- c("bin_mid", "cpg_count_in_bin", "assay_group", "meth_perc")

# df_1 <- subset(df_1, cpg_count_in_bin >=20)
df_2 <- melt(subset(df_1, select=c(bin_mid, meth_perc, cpg_count_in_bin, assay_group)), id=c("bin_mid", "assay_group"))

meth_table_names <- c("meth_perc" = "Methylation %", "cpg_count_in_bin" = "CpG counts")
ggplot(df_2, aes(x=bin_mid, y=value, fill=assay_group))+geom_bar(stat='identity', position='dodge') + 
          theme_gray() +
          facet_grid(variable~., scales = "free_y", labeller = as_labeller(meth_table_names) ) +
          labs(x="Distance from Center (bp)", y="values", fill="Assay Type")+
          theme(strip.text = element_text(face="bold", size=9),
          strip.background = element_rect(fill="lightblue", colour="black",size=1)) 
         # scale_y_continuous("Value") + 
         # scale_x_continuous("Distance from Center (bp)", limits=c(-20, 20))

dev.off()


### Type 2 plot : Combining barplot and boxplot across 20bp bins:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/ZBTB33_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))

pdf(file.path(output_dir,"ZBTB33_cpg_methylation_plot_2.pdf"))

### Barplot for meth distribution:
plot1 <- ggplot(mm, aes(x=bin_mid, y=value, fill=variable))+geom_bar(stat='identity', position='dodge') + theme_gray() +
          ##geom_text(aes(label=cpg_count_in_bin), size=1.1, hjust=0.5) + 
          ggtitle("ZBTB33 binding sites meth profile") +
          scale_x_continuous("Distance from Center (bp)") +
          scale_y_continuous("Methylation percent",limits=c(0,1)) +
          labs(fill="Assay Type") +
          theme(
            plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
            #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
            #axis.title.y = element_text(angle=0),
            legend.title=element_text(face="bold")
          )
### Barplot for CpG count in bins:
# plot2 <- ggplot(mm, aes(x=bin_mid, y=cpg_count_in_bin,fill=variable))+
#           geom_bar(stat='identity') + theme_gray() +
#           scale_x_continuous("Distance from Center (bp)") +
#           scale_y_continuous("CpG counts in bin") +
#           labs(fill="Assay") +
#           theme(
#             plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
#             #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
#             #axis.title.y = element_text(angle=0),
#             legend.title=element_text(face="bold")
#           )

### Box plot meth distribution for each CpG:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/ZBTB33_chip_wgbs_motifs_intersect_20bp_each_cpgs.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc)), id=c("bin_mid"))

plot2 <- ggplot(mm, aes(x=as.factor(bin_mid), y=value, fill=variable)) + 
      geom_boxplot(outlier.shape = NA) + theme_gray()+
      theme(
        legend.title=element_text(face="bold"), 
        axis.title.x = element_text(size=10, hjust = 0.6, colour = "black"),
        axis.ticks.x =element_line(size=0.2, color="black")) +
        scale_x_discrete(name="Distance from Center (bp)", breaks=c(-19.5,-10.5,0.5,10.5,19.5), labels=c("-20", "-10", "0", "10", "20")) +
        scale_y_continuous(limits=c(0,1)) +
        labs(y="Methylation percent", fill="Assay Type")

grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
dev.off()



### Type 3 plot : Combining barplot and violin plot across 20bp bins:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/ZBTB33_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
meth_plot_file <- meth_plot_file[meth_plot_file$cpg_count_in_bin >= 20]
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc, cpg_count_in_bin)), id=c("bin_mid", "cpg_count_in_bin"))

pdf(file.path(output_dir,"ZBTB33_cpg_methylation_plot_3.pdf"))

### Barplot for meth distribution:
plot1 <- ggplot(mm, aes(x=bin_mid, y=value, fill=variable))+geom_bar(stat='identity', position='dodge') + theme_gray() +
          ##geom_text(aes(label=cpg_count_in_bin), size=1.1, hjust=0.5) + 
          ggtitle("ZBTB33 binding sites meth profile") +
          scale_x_continuous("Distance from Center (bp)") +
          scale_y_continuous("Methylation percent",limits=c(0,1)) +
          labs(fill="Assay Type") +
          theme(
            plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
            #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
            #axis.title.y = element_text(angle=0),
            legend.title=element_text(face="bold")
          )
### Barplot for CpG count in bins:
# plot2 <- ggplot(mm, aes(x=bin_mid, y=cpg_count_in_bin,fill=variable))+
#           geom_bar(stat='identity') + theme_gray() +
#           scale_x_continuous("Distance from Center (bp)") +
#           scale_y_continuous("CpG counts in bin") +
#           labs(fill="Assay") +
#           theme(
#             plot.title=element_text(size=14, face="bold", hjust = 0.5, colour = "black"),
#             #axis.title = element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
#             #axis.title.y = element_text(angle=0),
#             legend.title=element_text(face="bold")
#           )

### violin plot meth distribution for each CpG:
meth_plot_file = fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/files/ZBTB33_chip_wgbs_motifs_intersect_20bp_each_cpgs.bed")
meth_plot_file$bin_mid <- with(meth_plot_file, (bin_start + bin_end)/2)
mm <- melt(subset(meth_plot_file, select=c(bin_mid, wgbs_meth_perc, chip_meth_perc)), id=c("bin_mid"))

plot2 <- ggplot(mm, aes(x=as.factor(bin_mid), y=value, fill=variable)) + 
      geom_violin(trim=FALSE) + theme_gray()+
      theme(
        legend.title=element_text(face="bold"), 
        axis.title.x = element_text(size=10, hjust = 0.6, colour = "black"),
        axis.ticks.x =element_line(size=0.2, color="black")) +
        scale_x_discrete(name="Distance from Center (bp)", breaks=c(-19.5,-10.5,0.5,10.5,19.5), labels=c("-20", "-10", "0", "10", "20")) +
        scale_y_continuous(limits=c(0,1)) +
        labs(y="Methylation percent", fill="Assay Type")


grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
dev.off()












#      scale_fill_hue(limits = c("Fair","Very Good"))
### Manually change the factors level that you want it to show up using indices:
# keeps<-c(1,4,6,7) # Indices for levels you want to show
# ggplot(data=diamonds)+
#   geom_bar(aes(x=color,y=carat),stat="identity")+
#   scale_x_discrete(
#     breaks=levels(d$color)[keeps],
#     labels=table(d$color)[keeps])



# require(gridExtra)
# grid.arrange(plot1, plot2, ncol = 1, heights = c(2, 1))
# gb1 <- ggplot_build(plot1)
# gb2 <- ggplot_build(p2)
# n1 <- length(gb1$panel$ranges[[1]]$y.labels)
# n2 <- length(gb2$panel$ranges[[1]]$y.labels)
# gA <- ggplot_gtable(gb1)
# gB <- ggplot_gtable(gb2)
# g <- gtable:::rbind_gtable(gA, gB, "last")
# panels <- g$layout$t[grep("panel", g$layout$name)]
# g$heights[panels] <- list(unit(n1*5, "null"), unit(n2,"null")) # change 5 to other int
# grid.newpage()
# grid.draw(g)


# ggplot2::ggplot(meth_plot_file, aes(bin_mid, wgbs_meth_perc)) +
#   #ggplot2::geom_dot(aes(color=annotation)) +
#   ggplot2::geom_line(color="red") +
#   ggplot2::geom_point(color="red") +
#   ggplot2::geom_point(aes(bin_mid, chip_meth_perc),color="blue") +
#   ggplot2::geom_line(aes(bin_mid, chip_meth_perc),color="blue") +

#   #ggplot2::ylab("Percent Methylation") +
#   #ggplot2::xlab("Distance from Center (bp)") + 
#   labs(x="Distance from Center (bp)", y="Percent Methylation", color="Assay")+
#   ggplot2::scale_x_continuous(expand=c(0,0)) +
#   ggplot2::scale_y_continuous(expand=c(0,0.0001), limits = c(0,1)) +
#   ggplot2::theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         # legend.position="none",
#         legend.key = element_blank(),
#         axis.line = element_line(),
#         panel.background = element_blank(),
#         axis.text.x = element_text(color="black", size=12),
#         axis.title.x = element_text(color="black", size=14, vjust=-1.5),
#         axis.text.y = element_text(color="black", size=12),
#         axis.title.y = element_text(color="black", size=14, vjust=3),
#         plot.margin = grid::unit(c(1,1,1,1), "cm"),
#         panel.border=element_blank(),
#         axis.ticks=element_line(size=0.6, color="black"))




# for (each in file_list) {
#   binned_perc_meth_table <- fread(file, sep="\t", header= TRUE)
#   names(binned_perc_meth_table) <- c("bin_start", "bin_end", "meth", "unmeth", "Percent_meth", "annotation")
#   binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table, (bin_start + bin_end)/2)
#   this_plot <- ggplot2::ggplot(binned_perc_meth_table, aes(bin_mid, Percent_meth)) +
#   ggplot2::geom_line(aes(color=annotation)) +
#   ggplot2::ylab("Percent Methylation") +
#   ggplot2::xlab("Distance from Center (bp)") +
#   ggplot2::scale_x_continuous(expand=c(0,0)) +
#   ggplot2::scale_y_continuous(expand=c(0,0.0001), limits = c(0,1)) +
#   ggplot2::theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         # legend.position="none",
#         legend.key = element_blank(),
#         axis.line = element_line(),
#         panel.background = element_blank(),
#         axis.text.x = element_text(color="black", size=12),
#         axis.title.x = element_text(color="black", size=14, vjust=-1.5),
#         axis.text.y = element_text(color="black", size=12),
#         axis.title.y = element_text(color="black", size=14, vjust=3),
#         plot.margin = grid::unit(c(1,1,1,1), "cm"),
#         panel.border=element_blank(),
#         axis.ticks=element_line(size=0.6, color="black"))

# }







file_dir <- "~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect" 
file_list <- Sys.glob(file.path(file_dir, "*cpg*.txt"))

### Read the file:
meth_table <- fread("~/Dropbox/encode_3/chip_bs_wgbs_analysis/chipbs_tf_intersect/ZBTB33_cpg_grouped_chipbs_wgbs_ratio.txt", sep="\t", header= TRUE)
#meth_table[meth_table$chip_coverage >100]

meth_table$coverage_anno <- ""
meth_table[meth_table$chip_coverage < 5]$coverage_anno <-  "<5"
meth_table[meth_table$chip_coverage >= 5 & meth_table$chip_coverage < 10]$coverage_anno <- "5-10"
meth_table[meth_table$chip_coverage >= 10 & meth_table$chip_coverage < 20]$coverage_anno <- "10-20"
meth_table[meth_table$chip_coverage >= 20 & meth_table$chip_coverage < 50]$coverage_anno <- "20-50"
meth_table[meth_table$chip_coverage >= 50 & meth_table$chip_coverage < 100]$coverage_anno <- "50-100"
meth_table[meth_table$chip_coverage >= 100 ]$coveno <- ">100"

#color_manual <-  c("red", "yellow", "black", "green", "blue", "darkblue")
pdf(file.path(output_dir,"JUND_cpg_methylation_profile.pdf"))

meth_table$coverage_anno <- factor(meth_table$coverage_anno, levels=c("<5","5-10","10-20","20-50","50-100",">100"))
ggplot(meth_table, aes(wgbs_meth_perc,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
#ggplot(meth_table, aesmeth_perc,chip_wgbs_ratio_norm_log, color=coverage_anno)) + 
  geom_point(alpha = 0.5) +
  geom_density2d(colour="black", h=0.3) + 
  labs(x="WGBS methylation fraction", y="log10 (ChIPBS:WGBS normalized ratio)", color="Coverage") + 
  theme_bw() + 
  ggtitle("ZBTB33 binding sites CpG methylation distribution") + 
  theme(
    axis.text.y = element_text(size=8),
    plot.title=element_text(size=14, face="bold", hjust = 0.6, colour = "black"),
    legend.title=element_text(face="bold")
  ) +
  scale_y_continuous(limits=c(-3,3))
#guides(fill=guide_legend(title="TF Category"))
#ggMarginal(this_plot + theme_gray(), col = "darkblue")  #scale_color_manual(values=color_manual)

dev.off()

















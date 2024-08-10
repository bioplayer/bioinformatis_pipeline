library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(qqman)

setwd("E:/data_analyse/GWAS")

gwas_data <- fread("GCST90081522_buildGRCh38.tsv.gz")


chrom = 5
jade2_start <- 134525670
jade2_end <- 134583227
jade2_region_start <- jade2_start - 30000
jade2_region_end <- jade2_end + 30000
tss_region_start <- jade2_start - 3000
tss_region_end <- jade2_start + 3000

annotation <- function(snps) {
  gr <- GRanges(
    seqnames = Rle(snps$chromosome),
    ranges = IRanges(start = snps$base_pair_location, end = snps$base_pair_location)
  )
  peakAnno <- annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db")
  peakAnno_df <- as.data.frame(peakAnno)
  return(peakAnno_df)
}





# JADE2的曼哈顿图
jade2_snps <- gwas_data %>%
  filter(chromosome == 5 & base_pair_location >= jade2_region_start & base_pair_location <= jade2_region_end)
jade2_snps$chromosome <- paste0("chr", jade2_snps$chromosome)


jade2_snps <- jade2_snps %>%
  mutate(logp = -log10(p_value))


ggplot(jade2_snps, aes(x = base_pair_location, y = logp)) +
  geom_point(color = "blue") +
  theme_minimal() +
  labs(title = "Manhattan Plot for SNPs in JADE2 +/- 30kb Region",
       x = "Genomic Position",
       y = "-log10(p-value)") +
  theme(legend.position = "none") +
  geom_vline(xintercept = c(jade2_start, jade2_end), linetype = "dashed", color = "red") +
  annotate("text", x = jade2_start, y = max(jade2_snps$logp), label = "JADE2 Start", vjust = -1, color = "red") +
  annotate("text", x = jade2_end, y = max(jade2_snps$logp), label = "JADE2 End", vjust = -1, color = "red")


# 绘制5号染色体的曼哈顿图
chr5_snps <- gwas_data %>%
  filter(chromosome == 5)
chr5_snps$chromosome <- paste0("chr", chr5_snps$chromosome)
aa <- annotation(chr5_snps)
annotated_snps <- merge(chr5_snps, aa, by.x = c("chromosome", "base_pair_location"), by.y = c("seqnames", "start"), all.x = TRUE)

annotated_snps$logp <- -log10(annotated_snps$p_value)
annotated_snps$logp[is.infinite(annotated_snps$logp)] <- NA

annotated_snps$chromosome <- gsub("chr", "", annotated_snps$chromosome)
annotated_snps$chromosome[annotated_snps$chromosome == "X"] <- "23"
annotated_snps$chromosome[annotated_snps$chromosome == "Y"] <- "24"
annotated_snps$chromosome[annotated_snps$chromosome == "MT"] <- "25"
annotated_snps$chromosome <- as.numeric(annotated_snps$chromosome)

significant <- annotated_snps %>%
  filter(logp > 5)

# 对每个基因选择-log10(p)值最大的SNP
max_snps <- significant %>%
  group_by(SYMBOL) %>%
  filter(logp == max(logp)) %>%
  ungroup()


pdf("manhattan_plot_chr5.pdf", width=30, heigh=30)
manhattan(annotated_snps, chr = "chromosome", bp = "base_pair_location", snp = "SYMBOL", p = "p_value", 
          main = "Manhattan Plot for Chromosome 5")


text(max_snps$base_pair_location, max_snps$logp, labels = max_snps$SYMBOL, pos = 3, cex = 0.8, col = "red")

# 关闭PDF设备
dev.off()

###############################
# 全基因组曼哈顿图
gwas_data$chromosome <- as.character(gwas_data$chromosome)
gwas_data <- gwas_data %>%
  filter(!is.na(chromosome)) %>%
  mutate(chromosome = case_when(
    chromosome == "X" ~ "23",
    chromosome == "Y" ~ "24",
    TRUE ~ chromosome
  )) %>%
  mutate(chromosome = as.numeric(chromosome))

# 计算-log10(p)值
gwas_data <- gwas_data %>%
  mutate(logp = -log10(p_value))

# 全基因组曼哈顿图
pdf("manhattan_plot_genomewide.pdf", width = 30, height = 15)
manhattan(gwas_data, chr = "chromosome", bp = "base_pair_location", snp = "rs_id", p = "p_value",
          main = "Manhattan Plot for Genome-wide", ylim = c(0, max(gwas_data$logp, na.rm = TRUE) + 2))
dev.off()


manhattan(gwasResults, 
          main = "Manhattan Plot", #设置主标题
          ylim = c(0, 10), #设置y轴范围
          cex = 0.6, #设置点的大小
          cex.axis = 0.9, #设置坐标轴字体大小
          col = c("blue4", "orange3","red"), #设置散点的颜色
          suggestiveline = F, genomewideline = F, #remove the suggestive and genome-wide significance lines
          chrlabs = c(paste0("chr",c(1:20)),"P","Q") #设置x轴染色体标签名
)



# 确保JADE2的区间内的所有SNP存在于数据集中
jade2_start <- 134525670
jade2_end <- 134583227

jade2_snps <- annotated_snps %>%
  filter(base_pair_location >= jade2_start & base_pair_location <= jade2_end)

if (nrow(jade2_snps) > 0) {
  significant <- annotated_snps %>%
    filter(logp > 5)
  significant <- bind_rows(significant, jade2_snps)
} else {
  significant <- annotated_snps %>%
    filter(logp > 5)
}

# 对每个基因选择-log10(p)值最大的SNP
max_snps <- significant %>%
  group_by(SYMBOL) %>%
  filter(logp == max(logp)) %>%
  ungroup()

# 找出JADE2区间内-log10(p)值最大的SNP
max_jade2_snp <- jade2_snps %>%
  filter(logp == max(logp)) %>%
  slice(1)  # 如果有多个相同的最大值，选择第一个

print("start plot")
# 绘制5号染色体曼哈顿图
pdf("manhattan_plot_chr5.pdf", width = 30, height = 15)
manhattan(annotated_snps, chr = "chromosome", bp = "base_pair_location", snp = "SYMBOL", p = "p_value",
          main = "Manhattan Plot for Chromosome 5", ylim = c(0, max(annotated_snps$logp, na.rm = TRUE) + 2),
          cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3", "red"),
          suggestiveline = FALSE, genomewideline = FALSE,
          chrlabs = c("chr5"))

text(max_snps$base_pair_location, max_snps$logp, labels = max_snps$SYMBOL, pos = 3, cex = 0.8, col = "red")

# 标注JADE2区间内的-log10(p)值最大的SNP
if (nrow(max_jade2_snp) > 0) {
  text(max_jade2_snp$base_pair_location, max_jade2_snp$logp, labels = "JADE2", pos = 3, cex = 0.8, col = "blue")
}

# 关闭PDF设备
dev.off()

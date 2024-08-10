setwd("E:/data_analyse/TCGA")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(ggplot2)
library(dplyr)

library(survival)
library(survminer)
projects <- c("TCGA-BRCA", "TCGA-SARC", "TCGA-ACC", "TCGA-UCEC", "TCGA-KIRC", "TCGA-LAML", "TCGA-SKCM",
              "TCGA-PAAD", "TCGA-TGCT", "TCGA-CESC", "TCGA-ESCA", "TCGA-THCA", "TCGA-LIHC", "TCGA-PRAD",
              "TCGA-READ", "TCGA-OV", "TCGA-UVM", "TCGA-BLCA", "TCGA-CHOL", "TCGA-GBM", "TCGA-UCS", 
              "TCGA-PCPG", "TCGA-MESO", "TCGA-DLBC", "TCGA-COAD", "TCGA-STAD", "TCGA-KIRP", "TCGA-THYM",
              "TCGA-KICH", "TCGA-LGG", "TCGA-LUSC", "TCGA-LUAD","TCGA-HNSC")

projects_test <- c("TCGA-BRCA", "TCGA-SARC", "TCGA-LAML", "TCGA-TGCT")
query <- GDCquery(
  project = projects_test,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
# GDCdownload(query)
data <- GDCprepare(query)
metadata <- colData(data)
sample_type_distribution <- table(metadata$project_id, metadata$sample_type)
print(sample_type_distribution)
write.csv(sample_type_distribution, "sample_type_distribution.csv", row.names = T)
with_normal_samples <- rownames(sample_type_distribution)[rowSums(sample_type_distribution[, "Solid Tissue Normal", drop=FALSE]) > 0]

query <- GDCquery(
  project = with_normal_samples,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
data <- GDCprepare(query)

# 2. 处理log数据和元数据
log_data_list <- list()
metadata_list <- list()

for (i in 1:length(projects_test)) {
  project <- projects_test[i]
  project_indices <- which(colData(data)$project_id == project)  # 获取与当前项目对应的列索引
  project_data <- assays(data)[[1]][, project_indices]
  
  log_data <- log2(project_data + 1)
  log_data_list[[project]] <- log_data
  
  # 创建metadata项目
  metadata_project <- data.frame(
    sample = colnames(log_data),
    cancer_type = rep(project, ncol(log_data)),
    sample_type = metadata$sample_type[project_indices]  # 根据索引获取sample_type
  )
  metadata_list[[project]] <- metadata_project
}

# 合并log数据和元数据
combined_log_data <- do.call(cbind, log_data_list)
combined_metadata <- do.call(rbind, metadata_list)

# 3. 去除版本号并过滤基因ID
gene_ids <- rownames(combined_log_data)
gene_ids_no_version <- sub("\\..*", "", gene_ids)

annotations <- read.csv("basic_human.txt", header = TRUE, sep = ";")
annotations$gene_id_no_version <- sub("\\..*", "", annotations$gene_id)

gene_symbols <- annotations %>%
  filter(gene_id_no_version %in% gene_ids_no_version)

rownames(combined_log_data) <- gene_ids_no_version
combined_log_data <- combined_log_data[gene_symbols$gene_id_no_version, ]
rownames(combined_log_data) <- gene_symbols$gene_name

# 4. 提取NCOA3基因表达量并合并数据
gene_expression <- combined_log_data["NCOA3", ]
gene_expression_df <- data.frame(expression = gene_expression)
gene_expression_df$sample <- rownames(gene_expression_df)

# 合并表达数据和元数据
merged_data <- merge(gene_expression_df, combined_metadata, by = "sample")

write.csv(merged_data, "merged_data.csv", row.names = T)
write.csv(combined_log_data, "combined_log_data.csv", row.names = T)
write.csv(combined_metadata, "combined_metadata.csv", row.names = T)

# 5. 处理数据和绘制箱线图
merged_data <- merged_data %>%
  filter(sample_type %in% c("Primary Tumor", "Metastatic", "Solid Tissue Normal"))

merged_data <- merged_data %>%
  group_by(cancer_type) %>%
  mutate(median_expression = median(expression)) %>%
  ungroup() %>%
  arrange(desc(median_expression))

# 将cancer_type转换为因子并按照中位数排序
merged_data$cancer_type <- factor(merged_data$cancer_type, levels = unique(merged_data$cancer_type))

# 绘制NCOA3基因在不同癌症类型中的表达箱线图
p1 <- ggplot(merged_data, aes(x = cancer_type, y = expression, fill = sample_type)) +
  geom_boxplot() +
  labs(title = "Expression of NCOA3 across different cancer types",
       x = "Cancer Type",
       y = "Expression (log2(FPKM + 1))") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic()
ggsave("Expression.pdf", plot = p1, width = 30, height = 15, units = "in")


#############################################################################################################
# 1. 计算正常样本的表达中位数
normal_medians <- merged_data %>%
  filter(sample_type == "Solid Tissue Normal") %>%
  group_by(cancer_type) %>%
  summarise(normal_median_expression = median(expression))

# 2. 合并正常样本中位数与所有样本数据
merged_data <- merge(merged_data, normal_medians, by = "cancer_type", all.x = TRUE)

# 3. 计算差异表达，并标记高表达样本
merged_data <- merged_data %>%
  mutate(expression_diff = expression - normal_median_expression,
         high_expression = ifelse(expression_diff > 0, "High", "Low"))

# 4. 计算高表达的样本数和比例，并按照高表达比例排序
high_expression_proportion <- merged_data %>%
  group_by(cancer_type) %>%
  summarise(
    high_expression_count = sum(high_expression == "High"),
    total_count = n(),
    high_expression_ratio = high_expression_count / total_count
  ) %>%
  arrange(desc(high_expression_ratio))

# 5. 绘制高表达比例柱状图
p2 <- ggplot(high_expression_proportion, aes(x = reorder(cancer_type, -high_expression_ratio), y = total_count)) +
  geom_bar(stat = "identity", aes(fill = "Total"), show.legend = FALSE, color = "black", width = 0.7) +
  geom_bar(data = high_expression_proportion %>%
             mutate(high_proportion = total_count * high_expression_ratio),
           aes(x = reorder(cancer_type, -high_expression_ratio), y = high_proportion, fill = "High Expression"), 
           stat = "identity", show.legend = FALSE, color = "black", width = 0.7) +
  geom_text(aes(label = paste0("N=", total_count, "\n(", scales::percent(high_expression_ratio, accuracy = 1), ")")),
            vjust = -0.5, color = "black", size = 3.5) +
  labs(title = "High Expression of NCOA3 Relative to Normal Samples in Different Cancer Types",
       x = "Cancer Type",
       y = "Number of Samples") +
  scale_fill_manual(values = c("Total" = "#D3D3D3", "High Expression" = "#4682B4")) +  # 调整颜色
  theme_minimal(base_size = 15) +  # 使用更简洁的主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18))  +# 添加网格线
  theme(panel.grid=element_blank())

# 6. 保存绘图为PDF文件
ggsave("plot.pdf", plot = p2, width = 30, height = 15, units = "in")

##############################################################################################################
merged_data <- read.csv("merged_data.csv")
projects <- c("TCGA-BRCA", "TCGA-SARC", "TCGA-UCEC", "TCGA-KIRC", "TCGA-SKCM",
              "TCGA-PAAD", "TCGA-CESC", "TCGA-ESCA", "TCGA-THCA", "TCGA-LIHC", "TCGA-PRAD",
              "TCGA-READ", "TCGA-BLCA", "TCGA-CHOL", "TCGA-GBM",
              "TCGA-PCPG", "TCGA-COAD", "TCGA-STAD", "TCGA-KIRP", "TCGA-THYM",
              "TCGA-KICH", "TCGA-LUSC", "TCGA-LUAD","TCGA-HNSC")
clin_query <- GDCquery(project = projects,
                       data.category = "Clinical",
                       data.type = "Clinical Supplement",
                       data.format = "BCR XML")

#follow_up_data <- GDCprepare_clinic(clin_query, clinical.info = "follow_up")
# 假设 merged_data 已包含基因表达量数据

# 提取病人ID
follow_up_data <- follow_up_data %>%
  mutate(patient_id = substr(bcr_followup_barcode, 1, 12))

# 提取样本ID
merged_data <- merged_data %>%
  mutate(patient_id = substr(sample, 1, 12))

# 提取独特的临床数据
clinical_data <- follow_up_data %>%
  distinct(patient_id, .keep_all = TRUE)

# 提取独特的表达数据
expression_data <- merged_data %>%
  distinct(patient_id, .keep_all = TRUE)

# 合并临床和表达数据
merged_data <- expression_data %>%
  left_join(clinical_data, by = "patient_id")

# 进一步处理数据，例如去除 NA 值
final_data <- merged_data %>%
  filter(!is.na(days_to_death) & !is.na(vital_status))

# 过滤掉缺失值
filtered_data <- final_data %>%
  filter(!is.na(days_to_death) & !is.na(vital_status) & !is.na(cancer_type))

# 按癌症类型进行生存分析
for (cancer in projects) {
  cancer_data <- filtered_data %>%
    filter(cancer_type == cancer)
  
  if (nrow(cancer_data) > 0) {
    # 计算正常样本的中位表达量
    normal_samples <- cancer_data %>%
      filter(sample_type == "Solid Tissue Normal")
    
    if (nrow(normal_samples) > 0) {
      median_expression <- median(normal_samples$expression, na.rm = TRUE)
    } else {
      median_expression <- median(cancer_data$expression, na.rm = TRUE)
    }
    
    # 根据中位数将样本分为高表达组和低表达组
    cancer_data <- cancer_data %>%
      mutate(expression_group = ifelse(expression > median_expression, "High", "Low"))
    
    surv_object <- Surv(time = cancer_data$days_to_death, 
                        event = cancer_data$vital_status == "Dead")
    
    surv_fit <- survfit(surv_object ~ expression_group, data = cancer_data)
    
    p <- ggsurvplot(
      surv_fit,
      data = cancer_data,
      pval = TRUE,
      risk.table = TRUE,
      palette = "jco", # 使用适合的调色板
      title = paste("Survival Curves for NCOA3 Expression in", cancer),
      xlab = "Time (days)",
      ylab = "Survival Probability",
      legend.title = "Expression Group",
      ggtheme = theme_classic()
    )
    
    # 保存图像
    ggsave(filename = paste0("survival_curve_", cancer, ".png"), plot = p$plot)
  } else {
    message(paste("No data available for", cancer))
  }
}
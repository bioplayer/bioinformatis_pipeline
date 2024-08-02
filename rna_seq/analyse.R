
deal = function(path) {
  all_files = grep("fc.txt$", list.files(path), value = TRUE)
  
  # deal with first file
  first_file = read.table(all_files[1], skip = 1, header = TRUE)
  column_name_first = gsub("_fc.txt$", "", all_files[1])
  column_name_first = gsub("-", "_", column_name_first)
  
  # deal with other file
  merged_data = data.frame(first_file[, 1], first_file[6], first_file[, 7])
  colnames(merged_data) = c("gene", "length", column_name_first)
  for (file in all_files) {
    column_name = gsub("_fc.txt$", "", file)
    column_name = gsub("-", "_", column_name)
    data = read.table(file, skip = 1, header = TRUE)
    selected_columns = data.frame(data[, 1], data[, 7])
    colnames(selected_columns) = c("gene", column_name)
    merged_data = merge(merged_data, selected_columns, all = TRUE)
  }
  
  # get information of matrix
  gene_list = merged_data$gene
  gene_len = merged_data$length
  
  merged_data <- subset(merged_data, select = -c(gene, length))
  merged_data = as.data.frame(sapply(merged_data, function(x) as.numeric(as.character(x))))
  rownames(merged_data) = gene_list
  merged_data$gene = gene_list
  merged_data$length = gene_len
  
  return(merged_data)
}

fpkm = function(raw_matrix, count_matrix) {
  kb = raw_matrix$length/1000
  rpk = count_matrix/kb
  fpkm = t(t(rpk)/colSums(count_matrix) * 10^6)
  return(as.data.frame(fpkm))
}

diff_deseq = function(gene_list, raw_data, FC, pvalue, rep_num) {
  variable_name = sub("_1$", "", colnames(raw_data)[1])
  rownames(raw_data) = gene_list
  condition = factor(c(rep("test", rep_num), rep("control", rep_num)))
  colData = data.frame(row.names = colnames(raw_data), condition)
  dds = DESeqDataSetFromMatrix(countData = raw_data, colData = colData, design = ~ condition)
  deseq = DESeq(dds) 
  res = results(deseq)
  out = data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  out = na.omit(out)
  out$threshold <- factor(
    ifelse(out$padj < pvalue & abs(out$log2FoldChange) > FC,
           ifelse(out$log2FoldChange > FC, 'Up', 'Down'),
           'None'),
    levels = c('Up', 'Down', 'None')
  )
  message("\n", variable_name, " FC=", FC, " p-value=",  pvalue,
          "\nup gene:", sum(out$threshold == 'Up'), "\ndown gene:", sum(out$threshold == 'Down'))
  return(out)
}

Enrichment_plot = function(matrix, filename, species) {
  if (species == "mouse") {
    BM_dataset = "mmusculus_gene_ensembl"
    GO_dataset = "org.Mm.eg.db"
    KEGG_dataset = "mmu"
    message("mouse data analyse")
  } else if (species == "human") {
    BM_dataset = "hsapiens_gene_ensembl"
    GO_dataset = "org.Hs.eg.db"
    KEGG_dataset = "hsa"
    message("human data analyse")
  } else {
    stop("Error: Species is not included.")
  }
  
  if (nrow(matrix > 10)) {
    rich = data.frame(gene=matrix$gene, threshold=matrix$threshold)
    entrez = getBM(attributes = "entrezgene_id", filters = "ensembl_gene_id",
                   values = unlist(lapply(rich$gene, function(x) sub("\\..*", "", x))), 
                   mart = useMart("ensembl", BM_dataset))
    GO = enrichGO(entrez$entrezgene_id, OrgDb = GO_dataset, keyType = "ENTREZID", ont = "ALL")
    png(file = paste0(filename,'_GO.png'), width = 2500, height = 1500, units = "px", res = 300)
    print(dotplot(GO, showCategory=10, title=paste0("gene=",nrow(rich))))
    dev.off()
    KEGG = enrichKEGG(entrez$entrezgene_id,organism = KEGG_dataset)
    png(file = paste0(filename,'_KEGG.png'), width = 2500, height = 1500, units = "px", res = 300)
    print(dotplot(KEGG, showCategory=10, title=paste0("gene=",nrow(rich))))
    dev.off()
  } else {
    message("\nCannot enrich pathways because there is no gene")
  }
  
  rich_up = rich[rich$threshold == "Up", ]
  if (nrow(rich_up) > 10) {
    entrez_up = getBM(attributes = "entrezgene_id", filters = "ensembl_gene_id",
                      values = unlist(lapply(rich_up$gene, function(x) sub("\\..*", "", x))), 
                      mart = useMart("ensembl", BM_dataset))
    GO_up = enrichGO(entrez_up$entrezgene_id, OrgDb = GO_dataset, keyType = "ENTREZID", ont = "ALL")
    png(file = paste0(filename,'_up_GO.png'), width = 2500, height = 1500, units = "px", res = 300)
    print(dotplot(GO_up, showCategory=10, title=paste0("gene=",nrow(rich_up))))
    dev.off()
    KEGG_up = enrichKEGG(entrez_up$entrezgene_id,organism = KEGG_dataset)
    png(file = paste0(filename,'_up_KEGG.png'), width = 2500, height = 1500, units = "px", res = 300)
    print(dotplot(KEGG_up, showCategory=10, title=paste0("gene=",nrow(rich_up))))
    dev.off()
  } else {
    message("\nCannot enrich pathways because there is no up gene")
  }
  
  rich_down = rich[rich$threshold == "Down", ]
  if (nrow(rich_down) > 10) {
    entrez_down = getBM(attributes = "entrezgene_id", filters = "ensembl_gene_id",
                        values = unlist(lapply(rich_down$gene, function(x) sub("\\..*", "", x))), 
                        mart = useMart("ensembl", BM_dataset))
    GO_down = enrichGO(entrez_down$entrezgene_id, OrgDb = GO_dataset, keyType = "ENTREZID", ont = "ALL")
    png(file = paste0(filename,'_down_GO.png'), width = 2500, height = 1500, units = "px", res = 300)
    print(dotplot(GO_down, showCategory=10, title=paste0("gene=",nrow(rich_down))))
    dev.off()
    KEGG_down = enrichKEGG(entrez_down$entrezgene_id,organism = KEGG_dataset)
    print(png(file = paste0(filename,'_down_KEGG.png'), width = 2500, height = 1500, units = "px", res = 300))
    dotplot(KEGG_down, showCategory=10, title=paste0("gene=",nrow(rich_down)))
    dev.off()
  } else {
    message("\nCannot enrich pathways because there is no down gene")
  }
  
}

diff_analyse = function(test, count_matrix, gene_list, FC, pvalue, rep_num, species) {
  result_list = list()
  result_filter_list = list()
  
  for (i in test) {
    df_name = paste0(i,"_vs_control")
    col_names = c()
    for (j in 1:rep_num) {
      col_names = c(col_names, paste0(i, "_", j))
    }
    col_names = c(col_names, paste0("CTR", "_", 1:rep_num))
    df = data.frame(count_matrix[, col_names])
    diff = diff_deseq(gene_list, df, FC, pvalue, rep_num)
    diff_filter = diff[diff$threshold %in% c("Up", "Down"), ]
    diff_filter$gene = rownames(diff_filter)
    
    # draw
    if (nrow(diff_filter) > 10) {
      message("\nStart gene enrichment")
      Enrichment_plot(diff_filter, i, species)
    } else {
      message("\nCannot enrich pathways because there is no different gene")
    }
    result_list[[df_name]] = diff
    result_filter_list[[df_name]] = diff_filter
  }
  return(list(result_list, result_filter_list))
}

pca_plot = function(fpkm_matrix, rep_num, point_size=8, text_size=3) {
  data.pca = prcomp(t(fpkm_matrix))
  pca_data = as.data.frame(data.pca$x)
  pca_data$group = sub(paste0("_[1-", rep_num, "]$"), "", rownames(t(fpkm_matrix)))
  p = ggplot(data = pca_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = group), size = point_size) +  
    scale_color_brewer(palette = "Set1") + 
    geom_text_repel(aes(label = group), size = text_size, show.legend = FALSE) + 
    theme(legend.position = "none", panel.background = element_rect(fill = 'white', color = 'black'), 
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle("Sample Distribution Overview")
  ggsave("pca_plot.pdf", plot = p, width = 10, height = 10)
}

library(DESeq2)

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(VennDiagram)
library(pheatmap)
library(ggpubr)


library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GenomicFeatures)



##############

FC = 0.5
pvalue = 0.05
rep_num = 2
species = "mouse" #or human

path = getwd()
raw_data = deal(path)
count_matrix = subset(raw_data, select = -c(gene, length))

count_matrix = subset(count_matrix, rowSums(count_matrix != 0) > 0)
fpkm_matrix = fpkm(raw_data, count_matrix)
gene_list = rownames(count_matrix)

pca_plot(log2(fpkm_matrix+1), rep_num)

colnames(count_matrix)  
test = c("ASA_0.5mM", "ASA_1mM", "CSA_0.5mM", "CSA_1mM", "NaAc_0.5mM", "NaAc_1mM", "NaCr_0.5mM", "NaCr_1mM")
control = c("CTR")

result = diff_analyse(test, count_matrix, gene_list, FC, pvalue, rep_num, species)



heatmap_gene = c()


XXX_1





gene_list = getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                filters = "ensembl_gene_id",
                values = unlist(lapply(rich$V1, function(x) sub("\\..*", "", x))),
                mart = useMart("ensembl", "mmusculus_gene_ensembl"))
gene_list = gene_list$external_gene_name

gene_list = getBM(attributes = c("entrezgene_id"),
                  filters = "ensembl_gene_id",
                  values = unlist(lapply(rich$V1, function(x) sub("\\..*", "", x))),
                  mart = useMart("ensembl", "mmusculus_gene_ensembl"))



##############################

venn.plot <- venn.diagram(
  x = list(
    CSA_1 = A0.5$gene,
    ASA_0.5= C1$gene
  ),
  filename = "out.png",
  lwd = 0,
  fill = c("cornflowerblue", "darkorchid1"),
  alpha = 0.75,
  label.col = "black",
  output_file = NULL,
  scaled = F
)

venn.plot <- venn.diagram(
  x = list(
    Nacr_1 = C1$gene,
    CSA_1 = NC1$gene
  ),
  filename = "out.png",
  lwd = 0,
  fill = c("cornflowerblue", "darkorchid1"),
  alpha = 0.75,
  label.col = "black",
  output_file = NULL,
  scaled = F
)

venn.plot <- venn.diagram(
  x = list(
    Nacr_1 = C1$gene,
    CSA_1 = NC1$gene
  ),
  filename = "out.png",
  lwd = 0,
  fill = c("cornflowerblue", "darkorchid1"),
  alpha = 0.75,
  label.col = "black",
  output_file = NULL,
  scaled = F
)
####################

heatmap_gene = unique(c(A0.5$gene, C0.5$gene, C1$gene, NC1$gene))
origin = fpkm_matrix
origin$gene = rownames(origin)
origin <- origin[origin$gene %in% heatmap_gene, ]
heatmap_data = data.frame(control_1=origin$CTR_1, control_2=origin$CTR_2,
                CSA_0.5_1=origin$CSA_0.5mM_1, CSA_0.5_2=origin$CSA_0.5mM_2,
                CSA_1_1=origin$CSA_1mM_1, CSA_1_2=origin$CSA_1mM_2,
                ASA_0.5_1=origin$ASA_0.5mM_1, ASA_0.5_2=origin$ASA_0.5mM_2,
                NaCr_1_1=origin$NaCr_1mM_1, NaCr_1_2=origin$NaAc_1mM_2
                )
rownames(heatmap_data) = origin$gene
pheatmap(heatmap_data, cluster_rows = FALSE, show_rownames = FALSE, show_colnames = TRUE, scale = "row")

out = log2(heatmap_data+1)
annotation_col = data.frame(class = c("Control", "Control", "CSA_0.5", "CSA_0.5", "CSA_1", "CSA_1", "ASA_0.5", "ASA_0.5", 
                                      "NaCr_1", "NaCr_1"))
row.names(annotation_col) <- colnames(heatmap_data)
pheatmap(out, cluster_rows = TRUE, show_rownames = FALSE, cluster_cols = FALSE, 
         scale = "row", treeheight_row = 0, annotation_col = annotation_col,annotation_legend = FALSE,
         main = "Heatmap of sample")


###################
N_plus_C = intersect(C1$gene, NC1$gene)
A_plus_C = intersect(A0.5$gene, C1$gene)

new_dge = as.data.frame(dge)
new_dge$gene = rownames(dge)
select_dge <- new_dge[new_dge$gene %in% N_plus_C, ]
select_dge = log2(select_dge[,-ncol(select_dge)]+1)
ggscatter(select_dge, x = "NaCr_1mM_1", y = "CSA_1mM_1", 
          add = "reg.line",conf.int = TRUE, 
          fill = "lightgray",
          add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
          cor.coef = T,
          cor.coeff.args = list(),
          cor.method = "pearson")
select_dge$NaCr_1mM_1
# ==================== Mantel Test相关性热图 ====================

# 清理工作环境
rm(list = ls())

# 安装并加载相关包
p_list <- c("ggplot2", "dplyr", "vegan", "RColorBrewer", "reshape2")
for(p in p_list){
  if (!requireNamespace(p, quietly = TRUE)){
    install.packages(p)
  }
  library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
}

# 安装并加载linkET
if(!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if(!requireNamespace("linkET", quietly = TRUE)){
  devtools::install_github("Hy4m/linkET", force = TRUE)
}
library(linkET)

# 设置工作目录
setwd("/Users/catherine/Downloads/董西洋铁-砷-汞等金属循环文献-冷泉/角鲨烯基于 contig 寻找的功能基因/")

cat("=== Mantel Test相关性热图分析 ===\n\n")

# ==================== 1. 数据准备 ====================

# 1.1 读取基因丰度数据（作为物种数据）
cat("1. 读取基因丰度数据...\n")
gene_abundance <- read.csv("样品与角鲨烯功能基因丰度 1.csv", 
                           header = TRUE, 
                           stringsAsFactors = FALSE,
                           check.names = FALSE)

# 设置行名并提取基因数据
rownames(gene_abundance) <- gene_abundance[, 1]
gene_data <- gene_abundance[, -1]

# 转换为数值型
for (col in colnames(gene_data)) {
  gene_data[[col]] <- as.numeric(as.character(gene_data[[col]]))
}

# 转置数据（行=样本，列=基因）
gene_data_t <- t(gene_data)

cat("基因数据维度:", dim(gene_data_t), " (样本×基因)\n")
print(head(gene_data_t[, 1:3]))

# 1.2 创建示例环境因子数据
# 在实际分析中，您应该有真实的环境因子数据
cat("\n2. 创建示例环境因子数据...\n")

# 创建示例环境因子
set.seed(123)
sample_names <- rownames(gene_data_t)
env_factors <- data.frame(
  row.names = sample_names,
  # 创建几个示例环境因子
  Temperature = rnorm(length(sample_names), mean = 20, sd = 5),
  pH = rnorm(length(sample_names), mean = 7, sd = 0.5),
  Salinity = rnorm(length(sample_names), mean = 30, sd = 5),
  Depth = runif(length(sample_names), 0, 100),
  Organic_Matter = runif(length(sample_names), 0.1, 5),
  Nitrate = runif(length(sample_names), 0.1, 10),
  Phosphate = runif(length(sample_names), 0.01, 1),
  Oxygen = runif(length(sample_names), 2, 10)
)

cat("环境因子数据维度:", dim(env_factors), "\n")
print(head(env_factors))

# ==================== 2. 执行Mantel检验 ====================

cat("\n3. 执行Mantel检验...\n")

# 对基因进行分组（根据角鲨烯合成途径）
gene_groups <- list(
  "Squalene_Synthases" = c("Squalene_synthase", "FPP_synthase_euk", "FPP_synthase_prok"),
  "Squalene_Modification" = c("Squalene_monooxygenase", "OSC_lanosterol_cyclase"),
  "HMG-CoA_Pathway" = c("HMG-CoA_synthase", "HMG-CoA_reductase"),
  "MEP_Pathway" = c("IspG_MEP_pathway")
)

# 筛选实际存在的基因
gene_groups_filtered <- list()
for (group_name in names(gene_groups)) {
  existing_genes <- gene_groups[[group_name]][gene_groups[[group_name]] %in% colnames(gene_data_t)]
  if (length(existing_genes) > 0) {
    gene_groups_filtered[[group_name]] <- existing_genes
  }
}

# 获取基因在转置数据中的列索引
spec_select_list <- list()
for (group_name in names(gene_groups_filtered)) {
  group_genes <- gene_groups_filtered[[group_name]]
  col_indices <- which(colnames(gene_data_t) %in% group_genes)
  if (length(col_indices) > 0) {
    spec_select_list[[group_name]] <- col_indices
  }
}

cat("基因分组:\n")
for (group_name in names(spec_select_list)) {
  cat("  ", group_name, ": ", length(spec_select_list[[group_name]]), "个基因\n", sep = "")
}

# 执行Mantel检验
mantel_result <- mantel_test(gene_data_t, env_factors,
                             spec_select = spec_select_list,
                             mantel_fun = 'mantel',
                             spec_dist_method = "bray",
                             env_dist_method = "euclidean")

cat("\nMantel检验完成，共", nrow(mantel_result), "个检验\n")

# 对Mantel的r和p值重新赋值（设置绘图标签）
mantel_plot_data <- mantel_result %>%
  mutate(
    r_cat = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf),
                labels = c("<0.25", "0.25-0.5", ">=0.5")),
    p_cat = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">=0.05"))
  )

write.csv(mantel_plot_data, "mantel_test_results.csv", row.names = FALSE)
cat("Mantel检验结果已保存到: mantel_test_results.csv\n")

# ==================== 3. 计算环境因子相关性 ====================

cat("\n4. 计算环境因子相关性...\n")

# 计算环境因子之间的Pearson相关性
env_corr <- correlate(env_factors)

# 保存相关性结果
env_corr_table <- env_corr %>% as_md_tbl()
write.csv(env_corr_table, "environment_factors_correlation.csv", row.names = TRUE)
cat("环境因子相关性矩阵已保存到: environment_factors_correlation.csv\n")

# ==================== 4. 绘制Mantel Test相关性热图 ====================

cat("\n5. 绘制Mantel Test相关性热图...\n")

# 4.1 模板1：热图在左下角，连线在右上角
cat("\n绘制模板1: 热图在左下角...\n")

# 自定义颜色
color_pal <- function(n) {
  if (n == 3) return(c("#4DAF4A", "#377EB8", "#E41A1C"))
  if (n == 4) return(c("#4DAF4A", "#377EB8", "#984EA3", "#E41A1C"))
  return(brewer.pal(n, "Set2"))
}

p1 <- qcorrplot(env_corr, 
                type = "lower", 
                diag = FALSE) +
  geom_square() +
  geom_mark(only_mark = TRUE) +
  geom_couple(
    aes(colour = p_cat, size = r_cat),
    data = mantel_plot_data,
    curvature = 0.1
  ) +
  geom_diag_label(
    mapping = aes(y = .y + 0.05),
    hjust = 0.15
  ) +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(4)) +
  guides(
    size = guide_legend(title = "Mantel's r",
                        override.aes = list(colour = "grey35"),
                        order = 2),
    colour = guide_legend(title = "Mantel's p",
                          override.aes = list(size = 3),
                          order = 1),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  ) +
  theme(
    axis.text.y = element_blank()
  ) +
  labs(title = "角鲨烯基因-环境因子Mantel Test分析",
       subtitle = "热图: 环境因子间相关性 | 连线: 基因-环境因子关联")

ggsave("mantel_heatmap_template1.pdf", p1, width = 10, height = 8)
ggsave("mantel_heatmap_template1.png", p1, width = 10, height = 8, dpi = 300)
cat("模板1已保存: mantel_heatmap_template1.pdf/png\n")

# 4.2 模板2：热图在右上角，连线在左下角
cat("\n绘制模板2: 热图在右上角...\n")

p2 <- qcorrplot(env_corr, 
                type = "upper", 
                diag = FALSE) +
  geom_square() +
  geom_mark(only_mark = TRUE) +
  geom_couple(
    aes(colour = p_cat, size = r_cat),
    data = mantel_plot_data,
    curvature = -0.1
  ) +
  geom_diag_label(
    mapping = aes(y = .y + 0.05),
    hjust = 0.05
  ) +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(4)) +
  guides(
    size = guide_legend(title = "Mantel's r",
                        override.aes = list(colour = "grey35"),
                        order = 2),
    colour = guide_legend(title = "Mantel's p",
                          override.aes = list(size = 3),
                          order = 1),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  ) +
  theme(
    axis.text.y = element_blank()
  ) +
  labs(title = "角鲨烯基因-环境因子Mantel Test分析",
       subtitle = "热图在右上角 | 连线在左下角")

ggsave("mantel_heatmap_template2.pdf", p2, width = 10, height = 8)
ggsave("mantel_heatmap_template2.png", p2, width = 10, height = 8, dpi = 300)
cat("模板2已保存: mantel_heatmap_template2.pdf/png\n")

# 4.3 模板3：圆形图
cat("\n绘制模板3: 圆形图...\n")

p3 <- qcorrplot(env_corr, 
                type = "upper", 
                diag = FALSE) +
  geom_circle2() +
  geom_mark(only_mark = TRUE) +
  geom_couple(
    aes(colour = p_cat, size = r_cat),
    data = mantel_plot_data,
    curvature = -0.1
  ) +
  geom_diag_label(
    mapping = aes(y = .y + 0.05),
    hjust = 0.05
  ) +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(4)) +
  guides(
    size = guide_legend(title = "Mantel's r",
                        override.aes = list(colour = "grey35"),
                        order = 2),
    colour = guide_legend(title = "Mantel's p",
                          override.aes = list(size = 3),
                          order = 1),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  ) +
  theme(
    axis.text.y = element_blank()
  ) +
  labs(title = "角鲨烯基因-环境因子Mantel Test分析（圆形图）")

ggsave("mantel_heatmap_template3.pdf", p3, width = 10, height = 8)
cat("模板3已保存: mantel_heatmap_template3.pdf\n")

# 4.4 模板4：饼图
cat("\n绘制模板4: 饼图...\n")

p4 <- qcorrplot(env_corr, 
                type = "upper", 
                diag = FALSE) +
  geom_pie2() +
  geom_mark(only_mark = TRUE) +
  geom_couple(
    aes(colour = p_cat, size = r_cat),
    data = mantel_plot_data,
    curvature = -0.1
  ) +
  geom_diag_label(
    mapping = aes(y = .y + 0.05),
    hjust = 0.05
  ) +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(4)) +
  guides(
    size = guide_legend(title = "Mantel's r",
                        override.aes = list(colour = "grey35"),
                        order = 2),
    colour = guide_legend(title = "Mantel's p",
                          override.aes = list(size = 3),
                          order = 1),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  ) +
  theme(
    axis.text.y = element_blank()
  ) +
  labs(title = "角鲨烯基因-环境因子Mantel Test分析（饼图）")

ggsave("mantel_heatmap_template4.pdf", p4, width = 10, height = 8)
cat("模板4已保存: mantel_heatmap_template4.pdf\n")

# ==================== 5. 创建基因-菌门版本的Mantel热图 ====================

cat("\n6. 创建基因-菌门版本的Mantel热图...\n")

# 5.1 准备菌门数据（简化示例）
# 在实际分析中，您需要从分类文件中提取真实的菌门数据
set.seed(123)
phyla <- c("Proteobacteria", "Actinobacteria", "Bacteroidetes", 
           "Firmicutes", "Acidobacteria", "Planctomycetes")

# 创建菌门丰度矩阵（示例数据）
phylum_matrix <- matrix(runif(length(sample_names) * length(phyla), 0, 100),
                        nrow = length(sample_names),
                        ncol = length(phyla),
                        dimnames = list(sample_names, phyla))

phylum_df <- as.data.frame(phylum_matrix)

# 5.2 执行Mantel检验（基因-菌门）
cat("执行基因-菌门Mantel检验...\n")

mantel_phylum_result <- mantel_test(gene_data_t, phylum_df,
                                    spec_select = spec_select_list,
                                    mantel_fun = 'mantel',
                                    spec_dist_method = "bray",
                                    env_dist_method = "bray")

mantel_phylum_plot <- mantel_phylum_result %>%
  mutate(
    r_cat = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf),
                labels = c("<0.25", "0.25-0.5", ">=0.5")),
    p_cat = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">=0.05"))
  )

# 5.3 绘制基因-菌门Mantel热图
p_phylum <- qcorrplot(correlate(phylum_df), 
                      type = "lower", 
                      diag = FALSE) +
  geom_square() +
  geom_mark(only_mark = TRUE) +
  geom_couple(
    aes(colour = p_cat, size = r_cat),
    data = mantel_phylum_plot,
    curvature = 0.1
  ) +
  geom_diag_label(
    mapping = aes(y = .y + 0.05),
    hjust = 0.15
  ) +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(4)) +
  guides(
    size = guide_legend(title = "Mantel's r",
                        override.aes = list(colour = "grey35"),
                        order = 2),
    colour = guide_legend(title = "Mantel's p",
                          override.aes = list(size = 3),
                          order = 1),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  ) +
  theme(
    axis.text.y = element_blank()
  ) +
  labs(title = "角鲨烯基因-菌门Mantel Test分析",
       subtitle = "热图: 菌门间相关性 | 连线: 基因-菌门关联")

ggsave("gene_phylum_mantel_heatmap.pdf", p_phylum, width = 10, height = 8)
ggsave("gene_phylum_mantel_heatmap.png", p_phylum, width = 10, height = 8, dpi = 300)
cat("基因-菌门Mantel热图已保存: gene_phylum_mantel_heatmap.pdf/png\n")

# ==================== 6. 生成分析报告 ====================

cat("\n7. 生成分析报告...\n")

sink("mantel_test_analysis_report.txt")

cat(rep("=", 70), "\n", sep = "")
cat("角鲨烯基因Mantel Test分析报告\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("分析时间:", format(Sys.time(), "%Y年%m月%d日 %H:%M:%S"), "\n\n")

cat("一、分析概述\n")
cat(rep("-", 70), "\n", sep = "")
cat("分析类型: Mantel Test相关性热图分析\n")
cat("基因数量:", ncol(gene_data_t), "\n")
cat("样本数量:", nrow(gene_data_t), "\n")
cat("环境因子数量:", ncol(env_factors), "\n")
cat("菌门数量:", ncol(phylum_df), "\n")
cat("基因分组:", length(spec_select_list), "组\n\n")

cat("二、基因分组信息\n")
cat(rep("-", 70), "\n", sep = "")
for (group_name in names(spec_select_list)) {
  genes_in_group <- colnames(gene_data_t)[spec_select_list[[group_name]]]
  cat(group_name, " (", length(genes_in_group), "个基因):\n", sep = "")
  cat("  ", paste(genes_in_group, collapse = ", "), "\n\n", sep = "")
}

cat("三、Mantel检验结果摘要\n")
cat(rep("-", 70), "\n", sep = "")
cat("基因-环境因子Mantel检验:\n")
cat("  总检验数:", nrow(mantel_plot_data), "\n")

# 统计显著性
sig_counts <- mantel_plot_data %>%
  group_by(p_cat) %>%
  summarise(Count = n())

for (i in 1:nrow(sig_counts)) {
  cat("  ", as.character(sig_counts$p_cat[i]), ": ", sig_counts$Count[i], "个检验\n", sep = "")
}

cat("\n基因-菌门Mantel检验:\n")
cat("  总检验数:", nrow(mantel_phylum_plot), "\n")

sig_counts_phylum <- mantel_phylum_plot %>%
  group_by(p_cat) %>%
  summarise(Count = n())

for (i in 1:nrow(sig_counts_phylum)) {
  cat("  ", as.character(sig_counts_phylum$p_cat[i]), ": ", 
      sig_counts_phylum$Count[i], "个检验\n", sep = "")
}

cat("\n四、最强的Mantel关联\n")
cat(rep("-", 70), "\n", sep = "")
cat("基因-环境因子 (前5强):\n")
top_env <- mantel_plot_data %>%
  arrange(desc(r)) %>%
  head(5)

for (i in 1:nrow(top_env)) {
  cat("  ", top_env$spec[i], " - 环境因子: r = ", round(top_env$r[i], 3),
      ", p = ", top_env$p[i], "\n", sep = "")
}

cat("\n基因-菌门 (前5强):\n")
top_phylum <- mantel_phylum_plot %>%
  arrange(desc(r)) %>%
  head(5)

for (i in 1:nrow(top_phylum)) {
  cat("  ", top_phylum$spec[i], " - ", top_phylum$env[i], 
      ": r = ", round(top_phylum$r[i], 3),
      ", p = ", top_phylum$p[i], "\n", sep = "")
}

cat("\n五、生成文件列表\n")
cat(rep("-", 70), "\n", sep = "")
cat("1. 主要图形文件:\n")
cat("   - mantel_heatmap_template1.pdf/png: 模板1 (热图左下角)\n")
cat("   - mantel_heatmap_template2.pdf/png: 模板2 (热图右上角)\n")
cat("   - mantel_heatmap_template3.pdf: 模板3 (圆形图)\n")
cat("   - mantel_heatmap_template4.pdf: 模板4 (饼图)\n")
cat("   - gene_phylum_mantel_heatmap.pdf/png: 基因-菌门版本\n\n")

cat("2. 数据文件:\n")
cat("   - mantel_test_results.csv: Mantel检验详细结果\n")
cat("   - environment_factors_correlation.csv: 环境因子相关性\n\n")

cat("3. 报告文件:\n")
cat("   - mantel_test_analysis_report.txt: 本报告\n\n")

cat("六、图形说明\n")
cat(rep("-", 70), "\n", sep = "")
cat("1. Mantel Test相关性热图包含两部分:\n")
cat("   - 左下/右上角的热图: 显示环境因子/菌门之间的Pearson相关性\n")
cat("   - 连线: 连接基因分组与环境因子/菌门，显示Mantel检验结果\n\n")

cat("2. 连线特征:\n")
cat("   - 线条粗细: 表示Mantel's r值大小\n")
cat("   - 线条颜色: 表示Mantel's p值显著性\n")
cat("   - 连线弯曲: 美观布局，无实际意义\n\n")

cat("3. 热图特征:\n")
cat("   - 颜色: 红色表示正相关，蓝色表示负相关\n")
cat("   - 显著性标记: *表示p<0.05, **表示p<0.01, ***表示p<0.001\n\n")

cat("七、分析结论\n")
cat(rep("-", 70), "\n", sep = "")
cat("1. 成功创建了标准的Mantel Test相关性热图\n")
cat("2. 提供了多种图形模板供选择\n")
cat("3. 分析了基因与环境因子、菌门之间的关联\n")
cat("4. 所有结果已保存为可编辑的PDF和CSV文件\n")

sink()

cat("\n8. 分析完成！\n")
cat("请查看以下主要文件:\n")
cat("1. mantel_heatmap_template1.pdf - 标准Mantel热图 (热图在左下角)\n")
cat("2. mantel_heatmap_template2.pdf - 热图在右上角的版本\n")
cat("3. gene_phylum_mantel_heatmap.pdf - 基因-菌门版本\n")
cat("4. mantel_test_analysis_report.txt - 详细分析报告\n\n")

cat("提示: 您可以使用这些模板，将示例数据替换为您的真实菌门数据。\n")
cat("只需将代码中的示例菌门数据替换为从'KO-MAG所有 BIN 所持有的功能基因 3.csv'\n")
cat("和'MAG 所有的 bin 的门纲目科属分类注释 2.csv'中提取的真实菌门丰度数据即可。\n")
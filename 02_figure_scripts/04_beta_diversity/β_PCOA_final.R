# ===== β多样性分析 - 修复维度匹配问题 ===============================
setwd("/Users/catherine/Downloads/metaphlan 作图/class 级别作图-柱状图")
# 加载包
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel) # 用于避免标签重叠

# 1. 转换文件编码
if(file.exists("class34_utf8.txt")) file.remove("class34_utf8.txt")
if(file.exists("class_group_utf8.txt")) file.remove("class_group_utf8.txt")
system("iconv -f UTF-16LE -t UTF-8 class34.txt > class34_utf8.txt")
system("iconv -f UTF-16LE -t UTF-8 class_group.txt > class_group_utf8.txt")

# 2. 读取转换后的文件
tax_tab <- read.delim("class34_utf8.txt", check.names = FALSE, stringsAsFactors = FALSE)
group_tab <- read.delim("class_group_utf8.txt", check.names = FALSE, stringsAsFactors = FALSE)
cat("物种表维度:", dim(tax_tab), "\n")
cat("分组表维度:", dim(group_tab), "\n")

# 3. 修复样本名问题 - 关键修复！
cat("\n=== 修复样本名问题 ===\n")
# 方法：手动处理数据转换，避免R自动修改列名
taxon_names <- tax_tab[[1]]  # 物种名（第一列）
# 提取样本名（列名），并修复R自动转换的点号问题
sample_names <- names(tax_tab)[-1]  # 排除第一列（物种名）
# 检查并移除"sum"样本
if("sum" %in% sample_names) {
  sample_names <- sample_names[sample_names != "sum"]
  tax_tab <- tax_tab[, c(1, which(names(tax_tab) %in% sample_names)), drop = FALSE]
}
cat("原始样本名数量:", length(sample_names), "\n")
cat("原始样本名 (前10个):", head(sample_names, 10), "\n")

# 提取数值数据
numeric_data <- tax_tab[, -1, drop = FALSE]
cat("数值数据维度:", dim(numeric_data), "\n")

# 确保数值转换
numeric_data <- as.data.frame(lapply(numeric_data, function(x) {
  x <- as.numeric(as.character(x))
  x[is.na(x)] <- 0
  return(x)
}))
cat("转换后数值数据维度:", dim(numeric_data), "\n")

# 手动创建矩阵 - 修复维度问题
tax_matrix <- as.matrix(t(numeric_data))  # 转置：样本在行，物种在列
# 设置正确的行名和列名
rownames(tax_matrix) <- sample_names  # 样本在行
colnames(tax_matrix) <- taxon_names   # 物种在列
cat("转换后矩阵维度:", dim(tax_matrix), "\n")
cat("样本数:", nrow(tax_matrix), "物种数:", ncol(tax_matrix), "\n")
cat("矩阵行名 (样本名, 前10个):", head(rownames(tax_matrix), 10), "\n")
cat("矩阵列名 (物种名, 前5个):", head(colnames(tax_matrix), 5), "\n")

# 4. β多样性分析
cat("\n开始β多样性分析...\n")
set.seed(123)
# 数据清洗
tax_matrix_clean <- tax_matrix[rowSums(tax_matrix, na.rm = TRUE) > 0, ]
tax_matrix_clean <- tax_matrix_clean[, colSums(tax_matrix_clean, na.rm = TRUE) > 0]
cat("清洗后矩阵维度:", dim(tax_matrix_clean), "\n")

# 计算Bray-Curtis距离
cat("计算Bray-Curtis距离...\n")
bc_dist <- vegdist(tax_matrix_clean, method = "bray", na.rm = FALSE)
# 检查距离矩阵的样本名
cat("距离矩阵样本名 (前5个):", head(labels(bc_dist), 5), "\n")

# NMDS分析
cat("进行NMDS排序...\n")
nmds <- metaMDS(bc_dist, k = 2, trymax = 100, autotransform = FALSE)
stress <- round(nmds$stress, 3)
cat("NMDS stress值:", stress, "\n")

# 检查NMDS样本名
cat("NMDS样本名 (前5个):", head(rownames(nmds$points), 5), "\n")

# 5. 提取坐标并合并分组
cat("处理坐标和分组数据...\n")
# 确保分组表的列名正确
names(group_tab) <- c("sample-id", "group")
# 现在样本名应该匹配了！
nmds_coord <- as.data.frame(nmds$points) %>%
  setNames(c("NMDS1", "NMDS2")) %>%
  tibble::rownames_to_column("sample") %>%
  left_join(group_tab, by = c("sample" = "sample-id")) %>%
  filter(!is.na(group)) %>%
  mutate(group = factor(group))
cat("有效样本数:", nrow(nmds_coord), "\n")
cat("分组情况:\n")
print(table(nmds_coord$group))

# 6. PERMANOVA分析
if(length(unique(nmds_coord$group)) >= 2) {
  cat("进行PERMANOVA分析...\n")
  adonis_mod <- adonis2(bc_dist ~ group, data = nmds_coord, method = "bray")
  p_manova  <- adonis_mod$Pr[1]
  r2_manova <- adonis_mod$R2[1]
  cat("PERMANOVA - p值:", p_manova, "R²:", r2_manova, "\n")
} else {
  cat("分组数不足，跳过PERMANOVA分析\n")
  p_manova <- NA
  r2_manova <- NA
}

# 7. PCOA分析 ===================================================
cat("\n=== 开始PCOA分析 ===\n")
# 进行PCOA分析
pcoa_result <- cmdscale(bc_dist, k = 2, eig = TRUE, add = TRUE)

# 计算每个主坐标的解释度
eigenvalues <- pcoa_result$eig
variance_explained <- round(eigenvalues[1:2] / sum(eigenvalues[eigenvalues > 0]) * 100, 2)
cat("PCOA轴1解释度:", variance_explained[1], "%\n")
cat("PCOA轴2解释度:", variance_explained[2], "%\n")

# 提取PCOA坐标
pcoa_coord <- as.data.frame(pcoa_result$points) %>%
  setNames(c("PCOA1", "PCOA2")) %>%
  tibble::rownames_to_column("sample") %>%
  left_join(group_tab, by = c("sample" = "sample-id")) %>%
  filter(!is.na(group)) %>%
  mutate(group = factor(group))

# 8. 绘图 =======================================================
cat("生成图形...\n")
# 颜色设置
group_col <- c(IS = "#00FF7F", AS = "#FF1493", ES = "#FFA500", NS = "#00FFFF")

# 创建统计标签
if(!is.na(p_manova)) {
  stat_lab <- sprintf("PERMANOVA p = %.3g\nR² = %.3f", p_manova, r2_manova)
  pcoa_stat_lab <- sprintf("PERMANOVA p = %.3g\nR² = %.3f\nPCOA1: %.1f%%\nPCOA2: %.1f%%", 
                           p_manova, r2_manova, variance_explained[1], variance_explained[2])
} else {
  stat_lab <- sprintf("stress = %.3f", stress)
  pcoa_stat_lab <- sprintf("PCOA1: %.1f%%\nPCOA2: %.1f%%", variance_explained[1], variance_explained[2])
}

# NMDS1箱线图
p_nmds1 <- ggplot(nmds_coord, aes(x = group, y = NMDS1, fill = group)) +
  geom_boxplot(width = 0.5, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = group_col) +
  labs(title = "NMDS1 Coordinate (Bray-Curtis)",
       x = "Group", y = "NMDS1") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# 添加统计标注
p_nmds1 <- p_nmds1 + 
  annotate("text", x = Inf, y = Inf, hjust = 1.2, vjust = 1.5,
           label = stat_lab, size = 4, colour = "black")

# 保存图形
ggsave("beta-NMDS1-box-Bayes.pdf", p_nmds1, width = 8, height = 6)
cat("图形已保存为: beta-NMDS1-box-Bayes.pdf\n")

# NMDS散点图
p_nmds_scatter <- ggplot(nmds_coord, aes(x = NMDS1, y = NMDS2, color = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = group_col) +
  stat_ellipse(level = 0.68) +
  labs(title = "NMDS Plot (Bray-Curtis)",
       x = paste("NMDS1 (stress =", stress, ")"),
       y = "NMDS2",
       color = "Group") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("beta-NMDS-scatter.pdf", p_nmds_scatter, width = 10, height = 8)
cat("散点图已保存为: beta-NMDS-scatter.pdf\n")

# PCOA散点图 ===================================================
p_pcoa <- ggplot(pcoa_coord, aes(x = PCOA1, y = PCOA2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = group_col) +
  stat_ellipse(level = 0.68, linetype = 2, size = 0.8) +
  labs(title = "PCOA Plot (Bray-Curtis)",
       x = paste("PCOA1 (", variance_explained[1], "%)"),
       y = paste("PCOA2 (", variance_explained[2], "%)"),
       color = "Group") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right") +
  # 添加统计信息
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = pcoa_stat_lab, size = 4, colour = "black")

# 保存PCOA图形
ggsave("beta-PCOA-scatter.pdf", p_pcoa, width = 10, height = 8)
cat("PCOA散点图已保存为: beta-PCOA-scatter.pdf\n")

# 9. 完成信息
cat("\n=== 分析完成！ ===\n")
cat("NMDS stress值:", stress, "(<0.2表示良好拟合)\n")
if(!is.na(p_manova)) {
  cat("PERMANOVA p值:", p_manova, "\n")
  cat("PERMANOVA R²:", r2_manova, "\n")
}
cat("PCOA解释度 - 轴1:", variance_explained[1], "%, 轴2:", variance_explained[2], "%\n")
cat("有效样本数:", nrow(nmds_coord), "\n")
cat("分组分布:\n")
print(table(nmds_coord$group))
cat("生成的PDF文件已保存到工作目录\n")

# 显示图形预览
print(p_nmds1)
print(p_nmds_scatter)
print(p_pcoa)

cat("\n🎉 恭喜！β多样性分析（包含PCOA）成功完成！\n")


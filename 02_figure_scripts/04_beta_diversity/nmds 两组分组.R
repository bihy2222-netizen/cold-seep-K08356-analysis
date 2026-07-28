# ============================================
# 微生物群落NMDS和PCoA分析脚本
# 简洁版：无灰色网格线，黑色边框，正方形，标注统计值
# ============================================

# 加载必要的包
library(vegan)
library(ggplot2)
library(patchwork)

# 设置工作目录
setwd("/Users/catherine/Downloads/NMDS/")

# ============================================
# 1. 数据读取与处理
# ============================================
cat("=== 读取与处理数据 ===\n")

# 清理文件并读取数据
clean_csv_file <- function(input_file) {
  # 读取原始文件
  lines <- readLines(input_file, warn = FALSE)
  lines_clean <- gsub("\x00", "", lines)
  lines_clean <- lines_clean[lines_clean != ""]
  
  temp_file <- "temp_species_cleaned.csv"
  writeLines(lines_clean, temp_file)
  
  # 读取清理后的文件
  species_raw <- read.csv(temp_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  
  # 处理重复行名
  species_names <- species_raw[,1]
  if(any(duplicated(species_names))) {
    species_names <- make.unique(species_names, sep = "_dup")
    species_raw[,1] <- species_names
  }
  
  # 设置行名
  rownames(species_raw) <- species_raw[,1]
  species_df <- species_raw[,-1, drop = FALSE]
  
  # 删除临时文件
  unlink(temp_file)
  
  return(species_df)
}

# 读取数据
species_df <- clean_csv_file("speciesallname.csv")
cat("数据维度:", dim(species_df), "\n")

# ============================================
# 2. 分组信息设置
# ============================================
cat("\n=== 设置分组信息 ===\n")

# 冷泉组样本
cold_seep_samples <- c(
  "SY365BB-0-4", "SY365BB-4-8", "SY365BB-8-12",
  "SY366YB-0-4", "SY366YB-4-8", "SY366YB-8-12",
  "SY366YW-0-4", "SY366YW-4-8", "SY366YW-8-12",
  "SY368YW-0-4", "SY368YW-4-8", "SY368YW-8-12",
  "SY456YB-0-4", "SY456YB-4-8", "SY456YB-8-12",
  "SY457BB-0-4", "SY457BB-4-8", "SY457BB-8-12",
  "SY459WG-0-4", "SY459WG-4-8", "SY459WG-8-12",
  "S1_0-3", "S1_6-9", "S1_9-12",
  "S2_0-3", "S2_3-6", "S2_12-15",
  "S4_12-15", "S4_9-12",
  "SQ_58_0-4", "SQ_58_4-8", "SQ_58_-8-12",
  "SQ_81_0-4", "SQ_81_4-8", "SQ_81_-8-12",
  "C1_0-6", "C1_12-18", "C1_6-12",
  "C2_0-6", "C2_12-18", "C2_6-12",
  "C3_0-6", "C3_12-18", "C3_6-12",
  "S13_0-2", "S14_4-6", "S15_8-10",
  "ES_2_0-6"
)

# 非冷泉组样本
non_seep_samples <- c(
  "NS_0-6",
  "S3_0-3", "S3_6-9", "S3_9-12",
  "R2111_N300_0-10", "R2111_N500_0-10", 
  "R2111_S300_0-10", "R2111_S500_0-10"
)

# 创建分组数据框
group_data <- data.frame(
  sample.id = c(cold_seep_samples, non_seep_samples),
  group = c(
    rep("Cold_Seep", length(cold_seep_samples)),
    rep("Non_seep", length(non_seep_samples))
  )
)

# 转置数据
species_df_t <- as.data.frame(t(species_df))

# 匹配样本
common_samples <- intersect(rownames(species_df_t), group_data$sample.id)
species_df_t <- species_df_t[common_samples, ]
group_data <- group_data[group_data$sample.id %in% common_samples, ]
group_data <- group_data[match(rownames(species_df_t), group_data$sample.id), ]

# 设置分组信息
group_col <- c(Cold_Seep = "#90EE90", Non_seep = "#87CEEB")
group_order <- c("Cold_Seep", "Non_seep")

cat("冷泉组样本数:", sum(group_data$group == "Cold_Seep"), "\n")
cat("非冷泉组样本数:", sum(group_data$group == "Non_seep"), "\n")
cat("总样本数:", nrow(group_data), "\n")

# ============================================
# 3. 自定义主题（白色背景，黑色边框，正方形）
# ============================================
custom_theme_square <- function(base_size = 12) {
  theme(
    # 图形边框设置
    plot.background = element_rect(fill = "white", color = "black", linewidth = 1.5),
    
    # 移除所有网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # 设置坐标轴线
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    
    # 坐标轴文字
    axis.text = element_text(color = "black", size = base_size * 0.9),
    axis.title = element_text(color = "black", size = base_size, face = "bold"),
    
    # 图例设置
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_text(face = "bold", size = base_size * 0.9),
    legend.text = element_text(size = base_size * 0.8),
    
    # 标题设置
    plot.title = element_text(hjust = 0.5, size = base_size + 2, face = "bold", 
                              margin = margin(b = 15)),
    
    # 正方形比例
    aspect.ratio = 1,
    
    # 边距
    plot.margin = margin(20, 20, 20, 20),
    
    # 移除背景和边框
    panel.background = element_rect(fill = "white"),
    panel.border = element_blank()
  )
}

# ============================================
# 4. NMDS分析
# ============================================
cat("\n=== 进行NMDS分析 ===\n")

# 计算距离矩阵
bray_dist <- vegdist(species_df_t, method = "bray")

# 执行NMDS
set.seed(123)
nmds_result <- metaMDS(bray_dist, k = 2, trymax = 100)
cat("NMDS应力值 (stress):", nmds_result$stress, "\n")

# 提取坐标
nmds_points <- as.data.frame(nmds_result$points)
colnames(nmds_points) <- c("NMDS1", "NMDS2")
nmds_points$Group <- factor(group_data$group, levels = group_order)

# ============================================
# 5. PCoA分析
# ============================================
cat("\n=== 进行PCoA分析 ===\n")

pcoa_result <- cmdscale(bray_dist, k = 2, eig = TRUE)
eig <- pcoa_result$eig
variance_explained <- eig[1:2] / sum(eig[eig > 0]) * 100
cat("PCoA轴1解释方差:", round(variance_explained[1], 2), "%\n")
cat("PCoA轴2解释方差:", round(variance_explained[2], 2), "%\n")

pcoa_points <- as.data.frame(pcoa_result$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points$Group <- factor(group_data$group, levels = group_order)

# ============================================
# 6. 统计检验
# ============================================
cat("\n=== 进行PERMANOVA统计检验 ===\n")
permanova_result <- adonis2(bray_dist ~ group_data$group, permutations = 999)
cat("PERMANOVA p值:", permanova_result$`Pr(>F)`[1], "\n")

# 准备统计标注
stress_label <- paste("Stress =", round(nmds_result$stress, 3))
pval_label <- paste("p =", sprintf("%.4f", permanova_result$`Pr(>F)`[1]))
pcoa_label1 <- paste("PERMANOVA\np =", sprintf("%.4f", permanova_result$`Pr(>F)`[1]))

# ============================================
# 7. 绘制NMDS图（简洁版）
# ============================================
cat("绘制NMDS图...\n")

nmds_plot <- ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Group, fill = Group)) +
  # 添加椭圆
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, linewidth = 0.5) +
  
  # 添加散点
  geom_point(size = 3, shape = 16, alpha = 0.8) +
  
  # 设置颜色
  scale_color_manual(
    name = "Group",
    values = group_col,
    labels = c("Cold_Seep" = "Cold Seep", "Non_seep" = "Non-seep")
  ) +
  scale_fill_manual(
    name = "Group",
    values = group_col,
    labels = c("Cold_Seep" = "Cold Seep", "Non_seep" = "Non-seep")
  ) +
  
  # 添加统计标注
  annotate("text", x = -Inf, y = Inf, 
           hjust = -0.1, vjust = 1.5,
           label = stress_label, 
           size = 4.5, color = "black", fontface = "bold") +
  
  annotate("text", x = -Inf, y = Inf, 
           hjust = -0.1, vjust = 3.0,
           label = pval_label, 
           size = 4.5, color = "black", fontface = "bold") +
  
  # 标签和标题
  labs(
    title = "NMDS Analysis",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  
  # 应用自定义主题
  theme_minimal() +
  custom_theme_square()

print(nmds_plot)

# 保存NMDS图
ggsave("NMDS_plot_clean.png", plot = nmds_plot, width = 8, height = 8, dpi = 300)
cat("NMDS图已保存: NMDS_plot_clean.png\n")

# ============================================
# 8. 绘制PCoA图（简洁版）
# ============================================
cat("绘制PCoA图...\n")

pcoa_plot <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = Group, fill = Group)) +
  # 添加椭圆
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, linewidth = 0.5) +
  
  # 添加散点
  geom_point(size = 3, shape = 16, alpha = 0.8) +
  
  # 设置颜色
  scale_color_manual(
    name = "Group",
    values = group_col,
    labels = c("Cold_Seep" = "Cold Seep", "Non_seep" = "Non-seep")
  ) +
  scale_fill_manual(
    name = "Group",
    values = group_col,
    labels = c("Cold_Seep" = "Cold Seep", "Non_seep" = "Non-seep")
  ) +
  
  # 添加统计标注
  annotate("text", x = -Inf, y = Inf, 
           hjust = -0.1, vjust = 1.5,
           label = pcoa_label1, 
           size = 4.5, color = "black", fontface = "bold") +
  
  # 标签和标题
  labs(
    title = "PCoA Analysis",
    x = paste("PCoA1 (", round(variance_explained[1], 1), "%)", sep = ""),
    y = paste("PCoA2 (", round(variance_explained[2], 1), "%)", sep = "")
  ) +
  
  # 应用自定义主题
  theme_minimal() +
  custom_theme_square()

print(pcoa_plot)

# 保存PCoA图
ggsave("PCoA_plot_clean.png", plot = pcoa_plot, width = 8, height = 8, dpi = 300)
cat("PCoA图已保存: PCoA_plot_clean.png\n")

# ============================================
# 9. 组合图形
# ============================================
cat("\n=== 创建组合图形 ===\n")

# 移除图例创建更简洁的版本
nmds_simple <- nmds_plot + theme(legend.position = "none")
pcoa_simple <- pcoa_plot + theme(legend.position = "none")

# 提取图例
legend_plot <- ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Group)) +
  geom_point() +
  scale_color_manual(
    name = "Group",
    values = group_col,
    labels = c("Cold_Seep" = "Cold Seep", "Non_seep" = "Non-seep")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold")
  )

# 使用patchwork组合
combined_plot <- (nmds_simple | pcoa_simple) / 
  guide_area() +
  plot_layout(
    heights = c(10, 1),
    guides = "collect"
  ) &
  theme(legend.position = "bottom")

print(combined_plot)

# 保存组合图形
ggsave("combined_plot_clean.png", plot = combined_plot, width = 16, height = 10, dpi = 300)
cat("组合图形已保存: combined_plot_clean.png\n")

# ============================================
# 10. 保存数据结果
# ============================================
# 保存坐标数据
nmds_points_output <- nmds_points
nmds_points_output$Sample <- rownames(species_df_t)
write.csv(nmds_points_output, "NMDS_coordinates_clean.csv", row.names = FALSE)

pcoa_points_output <- pcoa_points
pcoa_points_output$Sample <- rownames(species_df_t)
write.csv(pcoa_points_output, "PCoA_coordinates_clean.csv", row.names = FALSE)

# 保存统计结果
permanova_df <- as.data.frame(permanova_result)
write.csv(permanova_df, "PERMANOVA_results_clean.csv")

cat("\n=== 分析完成 ===\n")
cat("生成的文件:\n")
cat("1. NMDS坐标: NMDS_coordinates_clean.csv\n")
cat("2. PCoA坐标: PCoA_coordinates_clean.csv\n")
cat("3. PERMANOVA结果: PERMANOVA_results_clean.csv\n")
cat("4. NMDS图: NMDS_plot_clean.png\n")
cat("5. PCoA图: PCoA_plot_clean.png\n")
cat("6. 组合图: combined_plot_clean.png\n")

cat("\n图形特点:\n")
cat("✓ 白色背景，无灰色网格线\n")
cat("✓ 黑色外边框\n")
cat("✓ 正方形比例\n")
cat("✓ 标注stress值和p-value值\n")
cat("✓ 不显示样品名称\n")
cat("✓ 简洁美观的设计\n")
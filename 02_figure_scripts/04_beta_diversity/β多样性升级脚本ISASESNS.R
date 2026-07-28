# 加载必要的包
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)  # 添加ggpubr包用于统计检验和条形图
library(ggsignif) # 用于添加显著性标注
# 设置工作目录
setwd("/Users/catherine/Downloads/NMDS/")

# 读取物种丰度表
species_df <- read.csv("speciesall.csv", header = TRUE, row.names = 1, check.names = FALSE)

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
# 3. 自定义主题（新罗马字体，白色背景，黑色边框，正方形）
# ============================================
custom_theme_square_tnr <- function(base_size = 12) {
  theme_classic(base_family = font_family) +
    theme(
      # 图形边框设置（外边框）
      plot.background = element_rect(fill = "white", color = "black", linewidth = 1.5),
      
      # 移除所有网格线
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # 设置坐标轴线
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.15, "cm"),
      
      # 坐标轴文字（使用新罗马字体）
      axis.text = element_text(color = "black", size = base_size, family = font_family),
      axis.title = element_text(color = "black", size = base_size + 1, face = "bold", family = font_family),
      
      # 图例设置（放置在图形内部右上角）
      legend.position = c(0.95, 0.95),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.title = element_text(face = "bold", size = base_size, family = font_family),
      legend.text = element_text(size = base_size - 1, family = font_family),
      legend.key = element_rect(fill = "white", color = NA),
      legend.key.size = unit(0.8, "cm"),
      legend.margin = margin(5, 5, 5, 5),
      
      # 标题设置
      plot.title = element_text(hjust = 0.5, size = base_size + 2, face = "bold", 
                                family = font_family, margin = margin(b = 15)),
      
      # 正方形比例
      aspect.ratio = 1,
      
      # 边距（确保所有图形一致）
      plot.margin = margin(25, 25, 25, 25),
      
      # 坐标刻度线设置
      panel.background = element_rect(fill = "white")
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
permanova_p <- permanova_result$`Pr(>F)`[1]
cat("PERMANOVA p值:", permanova_p, "\n")

# 准备统计标注
stress_label <- paste("Stress =", sprintf("%.3f", nmds_result$stress))
pval_label <- paste("PERMANOVA\np =", ifelse(permanova_p < 0.001, "< 0.001", sprintf("%.3f", permanova_p)))

# ============================================
# 7. 统一坐标轴刻度函数
# ============================================
# 获取统一的坐标轴范围
get_unified_limits <- function(points1, points2, padding = 0.1) {
  # 获取NMDS和PCoA的范围
  nmds_range_x <- range(points1[,1])
  nmds_range_y <- range(points1[,2])
  pcoa_range_x <- range(points2[,1])
  pcoa_range_y <- range(points2[,2])
  
  # 获取最大范围
  unified_x <- range(c(nmds_range_x, pcoa_range_x))
  unified_y <- range(c(nmds_range_y, pcoa_range_y))
  
  # 添加padding
  x_range <- diff(unified_x) * padding
  y_range <- diff(unified_y) * padding
  
  unified_x[1] <- unified_x[1] - x_range
  unified_x[2] <- unified_x[2] + x_range
  unified_y[1] <- unified_y[1] - y_range
  unified_y[2] <- unified_y[2] + y_range
  
  return(list(x = unified_x, y = unified_y))
}

# 获取统一的坐标轴范围
unified_limits <- get_unified_limits(nmds_points[,1:2], pcoa_points[,1:2])

# ============================================
# 8. 绘制NMDS图（简洁版）
# ============================================
cat("绘制NMDS图...\n")

nmds_plot <- ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Group, fill = Group)) +
  # 添加椭圆
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, linewidth = 0.8) +
  
  # 添加散点
  geom_point(size = 3.5, shape = 16, alpha = 0.8) +
  
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
  
  # 统一坐标轴范围
  coord_cartesian(xlim = unified_limits$x, ylim = unified_limits$y) +
  
  # 添加统计标注（右上角）
  annotate("text", x = Inf, y = Inf, 
           hjust = 1.1, vjust = 1.5,
           label = stress_label, 
           size = 4, color = "black", family = font_family, fontface = "bold") +
  
  annotate("text", x = Inf, y = Inf, 
           hjust = 1.1, vjust = 2.8,
           label = pval_label, 
           size = 4, color = "black", family = font_family, fontface = "bold") +
  
  # 标签和标题
  labs(
    title = "NMDS Analysis",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  
  # 应用自定义主题
  custom_theme_square_tnr() +
  
  # 确保刻度线数量合适
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

print(nmds_plot)

# 保存NMDS图为PDF矢量图
ggsave("NMDS_plot_clean.pdf", plot = nmds_plot, width = 8, height = 8, device = cairo_pdf)
cat("NMDS图已保存: NMDS_plot_clean.pdf\n")

# ============================================
# 9. 绘制PCoA图（简洁版）
# ============================================
cat("绘制PCoA图...\n")

# 准备PCoA统计标注
pcoa_pval_label <- paste("PERMANOVA\np =", ifelse(permanova_p < 0.001, "< 0.001", sprintf("%.3f", permanova_p)))

pcoa_plot <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = Group, fill = Group)) +
  # 添加椭圆
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, linewidth = 0.8) +
  
  # 添加散点
  geom_point(size = 3.5, shape = 16, alpha = 0.8) +
  
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
  
  # 统一坐标轴范围
  coord_cartesian(xlim = unified_limits$x, ylim = unified_limits$y) +
  
  # 添加统计标注（右上角）
  annotate("text", x = Inf, y = Inf, 
           hjust = 1.1, vjust = 1.5,
           label = pcoa_pval_label, 
           size = 4, color = "black", family = font_family, fontface = "bold") +
  
  # 标签和标题
  labs(
    title = "PCoA Analysis",
    x = paste("PCoA1 (", round(variance_explained[1], 1), "%)", sep = ""),
    y = paste("PCoA2 (", round(variance_explained[2], 1), "%)", sep = "")
  ) +
  
  # 应用自定义主题
  custom_theme_square_tnr() +
  
  # 确保刻度线数量合适
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

print(pcoa_plot)

# 保存PCoA图为PDF矢量图
ggsave("PCoA_plot_clean.pdf", plot = pcoa_plot, width = 8, height = 8, device = cairo_pdf)
cat("PCoA图已保存: PCoA_plot_clean.pdf\n")

# ============================================
# 10. 组合图形
# ============================================
cat("\n=== 创建组合图形 ===\n")

# 调整图例位置为右下角（用于组合图）
custom_theme_combined <- custom_theme_square_tnr() +
  theme(
    legend.position = c(0.95, 0.05),  # 右下角
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = alpha("white", 0.8), color = "black", linewidth = 0.5)
  )

# 创建组合图形，移除各自的图例
nmds_combined <- nmds_plot + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14)
  ) +
  labs(title = "(A) NMDS Plot")

pcoa_combined <- pcoa_plot + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14)
  ) +
  labs(title = "(B) PCoA Plot")

# 创建共享图例
get_legend <- function(plot) {
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) {
    legend <- tmp$grobs[[leg]]
    return(legend)
  } else {
    return(NULL)
  }
}

# 提取图例
legend_plot <- ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(
    name = "Group",
    values = group_col,
    labels = c("Cold_Seep" = "Cold Seep", "Non_seep" = "Non-seep")
  ) +
  theme_minimal(base_family = font_family) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12, family = font_family),
    legend.text = element_text(size = 11, family = font_family),
    legend.key.size = unit(1, "cm")
  )

legend <- get_legend(legend_plot)

# 使用patchwork组合图形
if (!is.null(legend)) {
  combined_plot <- (nmds_combined + pcoa_combined) / 
    legend +
    plot_layout(heights = c(10, 1))
} else {
  combined_plot <- nmds_combined + pcoa_combined +
    plot_layout(ncol = 2)
}

# 添加整体标题
combined_plot <- combined_plot +
  plot_annotation(
    title = "Beta Diversity Analysis of Microbial Communities",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", family = font_family),
      plot.background = element_rect(fill = "white", color = "black", linewidth = 1.5)
    )
  )

print(combined_plot)

# 保存组合图形为PDF矢量图
ggsave("combined_plot_clean.pdf", plot = combined_plot, width = 16, height = 9, device = cairo_pdf)
cat("组合图形已保存: combined_plot_clean.pdf\n")

# ============================================
# 11. 保存数据结果
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

# ============================================
# 12. 生成分析报告
# ============================================
sink("analysis_report.txt", append = FALSE)
cat("微生物群落Beta多样性分析报告\n")
cat("=============================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("工作目录:", getwd(), "\n\n")

cat("样本信息:\n")
cat("- 冷泉组 (Cold Seep):", sum(group_data$group == "Cold_Seep"), "个样本\n")
cat("- 非冷泉组 (Non-seep):", sum(group_data$group == "Non_seep"), "个样本\n")
cat("- 总样本数:", nrow(group_data), "\n\n")

cat("分析结果:\n")
cat("- NMDS应力值 (Stress):", sprintf("%.3f", nmds_result$stress), "\n")
cat("- PCoA轴1解释方差:", round(variance_explained[1], 2), "%\n")
cat("- PCoA轴2解释方差:", round(variance_explained[2], 2), "%\n")
cat("- PERMANOVA p值:", ifelse(permanova_p < 0.001, "< 0.001", sprintf("%.3f", permanova_p)), "\n\n")

if (permanova_p < 0.05) {
  cat("结论: 冷泉组和非冷泉组微生物群落结构存在显著差异 (p ", 
      ifelse(permanova_p < 0.001, "< 0.001", paste("=", sprintf("%.3f", permanova_p))), ")\n", sep = "")
} else {
  cat("结论: 冷泉组和非冷泉组微生物群落结构无显著差异 (p =", sprintf("%.3f", permanova_p), ")\n")
}

cat("\n生成的文件:\n")
cat("1. NMDS坐标数据: NMDS_coordinates_clean.csv\n")
cat("2. PCoA坐标数据: PCoA_coordinates_clean.csv\n")
cat("3. PERMANOVA统计结果: PERMANOVA_results_clean.csv\n")
cat("4. NMDS图 (PDF): NMDS_plot_clean.pdf\n")
cat("5. PCoA图 (PDF): PCoA_plot_clean.pdf\n")
cat("6. 组合图 (PDF): combined_plot_clean.pdf\n")
cat("7. 分析报告: analysis_report.txt\n")
sink()

cat("\n=== 分析完成 ===\n")
cat("主要改进:\n")
cat("✓ 统一了NMDS和PCoA的坐标轴刻度\n")
cat("✓ 图例和统计标注放置在图形框内右上角\n")
cat("✓ 所有文字使用新罗马字体 (Times New Roman)\n")
cat("✓ 输出为PDF矢量图格式\n")
cat("✓ 黑色外边框，正方形比例\n")
cat("✓ 无灰色网格线，简洁美观\n")
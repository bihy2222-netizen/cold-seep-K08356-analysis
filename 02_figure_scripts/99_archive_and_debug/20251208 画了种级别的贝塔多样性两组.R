# ============================================
# 微生物群落NMDS和PCoA分析脚本
# ============================================

# 加载必要的包
library(vegan)
library(ggplot2)
library(ggrepel)  # 用于标签避免重叠
library(RColorBrewer)

# 设置工作目录
setwd("/Users/catherine/Downloads/师兄交作业系列 ddl/数据提交/")

# 读取物种丰度表
cat("=== 读取数据 ===\n")
species_df <- read.csv("speciesallname.csv", header = TRUE, row.names = 1, check.names = FALSE)

# 创建分组数据框（表层沉积物 vs 深层沉积物）
surface_samples <- c(
  "C1_0-6", "C2_0-6", "C3_0-6", "ES_2_0-6", "NS_0-6",
  "R2111_N300_0-10", "R2111_N500_0-10", "R2111_S300_0-10", "R2111_S500_0-10",
  "S1_0-3", "S13_0-2", "S2_0-3", "S2_3-6", "S3_0-3",
  "SQ_58_0-4", "SQ_81_0-4", "SY365BB-0-4", "SY366YB-0-4",
  "SY366YW-0-4", "SY368YW-0-4", "SY456YB-0-4", "SY457YW-0-4", "SY459YW-0-4"
)

deep_samples <- c(
  "C1_12-18", "C1_6-12", "C2_12-18", "C2_6-12", "C3_12-18", "C3_6-12",
  "S14_4-6", "S15_8-10", "S1_6-9", "S1_9-12", "S2_12-15",
  "S3_6-9", "S3_9-12", "S4_12-15", "S4_9-12",
  "SQ_58_4-8", "SQ_58_-8-12", "SQ_81_4-8", "SQ_81_-8-12",
  "SY365BB-4-8", "SY365BB-8-12", "SY366YB-4-8", "SY366YB-8-12",
  "SY366YW-4-8", "SY366YW-8-12", "SY368YW-4-8", "SY368YW-8-12",
  "SY456YB-4-8", "SY456YB-8-12", "SY457YW-4-8", "SY457YW-8-12",
  "SY459YW-4-8", "SY459YW-8-12"
)

# 创建分组数据框
group_data <- data.frame(
  sample.id = c(surface_samples, deep_samples),
  group = c(
    rep("surface_sediment", length(surface_samples)),
    rep("deep_sediment", length(deep_samples))
  )
)

cat("表层沉积物样本数:", length(surface_samples), "\n")
cat("深层沉积物样本数:", length(deep_samples), "\n")

# 转置物种数据，使样本为行，物种为列
species_df_t <- as.data.frame(t(species_df))

# 检查样本匹配
common_samples <- intersect(rownames(species_df_t), group_data$sample.id)
cat("共有样本数量:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  cat("样本名不匹配！尝试模糊匹配...\n")
  species_samples <- rownames(species_df_t)
  group_samples <- group_data$sample.id
  
  # 查找可能的匹配
  matched_indices <- sapply(species_samples, function(x) {
    matches <- agrep(x, group_samples, max.distance = 0.1, value = TRUE)
    if (length(matches) > 0) return(matches[1]) else return(NA)
  })
  
  matched_pairs <- na.omit(data.frame(
    species_sample = species_samples,
    group_sample = matched_indices
  ))
  
  if (nrow(matched_pairs) > 0) {
    cat("找到", nrow(matched_pairs), "个匹配\n")
    common_samples <- matched_pairs$species_sample
  } else {
    stop("无法匹配样本")
  }
}

# 过滤数据，只保留共有的样本
species_df_t <- species_df_t[common_samples, ]
group_data <- group_data[group_data$sample.id %in% common_samples, ]

# 确保样本顺序一致
group_data <- group_data[match(rownames(species_df_t), group_data$sample.id), ]

# 设置颜色和形状
group_col <- c(surface_sediment = "#00FF7F", deep_sediment = "#FF1493")
group_order <- c("surface_sediment", "deep_sediment")

# ============================================
# 1. NMDS分析（基于Bray-Curtis距离）
# ============================================
cat("\n=== 进行NMDS分析 ===\n")

# 计算Bray-Curtis距离矩阵
bray_dist <- vegdist(species_df_t, method = "bray")

# 执行NMDS分析
set.seed(123)  # 设置随机种子以确保结果可重复
nmds_result <- metaMDS(bray_dist, k = 2, trymax = 100)

cat("NMDS应力值 (stress):", nmds_result$stress, "\n")
cat("NMDS收敛:", nmds_result$converged, "\n")

# 提取NMDS坐标
nmds_points <- as.data.frame(nmds_result$points)
colnames(nmds_points) <- c("NMDS1", "NMDS2")
nmds_points$Sample <- rownames(nmds_points)
nmds_points$Group <- group_data$group
nmds_points$Group <- factor(nmds_points$Group, levels = group_order)

# 计算各组中心点
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Group, data = nmds_points, FUN = mean)

# 绘制NMDS图
cat("绘制NMDS图...\n")
nmds_plot <- ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Group, shape = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  # 添加样本标签（可选）
  geom_text_repel(aes(label = Sample), size = 2.5, max.overlaps = 20, show.legend = FALSE) +
  # 添加各组中心点
  geom_point(data = centroids, aes(x = NMDS1, y = NMDS2), 
             size = 5, shape = 17, alpha = 0.7, show.legend = FALSE) +
  # 添加椭圆（95%置信区间）
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.2, level = 0.95) +
  # 设置颜色和形状
  scale_color_manual(values = group_col, 
                     labels = c("surface_sediment" = "Surface Sediment", 
                                "deep_sediment" = "Deep Sediment")) +
  scale_fill_manual(values = group_col,
                    labels = c("surface_sediment" = "Surface Sediment", 
                               "deep_sediment" = "Deep Sediment")) +
  scale_shape_manual(values = c(16, 18),
                     labels = c("surface_sediment" = "Surface Sediment", 
                                "deep_sediment" = "Deep Sediment")) +
  # 主题和标签
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = paste("NMDS Plot of Microbial Communities (Stress =", 
                  round(nmds_result$stress, 3), ")"),
    x = "NMDS1",
    y = "NMDS2",
    caption = "Based on Bray-Curtis distance"
  )

print(nmds_plot)

# 保存NMDS图
ggsave("NMDS_plot.png", plot = nmds_plot, width = 10, height = 8, dpi = 300)
ggsave("NMDS_plot.pdf", plot = nmds_plot, width = 10, height = 8)
cat("NMDS图已保存: NMDS_plot.png, NMDS_plot.pdf\n")

# 保存NMDS坐标数据
write.csv(nmds_points, "NMDS_coordinates.csv", row.names = FALSE)
cat("NMDS坐标数据已保存: NMDS_coordinates.csv\n")

# ============================================
# 2. PCoA分析（基于Bray-Curtis距离）
# ============================================
cat("\n=== 进行PCoA分析 ===\n")

# 执行PCoA分析
pcoa_result <- cmdscale(bray_dist, k = 2, eig = TRUE)

# 计算解释的方差比例
eig <- pcoa_result$eig
variance_explained <- eig[1:2] / sum(eig[eig > 0]) * 100

cat("PCoA轴1解释方差:", round(variance_explained[1], 2), "%\n")
cat("PCoA轴2解释方差:", round(variance_explained[2], 2), "%\n")

# 提取PCoA坐标
pcoa_points <- as.data.frame(pcoa_result$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points$Sample <- rownames(pcoa_points)
pcoa_points$Group <- group_data$group
pcoa_points$Group <- factor(pcoa_points$Group, levels = group_order)

# 计算各组中心点
centroids_pcoa <- aggregate(cbind(PCoA1, PCoA2) ~ Group, data = pcoa_points, FUN = mean)

# 绘制PCoA图
cat("绘制PCoA图...\n")
pcoa_plot <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = Group, shape = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  # 添加样本标签（可选）
  geom_text_repel(aes(label = Sample), size = 2.5, max.overlaps = 20, show.legend = FALSE) +
  # 添加各组中心点
  geom_point(data = centroids_pcoa, aes(x = PCoA1, y = PCoA2), 
             size = 5, shape = 17, alpha = 0.7, show.legend = FALSE) +
  # 添加椭圆（95%置信区间）
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.2, level = 0.95) +
  # 设置颜色和形状
  scale_color_manual(values = group_col, 
                     labels = c("surface_sediment" = "Surface Sediment", 
                                "deep_sediment" = "Deep Sediment")) +
  scale_fill_manual(values = group_col,
                    labels = c("surface_sediment" = "Surface Sediment", 
                               "deep_sediment" = "Deep Sediment")) +
  scale_shape_manual(values = c(16, 18),
                     labels = c("surface_sediment" = "Surface Sediment", 
                                "deep_sediment" = "Deep Sediment")) +
  # 主题和标签
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = paste("PCoA Plot of Microbial Communities"),
    x = paste("PCoA1 (", round(variance_explained[1], 1), "%)", sep = ""),
    y = paste("PCoA2 (", round(variance_explained[2], 1), "%)", sep = ""),
    caption = "Based on Bray-Curtis distance"
  )

print(pcoa_plot)

# 保存PCoA图
ggsave("PCoA_plot.png", plot = pcoa_plot, width = 10, height = 8, dpi = 300)
ggsave("PCoA_plot.pdf", plot = pcoa_plot, width = 10, height = 8)
cat("PCoA图已保存: PCoA_plot.png, PCoA_plot.pdf\n")

# 保存PCoA坐标数据
write.csv(pcoa_points, "PCoA_coordinates.csv", row.names = FALSE)
cat("PCoA坐标数据已保存: PCoA_coordinates.csv\n")

# ============================================
# 3. 统计检验（PERMANOVA）
# ============================================
cat("\n=== 进行PERMANOVA统计检验 ===\n")

# 执行PERMANOVA分析
permanova_result <- adonis2(bray_dist ~ group_data$group, permutations = 999)
cat("PERMANOVA结果:\n")
print(permanova_result)

# 保存PERMANOVA结果
write.csv(as.data.frame(permanova_result), "PERMANOVA_results.csv")
cat("PERMANOVA结果已保存: PERMANOVA_results.csv\n")

# 进行成对比较（如果PERMANOVA显著）
if (permanova_result$`Pr(>F)`[1] < 0.05) {
  cat("\n进行成对PERMANOVA比较...\n")
  
  # 手动进行成对比较
  groups <- unique(group_data$group)
  pairwise_results <- data.frame()
  
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      group1 <- groups[i]
      group2 <- groups[j]
      
      # 选择这两个组的样本
      idx1 <- group_data$group == group1
      idx2 <- group_data$group == group2
      idx <- idx1 | idx2
      
      # 提取子集
      sub_dist <- as.dist(as.matrix(bray_dist)[idx, idx])
      sub_group <- group_data$group[idx]
      
      # 执行PERMANOVA
      pair_result <- adonis2(sub_dist ~ sub_group, permutations = 999)
      
      pairwise_results <- rbind(pairwise_results, data.frame(
        Comparison = paste(group1, "vs", group2),
        R2 = pair_result$R2[1],
        F_value = pair_result$F[1],
        p_value = pair_result$`Pr(>F)`[1]
      ))
    }
  }
  
  # 添加显著性标记
  pairwise_results$Significance <- ifelse(pairwise_results$p_value < 0.001, "***",
                                          ifelse(pairwise_results$p_value < 0.01, "**",
                                                 ifelse(pairwise_results$p_value < 0.05, "*", "ns")))
  
  cat("成对比较结果:\n")
  print(pairwise_results)
  
  # 保存成对比较结果
  write.csv(pairwise_results, "pairwise_PERMANOVA_results.csv", row.names = FALSE)
  cat("成对PERMANOVA结果已保存: pairwise_PERMANOVA_results.csv\n")
}

# ============================================
# 4. 组合图形（可选）
# ============================================
cat("\n=== 创建组合图形 ===\n")

# 创建一个包含NMDS和PCoA的组合图形
library(patchwork)

# 简化图形以便组合
nmds_simple <- ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(level = 0.95, alpha = 0.2) +
  scale_color_manual(values = group_col) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 10)
  ) +
  labs(
    title = paste("NMDS (Stress =", round(nmds_result$stress, 3), ")"),
    x = "NMDS1",
    y = "NMDS2"
  )

pcoa_simple <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(level = 0.95, alpha = 0.2) +
  scale_color_manual(values = group_col) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 10)
  ) +
  labs(
    title = paste("PCoA"),
    x = paste("PCoA1 (", round(variance_explained[1], 1), "%)", sep = ""),
    y = paste("PCoA2 (", round(variance_explained[2], 1), "%)", sep = "")
  )

# 组合图形
combined_plot <- nmds_simple + pcoa_simple +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom')

print(combined_plot)

# 保存组合图形
ggsave("combined_NMDS_PCoA.png", plot = combined_plot, width = 12, height = 6, dpi = 300)
cat("组合图形已保存: combined_NMDS_PCoA.png\n")

# ============================================
# 5. 生成分析报告
# ============================================
cat("\n=== 生成分析报告 ===\n")

sink("beta_diversity_analysis_report.txt")
cat("微生物群落Beta多样性分析报告\n")
cat("=============================\n\n")
cat("分析时间:", Sys.time(), "\n")
cat("工作目录:", getwd(), "\n")
cat("输入文件: speciesallname.csv\n\n")

cat("样本信息:\n")
cat("- 表层沉积物样本数:", sum(group_data$group == "surface_sediment"), "\n")
cat("- 深层沉积物样本数:", sum(group_data$group == "deep_sediment"), "\n")
cat("- 总样本数:", nrow(group_data), "\n\n")

cat("NMDS分析结果:\n")
cat("- NMDS应力值 (Stress):", nmds_result$stress, "\n")
cat("- 收敛状态:", nmds_result$converged, "\n")
cat("- 迭代次数:", nmds_result$iters, "\n\n")

cat("PCoA分析结果:\n")
cat("- PCoA轴1解释方差:", round(variance_explained[1], 2), "%\n")
cat("- PCoA轴2解释方差:", round(variance_explained[2], 2), "%\n")
cat("- 总解释方差:", round(sum(variance_explained), 2), "%\n\n")

cat("PERMANOVA统计检验结果:\n")
cat("R² =", permanova_result$R2[1], "\n")
cat("F值 =", permanova_result$F[1], "\n")
cat("P值 =", permanova_result$`Pr(>F)`[1], "\n")

if (permanova_result$`Pr(>F)`[1] < 0.05) {
  cat("结论: 组间微生物群落结构存在显著差异 (p < 0.05)\n")
} else {
  cat("结论: 组间微生物群落结构无显著差异\n")
}

cat("\n生成的文件:\n")
cat("1. NMDS坐标: NMDS_coordinates.csv\n")
cat("2. PCoA坐标: PCoA_coordinates.csv\n")
cat("3. PERMANOVA结果: PERMANOVA_results.csv\n")
if (exists("pairwise_results")) {
  cat("4. 成对比较结果: pairwise_PERMANOVA_results.csv\n")
}
cat("5. 图形文件:\n")
cat("   - NMDS_plot.png/pdf\n")
cat("   - PCoA_plot.png/pdf\n")
cat("   - combined_NMDS_PCoA.png\n")
cat("6. 分析报告: beta_diversity_analysis_report.txt\n")
sink()

cat("分析报告已保存: beta_diversity_analysis_report.txt\n")

# ============================================
# 6. 显示关键结果
# ============================================
cat("\n" , strrep("=", 60), "\n", sep = "")
cat("Beta多样性分析完成！\n")
cat(strrep("=", 60), "\n")

cat("\n关键结果:\n")
cat("1. NMDS应力值:", round(nmds_result$stress, 3), "\n")
cat("2. PCoA解释方差: 轴1 =", round(variance_explained[1], 1), 
    "%, 轴2 =", round(variance_explained[2], 1), "%\n")
cat("3. PERMANOVA P值:", permanova_result$`Pr(>F)`[1], "\n")

if (permanova_result$`Pr(>F)`[1] < 0.05) {
  cat("4. 结论: 表层和深层沉积物微生物群落结构存在显著差异！\n")
} else {
  cat("4. 结论: 表层和深层沉积物微生物群落结构无显著差异。\n")
}

cat("\n所有结果文件已保存在当前目录。\n")
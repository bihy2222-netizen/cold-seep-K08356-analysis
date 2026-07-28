# 加载必要的包
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# 读取物种丰度表
species_df <- read.csv("speciesall.csv", header = TRUE, row.names = 1, check.names = FALSE)

# 创建正确的分组数据框（基于您提供的信息）
group_data <- data.frame(
  sample.id = c(
    "C1_0-6", "C1_12-18", "C1_6-12", "C2_0-6", "C2_12-18", "C2_6-12",
    "C3_0-6", "C3_12-18", "C3_6-12", "ES_2_0-6", "NS_0-6", "R2111_N300_0-10",
    "R2111_N500_0-10", "R2111_S300_0-10", "R2111_S500_0-10", "S1_0-3", "S13_0-2",
    "S14_4-6", "S15_8-10", "S1_6-9", "S1_9-12", "S2_0-3", "S2_12-15", "S2_3-6",
    "S3_0-3", "S3_6-9", "S3_9-12", "S4_12-15", "S4_9-12", "SQ_58_0-4", "SQ_58_4-8",
    "SQ_58_-8-12", "SQ_81_0-4", "SQ_81_4-8", "SQ_81_-8-12", "SY365BB-0-4", "SY365BB-4-8",
    "SY365BB-8-12", "SY366YB-0-4", "SY366YB-4-8", "SY366YB-8-12", "SY366YW-0-4",
    "SY366YW-4-8", "SY366YW-8-12", "SY368YW-0-4", "SY368YW-4-8", "SY368YW-8-12",
    "SY456YB-0-4", "SY456YB-4-8", "SY456YB-8-12", "SY457YW-0-4", "SY457YW-4-8",
    "SY457YW-8-12", "SY459YW-0-4", "SY459YW-4-8", "SY459YW-8-12"
  ),
  group = c(
    "ES", "ES", "ES", "ES", "ES", "ES", "ES", "ES", "ES", "ES", "NS", "NS", "NS", "NS", "NS",
    "AS", "ES", "ES", "ES", "AS", "AS", "AS", "AS", "AS", "NS", "NS", "NS", "AS", "AS", "AS",
    "AS", "AS", "AS", "AS", "AS", "IS", "IS", "IS", "IS", "IS", "IS", "IS", "IS", "IS", "IS",
    "IS", "IS", "IS", "IS", "IS", "IS", "IS", "IS", "IS", "IS", "IS"
  )
)

# 保存正确的分组文件
write.csv(group_data, "group_corrected.csv", row.names = FALSE)
cat("已创建正确的分组文件: group_corrected.csv\n")

# 使用正确的分组数据
group_df <- group_data

# 转置物种数据，使样本为行，物种为列
species_df_t <- as.data.frame(t(species_df))

# 检查样本匹配
cat("物种表样本数量:", nrow(species_df_t), "\n")
cat("分组文件样本数量:", nrow(group_df), "\n")

common_samples <- intersect(rownames(species_df_t), group_df$sample.id)
cat("共有样本数量:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  cat("样本名不匹配！详细信息：\n")
  cat("物种表样本示例:", head(rownames(species_df_t), 3), "\n")
  cat("分组文件样本示例:", head(group_df$sample.id, 3), "\n")
  
  # 尝试自动匹配
  cat("尝试自动匹配样本...\n")
  species_samples <- rownames(species_df_t)
  group_samples <- group_df$sample.id
  
  # 检查是否有完全匹配的
  matched <- species_samples %in% group_samples
  cat("完全匹配的样本数量:", sum(matched), "\n")
  
  if (sum(matched) > 0) {
    common_samples <- species_samples[matched]
    cat("使用匹配的样本继续分析...\n")
  } else {
    stop("无法匹配样本，请检查样本名称")
  }
}

# 过滤数据，只保留共有的样本
species_df_t <- species_df_t[common_samples, ]
group_df <- group_df[group_df$sample.id %in% common_samples, ]

# 确保样本顺序一致
group_df <- group_df[match(rownames(species_df_t), group_df$sample.id), ]

# 显示分组信息
cat("所有组别:", unique(group_df$group), "\n")
cat("各组样本数量:\n")
print(table(group_df$group))

# 检查数据是否为数值型
cat("数据格式检查:\n")
cat("物种数据类:", class(species_df_t[,1]), "\n")

# 将数据转换为数值型（如果需要）
species_df_t_numeric <- as.data.frame(lapply(species_df_t, function(x) {
  if (!is.numeric(x)) {
    as.numeric(as.character(x))
  } else {
    x
  }
}))
rownames(species_df_t_numeric) <- rownames(species_df_t)

# 将物种数据中的零替换为一个极小值，避免对数计算错误
species_df_t_numeric[species_df_t_numeric == 0] <- 1e-10

# 计算 Alpha 多样性指数
cat("计算 Alpha 多样性指数...\n")

# 检查数据是否适合计算多样性
if (any(is.na(species_df_t_numeric))) {
  cat("发现NA值，用0替换...\n")
  species_df_t_numeric[is.na(species_df_t_numeric)] <- 1e-10
}

shannon <- diversity(species_df_t_numeric, index = "shannon")
simpson <- diversity(species_df_t_numeric, index = "simpson")
observed_species <- rowSums(species_df_t_numeric > 1e-10)

# 对于 Chao1，需要整数数据
species_df_int <- round(species_df_t_numeric * 1000)  # 放大后四舍五入
chao1 <- estimateR(species_df_int)[2, ]

# 构建 Alpha 多样性数据框
alpha_df <- data.frame(
  Sample = rownames(species_df_t_numeric),
  Group = group_df$group,
  Shannon = shannon,
  Simpson = simpson,
  Observed_Species = observed_species,
  Chao1 = chao1
)

# 输出 Alpha 多样性表格
write.csv(alpha_df, "alpha_diversity_results.csv", row.names = FALSE)
cat("Alpha 多样性结果已保存: alpha_diversity_results.csv\n")

# 设置组别顺序和颜色（按照您的要求）
group_order <- c("IS", "AS", "ES", "NS")
group_col <- c(IS = "#00FF7F", AS = "#FF1493", ES = "#FFA500", NS = "#00FFFF")

alpha_long <- alpha_df %>%
  pivot_longer(cols = c(Shannon, Simpson, Observed_Species, Chao1),
               names_to = "Index",
               values_to = "Value")

# 设置因子水平
alpha_long$Group <- factor(alpha_long$Group, levels = group_order)
alpha_df$Group <- factor(alpha_df$Group, levels = group_order)

# 绘制 Alpha 多样性箱线图
cat("绘制 Alpha 多样性箱线图...\n")
p <- ggplot(alpha_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  facet_wrap(~ Index, scales = "free_y", ncol = 2) +
  theme_minimal() +
  labs(title = "Alpha Diversity Across Groups (IS, AS, ES, NS)",
       y = "Diversity Index Value",
       x = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.text = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = group_col)

print(p)

# 保存图形
ggsave("alpha_diversity_plot.png", plot = p, width = 10, height = 8, dpi = 300)
cat("Alpha 多样性图形已保存: alpha_diversity_plot.png\n")

# 进行组间 Alpha 多样性差异检验（Kruskal-Wallis 检验）
cat("进行组间差异检验...\n")
alpha_test_results <- data.frame(
  Index = character(),
  Chi_squared = numeric(),
  df = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (idx in c("Shannon", "Simpson", "Observed_Species", "Chao1")) {
  test_result <- kruskal.test(alpha_df[[idx]] ~ alpha_df$Group)
  alpha_test_results <- rbind(alpha_test_results, data.frame(
    Index = idx,
    Chi_squared = round(test_result$statistic, 4),
    df = test_result$parameter,
    p_value = round(test_result$p.value, 6)
  ))
}

# 添加显著性标记
alpha_test_results$Significance <- ifelse(alpha_test_results$p_value < 0.001, "***",
                                          ifelse(alpha_test_results$p_value < 0.01, "**",
                                                 ifelse(alpha_test_results$p_value < 0.05, "*", "ns")))

# 输出检验结果
write.csv(alpha_test_results, "alpha_diversity_kruskal_test.csv", row.names = FALSE)
cat("Kruskal-Wallis 检验结果已保存: alpha_diversity_kruskal_test.csv\n")

# 如果有显著差异，进行事后检验（Dunn test）
if (any(alpha_test_results$p_value < 0.05)) {
  cat("进行事后检验...\n")
  if (!require(FSA, quietly = TRUE)) {
    install.packages("FSA")
    library(FSA)
  }
  
  posthoc_results <- list()
  
  for (idx in c("Shannon", "Simpson", "Observed_Species", "Chao1")) {
    if (alpha_test_results$p_value[alpha_test_results$Index == idx] < 0.05) {
      cat("对", idx, "指数进行 Dunn 事后检验...\n")
      dunn_result <- dunnTest(alpha_df[[idx]] ~ alpha_df$Group, method = "bh")
      posthoc_results[[idx]] <- dunn_result$res
      
      # 保存事后检验结果
      write.csv(dunn_result$res, paste0("posthoc_dunn_", idx, ".csv"), row.names = FALSE)
      cat("事后检验结果已保存: posthoc_dunn_", idx, ".csv\n", sep = "")
    }
  }
} else {
  cat("所有 Alpha 多样性指数在组间均无显著差异 (p > 0.05)\n")
}

# 输出各组 Alpha 多样性的描述性统计
desc_stats <- alpha_df %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    Shannon_mean = round(mean(Shannon, na.rm = TRUE), 3),
    Shannon_sd = round(sd(Shannon, na.rm = TRUE), 3),
    Simpson_mean = round(mean(Simpson, na.rm = TRUE), 3),
    Simpson_sd = round(sd(Simpson, na.rm = TRUE), 3),
    Observed_mean = round(mean(Observed_Species, na.rm = TRUE), 1),
    Observed_sd = round(sd(Observed_Species, na.rm = TRUE), 1),
    Chao1_mean = round(mean(Chao1, na.rm = TRUE), 1),
    Chao1_sd = round(sd(Chao1, na.rm = TRUE), 1)
  )

write.csv(desc_stats, "alpha_diversity_descriptive_stats.csv", row.names = FALSE)
cat("描述性统计已保存: alpha_diversity_descriptive_stats.csv\n")

# 显示汇总信息
cat("\n=== Alpha 多样性分析完成 ===\n")
cat("分析的组别:", paste(group_order, collapse = ", "), "\n")
cat("组别数量:", length(unique(alpha_df$Group)), "\n")
cat("总样本数量:", nrow(alpha_df), "\n")
cat("\n结果文件:\n")
cat("- Alpha 多样性结果: alpha_diversity_results.csv\n")
cat("- 描述性统计: alpha_diversity_descriptive_stats.csv\n")
cat("- Kruskal-Wallis 检验: alpha_diversity_kruskal_test.csv\n")
if (any(alpha_test_results$p_value < 0.05)) {
  cat("- 事后检验结果: posthoc_dunn_*.csv\n")
}
cat("- 图形文件: alpha_diversity_plot.png\n")
cat("- 分组文件: group_corrected.csv\n")

# 显示检验结果
cat("\nKruskal-Wallis 检验结果:\n")
print(alpha_test_results)

# 显示描述性统计
cat("\n各组 Alpha 多样性描述性统计:\n")
print(desc_stats)

# 显示各组样本数量
cat("\n各组样本数量:\n")
sample_counts <- table(alpha_df$Group)
print(sample_counts)

# 保存样本数量信息
sample_count_df <- data.frame(
  Group = names(sample_counts),
  Sample_Count = as.numeric(sample_counts)
)
write.csv(sample_count_df, "group_sample_counts.csv", row.names = FALSE)
cat("各组样本数量已保存: group_sample_counts.csv\n")
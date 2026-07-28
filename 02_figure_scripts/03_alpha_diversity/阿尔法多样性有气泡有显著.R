# 加载必要的包
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggsignif)

# 设置工作目录到phylum文件所在路径
setwd("/Users/catherine/Downloads/metaphlan 作图/class 级别作图-柱状图/phylum/")

# 首先检查文件内容
cat("检查phylum.txt文件内容...\n")
file_content <- readLines("phylum.txt", n = 10)
cat("文件前10行:\n")
print(file_content)

# 尝试不同的读取方法
cat("尝试读取phylum.txt文件...\n")

# 方法1：使用read.delim
tryCatch({
  phylum_df <- read.delim("phylum.txt", header = TRUE, row.names = 1, 
                          check.names = FALSE, sep = "\t")
  cat("方法1成功\n")
}, error = function(e) {
  cat("方法1失败:", e$message, "\n")
  
  # 方法2：使用read.table
  tryCatch({
    phylum_df <- read.table("phylum.txt", header = TRUE, row.names = 1,
                            check.names = FALSE, sep = "\t", fill = TRUE)
    cat("方法2成功\n")
  }, error = function(e) {
    cat("方法2失败:", e$message, "\n")
    
    # 方法3：手动处理
    cat("使用手动处理方法...\n")
    raw_data <- readLines("phylum.txt")
    # 清理数据
    raw_data <- gsub('"', '', raw_data)  # 移除引号
    raw_data <- raw_data[raw_data != ""]  # 移除空行
    
    # 提取列名（第一行）
    headers <- strsplit(raw_data[1], "\t")[[1]]
    
    # 创建数据框
    data_list <- list()
    for(i in 2:length(raw_data)) {
      row_data <- strsplit(raw_data[i], "\t")[[1]]
      if(length(row_data) == length(headers)) {
        row_name <- row_data[1]
        numeric_data <- as.numeric(row_data[-1])
        data_list[[row_name]] <- numeric_data
      }
    }
    
    phylum_df <- as.data.frame(do.call(rbind, data_list))
    colnames(phylum_df) <- headers[-1]
    cat("方法3成功\n")
  })
})

# 检查数据
cat("数据维度:", dim(phylum_df), "\n")
cat("数据类型:", class(phylum_df), "\n")
cat("前5行5列:\n")
print(phylum_df[1:5, 1:5])

# 创建分组数据框
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

# 保存分组文件
write.csv(group_data, "group_corrected_phylum.csv", row.names = FALSE)
cat("已创建分组文件: group_corrected_phylum.csv\n")

# 使用分组数据
group_df <- group_data

# 转置门水平数据，使样本为行，物种为列
phylum_df_t <- as.data.frame(t(phylum_df))

# 检查样本匹配
cat("门水平数据样本数量:", nrow(phylum_df_t), "\n")
cat("分组文件样本数量:", nrow(group_df), "\n")

common_samples <- intersect(rownames(phylum_df_t), group_df$sample.id)
cat("共有样本数量:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  cat("样本名不匹配！详细信息：\n")
  cat("门水平数据样本示例:", head(rownames(phylum_df_t), 3), "\n")
  cat("分组文件样本示例:", head(group_df$sample.id, 3), "\n")
  stop("无法匹配样本，请检查样本名称")
}

# 过滤数据，只保留共有的样本
phylum_df_t <- phylum_df_t[common_samples, ]
group_df <- group_df[group_df$sample.id %in% common_samples, ]

# 确保样本顺序一致
group_df <- group_df[match(rownames(phylum_df_t), group_df$sample.id), ]

# 显示分组信息
cat("所有组别:", unique(group_df$group), "\n")
cat("各组样本数量:\n")
print(table(group_df$group))

# 将数据转换为数值型
phylum_df_t_numeric <- as.data.frame(lapply(phylum_df_t, function(x) {
  if (!is.numeric(x)) {
    as.numeric(as.character(x))
  } else {
    x
  }
}))
rownames(phylum_df_t_numeric) <- rownames(phylum_df_t)

# 将数据中的零替换为一个极小值，避免对数计算错误
phylum_df_t_numeric[phylum_df_t_numeric == 0] <- 1e-10

# 计算 Alpha 多样性指数
cat("计算 Alpha 多样性指数...\n")

# 检查数据是否适合计算多样性
if (any(is.na(phylum_df_t_numeric))) {
  cat("发现NA值，用极小值替换...\n")
  phylum_df_t_numeric[is.na(phylum_df_t_numeric)] <- 1e-10
}

shannon <- diversity(phylum_df_t_numeric, index = "shannon")
simpson <- diversity(phylum_df_t_numeric, index = "simpson")
observed_species <- rowSums(phylum_df_t_numeric > 1e-10)

# 对于 Chao1，需要整数数据
phylum_df_int <- round(phylum_df_t_numeric * 1000)
chao1 <- estimateR(phylum_df_int)[2, ]

# 构建 Alpha 多样性数据框
alpha_df <- data.frame(
  Sample = rownames(phylum_df_t_numeric),
  Group = group_df$group,
  Shannon = shannon,
  Simpson = simpson,
  Observed_Species = observed_species,
  Chao1 = chao1
)

# 输出 Alpha 多样性表格
write.csv(alpha_df, "alpha_diversity_results_phylum.csv", row.names = FALSE)
cat("Alpha 多样性结果已保存: alpha_diversity_results_phylum.csv\n")

# 设置组别顺序和颜色
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
  labs(title = "Alpha Diversity Across Groups (IS, AS, ES, NS) - Phylum Level",
       y = "Diversity Index Value",
       x = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.text = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = group_col)

print(p)

# 保存图形
ggsave("alpha_diversity_plot_phylum.png", plot = p, width = 10, height = 8, dpi = 300)
cat("Alpha 多样性图形已保存: alpha_diversity_plot_phylum.png\n")

# 绘制Shannon指数条形图
cat("绘制Shannon指数条形图并进行统计检验...\n")

mycol <- group_col

shannon_plot <- ggbarplot(alpha_df, x = "Group", y = "Shannon", fill = "Group",
                          palette = mycol, legend = "none",
                          add = "mean_se", 
                          add.params = list(width = 0.3),
                          label = FALSE) +
  geom_jitter(aes(color = Group), width = 0.2, alpha = 0.7, size = 2) +
  scale_color_manual(values = mycol) +
  labs(title = "Shannon Diversity Index by Group - Phylum Level",
       x = "Group", y = "Shannon Diversity Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none")

# 进行成对Wilcoxon检验
cat("进行成对Wilcoxon检验...\n")
comparison_pairs <- list(
  c('IS', 'AS'), c('IS', 'ES'), c('IS', 'NS'),
  c('AS', 'ES'), c('AS', 'NS'), c('ES', 'NS')
)

pairwise_results <- data.frame()
for (pair in comparison_pairs) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  data1 <- alpha_df$Shannon[alpha_df$Group == group1]
  data2 <- alpha_df$Shannon[alpha_df$Group == group2]
  
  test_result <- wilcox.test(data1, data2)
  p_val <- test_result$p.value
  
  pairwise_results <- rbind(pairwise_results, 
                            data.frame(Group1 = group1, Group2 = group2, p_value = p_val))
}

# 添加显著性标注
y_max <- max(alpha_df$Shannon)
step <- y_max * 0.08
i <- 0

for (pair in comparison_pairs) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  p_val <- pairwise_results$p_value[pairwise_results$Group1 == group1 & 
                                      pairwise_results$Group2 == group2]
  
  sig_label <- ifelse(p_val < 0.001, "***",
                      ifelse(p_val < 0.01, "**",
                             ifelse(p_val < 0.05, "*", "ns")))
  
  if(p_val < 0.05) {
    i <- i + 1
    y_position <- y_max + step * i
    
    shannon_plot <- shannon_plot + 
      geom_signif(
        comparisons = list(c(group1, group2)),
        annotations = sig_label,
        y_position = y_position,
        tip_length = 0.01,
        vjust = 0.5,
        textsize = 4
      )
  }
}

# 添加统计检验结果
anova_result <- kruskal.test(Shannon ~ Group, data = alpha_df)
anova_p <- anova_result$p.value
anova_label <- ifelse(anova_p < 0.001, "Kruskal-Wallis: p < 0.001",
                      ifelse(anova_p < 0.01, paste0("Kruskal-Wallis: p = ", round(anova_p, 4)),
                             paste0("Kruskal-Wallis: p = ", round(anova_p, 3))))

shannon_plot <- shannon_plot + 
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = anova_label, size = 4, colour = "black")

print(shannon_plot)

# 保存Shannon指数条形图
ggsave("shannon_diversity_barplot_phylum.png", plot = shannon_plot, width = 8, height = 8, dpi = 300)
cat("Shannon指数条形图已保存: shannon_diversity_barplot_phylum.png\n")

# 进行统计检验
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
write.csv(alpha_test_results, "alpha_diversity_kruskal_test_phylum.csv", row.names = FALSE)
cat("Kruskal-Wallis 检验结果已保存: alpha_diversity_kruskal_test_phylum.csv\n")

# 保存成对比较结果
write.csv(pairwise_results, "shannon_pairwise_wilcox_phylum.csv", row.names = FALSE)
cat("Shannon指数成对比较结果已保存: shannon_pairwise_wilcox_phylum.csv\n")

# 显示汇总信息
cat("\n=== Alpha 多样性分析完成（门水平）===\n")
cat("分析的组别:", paste(group_order, collapse = ", "), "\n")
cat("组别数量:", length(unique(alpha_df$Group)), "\n")
cat("总样本数量:", nrow(alpha_df), "\n")
cat("\n结果文件:\n")
cat("- Alpha 多样性结果: alpha_diversity_results_phylum.csv\n")
cat("- Kruskal-Wallis 检验: alpha_diversity_kruskal_test_phylum.csv\n")
cat("- Shannon指数成对比较: shannon_pairwise_wilcox_phylum.csv\n")
cat("- Alpha多样性箱线图: alpha_diversity_plot_phylum.png\n")
cat("- Shannon指数条形图: shannon_diversity_barplot_phylum.png\n")
cat("- 分组文件: group_corrected_phylum.csv\n")

# 显示检验结果
cat("\nKruskal-Wallis 检验结果:\n")
print(alpha_test_results)

# 显示各组样本数量
cat("\n各组样本数量:\n")
sample_counts <- table(alpha_df$Group)
print(sample_counts)
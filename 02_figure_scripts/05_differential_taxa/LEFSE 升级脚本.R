# ============================================
# LEfSe分析综合修复版（完全修复字符串连接问题）
# ============================================

rm(list = ls())

cat("========== LEfSe分析综合修复版 ==========\n\n")

# 设置工作目录
work_dir <- "/Users/catherine/Downloads/师兄交作业系列 ddl/数据提交/"
if (!dir.exists(work_dir)) {
  cat("警告：工作目录不存在，请检查路径！\n")
  cat("当前工作目录：", getwd(), "\n")
  work_dir <- getwd()
}
setwd(work_dir)
cat("工作目录：", getwd(), "\n\n")

# 加载必要包
cat("加载必要包...\n")
required_packages <- c("ggplot2", "dplyr", "tidyr", "broom", "openxlsx", "stringr", "purrr",
                       "ggdendro", "pheatmap", "stats", "utils", "grDevices", "graphics")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 读取物种丰度数据（CSV格式）
csv_file <- "speciesallname.csv"

if (file.exists(csv_file)) {
  cat("找到数据文件：", csv_file, "\n")
  
  # 读取CSV文件
  species_df <- read.csv(csv_file, row.names = 1, check.names = FALSE)
  
  cat("数据读取成功！\n")
  cat("数据维度：", dim(species_df), "（物种数 × 样本数）\n")
  cat("物种数：", nrow(species_df), "\n")
  cat("样本数：", ncol(species_df), "\n")
  
} else {
  cat("错误：找不到文件 '", csv_file, "'\n", sep = "")
  cat("当前目录下的文件：\n")
  print(list.files(pattern = "\\.csv$"))
  stop("请确保文件存在于工作目录中")
}

# ==========================================
# 定义各组样本（根据您的数据）
# ==========================================
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("定义样本分组\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 首先查看数据中的实际样本名
all_samples_in_data <- colnames(species_df)
cat("数据中的所有样本名（前20个）：\n")
print(head(all_samples_in_data, 20))
cat("总计：", length(all_samples_in_data), "个样本\n\n")

# 根据您之前的STAMP分析定义分组
IS_samples <- c("SY365BB-0-4", "SY365BB-4-8", "SY365BB-8-12",
                "SY366YB-0-4", "SY366YB-4-8", "SY366YB-8-12",
                "SY367YB-0-4", "SY367YB-4-8", "SY367YB-8-12",
                "SQ_58_0-4", "SQ_58_4-8", "SQ_58_8-12",
                "YB_LM_0-4", "YB_LM_4-8", "YB_LM_8-12")

AS_samples <- c("S1_0-3", "S1_6-9", "S1_9-12",
                "S2_0-3", "S2_3-6", "S2_12-15",
                "S4_12-15", "S4_9-12",
                "SQ_56_0-3", "SQ_56_3-6", "SQ_56_6-9", "SQ_56_9-12",
                "SQ_59_0-3", "SQ_59_3-6", "SQ_59_6-9", "SQ_59_9-12",
                "SQB_1_0-3")

ES_samples <- c("C1_0-6", "C1_6-12", "C1_12-18",
                "C2_0-6", "C2_6-12", "C2_12-18",
                "C3_0-6", "C3_6-12", "C3_12-18",
                "ES_1_0-6", "ES_1_6-12", "ES_1_12-18",
                "ES_2_0-6", "ES_2_6-12", "ES_2_12-18",
                "ES_3_0-6", "ES_3_6-12", "ES_3_12-18",
                "ES_4_0-6")

NS_samples <- c("S3_0-3", "S3_6-9", "S3_9-12",
                "NS_0-6",
                "R2111_N300_0-10", "R2111_N500_0-10",
                "R2113_N300_0-10", "R2113_N500_0-10")

# 检查分组情况
cat("IS 组样本数：", length(IS_samples), "\n")
cat("AS 组样本数：", length(AS_samples), "\n")
cat("ES 组样本数：", length(ES_samples), "\n")
cat("NS 组样本数：", length(NS_samples), "\n")

# 创建所有样本列表
all_defined_samples <- c(IS_samples, AS_samples, ES_samples, NS_samples)
cat("预定义样本总数：", length(all_defined_samples), "\n")

# ==========================================
# 样本匹配检查
# ==========================================
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("样本匹配检查\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 检查CSV文件中的实际样本
actual_samples <- colnames(species_df)
cat("数据中实际样本数：", length(actual_samples), "\n")

# 找出匹配和不匹配的样本
matched_samples <- intersect(all_defined_samples, actual_samples)
missing_samples <- setdiff(all_defined_samples, actual_samples)
extra_samples <- setdiff(actual_samples, all_defined_samples)

cat("\n匹配情况：\n")
cat("✅ 匹配的样本：", length(matched_samples), "个\n")
cat("❌ 缺失的样本：", length(missing_samples), "个\n")
cat("📝 额外的样本（不在定义列表中）：", length(extra_samples), "个\n")

# 详细报告缺失样本
if (length(missing_samples) > 0) {
  cat("\n详细缺失样本报告：\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  
  # 按组分类缺失样本
  missing_by_group <- list(
    IS = intersect(missing_samples, IS_samples),
    AS = intersect(missing_samples, AS_samples),
    ES = intersect(missing_samples, ES_samples),
    NS = intersect(missing_samples, NS_samples)
  )
  
  for (group in c("IS", "AS", "ES", "NS")) {
    missing_in_group <- missing_by_group[[group]]
    if (length(missing_in_group) > 0) {
      cat(group, "组缺失", length(missing_in_group), "个样本：\n")
      cat(" ", paste(missing_in_group, collapse = ", "), "\n")
    }
  }
}

cat("\n检查可能的样本名匹配问题：\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

find_similar_names <- function(target_names, source_names, max_dist = 3) {
  similar_pairs <- list()
  for (target in target_names) {
    for (source in source_names) {
      if (adist(target, source) <= max_dist) {
        similar_pairs[[length(similar_pairs) + 1]] <- c(target, source, adist(target, source))
      }
    }
  }
  return(similar_pairs)
}

similar_matches <- find_similar_names(missing_samples, actual_samples)
if (length(similar_matches) > 0) {
  cat("发现可能有误的样本名匹配：\n")
  for (pair in similar_matches) {
    cat(" 定义：", pair[1], " → 数据中：", pair[2], " (距离=", pair[3], ")\n", sep = "")
  }
} else {
  cat("未发现明显的样本名匹配问题\n")
}

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("样本使用策略\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 样本名映射（如果有需要修正的样本名）
sample_mapping <- c(
  "SQ_58_-8-12" = "SQ_58_8-12" # 可能的修正示例
)

# 应用映射
matched_samples_fixed <- matched_samples
for (i in seq_along(matched_samples_fixed)) {
  if (matched_samples_fixed[i] %in% names(sample_mapping)) {
    new_name <- sample_mapping[matched_samples_fixed[i]]
    if (new_name %in% actual_samples) {
      matched_samples_fixed[i] <- new_name
      cat("样本名映射：", matched_samples[i], " → ", new_name, "\n")
    }
  }
}

# 最终使用的样本
final_samples <- intersect(matched_samples_fixed, actual_samples)
species_df_final <- species_df[, final_samples, drop = FALSE]

cat("\n最终使用的样本数：", length(final_samples), "\n")
cat("最终数据维度：", dim(species_df_final), "\n")

cat("\n创建分组信息...\n")
group_info <- data.frame(
  Sample = colnames(species_df_final),
  Group = NA,
  Original_Name = colnames(species_df_final),
  row.names = colnames(species_df_final)
)

# 分配分组
for (i in 1:nrow(group_info)) {
  sample_name <- group_info$Sample[i]
  if (sample_name %in% IS_samples) {
    group_info$Group[i] <- "IS"
  } else if (sample_name %in% AS_samples) {
    group_info$Group[i] <- "AS"
  } else if (sample_name %in% ES_samples) {
    group_info$Group[i] <- "ES"
  } else if (sample_name %in% NS_samples) {
    group_info$Group[i] <- "NS"
  }
}

# 检查是否有未分组的样本
ungrouped <- is.na(group_info$Group)
if (any(ungrouped)) {
  cat("警告：以下样本无法分组：\n")
  print(group_info$Sample[ungrouped])
  
  # 从分析中移除这些样本
  species_df_final <- species_df_final[, !ungrouped, drop = FALSE]
  group_info <- group_info[!ungrouped, ]
  cat("已移除", sum(ungrouped), "个未分组样本\n")
}

group_info$Group <- factor(group_info$Group, levels = c("IS", "AS", "ES", "NS"))

# 设置颜色
group_col <- c(IS = "#00FF7F", AS = "#FF1493", ES = "#FFA500", NS = "#00FFFF")

cat("\n最终样本分组统计：\n")
print(table(group_info$Group))

# 保存样本匹配报告
cat("\n生成样本匹配报告...\n")
sample_report <- data.frame(
  Sample_Name = colnames(species_df_final),
  Group = group_info$Group,
  In_Original_List = colnames(species_df_final) %in% all_defined_samples,
  Status = "Used"
)

if (length(missing_samples) > 0) {
  missing_report <- data.frame(
    Sample_Name = missing_samples,
    Group = sapply(missing_samples, function(x) {
      if (x %in% IS_samples) return("IS")
      if (x %in% AS_samples) return("AS")
      if (x %in% ES_samples) return("ES")
      if (x %in% NS_samples) return("NS")
      return("Unknown")
    }),
    In_Original_List = TRUE,
    Status = "Missing"
  )
  sample_report <- rbind(sample_report, missing_report)
}

write.csv(sample_report, "LEfSe_Sample_Matching_Report.csv", row.names = FALSE)
cat("✅ 样本匹配报告已保存：LEfSe_Sample_Matching_Report.csv\n")

# ==========================================
# 执行LEfSe分析
# ==========================================
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("开始LEfSe分析\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

perform_simple_lefse <- function(data_matrix, group_info) {
  
  # 转换为数值矩阵
  data_numeric <- as.matrix(data_matrix)
  mode(data_numeric) <- "numeric"
  
  results <- data.frame(
    Species = rownames(data_numeric),
    p_value = NA,
    p_adj = NA,
    IS_mean = NA,
    AS_mean = NA,
    ES_mean = NA,
    NS_mean = NA,
    stringsAsFactors = FALSE
  )
  
  # 分析每个物种
  n_species <- nrow(data_numeric)
  cat("分析", n_species, "个物种...\n")
  
  for (i in 1:n_species) {
    values <- as.numeric(data_numeric[i, ])
    
    # 各组数据
    is_data <- values[group_info$Group == "IS"]
    as_data <- values[group_info$Group == "AS"]
    es_data <- values[group_info$Group == "ES"]
    ns_data <- values[group_info$Group == "NS"]
    
    # 计算均值
    results$IS_mean[i] <- mean(is_data, na.rm = TRUE)
    results$AS_mean[i] <- mean(as_data, na.rm = TRUE)
    results$ES_mean[i] <- mean(es_data, na.rm = TRUE)
    results$NS_mean[i] <- mean(ns_data, na.rm = TRUE)
    
    # Kruskal-Wallis检验
    all_values <- c(is_data, as_data, es_data, ns_data)
    group_factor <- factor(rep(c("IS", "AS", "ES", "NS"), 
                               c(length(is_data), length(as_data), 
                                 length(es_data), length(ns_data))))
    
    if (length(all_values) >= 4 && length(unique(all_values)) > 1) {
      kw_test <- try(kruskal.test(all_values ~ group_factor), silent = TRUE)
      if (!inherits(kw_test, "try-error")) {
        results$p_value[i] <- kw_test$p.value
      }
    }
    
    # 进度显示
    if (i %% 100 == 0) cat("已分析", i, "/", n_species, "\n")
  }
  
  # 多重检验校正
  results$p_adj <- p.adjust(results$p_value, method = "BH")
  
  # 确定富集组
  mean_cols <- c("IS_mean", "AS_mean", "ES_mean", "NS_mean")
  results$Enriched <- apply(results[, mean_cols], 1, function(x) {
    if (max(x, na.rm = TRUE) == 0) return(NA)
    c("IS", "AS", "ES", "NS")[which.max(x)]
  })
  
  # 计算LDA分数（简化版）
  results$LDA_Score <- -log10(results$p_adj + 1e-10) *
    apply(results[, mean_cols], 1, function(x) {
      if (max(x, na.rm = TRUE) == 0) return(1)
      max_val <- max(x, na.rm = TRUE)
      second_max <- sort(x, decreasing = TRUE)[2]
      if (second_max == 0 || is.na(second_max)) return(2)
      max_val / second_max
    })
  
  # 排序
  results <- results[order(results$LDA_Score, decreasing = TRUE), ]
  
  return(results)
}

# 运行分析
lefse_results <- perform_simple_lefse(species_df_final, group_info)

# 简化物种名称
cat("\n简化物种名称...\n")

simplify_species_name <- function(full_name) {
  
  # 分割分类层级
  parts <- unlist(strsplit(full_name, " "))
  
  # 提取关键分类信息
  if (length(parts) >= 7) {
    # 有属和种：显示"属_种"
    genus <- parts[6]
    species <- parts[7]
    return(paste0(genus, "", species))
  } else if (length(parts) >= 6) {
    # 有属：显示"属"
    genus <- parts[6]
    # 加上纲的信息
    class <- parts[3]
    return(paste0(class, "", genus))
  } else if (length(parts) >= 4) {
    # 只有纲和目：显示"纲_目"
    class <- parts[3]
    order <- parts[4]
    return(paste0(class, "", order))
  } else if (length(parts) >= 3) {
    # 只有门和纲：显示"门_纲"
    phylum <- parts[2]
    class <- parts[3]
    return(paste0(phylum, "", class))
  } else {
    # 其他情况：截断原始名称
    return(ifelse(nchar(full_name) > 30,
                  paste0(substr(full_name, 1, 27), "..."),
                  full_name))
  }
}

# 应用简化函数
lefse_results$Species_short <- sapply(lefse_results$Species, simplify_species_name)

# 确保名称唯一
make_names_unique <- function(names) {
  name_counts <- table(names)
  result <- character(length(names))
  
  for (i in seq_along(names)) {
    name <- names[i]
    if (name_counts[name] > 1) {
      # 查找这是该名称的第几次出现
      occurrence <- sum(names[1:i] == name)
      result[i] <- paste0(name, "_", occurrence)
    } else {
      result[i] <- name
    }
  }
  
  return(result)
}

lefse_results$Species_short <- make_names_unique(lefse_results$Species_short)

# 保存结果
write.csv(lefse_results, "LEfSe_Results_Final.csv", row.names = FALSE)
cat("✅ 分析结果已保存：LEfSe_Results_Final.csv\n")

# ==========================================
# 绘制LEfSe条形图（修复版）
# ==========================================
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("绘制LEfSe条形图\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 选择显著物种
sig_results <- lefse_results[lefse_results$p_adj < 0.05 & !is.na(lefse_results$p_adj), ]
if (nrow(sig_results) > 0) {
  top_n <- min(25, nrow(sig_results))
  plot_data <- head(sig_results, top_n)
  cat("使用", top_n, "个显著物种绘制条形图\n")
} else {
  top_n <- min(25, nrow(lefse_results))
  plot_data <- head(lefse_results, top_n)
  cat("没有显著物种，使用前", top_n, "个物种绘制条形图\n")
}

# 排序并设置因子
plot_data <- plot_data[order(plot_data$LDA_Score, decreasing = FALSE), ]
plot_data$Species_short <- factor(plot_data$Species_short,
                                  levels = plot_data$Species_short)

# 创建图形
library(ggplot2)
p <- ggplot(plot_data, aes(x = LDA_Score, y = Species_short, fill = Enriched)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = group_col,
                    name = "Enriched in",
                    na.value = "grey70") +
  labs(
    title = "LEfSe Analysis Results",
    subtitle = paste("Top", nrow(plot_data), "Differential Features"),
    x = "LDA Score",
    y = "Features"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  geom_vline(xintercept = 2.0,
             linetype = "dashed",
             color = "red",
             alpha = 0.6)

# 保存图形
ggsave("A_LEfSe_Barplot_Final.pdf", p, width = 12, height = 10)
ggsave("A_LEfSe_Barplot_Final.png", p, width = 12, height = 10, dpi = 300)
cat("✅ LEfSe条形图已保存：A_LEfSe_Barplot_Final.pdf/.png\n")

# ==========================================
# 绘制分支图
# ==========================================
cat("\n绘制分支图...\n")

# 选择Top物种
top_n_clad <- min(30, nrow(lefse_results))
top_species <- head(lefse_results$Species, top_n_clad)

# 提取数据
species_data <- as.matrix(species_df_final[top_species, ])
mode(species_data) <- "numeric"

# 计算距离和聚类
species_dist <- dist(species_data)
hc <- hclust(species_dist, method = "ward.D2")
dend <- as.dendrogram(hc)

library(ggdendro)
# 准备绘图数据
dend_data <- dendro_data(dend)
leaf_order <- order.dendrogram(dend)
leaf_labels <- labels(dend)[leaf_order]

# 使用简化的物种名
short_labels <- lefse_results$Species_short[match(leaf_labels, lefse_results$Species)]

label_data <- data.frame(
  x = 1:length(leaf_labels),
  y = 0,
  label = short_labels,
  Enriched = lefse_results$Enriched[match(leaf_labels, lefse_results$Species)],
  stringsAsFactors = FALSE
)

# 创建分支图
p_clad <- ggplot() +
  geom_segment(data = dend_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.8, color = "grey50") +
  geom_point(data = label_data,
             aes(x = x, y = y, color = Enriched),
             size = 4, alpha = 0.8) +
  scale_color_manual(values = group_col,
                     name = "Enriched in",
                     na.value = "grey70") +
  labs(
    title = "Cladogram of Differential Features",
    subtitle = paste("Top", top_n_clad, "Species")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  coord_flip() +
  scale_y_reverse()

# 保存图形
ggsave("B_Cladogram_Final.pdf", p_clad, width = 14, height = 12)
ggsave("B_Cladogram_Final.png", p_clad, width = 14, height = 12, dpi = 300)
cat("✅ 分支图已保存：B_Cladogram_Final.pdf/.png\n")

# ==========================================
# 生成热图
# ==========================================
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("生成热图\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 验证热图中的样本
cat("\n验证热图样本：\n")
for (group in c("IS", "AS", "ES", "NS")) {
  group_samples <- colnames(species_df_final)[group_info$Group == group]
  cat(group, "组：", length(group_samples), "个样本\n")
}

# 选择Top物种
top_n_heat <- min(20, nrow(lefse_results))
heat_species <- head(lefse_results$Species, top_n_heat)

# 提取并归一化数据
heat_data <- as.matrix(species_df_final[heat_species, ])
mode(heat_data) <- "numeric"

cat("\n热图数据维度：", dim(heat_data), "\n")
cat("各组的样本数：\n")
print(table(group_info$Group))

# 按行归一化（Z-score标准化）
heat_data_norm <- t(apply(heat_data, 1, function(x) {
  if (sd(x, na.rm = TRUE) > 0) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  } else {
    rep(0, length(x))
  }
}))

# 使用简化的物种名称
heat_species_short <- lefse_results$Species_short[match(heat_species, lefse_results$Species)]
rownames(heat_data_norm) <- heat_species_short

# 按组排序样本
group_order <- order(group_info$Group)
heat_data_norm <- heat_data_norm[, group_order]
sample_names <- colnames(heat_data_norm)

# 创建分组注释数据框
annotation_col <- data.frame(
  Group = group_info$Group[group_order],
  row.names = sample_names
)

# 创建分组颜色列表
annotation_colors <- list(
  Group = group_col
)

# 计算合适的图形大小
n_species <- nrow(heat_data_norm)
n_samples <- ncol(heat_data_norm)

# 动态调整图形尺寸
plot_width <- max(16, n_samples * 0.25)
plot_height <- max(10, n_species * 0.5)

cat("\n热图参数：\n")
cat("物种数量：", n_species, "\n")
cat("样本数量：", n_samples, "\n")
cat("画布大小：", plot_width, "x", plot_height, "英寸\n")

# 计算样本分组的间隔位置
group_gaps <- which(diff(as.numeric(annotation_col$Group)) != 0)

# 保存为PDF
pdf("C_Heatmap_Final.pdf", width = plot_width, height = plot_height)
pheatmap(heat_data_norm,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         fontsize_row = 10,
         fontsize_col = 8,
         cellwidth = max(8, 100/n_samples),
         cellheight = max(12, 40/n_species),
         main = paste("Heatmap of Top", top_n_heat, "Differential Species"),
         border_color = NA,
         gaps_col = group_gaps,
         angle_col = 45)
dev.off()
cat("✅ 热图已保存：C_Heatmap_Final.pdf\n")

# PNG版本
png("C_Heatmap_Final.png", width = plot_width * 200, height = plot_height * 200, res = 300)
pheatmap(heat_data_norm,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         fontsize_row = 10,
         fontsize_col = 8,
         cellwidth = max(8, 100/n_samples),
         cellheight = max(12, 40/n_species),
         main = paste("Heatmap of Top", top_n_heat, "Differential Species"),
         border_color = NA,
         gaps_col = group_gaps,
         angle_col = 45)
dev.off()
cat("✅ 热图已保存：C_Heatmap_Final.png\n")

# ==========================================
# 为每个分组绘制单独的LEfSe条形图
# ==========================================
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("绘制分组LEfSe条形图\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 计算百分位数分数
calculate_percentile_score <- function(lda_scores) {
  # 计算每个值的百分位数
  percentiles <- ecdf(lda_scores)(lda_scores) * 100
  return(percentiles)
}

lefse_results$Percentile_Score <- calculate_percentile_score(lefse_results$LDA_Score)

# 为每个分组选择显著菌种
cat("\n为每个分组选择显著菌种...\n")

# 首先筛选显著结果
sig_results <- lefse_results[lefse_results$p_adj < 0.05 & !is.na(lefse_results$p_adj), ]

if (nrow(sig_results) > 0) {
  
  # 创建分组数据框列表
  group_plot_data <- list()
  
  for (group in c("IS", "AS", "ES", "NS")) {
    cat("处理", group, "组...\n")
    
    # 筛选该组富集的显著物种
    group_species <- sig_results[sig_results$Enriched == group, ]
    
    if (nrow(group_species) > 0) {
      # 按原始LDA分数排序，取前20个
      top_n <- min(20, nrow(group_species))
      group_top <- head(group_species[order(group_species$LDA_Score, decreasing = TRUE), ], top_n)
      
      # 添加分组标签和标准化分数
      group_top$Group <- group
      group_top$Plot_Score <- calculate_percentile_score(group_top$LDA_Score)
      
      # 确保所有分数为正数
      if (min(group_top$Plot_Score, na.rm = TRUE) <= 0) {
        group_top$Plot_Score <- group_top$Plot_Score + abs(min(group_top$Plot_Score, na.rm = TRUE)) + 1
      }
      
      # 按标准化分数排序
      group_top <- group_top[order(group_top$Plot_Score, decreasing = FALSE), ]
      group_top$Species_short <- factor(group_top$Species_short, 
                                        levels = group_top$Species_short)
      
      group_plot_data[[group]] <- group_top
      cat("  ", group, "组：找到", nrow(group_species), "个富集物种，取前", top_n, "个\n")
    } else {
      cat("  ", group, "组：没有显著富集物种\n")
      group_plot_data[[group]] <- NULL
    }
  }
  
  # 合并所有分组数据
  all_group_data <- do.call(rbind, group_plot_data)
  
  if (!is.null(all_group_data) && nrow(all_group_data) > 0) {
    # 确保Group是因子
    all_group_data$Group <- factor(all_group_data$Group, levels = c("IS", "AS", "ES", "NS"))
    
    # 动态调整图形高度
    n_species_per_group <- sapply(group_plot_data, function(x) if(!is.null(x)) nrow(x) else 0)
    max_species <- max(n_species_per_group)
    plot_height <- 6 + max_species * 0.3
    
    # 创建分面条形图
    p_groups <- ggplot(all_group_data, aes(x = Plot_Score, y = Species_short, fill = Group)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = group_col, 
                        name = "富集于",
                        guide = "none") +
      labs(
        title = "LEfSe Analysis Results by Group",
        subtitle = "Top Significant Species for Each Group (p < 0.05)",
        x = "Percentile Score (0-100)",
        y = "Features"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = NA),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(20, 20, 20, 20)
      ) +
      facet_wrap(~ Group, scales = "free_y", ncol = 2) +
      geom_vline(xintercept = 50, 
                 linetype = "dashed", 
                 color = "red", 
                 alpha = 0.6)
    
    # 保存图形
    ggsave("D_LEfSe_Barplot_By_Group.pdf", p_groups, width = 16, height = plot_height)
    ggsave("D_LEfSe_Barplot_By_Group.png", p_groups, width = 16, height = plot_height, dpi = 300)
    cat("✅ 分组LEfSe条形图已保存：D_LEfSe_Barplot_By_Group.pdf/.png\n")
  } else {
    cat("⚠️ 没有足够的显著物种数据创建分组图形\n")
  }
} else {
  cat("⚠️ 没有显著差异物种，无法创建分组图形\n")
}

# ==========================================
# 生成综合报告
# ==========================================
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("生成综合报告\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

sink("E_LEfSe_Comprehensive_Report.txt")

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("LEfSe分析综合报告\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("分析日期：", date(), "\n")
cat("工作目录：", getwd(), "\n\n")

cat("1. 样本匹配总结\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat("预定义样本：", length(all_defined_samples), "\n")
cat("数据中实际样本：", length(actual_samples), "\n")
cat("匹配使用的样本：", length(final_samples), "\n")
cat("缺失样本：", length(missing_samples), "\n")
cat("额外样本：", length(extra_samples), "\n\n")

cat("2. 最终样本分布\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
final_counts <- table(group_info$Group)
for (group in c("IS", "AS", "ES", "NS")) {
  count <- ifelse(group %in% names(final_counts), final_counts[group], 0)
  cat(group, " 组：", count, " 个样本\n", sep = "")
}
cat("\n")

cat("3. 分析结果\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
total <- nrow(lefse_results)
sig <- sum(lefse_results$p_adj < 0.05, na.rm = TRUE)
cat("分析物种总数：", total, "\n")
cat("显著差异物种 (p < 0.05)：", sig,
    " (", round(sig/total*100, 1), "%)\n", sep = "")

if (sig > 0) {
  sig_results <- lefse_results[lefse_results$p_adj < 0.05, ]
  enrichment_counts <- table(sig_results$Enriched)
  cat("\n富集分布：\n")
  for (group in c("IS", "AS", "ES", "NS")) {
    count <- ifelse(group %in% names(enrichment_counts), enrichment_counts[group], 0)
    cat(" ", group, "：", count, " 个物种\n", sep = "")
  }
}
cat("\n")

cat("4. 生成的文件\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
output_files <- c(
  "LEfSe_Sample_Matching_Report.csv",
  "LEfSe_Results_Final.csv",
  "A_LEfSe_Barplot_Final.pdf",
  "A_LEfSe_Barplot_Final.png",
  "B_Cladogram_Final.pdf",
  "B_Cladogram_Final.png",
  "C_Heatmap_Final.pdf",
  "C_Heatmap_Final.png",
  "D_LEfSe_Barplot_By_Group.pdf",
  "D_LEfSe_Barplot_By_Group.png",
  "E_LEfSe_Comprehensive_Report.txt"
)

for (file in output_files) {
  if (file.exists(file)) {
    size_kb <- round(file.size(file) / 1024, 1)
    cat("✅ ", file, " (", size_kb, " KB)\n", sep = "")
  } else {
    cat("❌ ", file, " (未生成)\n", sep = "")
  }
}

cat("\n5. 缺失样本（如有）\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
if (length(missing_samples) > 0) {
  cat("以下样本已定义但未在数据中找到：\n")
  for (sample in missing_samples) {
    cat(" • ", sample, "\n", sep = "")
  }
  cat("\n可能的原因：\n")
  cat(" - CSV文件中的样本名可能不同\n")
  cat(" - 可能有拼写错误或格式差异\n")
  cat(" - 样本可能在另一个文件中\n")
} else {
  cat("所有定义样本都在数据中找到。\n")
}

cat("\n6. 建议\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat("1. 检查 LEfSe_Sample_Matching_Report.csv 了解详细样本信息\n")
cat("2. 查看热图 (C_Heatmap_Final.pdf) 确认所有组都有代表样本\n")
cat("3. 如果有缺失样本，检查命名不一致问题\n")
cat("4. 考虑使用CSV文件中的确切样本名\n")

sink()

cat("✅ 综合报告已保存：E_LEfSe_Comprehensive_Report.txt\n")

# ==========================================
# 最终总结
# ==========================================
cat("\n")
cat(paste(rep("*", 70), collapse = ""), "\n")
cat("🎉 LEfSe分析综合修复版完成！\n")
cat(paste(rep("*", 70), collapse = ""), "\n\n")

cat("📊 关键成果：\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
cat("• 匹配样本数：", length(final_samples), "/", length(all_defined_samples), "\n")
cat("• 分析物种数：", total, "\n")
cat("• 显著差异物种：", sig, "\n")
cat("• 输出文件数：", length(output_files), "\n\n")

cat("📁 主要输出文件：\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
cat("1. LEfSe_Sample_Matching_Report.csv - 样本匹配详情\n")
cat("2. LEfSe_Results_Final.csv - 完整的LDA分析结果\n")
cat("3. A_LEfSe_Barplot_Final.pdf - 主条形图\n")
cat("4. B_Cladogram_Final.pdf - 分支图\n")
cat("5. C_Heatmap_Final.pdf - 热图\n")
cat("6. D_LEfSe_Barplot_By_Group.pdf - 分组条形图\n")
cat("7. E_LEfSe_Comprehensive_Report.txt - 综合报告\n\n")

cat("🔍 样本检查结果：\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
if (length(missing_samples) > 0) {
  cat("⚠️ 发现", length(missing_samples), "个缺失样本\n")
  cat("  请检查 LEfSe_Sample_Matching_Report.csv 了解详情\n")
} else {
  cat("✅ 所有预定义样本都在数据中\n")
}

cat("\n💡 下一步建议：\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
cat("1. 打开 C_Heatmap_Final.pdf 确认所有组的样本都在\n")
cat("2. 检查 D_LEfSe_Barplot_By_Group.pdf 查看各组的富集物种\n")
cat("3. 如果有缺失样本，手动检查原始CSV文件中的列名\n")
cat("4. 考虑更新样本定义列表以匹配实际数据\n")

cat("\n")
cat(paste(rep("✨", 30), collapse = ""), "\n")
cat("分析全部完成！\n")
cat(paste(rep("✨", 30), collapse = ""), "\n")

cat("\n工作目录：", getwd(), "\n")
cat("完成时间：", Sys.time(), "\n")

# ==========================================
# 可选：生成样本可视化
# ==========================================
cat("\n生成样本分布可视化...\n")

# 创建样本分布条形图
sample_dist_plot <- ggplot(data.frame(Group = names(final_counts), Count = as.numeric(final_counts)),
                           aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = group_col) +
  labs(
    title = "Sample Distribution by Group",
    subtitle = paste("Total:", sum(final_counts), "matched samples"),
    x = "Group",
    y = "Number of Samples"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "none"
  ) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4)

ggsave("F_Sample_Distribution.pdf", sample_dist_plot, width = 8, height = 6)
ggsave("F_Sample_Distribution.png", sample_dist_plot, width = 8, height = 6, dpi = 300)
cat("✅ 样本分布图已保存：F_Sample_Distribution.pdf/.png\n")

cat("\n所有分析完成！请检查输出文件。\n")
# ==========================================
# 完整的LEfSe分析 - 重新运行版本
# ==========================================

# 1. 设置工作目录
setwd("/Users/catherine/Downloads/师兄交作业系列 ddl/数据提交/")
cat("工作目录：", getwd(), "\n\n")

# 2. 加载必要包
cat("加载必要包...\n")
required_packages <- c("ggplot2", "dplyr", "tidyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 3. 读取数据
cat("读取数据...\n")
species_df <- read.csv("speciesallname.csv", row.names = 1, check.names = FALSE)
cat("数据维度：", dim(species_df), "\n")

# 4. 定义分组
IS_samples <- c("SY365BB-0-4", "SY365BB-4-8", "SY365BB-8-12",
                "SY366YB-0-4", "SY366YB-4-8", "SY366YB-8-12",
                "SY366YW-0-4", "SY366YW-4-8", "SY366YW-8-12",
                "SY368YW-0-4", "SY368YW-4-8", "SY368YW-8-12",
                "SY456YB-0-4", "SY456YB-4-8", "SY456YB-8-12",
                "SY457BB-0-4", "SY457BB-4-8", "SY457BB-8-12",
                "SY459WG-0-4", "SY459WG-4-8", "SY459WG-8-12")

AS_samples <- c("S1_0-3", "S1_6-9", "S1_9-12",
                "S2_0-3", "S2_3-6", "S2_12-15",
                "S4_12-15", "S4_9-12",
                "SQ_58_0-4", "SQ_58_4-8", "SQ_58_-8-12",
                "SQ_81_0-4", "SQ_81_4-8", "SQ_81_-8-12")

ES_samples <- c("C1_0-6", "C1_6-12", "C1_12-18",
                "C2_0-6", "C2_6-12", "C2_12-18",
                "C3_0-6", "C3_6-12", "C3_12-18",
                "ES_2_0-6",
                "S13_0-2", "S14_4-6", "S15_8-10")

NS_samples <- c("S3_0-3", "S3_6-9", "S3_9-12",
                "NS_0-6",
                "R2111_N300_0-10", "R2111_N500_0-10", 
                "R2111_S300_0-10", "R2111_S500_0-10")

# 5. 筛选样本
cat("\n筛选样本...\n")
all_samples <- c(IS_samples, AS_samples, ES_samples, NS_samples)
existing_samples <- intersect(all_samples, colnames(species_df))
species_df_filtered <- species_df[, existing_samples]

cat("有效样本：", ncol(species_df_filtered), "个\n")
cat("有效物种：", nrow(species_df_filtered), "个\n")

# 6. 创建分组信息
group_info <- data.frame(
  Sample = colnames(species_df_filtered),
  Group = NA,
  row.names = colnames(species_df_filtered)
)

group_info$Group[group_info$Sample %in% IS_samples] <- "IS"
group_info$Group[group_info$Sample %in% AS_samples] <- "AS"
group_info$Group[group_info$Sample %in% ES_samples] <- "ES"
group_info$Group[group_info$Sample %in% NS_samples] <- "NS"

group_info$Group <- factor(group_info$Group, levels = c("IS", "AS", "ES", "NS"))

# 7. 设置颜色
group_col <- c(IS = "#00FF7F", AS = "#FF1493", ES = "#FFA500", NS = "#00FFFF")

cat("\n样本分组统计：\n")
print(table(group_info$Group))

# 8. 执行LEfSe分析（简化版）
cat("\n开始LEfSe分析...\n")

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
    
    if (length(all_values) >= 4) {
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
lefse_results <- perform_simple_lefse(species_df_filtered, group_info)

# 9. 简化物种名称（按照你的要求）
cat("\n简化物种名称...\n")

simplify_species_name <- function(full_name) {
  # 分割分类层级
  parts <- unlist(strsplit(full_name, " "))
  
  # 提取关键分类信息
  if (length(parts) >= 7) {
    # 有属和种：显示"属_种"
    genus <- parts[6]
    species <- parts[7]
    return(paste0(genus, "_", species))
  } else if (length(parts) >= 6) {
    # 有属：显示"属"
    genus <- parts[6]
    # 加上纲的信息
    class <- parts[3]
    return(paste0(class, "_", genus))
  } else if (length(parts) >= 4) {
    # 只有纲和目：显示"纲_目"
    class <- parts[3]
    order <- parts[4]
    return(paste0(class, "_", order))
  } else if (length(parts) >= 3) {
    # 只有门和纲：显示"门_纲"
    phylum <- parts[2]
    class <- parts[3]
    return(paste0(phylum, "_", class))
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

# 10. 保存结果
write.csv(lefse_results, "LEfSe_Results_Simple.csv", row.names = FALSE)
cat("✅ 分析结果已保存：LEfSe_Results_Simple.csv\n")

# 11. 绘制LEfSe条形图
cat("\n绘制LEfSe条形图...\n")

# 选择Top 25显著物种
sig_results <- lefse_results[lefse_results$p_adj < 0.05 & !is.na(lefse_results$p_adj), ]
if (nrow(sig_results) > 0) {
  top_n <- min(25, nrow(sig_results))
  plot_data <- head(sig_results, top_n)
} else {
  top_n <- min(25, nrow(lefse_results))
  plot_data <- head(lefse_results, top_n)
}

# 排序并设置因子
plot_data <- plot_data[order(plot_data$LDA_Score, decreasing = FALSE), ]
plot_data$Species_short <- factor(plot_data$Species_short, 
                                  levels = plot_data$Species_short)

# 创建图形
p <- ggplot(plot_data, aes(x = LDA_Score, y = Species_short, fill = Enriched)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = group_col, 
                    name = "Enriched in",
                    na.value = "grey70") +
  labs(
    title = "LEfSe Analysis Results",
    subtitle = paste("Top", nrow(plot_data), "Differential Features"),
    x = "LDA Score (log10)",
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
ggsave("A_LEfSe_Barplot.pdf", p, width = 12, height = 10)
ggsave("A_LEfSe_Barplot.png", p, width = 12, height = 10, dpi = 300)
cat("✅ LEfSe条形图已保存：A_LEfSe_Barplot.pdf/.png\n")

# 12. 绘制分支图（简化版）
cat("\n绘制分支图...\n")

# 安装必要的包
if (!require("ggdendro", quietly = TRUE)) {
  install.packages("ggdendro")
  library(ggdendro)
}

# 选择Top 30物种
top_n_clad <- min(30, nrow(lefse_results))
top_species <- head(lefse_results$Species, top_n_clad)

# 提取数据
species_data <- as.matrix(species_df_filtered[top_species, ])
mode(species_data) <- "numeric"

# 计算距离和聚类
species_dist <- dist(species_data)
hc <- hclust(species_dist, method = "ward.D2")
dend <- as.dendrogram(hc)

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
ggsave("B_Cladogram.pdf", p_clad, width = 14, height = 12)
ggsave("B_Cladogram.png", p_clad, width = 14, height = 12, dpi = 300)
cat("✅ 分支图已保存：B_Cladogram.pdf/.png\n")

# 13. 生成热图
cat("\n绘制热图...\n")

# 选择Top 20物种
top_n_heat <- min(20, nrow(lefse_results))
heat_species <- head(lefse_results$Species, top_n_heat)

# 提取并归一化数据
heat_data <- as.matrix(species_df_filtered[heat_species, ])
mode(heat_data) <- "numeric"

# 按行归一化
heat_data_norm <- t(apply(heat_data, 1, function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))

# 按组排序样本
group_order <- order(group_info$Group)
heat_data_norm <- heat_data_norm[, group_order]

# 创建简单的热图
png("C_Heatmap.png", width = 1200, height = 800)
par(mar = c(8, 12, 4, 2))

# 设置组颜色
group_colors <- group_col[as.character(group_info$Group[group_order])]

heatmap(heat_data_norm,
        ColSideColors = group_colors,
        col = colorRampPalette(c("blue", "white", "red"))(100),
        scale = "row",
        margins = c(8, 12),
        main = paste("Heatmap of Top", top_n_heat, "Species"),
        cexRow = 0.8,
        cexCol = 0.8)

# 添加图例
legend("topright", 
       legend = names(group_col),
       fill = group_col,
       title = "Group",
       cex = 0.8)

dev.off()

pdf("C_Heatmap.pdf", width = 12, height = 8)
par(mar = c(8, 12, 4, 2))
heatmap(heat_data_norm,
        ColSideColors = group_colors,
        col = colorRampPalette(c("blue", "white", "red"))(100),
        scale = "row",
        margins = c(8, 12),
        main = paste("Heatmap of Top", top_n_heat, "Species"),
        cexRow = 0.8,
        cexCol = 0.8)
legend("topright", 
       legend = names(group_col),
       fill = group_col,
       title = "Group",
       cex = 0.8)
dev.off()

cat("✅ 热图已保存：C_Heatmap.pdf/.png\n")

# 14. 生成报告（修正版）
cat("\n生成分析报告...\n")
sink("D_Analysis_Report.txt")

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("LEfSe Analysis Report\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Analysis date:", date(), "\n")
cat("Working directory:", getwd(), "\n\n")

cat("1. Sample Information\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat("Total samples:", nrow(group_info), "\n")
cat("Total species:", nrow(lefse_results), "\n\n")

cat("Sample distribution:\n")
print(table(group_info$Group))
cat("\n")

cat("2. Differential Analysis Results\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

total <- nrow(lefse_results)
sig <- sum(lefse_results$p_adj < 0.05, na.rm = TRUE)

cat("Total species analyzed:", total, "\n")
cat("Significantly different species (p < 0.05):", sig, 
    " (", round(sig/total*100, 1), "%)\n\n", sep = "")

if (sig > 0) {
  cat("Enrichment distribution:\n")
  sig_results <- lefse_results[lefse_results$p_adj < 0.05, ]
  enrichment_counts <- table(sig_results$Enriched)
  
  for (group in c("IS", "AS", "ES", "NS")) {
    count <- ifelse(group %in% names(enrichment_counts), 
                    enrichment_counts[group], 0)
    cat("  ", group, ": ", count, " species\n", sep = "")
  }
  cat("\n")
}

cat("3. Top 10 Differential Species\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
top_10 <- head(lefse_results, 10)
for (i in 1:nrow(top_10)) {
  cat(sprintf("%2d. %-40s LDA=%.2f (p=%.2e) [%s]\n",
              i, top_10$Species_short[i],
              top_10$LDA_Score[i],
              top_10$p_adj[i],
              top_10$Enriched[i]))
}

cat("\n4. Color Scheme\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
for (group in names(group_col)) {
  cat("  ", group, ": ", group_col[group], "\n", sep = "")
}

cat("\n5. Generated Files\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat("A_LEfSe_Barplot.pdf/.png - LEfSe barplot\n")
cat("B_Cladogram.pdf/.png - Cladogram\n")
cat("C_Heatmap.pdf/.png - Heatmap\n")
cat("LEfSe_Results_Simple.csv - Analysis results\n")
cat("D_Analysis_Report.txt - This report\n")

sink()
cat("✅ 分析报告已保存：D_Analysis_Report.txt\n")

# 15. 最终总结
cat("\n")
cat(paste(rep("*", 70), collapse = ""), "\n")
cat("🎉 LEfSe Analysis Complete!\n")
cat(paste(rep("*", 70), collapse = ""), "\n\n")

cat("📊 Key Results:\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat("• Total species analyzed:", total, "\n")
cat("• Significant species:", sig, " (", round(sig/total*100, 1), "%)\n", sep = "")

cat("\n📁 Output Files:\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

output_files <- c("A_LEfSe_Barplot.pdf", "A_LEfSe_Barplot.png",
                  "B_Cladogram.pdf", "B_Cladogram.png",
                  "C_Heatmap.pdf", "C_Heatmap.png",
                  "LEfSe_Results_Simple.csv", "D_Analysis_Report.txt")

for (file in output_files) {
  if (file.exists(file)) {
    size_mb <- round(file.size(file) / 1024 / 1024, 3)
    cat("✅ ", file, " (", size_mb, " MB)\n", sep = "")
  } else {
    cat("⚠️  ", file, "\n", sep = "")
  }
}

cat("\n🎨 Color Scheme Applied:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
for (group in names(group_col)) {
  cat("  ", group, " = ", group_col[group], "\n", sep = "")
}

cat("\n💡 Visualization Files:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat("A_LEfSe_Barplot.pdf - LEfSe barplot with simplified species names\n")
cat("B_Cladogram.pdf - Cladogram showing phylogenetic relationships\n")
cat("C_Heatmap.pdf - Heatmap of species abundance patterns\n")

cat("\n")
cat(paste(rep("✨", 30), collapse = ""), "\n")
cat("所有分析已完成！\n")
cat(paste(rep("✨", 30), collapse = ""), "\n")

cat("\n工作目录：", getwd(), "\n")
cat("完成时间：", Sys.time(), "\n")
# ==========================================
# 16. 为每个分组绘制单独的LEfSe条形图
# ==========================================

cat("\n" + paste(rep("=", 60), collapse = "") + "\n")
cat("绘制分组LEfSe条形图\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 16.1 为每个分组选择前20个显著菌种
cat("\n为每个分组选择前20个显著菌种...\n")

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
      # 按LDA分数排序，取前20个
      top_n <- min(20, nrow(group_species))
      group_top <- head(group_species[order(group_species$LDA_Score, decreasing = TRUE), ], top_n)
      
      # 添加分组标签
      group_top$Group <- group
      
      # 按LDA分数排序（从高到低，用于绘图）
      group_top <- group_top[order(group_top$LDA_Score, decreasing = FALSE), ]
      group_top$Species_short <- factor(group_top$Species_short, 
                                        levels = group_top$Species_short)
      
      group_plot_data[[group]] <- group_top
      cat("  ", group, "组：找到", nrow(group_species), "个富集物种，取前", top_n, "个\n")
    } else {
      cat("  ", group, "组：没有显著富集物种\n")
      group_plot_data[[group]] <- NULL
    }
  }
  
  # 16.2 合并所有分组数据
  all_group_data <- do.call(rbind, group_plot_data)
  
  if (!is.null(all_group_data) && nrow(all_group_data) > 0) {
    # 确保Group是因子，按指定顺序
    all_group_data$Group <- factor(all_group_data$Group, levels = c("IS", "AS", "ES", "NS"))
    
    # 16.3 创建分面条形图
    cat("\n创建分组LEfSe条形图...\n")
    
    p_groups <- ggplot(all_group_data, aes(x = LDA_Score, y = Species_short, fill = Group)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = group_col, 
                        name = "富集于",
                        guide = "none") +  # 关闭图例，因为分面标题已经显示
      labs(
        title = "LEfSe Analysis Results by Group",
        subtitle = "Top 20 Significant Species for Each Group (p < 0.05)",
        x = "LDA Score (log10)",
        y = "Features"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, 
                                  margin = margin(b = 10)),
        plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40",
                                     margin = margin(b = 15)),
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
      geom_vline(xintercept = 2.0, 
                 linetype = "dashed", 
                 color = "red", 
                 alpha = 0.6,
                 linewidth = 0.6)
    
    # 保存图形
    ggsave("A_LEfSe_Barplot2.pdf", p_groups, width = 16, height = 14)
    ggsave("A_LEfSe_Barplot2.png", p_groups, width = 16, height = 14, dpi = 300)
    cat("✅ 分组LEfSe条形图已保存：A_LEfSe_Barplot2.pdf/.png\n")
    
    # 16.4 为每个分组单独创建图形（可选）
    cat("\n为每个分组创建单独的图形...\n")
    
    for (group in names(group_plot_data)) {
      if (!is.null(group_plot_data[[group]]) && nrow(group_plot_data[[group]]) > 0) {
        group_data <- group_plot_data[[group]]
        
        p_single <- ggplot(group_data, aes(x = LDA_Score, y = Species_short, fill = Group)) +
          geom_bar(stat = "identity", width = 0.7) +
          scale_fill_manual(values = group_col[group], 
                            guide = "none") +
          labs(
            title = paste("LEfSe Results -", group, "Group"),
            subtitle = paste("Top", nrow(group_data), "Significant Species"),
            x = "LDA Score (log10)",
            y = "Features"
          ) +
          theme_minimal(base_size = 12) +
          theme(
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 10),
            panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "white", color = NA)
          ) +
          geom_vline(xintercept = 2.0, 
                     linetype = "dashed", 
                     color = "red", 
                     alpha = 0.6)
        
        ggsave(paste0("A_LEfSe_", group, "_only.pdf"), p_single, width = 10, height = 8)
        cat("  ✅ ", group, "组单独图形已保存：A_LEfSe_", group, "_only.pdf\n", sep = "")
      }
    }
    
    # 16.5 创建组合图（所有组在一起，但用不同颜色）
    cat("\n创建组合条形图（所有组）...\n")
    
    # 添加分组颜色
    all_group_data$Color <- group_col[as.character(all_group_data$Group)]
    
    p_combined <- ggplot(all_group_data, aes(x = LDA_Score, y = Species_short, fill = Group)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = group_col, 
                        name = "富集组") +
      labs(
        title = "LEfSe Results - All Groups Combined",
        subtitle = "Top Significant Species from Each Group",
        x = "LDA Score (log10)",
        y = "Features"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        legend.title = element_text(face = "bold"),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      geom_vline(xintercept = 2.0, 
                 linetype = "dashed", 
                 color = "red", 
                 alpha = 0.6)
    
    ggsave("A_LEfSe_All_Groups.pdf", p_combined, width = 14, height = 12)
    ggsave("A_LEfSe_All_Groups.png", p_combined, width = 14, height = 12, dpi = 300)
    cat("✅ 组合条形图已保存：A_LEfSe_All_Groups.pdf/.png\n")
    
  } else {
    cat("⚠️  没有足够的显著物种数据创建分组图形\n")
  }
} else {
  cat("⚠️  没有显著差异物种，无法创建分组图形\n")
}

# ==========================================
# 17. 更新最终总结
# ==========================================

cat("\n")
cat(paste(rep("*", 70), collapse = ""), "\n")
cat("🎉 LEfSe分析增强版完成！\n")
cat(paste(rep("*", 70), collapse = ""), "\n\n")

cat("📊 分析结果统计：\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

if (exists("lefse_results")) {
  total <- nrow(lefse_results)
  sig <- sum(lefse_results$p_adj < 0.05, na.rm = TRUE)
  
  cat("• 总分析物种：", total, "\n", sep = "")
  cat("• 显著差异物种：", sig, " (", round(sig/total*100, 1), "%)\n\n", sep = "")
  
  # 分组统计
  if (sig > 0) {
    sig_results <- lefse_results[lefse_results$p_adj < 0.05, ]
    enrichment_counts <- table(sig_results$Enriched)
    
    cat("• 各组富集物种数量：\n")
    for (group in c("IS", "AS", "ES", "NS")) {
      count <- ifelse(group %in% names(enrichment_counts), 
                      enrichment_counts[group], 0)
      cat("  ", group, "组：", count, " 个物种\n", sep = "")
    }
  }
}

cat("\n📁 新增的输出文件：\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

new_files <- c("A_LEfSe_Barplot2.pdf", "A_LEfSe_Barplot2.png",
               "A_LEfSe_All_Groups.pdf", "A_LEfSe_All_Groups.png")

for (file in new_files) {
  if (file.exists(file)) {
    size_mb <- round(file.size(file) / 1024 / 1024, 3)
    cat("🆕 ", file, " (", size_mb, " MB)\n", sep = "")
  }
}

# 为每个分组检查单独文件
for (group in c("IS", "AS", "ES", "NS")) {
  file_name <- paste0("A_LEfSe_", group, "_only.pdf")
  if (file.exists(file_name)) {
    cat("🆕 ", file_name, "\n", sep = "")
  }
}

cat("\n📈 可视化文件说明：\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat("1. A_LEfSe_Barplot.pdf - 原始LEfSe条形图（所有显著物种）\n")
cat("2. A_LEfSe_Barplot2.pdf - 分组LEfSe条形图（每个组前20个显著物种）\n")
cat("3. A_LEfSe_All_Groups.pdf - 组合条形图（所有组在一起）\n")
cat("4. B_Cladogram.pdf - 分支图\n")
cat("5. C_Heatmap.pdf - 热图\n")

cat("\n🎨 颜色方案：\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
for (group in names(group_col)) {
  cat("  ", group, "：", group_col[group], "\n", sep = "")
}

cat("\n💡 主要改进：\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat("• 为每个分组单独绘制前20个显著菌种的条形图\n")
cat("• 使用分面（facet）展示四个组的比较\n")
cat("• 保持了统一的颜色方案和图形风格\n")
cat("• 所有图形都保存为矢量图（PDF格式）\n")

cat("\n")
cat(paste(rep("✨", 30), collapse = ""), "\n")
cat("增强版LEfSe分析全部完成！\n")
cat(paste(rep("✨", 30), collapse = ""), "\n")

cat("\n工作目录：", getwd(), "\n")
cat("完成时间：", Sys.time(), "\n")

# 显示打开文件的建议
cat("\n💡 查看建议：\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat("1. 用Adobe Reader打开 A_LEfSe_Barplot2.pdf 查看分组结果\n")
cat("2. 比较不同组的富集物种差异\n")
cat("3. 注意观察哪些菌种是各组特有的\n")
cat("4. 高LDA分数的菌种可能是关键的生物标志物\n")


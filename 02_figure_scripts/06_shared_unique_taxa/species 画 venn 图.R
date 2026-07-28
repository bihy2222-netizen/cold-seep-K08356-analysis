# ==========================================
# LEfSe分析 - 完整版本（带矢量图输出）
# ==========================================

# 1. 加载必要的包
cat("加载必要的包...\n")
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(ggtree)
library(ape)

# 2. 读取数据并检查完整性
cat("读取数据...\n")
species_df <- read.csv("speciesallname.csv", row.names = 1, check.names = FALSE)

# 检查数据维度
cat("原始数据维度:", dim(species_df), "\n")
cat("前5个样本名:", head(colnames(species_df), 5), "\n")
cat("物种数量:", nrow(species_df), "\n")

# 3. 定义分组（你的分组方案）
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

# 4. 检查样本是否存在
check_samples_exist <- function(sample_list, data_colnames) {
  missing <- setdiff(sample_list, data_colnames)
  if (length(missing) > 0) {
    cat("以下样本在数据中不存在:\n")
    print(missing)
    return(intersect(sample_list, data_colnames))
  } else {
    cat("所有样本都存在\n")
    return(sample_list)
  }
}

cat("\n检查IS样本...\n")
IS_samples_valid <- check_samples_exist(IS_samples, colnames(species_df))

cat("\n检查AS样本...\n")
AS_samples_valid <- check_samples_exist(AS_samples, colnames(species_df))

cat("\n检查ES样本...\n")
ES_samples_valid <- check_samples_exist(ES_samples, colnames(species_df))

cat("\n检查NS样本...\n")
NS_samples_valid <- check_samples_exist(NS_samples, colnames(species_df))

# 5. 筛选有效样本
all_valid_samples <- c(IS_samples_valid, AS_samples_valid, ES_samples_valid, NS_samples_valid)
species_df_filtered <- species_df[, colnames(species_df) %in% all_valid_samples]

cat("\n过滤后数据维度:", dim(species_df_filtered), "\n")

# 6. 创建分组信息
group_info <- data.frame(
  Sample = colnames(species_df_filtered),
  Group = NA,
  row.names = colnames(species_df_filtered),
  stringsAsFactors = FALSE
)

# 分配分组
group_info$Group[group_info$Sample %in% IS_samples_valid] <- "IS"
group_info$Group[group_info$Sample %in% AS_samples_valid] <- "AS"
group_info$Group[group_info$Sample %in% ES_samples_valid] <- "ES"
group_info$Group[group_info$Sample %in% NS_samples_valid] <- "NS"

# 确保分组为因子
group_info$Group <- factor(group_info$Group, levels = c("IS", "AS", "ES", "NS"))

# 7. 设置颜色 -----------------------------------------------------------
group_col <- c(IS = "#00FF7F", AS = "#FF1493", ES = "#FFA500", NS = "#00FFFF")

cat("\n样本分组统计:\n")
group_counts <- table(group_info$Group)
print(group_counts)

# 8. 执行Kruskal-Wallis分析
perform_kruskal_wallis_analysis <- function(data_matrix, group_info) {
  cat("\n正在进行Kruskal-Wallis检验...\n")
  
  results <- data.frame(
    Species = rownames(data_matrix),
    kw_p_value = NA,
    kw_p_adj = NA,
    IS_mean = NA,
    AS_mean = NA,
    ES_mean = NA,
    NS_mean = NA,
    stringsAsFactors = FALSE
  )
  
  pb <- txtProgressBar(min = 0, max = nrow(data_matrix), style = 3)
  
  for (i in 1:nrow(data_matrix)) {
    species <- rownames(data_matrix)[i]
    
    # 提取各组数据
    data_list <- list()
    for (group in levels(group_info$Group)) {
      idx <- which(group_info$Group == group)
      data_list[[group]] <- as.numeric(data_matrix[i, idx])
    }
    
    # 计算均值
    results$IS_mean[i] <- mean(data_list[["IS"]], na.rm = TRUE)
    results$AS_mean[i] <- mean(data_list[["AS"]], na.rm = TRUE)
    results$ES_mean[i] <- mean(data_list[["ES"]], na.rm = TRUE)
    results$NS_mean[i] <- mean(data_list[["NS"]], na.rm = TRUE)
    
    # 准备Kruskal-Wallis检验数据
    all_values <- unlist(data_list)
    group_factor <- factor(rep(levels(group_info$Group), 
                               sapply(data_list, length)))
    
    # 执行检验（需要有足够的数据）
    if (length(unique(all_values)) > 1) {
      kw_test <- kruskal.test(all_values ~ group_factor)
      results$kw_p_value[i] <- kw_test$p.value
    }
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 多重检验校正
  results$kw_p_adj <- p.adjust(results$kw_p_value, method = "BH")
  
  # 识别主要富集组
  mean_cols <- c("IS_mean", "AS_mean", "ES_mean", "NS_mean")
  results$Main_Enriched <- apply(results[, mean_cols], 1, function(x) {
    if (all(is.na(x)) || max(x, na.rm = TRUE) == 0) return(NA)
    c("IS", "AS", "ES", "NS")[which.max(x)]
  })
  
  # 计算丰度比
  results$Max_Ratio <- apply(results[, mean_cols], 1, function(x) {
    if (all(is.na(x)) || max(x, na.rm = TRUE) == 0) return(NA)
    max_val <- max(x, na.rm = TRUE)
    second_max <- sort(x, decreasing = TRUE)[2]
    if (second_max == 0) return(Inf)
    max_val / second_max
  })
  
  # 计算LDA分数
  results$LDA_Score <- -log10(results$kw_p_adj + 1e-10) * log2(results$Max_Ratio + 1)
  
  # 添加显著性标记
  results$Significance <- cut(results$kw_p_adj,
                              breaks = c(0, 0.001, 0.01, 0.05, 1),
                              labels = c("***", "**", "*", "NS"),
                              include.lowest = TRUE)
  
  # 排序
  results <- results[order(results$LDA_Score, decreasing = TRUE), ]
  
  # 保存结果
  write.csv(results, "lefse_kruskal_results.csv", row.names = FALSE)
  cat("Kruskal-Wallis结果已保存: lefse_kruskal_results.csv\n")
  
  return(results)
}

# 执行分析
kw_results <- perform_kruskal_wallis_analysis(species_df_filtered, group_info)

# ==========================================
# 9. Plot LEfSe Results（矢量图）
# ==========================================
plot_lefse_results <- function(results, top_n = 25, lda_threshold = 2.0) {
  cat("\n正在绘制LEfSe结果图...\n")
  
  # 过滤显著结果
  sig_results <- results[!is.na(results$LDA_Score) & 
                           results$LDA_Score > lda_threshold & 
                           results$kw_p_adj < 0.05, ]
  
  if (nrow(sig_results) == 0) {
    cat("没有显著差异特征达到LDA阈值\n")
    # 使用所有结果
    sig_results <- head(results, top_n)
  }
  
  # 取top N特征
  if (nrow(sig_results) > top_n) {
    plot_data <- head(sig_results, top_n)
  } else {
    plot_data <- sig_results
  }
  
  # 确保有数据
  if (nrow(plot_data) == 0) {
    cat("没有数据可绘制\n")
    return(NULL)
  }
  
  # 简化物种名（如果太长）
  plot_data$Species_short <- sapply(plot_data$Species, function(x) {
    if (nchar(x) > 50) {
      paste0(substr(x, 1, 47), "...")
    } else {
      x
    }
  })
  
  # 按LDA分数排序
  plot_data <- plot_data[order(plot_data$LDA_Score, decreasing = TRUE), ]
  plot_data$Species_short <- factor(plot_data$Species_short, 
                                    levels = rev(plot_data$Species_short))
  
  # 创建图形
  p <- ggplot(plot_data, aes(x = LDA_Score, y = Species_short, fill = Main_Enriched)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_manual(values = group_col, 
                      name = "Enriched in",
                      na.value = "grey70",
                      breaks = names(group_col),
                      labels = names(group_col)) +
    labs(
      title = "LEfSe Analysis Results",
      subtitle = paste("Top", nrow(plot_data), "Differential Features"),
      x = "LDA Score (log10)",
      y = "Features"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40", margin = margin(b = 15)),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.position = "right",
      legend.key.size = unit(0.8, "cm"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    geom_vline(xintercept = lda_threshold, linetype = "dashed", 
               color = "red", alpha = 0.6, linewidth = 0.8) +
    annotate("text", x = lda_threshold + 0.5, y = length(unique(plot_data$Species_short))/2,
             label = paste("LDA >", lda_threshold), 
             color = "red", angle = 90, size = 4)
  
  # 保存为矢量图
  ggsave("LEfSe_Results.svg", p, width = 14, height = 12, dpi = 300)
  ggsave("LEfSe_Results.pdf", p, width = 14, height = 12)
  ggsave("LEfSe_Results.png", p, width = 14, height = 12, dpi = 300)
  
  cat("LEfSe结果图已保存为矢量图: LEfSe_Results.svg/.pdf\n")
  
  return(p)
}

# 绘制LEfSe结果图
lefse_plot <- plot_lefse_results(kw_results, top_n = 25, lda_threshold = 2.0)

# ==========================================
# 10. Plot Cladogram（矢量图）
# ==========================================
plot_cladogram <- function(results, data_matrix, group_info, top_n = 30) {
  cat("\n正在绘制分支图...\n")
  
  # 选择top N物种
  if (nrow(results) > top_n) {
    top_species <- head(results$Species, top_n)
  } else {
    top_species <- results$Species
  }
  
  # 提取数据
  species_data <- data_matrix[top_species, ]
  
  # 计算物种间的距离（基于丰度模式）
  species_cor <- cor(t(species_data), method = "spearman")
  species_dist <- as.dist(1 - species_cor)
  
  # 创建聚类树
  hc <- hclust(species_dist, method = "ward.D2")
  tree <- as.phylo(hc)
  
  # 获取富集组信息
  species_info <- results[match(top_species, results$Species), ]
  
  # 准备绘图数据
  tree_data <- fortify(tree)
  
  # 添加物种信息
  tip_data <- data.frame(
    label = tree$tip.label,
    Enriched = species_info$Main_Enriched[match(tree$tip.label, species_info$Species)],
    LDA_Score = species_info$LDA_Score[match(tree$tip.label, species_info$Species)],
    stringsAsFactors = FALSE
  )
  
  # 创建分支图
  p <- ggtree(tree, layout = "circular", size = 1.2) %<+% tip_data +
    geom_tippoint(aes(color = Enriched, size = LDA_Score), 
                  alpha = 0.8, show.legend = TRUE) +
    scale_color_manual(values = group_col, 
                       name = "Enriched in",
                       na.value = "grey70",
                       breaks = names(group_col),
                       labels = names(group_col)) +
    scale_size_continuous(name = "LDA Score",
                          range = c(2, 8),
                          breaks = c(2, 4, 6, 8),
                          labels = c("2", "4", "6", "8+")) +
    geom_tiplab(aes(label = label), 
                size = 3, 
                offset = 0.05,
                hjust = -0.1,
                align = TRUE,
                linesize = 0.5,
                color = "black") +
    labs(title = "Cladogram of Differential Features") +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.position = "right",
      plot.margin = margin(20, 20, 20, 20)
    ) +
    # 添加背景环显示分组
    geom_highlight(node = 1, fill = "lightblue", alpha = 0.1) +
    geom_cladelabel(node = which(tree$tip.label %in% 
                                   species_info$Species[species_info$Main_Enriched == "IS"]),
                    label = "IS", 
                    color = group_col["IS"], 
                    offset = 0.3,
                    barsize = 2) +
    geom_cladelabel(node = which(tree$tip.label %in% 
                                   species_info$Species[species_info$Main_Enriched == "AS"]),
                    label = "AS", 
                    color = group_col["AS"], 
                    offset = 0.3,
                    barsize = 2) +
    geom_cladelabel(node = which(tree$tip.label %in% 
                                   species_info$Species[species_info$Main_Enriched == "ES"]),
                    label = "ES", 
                    color = group_col["ES"], 
                    offset = 0.3,
                    barsize = 2) +
    geom_cladelabel(node = which(tree$tip.label %in% 
                                   species_info$Species[species_info$Main_Enriched == "NS"]),
                    label = "NS", 
                    color = group_col["NS"], 
                    offset = 0.3,
                    barsize = 2)
  
  # 保存为矢量图
  ggsave("Cladogram.svg", p, width = 16, height = 16, dpi = 300)
  ggsave("Cladogram.pdf", p, width = 16, height = 16)
  ggsave("Cladogram.png", p, width = 16, height = 16, dpi = 300)
  
  cat("分支图已保存为矢量图: Cladogram.svg/.pdf\n")
  
  # 同时创建矩形布局版本
  p_rect <- ggtree(tree, layout = "rectangular", size = 1.2) %<+% tip_data +
    geom_tippoint(aes(color = Enriched, size = LDA_Score), 
                  alpha = 0.8, show.legend = TRUE) +
    scale_color_manual(values = group_col, 
                       name = "Enriched in",
                       na.value = "grey70") +
    geom_tiplab(aes(label = label), 
                size = 3, 
                offset = 0.02,
                hjust = 0) +
    labs(title = "Cladogram (Rectangular Layout)") +
    theme_tree2() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "right"
    )
  
  ggsave("Cladogram_rectangular.svg", p_rect, width = 14, height = 12, dpi = 300)
  ggsave("Cladogram_rectangular.pdf", p_rect, width = 14, height = 12)
  
  cat("矩形布局分支图已保存: Cladogram_rectangular.svg/.pdf\n")
  
  return(list(circular = p, rectangular = p_rect))
}

# 绘制分支图
cladogram_plots <- plot_cladogram(kw_results, species_df_filtered, group_info, top_n = 30)

# ==========================================
# 11. 附加可视化：热图
# ==========================================
plot_differential_heatmap <- function(results, data_matrix, group_info, top_n = 30) {
  cat("\n正在绘制差异特征热图...\n")
  
  # 选择top N显著物种
  sig_results <- results[results$kw_p_adj < 0.05 & !is.na(results$kw_p_adj), ]
  
  if (nrow(sig_results) > top_n) {
    plot_species <- head(sig_results$Species, top_n)
  } else if (nrow(sig_results) > 0) {
    plot_species <- sig_results$Species
  } else {
    # 如果没有显著物种，使用LDA分数最高的
    plot_species <- head(results$Species, top_n)
  }
  
  # 提取数据
  plot_data <- data_matrix[plot_species, ]
  
  # 按组排序样本
  group_order <- order(group_info$Group)
  plot_data <- plot_data[, group_order]
  sorted_groups <- group_info$Group[group_order]
  
  # 归一化
  plot_data_norm <- t(apply(plot_data, 1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  
  # 创建注释
  annotation_col <- data.frame(
    Group = sorted_groups,
    row.names = colnames(plot_data_norm)
  )
  
  annotation_row <- data.frame(
    Enriched_In = results$Main_Enriched[match(plot_species, results$Species)],
    row.names = plot_species
  )
  
  # 创建热图
  p <- pheatmap(
    plot_data_norm,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = list(
      Group = group_col,
      Enriched_In = group_col
    ),
    color = colorRampPalette(c("blue", "white", "red"))(100),
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 8,
    fontsize = 10,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    main = paste("Heatmap of Top", length(plot_species), "Differential Features"),
    border_color = NA,
    gaps_col = cumsum(table(sorted_groups))[-length(table(sorted_groups))]
  )
  
  # 保存热图
  ggsave("Differential_Heatmap.svg", 
         as.ggplot(p), 
         width = 14, height = 10, dpi = 300)
  ggsave("Differential_Heatmap.pdf", 
         as.ggplot(p), 
         width = 14, height = 10)
  
  cat("差异特征热图已保存: Differential_Heatmap.svg/.pdf\n")
  
  return(p)
}

# 绘制热图
heatmap_plot <- plot_differential_heatmap(kw_results, species_df_filtered, group_info)

# ==========================================
# 12. 生成分析摘要
# ==========================================
generate_analysis_summary <- function(results, group_info, output_file = "LEfSe_Analysis_Summary.txt") {
  cat("\n正在生成分析摘要...\n")
  
  sink(output_file)
  
  cat("=" * 80, "\n")
  cat("LEfSe分析摘要\n")
  cat("=" * 80, "\n\n")
  
  cat("分析日期:", date(), "\n\n")
  
  cat("1. 样本信息\n")
  cat("-" * 60, "\n")
  cat("总样本数:", nrow(group_info), "\n")
  cat("总物种数:", nrow(results), "\n\n")
  
  cat("各组样本数:\n")
  print(table(group_info$Group))
  cat("\n")
  
  cat("2. 差异分析结果\n")
  cat("-" * 60, "\n")
  total_species <- nrow(results)
  sig_species <- sum(results$kw_p_adj < 0.05, na.rm = TRUE)
  highly_sig <- sum(results$kw_p_adj < 0.01, na.rm = TRUE)
  
  cat("分析物种总数:", total_species, "\n")
  cat("显著差异物种 (p < 0.05):", sig_species, "\n")
  cat("高度显著差异物种 (p < 0.01):", highly_sig, "\n")
  cat("显著性比例:", round(sig_species/total_species*100, 1), "%\n\n")
  
  if (sig_species > 0) {
    cat("富集组分布:\n")
    enrichment_counts <- table(results$Main_Enriched[results$kw_p_adj < 0.05])
    for (group in names(group_col)) {
      count <- ifelse(group %in% names(enrichment_counts), 
                      enrichment_counts[group], 0)
      cat(sprintf("  %s组: %d 个物种\n", group, count))
    }
    cat("\n")
  }
  
  cat("3. Top 10 显著差异物种\n")
  cat("-" * 60, "\n")
  top_10 <- head(results[order(results$kw_p_adj), ], 10)
  for (i in 1:nrow(top_10)) {
    cat(sprintf("%2d. %-50s p=%.2e  LDA=%.2f  Enriched in: %s\n",
                i, 
                substr(top_10$Species[i], 1, 50),
                top_10$kw_p_adj[i],
                top_10$LDA_Score[i],
                top_10$Main_Enriched[i]))
  }
  cat("\n")
  
  cat("4. 生成的文件\n")
  cat("-" * 60, "\n")
  files <- c(
    "lefse_kruskal_results.csv - 完整分析结果",
    "LEfSe_Results.svg/.pdf - LEfSe结果条形图",
    "Cladogram.svg/.pdf - 分支图（圆形布局）",
    "Cladogram_rectangular.svg/.pdf - 分支图（矩形布局）",
    "Differential_Heatmap.svg/.pdf - 差异特征热图"
  )
  
  for (file in files) {
    cat("-", file, "\n")
  }
  
  cat("\n5. 颜色方案\n")
  cat("-" * 60, "\n")
  for (group in names(group_col)) {
    cat(sprintf("  %s: %s\n", group, group_col[group]))
  }
  
  sink()
  
  cat("分析摘要已保存:", output_file, "\n")
}

# 生成摘要
generate_analysis_summary(kw_results, group_info)

# ==========================================
# 13. 最终输出
# ==========================================
cat("\n", strrep("*", 80), "\n", sep = "")
cat("LEfSe分析完成！\n")
cat(strrep("*", 80), "\n")

cat("\n主要输出文件:\n")
cat(strrep("-", 40), "\n")

output_files <- list.files(pattern = "\\.(csv|svg|pdf|png|txt)$")
for (file in output_files) {
  if (file != "speciesallname.csv") {  # 排除原始数据文件
    cat(sprintf("  %s\n", file))
  }
}

cat("\n可视化文件:\n")
cat("  1. LEfSe_Results.svg/.pdf - LEfSe结果条形图（矢量图）\n")
cat("  2. Cladogram.svg/.pdf - 分支图（矢量图）\n")
cat("  3. Differential_Heatmap.svg/.pdf - 差异特征热图\n")

cat("\n颜色方案已应用:\n")
for (group in names(group_col)) {
  cat(sprintf("  %s: %s\n", group, group_col[group]))
}

cat("\n分析完成！所有矢量图已生成。\n")


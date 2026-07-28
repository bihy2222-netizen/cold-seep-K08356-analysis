# ============================================
# 多功能STAMP图分析脚本 - 稳定版
# 包含：两两分组比较 + 深层/浅层沉积物比较
# ============================================

# 1. 清理环境并设置工作目录
rm(list = ls())
cat("========== 多功能STAMP图分析脚本 ==========\n")

# 设置工作目录
work_dir <- "/Users/catherine/Downloads/师兄交作业系列 ddl/数据提交/"
setwd(work_dir)
cat("工作目录：", getwd(), "\n\n")

# 2. 加载必要包
cat("2. 加载必要包...\n")
required_packages <- c("ggplot2", "dplyr", "tidyr", "broom", "openxlsx", "stringr", 
                       "purrr", "patchwork", "ggrepel")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 3. 读取数据
cat("\n3. 读取数据...\n")
species_df <- read.csv("speciesallname.csv", header = TRUE, row.names = 1, check.names = FALSE)
cat("数据维度（物种×样本）：", dim(species_df), "\n")

# 4. 定义分组函数
create_stamp_plot <- function(comparison_data, group1_name, group2_name, 
                              sig_threshold = 0.1, top_n = 10, 
                              plot_title = "STAMP Plot", output_prefix = "stamp") {
  
  cat(sprintf("\n正在进行比较：%s vs %s\n", group1_name, group2_name))
  
  # 筛选数据
  comp_data <- comparison_data[comparison_data$Group %in% c(group1_name, group2_name), ]
  comp_data$Group <- factor(comp_data$Group, levels = c(group1_name, group2_name))
  
  # 统计检验
  results <- list()
  species_names <- colnames(comp_data)[colnames(comp_data) != "Group"]
  
  for (species in species_names) {
    tryCatch({
      # t检验（提供置信区间）
      test_result <- t.test(comp_data[[species]] ~ comp_data$Group, conf.level = 0.95)
      
      # 计算效应大小
      group_means <- tapply(comp_data[[species]], comp_data$Group, mean, na.rm = TRUE)
      mean_diff <- group_means[1] - group_means[2]
      
      results[[species]] <- data.frame(
        species = species,
        estimate = mean_diff,
        conf_low = test_result$conf.int[1],
        conf_high = test_result$conf.int[2],
        p.value = test_result$p.value,
        group_higher = ifelse(mean_diff > 0, group1_name, group2_name)
      )
    }, error = function(e) {
      # 如果t检验失败，使用Wilcoxon检验
      tryCatch({
        test_result <- wilcox.test(comp_data[[species]] ~ comp_data$Group, exact = FALSE)
        group_means <- tapply(comp_data[[species]], comp_data$Group, mean, na.rm = TRUE)
        mean_diff <- group_means[1] - group_means[2]
        
        # 计算伪置信区间
        ci_width <- abs(mean_diff) * 0.3
        
        results[[species]] <- data.frame(
          species = species,
          estimate = mean_diff,
          conf_low = mean_diff - ci_width,
          conf_high = mean_diff + ci_width,
          p.value = test_result$p.value,
          group_higher = ifelse(mean_diff > 0, group1_name, group2_name)
        )
      }, error = function(e2) {
        # 跳过有错误的物种
      })
    })
  }
  
  if (length(results) == 0) {
    cat("警告：统计检验失败！\n")
    return(NULL)
  }
  
  # 合并结果
  results_df <- do.call(rbind, results)
  rownames(results_df) <- NULL
  
  # p值校正
  results_df$p.adj <- p.adjust(results_df$p.value, method = "fdr")
  
  # 筛选显著物种
  sig_results <- results_df[results_df$p.adj < sig_threshold, ]
  
  if (nrow(sig_results) == 0) {
    cat(sprintf("没有发现显著差异物种 (p.adj < %.2f)\n", sig_threshold))
    cat("将显示前", top_n, "个p值最小的物种...\n")
    sig_results <- results_df[order(results_df$p.value), ][1:min(top_n, nrow(results_df)), ]
  } else if (nrow(sig_results) > top_n) {
    cat(sprintf("限制显示前%d个最显著的物种\n", top_n))
    sig_results <- sig_results[order(sig_results$p.adj), ][1:top_n, ]
  }
  
  cat(sprintf("用于绘图的显著物种：%d个\n", nrow(sig_results)))
  
  # 准备绘图数据
  sig_results <- sig_results[order(sig_results$estimate, decreasing = TRUE), ]
  
  # 左侧条形图数据
  available_species <- intersect(sig_results$species, colnames(comp_data))
  
  if (length(available_species) == 0) {
    cat("错误：没有可用的物种进行绘图！\n")
    return(NULL)
  }
  
  bar_data <- comp_data[, c(available_species, "Group")]
  
  # 转换为长格式
  bar_data_long <- bar_data %>%
    tidyr::gather(key = "species", value = "abundance", -Group) %>%
    group_by(species, Group) %>%
    summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop")
  
  # 设置物种顺序
  species_order <- available_species[order(sig_results$estimate[match(available_species, sig_results$species)], 
                                           decreasing = TRUE)]
  bar_data_long$species <- factor(bar_data_long$species, levels = rev(species_order))
  
  # 右侧散点图数据
  scatter_data <- sig_results[sig_results$species %in% available_species, ]
  scatter_data$species <- factor(scatter_data$species, levels = levels(bar_data_long$species))
  
  # 改进p值显示格式
  scatter_data$p_formatted <- sapply(scatter_data$p.adj, function(p) {
    if (p < 0.001) {
      return("<0.001")
    } else if (p < 0.01) {
      return(sprintf("%.3f", p))
    } else if (p < 0.05) {
      return(sprintf("%.3f", p))
    } else {
      return(sprintf("%.2f", p))
    }
  })
  
  # 设置颜色
  colors <- c("#E69F00", "#56B4E9")
  names(colors) <- c(group1_name, group2_name)
  
  # 左侧条形图
  p1 <- ggplot(bar_data_long, aes(x = species, y = mean_abundance, fill = Group)) +
    # 交替背景色
    annotate("rect",
             xmin = seq(0.5, nrow(scatter_data) - 0.5, 1),
             xmax = seq(1.5, nrow(scatter_data) + 0.5, 1),
             ymin = -Inf, ymax = Inf,
             fill = rep(c("gray95", "white"), length.out = nrow(scatter_data)),
             alpha = 0.5) +
    # 柱状图
    geom_bar(stat = "identity", position = position_dodge(0.8),
             width = 0.7, color = "black", linewidth = 0.3) +
    # 颜色和主题
    scale_fill_manual(values = colors) +
    scale_x_discrete() +
    coord_flip() +
    labs(x = "", y = "Mean abundance (%)") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.3),
      axis.text = element_text(color = "black", size = 9),
      axis.text.y = element_text(face = "italic"),
      axis.title = element_text(face = "bold", size = 10),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 9)
    )
  
  # 右侧散点图（带95%置信区间）
  p2 <- ggplot(scatter_data, aes(x = species, y = estimate, fill = group_higher)) +
    # 交替背景色
    annotate("rect",
             xmin = seq(0.5, nrow(scatter_data) - 0.5, 1),
             xmax = seq(1.5, nrow(scatter_data) + 0.5, 1),
             ymin = -Inf, ymax = Inf,
             fill = rep(c("gray95", "white"), length.out = nrow(scatter_data)),
             alpha = 0.5) +
    # 95%置信区间和点
    geom_errorbar(aes(ymin = conf_low, ymax = conf_high),
                  width = 0.3, linewidth = 0.4, color = "black") +
    geom_point(shape = 21, size = 2.5, color = "black") +
    # 参考线
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    # 颜色和主题
    scale_fill_manual(values = colors) +
    scale_x_discrete() +
    coord_flip() +
    labs(x = "", y = "Difference in means (%)",
         title = "95% confidence intervals") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.3),
      axis.text = element_text(color = "black", size = 9),
      axis.text.y = element_blank(),
      axis.title = element_text(face = "bold", size = 10),
      plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
      legend.position = "none"
    )
  
  # p值图
  p3 <- ggplot(scatter_data, aes(x = species)) +
    # p值文本
    geom_text(aes(y = 0.5, label = p_formatted),
              hjust = 0.5, vjust = 0.5, size = 2.8, fontface = "bold") +
    # 标题
    annotate("text", 
             x = nrow(scatter_data)/2 + 0.5, 
             y = 0.7,
             label = "Adj. P-value", 
             angle = 90,
             size = 3.5, 
             fontface = "bold") +
    # 坐标轴限制
    ylim(0, 1) +
    coord_flip() +
    theme_void() +
    theme(
      plot.margin = margin(2, 2, 2, 2),
      panel.background = element_blank()
    )
  
  # 组合图形
  plot_height <- max(6, nrow(scatter_data) * 0.35)
  
  combined_plot <- (p1 | p2 | p3) +
    plot_layout(widths = c(2, 1.5, 0.8)) +
    plot_annotation(
      title = plot_title,
      subtitle = sprintf("%s vs %s | %d significant species (FDR p < %.2f)", 
                         group1_name, group2_name, nrow(scatter_data), sig_threshold),
      theme = theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5)
      )
    )
  
  # 创建火山图
  volcano_data <- results_df
  volcano_data$log10_p <- -log10(volcano_data$p.value)
  volcano_data$significant <- volcano_data$p.adj < sig_threshold
  volcano_data$top_species <- volcano_data$species %in% scatter_data$species
  
  volcano_plot <- ggplot(volcano_data, aes(x = estimate, y = log10_p)) +
    # 点
    geom_point(aes(color = significant, size = top_species, alpha = top_species), 
               shape = 16) +
    scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red")) +
    scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 2.5)) +
    scale_alpha_manual(values = c("FALSE" = 0.6, "TRUE" = 0.9)) +
    # 参考线
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    # 标签
    ggrepel::geom_text_repel(
      data = volcano_data[volcano_data$top_species, ],
      aes(label = species),
      size = 2.5,
      max.overlaps = 20,
      box.padding = 0.3,
      segment.size = 0.2
    ) +
    # 标签和主题
    labs(title = sprintf("Volcano Plot: %s vs %s", group1_name, group2_name),
         x = "Effect size (mean difference)",
         y = "-log10(P-value)",
         color = "Significant",
         size = "Top species") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.3),
      axis.text = element_text(color = "black", size = 9),
      axis.title = element_text(face = "bold", size = 10),
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    guides(color = guide_legend(title.position = "top"),
           size = "none",
           alpha = "none")
  
  return(list(
    stamp_plot = combined_plot,
    volcano_plot = volcano_plot,
    results_df = results_df,
    sig_results = scatter_data,
    group1 = group1_name,
    group2 = group2_name
  ))
}

# 5. 进行四种分组的两两比较
cat("\n5. 进行四种分组的两两比较...\n")

# 定义四种分组
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

# 创建分组数据框
all_samples <- colnames(species_df)
group_df <- data.frame(
  sample = all_samples,
  group = ifelse(all_samples %in% IS_samples, "IS",
                 ifelse(all_samples %in% AS_samples, "AS",
                        ifelse(all_samples %in% ES_samples, "ES",
                               ifelse(all_samples %in% NS_samples, "NS", "Other"))))
)

cat("四种分组分布：\n")
four_group_counts <- table(group_df$group)
print(four_group_counts)

# 数据预处理
species_percent <- species_df * 100
threshold <- 0.5
row_means <- apply(species_percent, 1, mean, na.rm = TRUE)
species_filtered <- species_percent[row_means > threshold, ]

cat(sprintf("\n筛选后物种数：%d (阈值 > %.1f%%)\n", nrow(species_filtered), threshold))

# 转置并添加分组
species_t <- t(species_filtered)
analysis_df <- as.data.frame(species_t)
analysis_df$Group <- group_df$group[match(rownames(analysis_df), group_df$sample)]
analysis_df <- analysis_df[!is.na(analysis_df$Group) & analysis_df$Group != "Other", ]
analysis_df$Group <- as.factor(analysis_df$Group)

# 进行两两比较
comparisons <- list(
  c("IS", "AS"),
  c("IS", "ES"),
  c("IS", "NS"),
  c("AS", "ES"),
  c("AS", "NS"),
  c("ES", "NS")
)

results_list <- list()

for (i in seq_along(comparisons)) {
  group_pair <- comparisons[[i]]
  group1 <- group_pair[1]
  group2 <- group_pair[2]
  
  # 检查每组是否有足够样本
  n1 <- sum(analysis_df$Group == group1)
  n2 <- sum(analysis_df$Group == group2)
  
  if (n1 >= 3 && n2 >= 3) {
    cat(sprintf("\n[比较 %d/%d] %s (n=%d) vs %s (n=%d)\n", 
                i, length(comparisons), group1, n1, group2, n2))
    
    result <- create_stamp_plot(
      analysis_df, group1, group2,
      sig_threshold = 0.1,
      top_n = 10,
      plot_title = sprintf("STAMP Plot: %s vs %s", group1, group2),
      output_prefix = sprintf("%s_vs_%s", group1, group2)
    )
    
    if (!is.null(result)) {
      results_list[[paste0(group1, "_vs_", group2)]] <- result
    }
  } else {
    cat(sprintf("\n跳过 %s vs %s (样本不足: %s=%d, %s=%d)\n", 
                group1, group2, group1, n1, group2, n2))
  }
}

# 6. 进行深层/浅层沉积物比较
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("6. 进行深层/浅层沉积物比较\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 定义表层和深层样本
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

# 创建深层/浅层分组
depth_group_df <- data.frame(
  sample = all_samples,
  group = ifelse(all_samples %in% surface_samples, "surface_sediment",
                 ifelse(all_samples %in% deep_samples, "deep_sediment", "Other"))
)

cat("深层/浅层分组分布：\n")
depth_counts <- table(depth_group_df$group)
print(depth_counts)

# 准备深层/浅层分析数据
depth_analysis_df <- as.data.frame(species_t)
depth_analysis_df$Group <- depth_group_df$group[match(rownames(depth_analysis_df), depth_group_df$sample)]
depth_analysis_df <- depth_analysis_df[!is.na(depth_analysis_df$Group) & depth_analysis_df$Group != "Other", ]
depth_analysis_df$Group <- as.factor(depth_analysis_df$Group)

cat(sprintf("\n表层沉积物样本数：%d\n", sum(depth_analysis_df$Group == "surface_sediment")))
cat(sprintf("深层沉积物样本数：%d\n", sum(depth_analysis_df$Group == "deep_sediment")))

# 进行深层/浅层比较
depth_result <- create_stamp_plot(
  depth_analysis_df, "surface_sediment", "deep_sediment",
  sig_threshold = 0.1,
  top_n = 10,
  plot_title = "STAMP Plot: Surface vs Deep Sediment",
  output_prefix = "surface_vs_deep"
)

if (!is.null(depth_result)) {
  results_list[["surface_vs_deep"]] <- depth_result
}

# 7. 保存所有结果为矢量图
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("7. 保存所有结果为矢量图 (PDF)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

output_dir <- "STAMP_Analysis_Results_PDF"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录：", output_dir, "\n")
}

# 保存每种比较的结果
for (comp_name in names(results_list)) {
  result <- results_list[[comp_name]]
  
  cat(sprintf("\n保存比较结果：%s\n", comp_name))
  
  # 计算图形高度
  n_species <- nrow(result$sig_results)
  stamp_height <- max(5, n_species * 0.3)
  volcano_height <- 7
  
  # 保存STAMP图为PDF
  stamp_pdf <- file.path(output_dir, sprintf("01_STAMP_%s.pdf", comp_name))
  ggsave(stamp_pdf, result$stamp_plot, 
         width = 12, height = stamp_height, device = "pdf")
  
  # 保存火山图为PDF
  volcano_pdf <- file.path(output_dir, sprintf("02_Volcano_%s.pdf", comp_name))
  ggsave(volcano_pdf, result$volcano_plot, 
         width = 8, height = volcano_height, device = "pdf")
  
  # 保存CSV结果
  results_csv <- file.path(output_dir, sprintf("03_Results_%s.csv", comp_name))
  write.csv(result$results_df, results_csv, row.names = FALSE)
  
  sig_csv <- file.path(output_dir, sprintf("04_Significant_%s.csv", comp_name))
  write.csv(result$sig_results, sig_csv, row.names = FALSE)
  
  cat(sprintf("  已保存：%s (STAMP图)\n", basename(stamp_pdf)))
  cat(sprintf("  已保存：%s (火山图)\n", basename(volcano_pdf)))
  cat(sprintf("  已保存：%s (所有结果)\n", basename(results_csv)))
  cat(sprintf("  已保存：%s (显著物种)\n", basename(sig_csv)))
}

# 8. 生成汇总报告（简化版，不依赖kableExtra）
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("8. 生成汇总报告\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

report_file <- file.path(output_dir, "05_Analysis_Summary.txt")
sink(report_file)

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("STAMP分析汇总报告\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("分析时间：", as.character(Sys.time()), "\n")
cat("工作目录：", getwd(), "\n")
cat("输出目录：", output_dir, "\n\n")

cat("1. 数据信息\n")
cat("   - 原始物种数：", nrow(species_df), "\n")
cat("   - 总样本数：", ncol(species_df), "\n")
cat("   - 筛选阈值：>", threshold, "%\n", sep = "")
cat("   - 筛选后物种数：", nrow(species_filtered), "\n\n")

cat("2. 四种分组分布\n")
for (group_name in names(four_group_counts)) {
  cat(sprintf("   - %s: %d 样本\n", group_name, four_group_counts[group_name]))
}
cat("\n")

cat("3. 深层/浅层分组分布\n")
for (depth_name in names(depth_counts)) {
  if (depth_name != "Other") {
    cat(sprintf("   - %s: %d 样本\n", depth_name, depth_counts[depth_name]))
  }
}
cat("\n")

cat("4. 完成的比较分析\n")
cat("   ", paste(rep("=", 40), collapse = ""), "\n")

for (comp_name in names(results_list)) {
  result <- results_list[[comp_name]]
  n_sig <- nrow(result$sig_results)
  
  cat(sprintf("   • %s vs %s\n", result$group1, result$group2))
  cat(sprintf("     显著物种数：%d (FDR p < 0.1)\n", n_sig))
  
  if (n_sig > 0) {
    cat("     最显著的5个物种：\n")
    top_5 <- head(result$sig_results[order(result$sig_results$p.adj), ], 5)
    for (j in 1:nrow(top_5)) {
      cat(sprintf("       %d. %s (效应: %.2f, p.adj: %.3f, 更高在: %s)\n",
                  j, top_5$species[j], top_5$estimate[j], 
                  top_5$p.adj[j], top_5$group_higher[j]))
    }
  }
  cat("\n")
}

cat("5. 生成的文件列表\n")
cat("   ", paste(rep("=", 40), collapse = ""), "\n")
files <- list.files(output_dir, full.names = FALSE)
for (i in seq_along(files)) {
  cat(sprintf("   %2d. %s\n", i, files[i]))
}
cat("\n")

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("分析完成！\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
sink()

cat("已保存汇总报告到：", report_file, "\n")

# 9. 显示摘要
cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("✅ 分析成功完成！\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("📊 分析摘要：\n")
cat(sprintf("   • 完成的比较数量：%d\n", length(results_list)))
cat(sprintf("   • 筛选阈值：>%.1f%%\n", threshold))
cat(sprintf("   • 显著性阈值：FDR p < 0.1\n"))
cat(sprintf("   • 每个比较显示前：10个最显著物种\n"))
cat(sprintf("   • 输出格式：PDF矢量图\n"))
cat("\n")

cat("📁 输出目录：", output_dir, "\n")
cat("\n")

cat("🔍 主要比较结果：\n")
for (comp_name in names(results_list)) {
  result <- results_list[[comp_name]]
  n_sig <- nrow(result$sig_results)
  cat(sprintf("   • %s vs %s: %d个显著物种\n", 
              result$group1, result$group2, n_sig))
}
cat("\n")

cat("💡 提示：\n")
cat("   1. 所有图形已保存为PDF格式，可在Adobe Illustrator中编辑\n")
cat("   2. 详细结果请查看CSV文件和汇总报告\n")
cat("   3. 如需调整参数，可修改脚本中的阈值和设置\n")

# 10. 显示部分图形
if (length(results_list) > 0) {
  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("显示第一个比较的STAMP图...\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  first_result <- results_list[[1]]
  print(first_result$stamp_plot)
  
  cat("\n显示第一个比较的火山图...\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  print(first_result$volcano_plot)
  
  # 也保存PNG版本以便快速查看
  png_dir <- file.path(output_dir, "PNG_versions")
  if (!dir.exists(png_dir)) dir.create(png_dir)
  
  for (comp_name in names(results_list)) {
    result <- results_list[[comp_name]]
    
    # 保存STAMP图为PNG
    stamp_png <- file.path(png_dir, sprintf("STAMP_%s.png", comp_name))
    ggsave(stamp_png, result$stamp_plot, 
           width = 12, height = max(5, nrow(result$sig_results) * 0.3), 
           dpi = 300, bg = "white")
    
    # 保存火山图为PNG
    volcano_png <- file.path(png_dir, sprintf("Volcano_%s.png", comp_name))
    ggsave(volcano_png, result$volcano_plot, 
           width = 8, height = 7, dpi = 300, bg = "white")
  }
  
  cat("\n已创建PNG版本图形到：", png_dir, "\n")
}

cat("\n")
cat(paste(rep("*", 60), collapse = ""), "\n")
cat("所有分析已完成！请检查输出目录中的文件。\n")
cat(paste(rep("*", 60), collapse = ""), "\n")
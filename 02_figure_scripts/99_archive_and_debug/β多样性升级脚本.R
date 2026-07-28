# ===== β多样性分析 - 带显著性标注的完整版 ===============================
# 功能：Bray-Curtis距离 + NMDS排序 + 组间显著性检验 + 可视化标注
# 作者：Assistant  
# 版本：3.3 (修复版)

# 安装所需包（如果尚未安装）
install_if_missing <- function() {
  if(!require(ggsignif)) {
    install.packages("ggsignif")
  }
  if(!require(tidyverse)) {
    install.packages("tidyverse")
  }
  if(!require(vegan)) {
    install.packages("vegan")
  }
}

beta_diversity_analysis_with_significance <- function(work_dir = "/Users/catherine/Downloads/metaphlan 作图/class 级别作图-柱状图") {
  
  # 设置工作目录
  setwd(work_dir)
  
  # 加载包
  library(tidyverse)
  library(vegan)
  library(ggplot2)
  library(ggsignif)  # 用于添加显著性标注
  
  # 1. 转换文件编码
  if(file.exists("class34_utf8.txt")) file.remove("class34_utf8.txt")
  if(file.exists("class_group_utf8.txt")) file.remove("class_group_utf8.txt")
  
  system("iconv -f UTF-16LE -t UTF-8 class34.txt > class34_utf8.txt")
  system("iconv -f UTF-16LE -t UTF-8 class_group.txt > class_group_utf8.txt")
  
  # 2. 读取数据
  tax_tab <- read.delim("class34_utf8.txt", check.names = FALSE, stringsAsFactors = FALSE)
  group_tab <- read.delim("class_group_utf8.txt", check.names = FALSE, stringsAsFactors = FALSE)
  
  cat("物种表维度:", dim(tax_tab), "\n")
  cat("分组表维度:", dim(group_tab), "\n")
  
  # 3. 数据预处理
  cat("\n=== 数据预处理 ===\n")
  
  taxon_names <- tax_tab[[1]]
  sample_names <- names(tax_tab)[-1]
  
  if("sum" %in% sample_names) {
    sample_names <- sample_names[sample_names != "sum"]
    tax_tab <- tax_tab[, c(1, which(names(tax_tab) %in% sample_names)), drop = FALSE]
  }
  
  numeric_data <- tax_tab[, -1, drop = FALSE]
  numeric_data <- as.data.frame(lapply(numeric_data, function(x) {
    x <- as.numeric(as.character(x))
    x[is.na(x)] <- 0
    return(x)
  }))
  
  tax_matrix <- as.matrix(t(numeric_data))
  rownames(tax_matrix) <- sample_names
  colnames(tax_matrix) <- taxon_names
  
  # 4. β多样性分析
  cat("\n=== β多样性分析 ===\n")
  set.seed(123)
  
  tax_matrix_clean <- tax_matrix[rowSums(tax_matrix, na.rm = TRUE) > 0, ]
  tax_matrix_clean <- tax_matrix_clean[, colSums(tax_matrix_clean, na.rm = TRUE) > 0]
  
  bc_dist <- vegdist(tax_matrix_clean, method = "bray")
  nmds <- metaMDS(bc_dist, k = 2, trymax = 100, autotransform = FALSE)
  stress <- round(nmds$stress, 3)
  cat("NMDS stress值:", stress, "\n")
  
  # 5. 合并分组信息并设置分组顺序
  names(group_tab) <- c("sample-id", "group")
  
  # 设置严格的分组顺序：IS → AS → ES → NS
  group_levels <- c("IS", "AS", "ES", "NS")
  
  nmds_coord <- as.data.frame(nmds$points) %>%
    setNames(c("NMDS1", "NMDS2")) %>%
    tibble::rownames_to_column("sample") %>%
    left_join(group_tab, by = c("sample" = "sample-id")) %>%
    filter(!is.na(group)) %>%
    mutate(group = factor(group, levels = group_levels))  # 设置因子顺序
  
  cat("有效样本数:", nrow(nmds_coord), "\n")
  cat("分组分布（严格按照IS→AS→ES→NS顺序）:\n")
  print(table(nmds_coord$group))
  
  # 6. 统计检验 - 整体PERMANOVA
  valid_samples <- nmds_coord$sample
  bc_dist_subset <- as.dist(as.matrix(bc_dist)[valid_samples, valid_samples])
  
  if(length(unique(nmds_coord$group)) >= 2) {
    adonis_mod <- adonis2(bc_dist_subset ~ group, data = nmds_coord, method = "bray")
    p_manova <- adonis_mod$Pr[1]
    r2_manova <- adonis_mod$R2[1]
    cat("整体PERMANOVA - p值:", p_manova, "R²:", r2_manova, "\n")
  } else {
    p_manova <- NA
    r2_manova <- NA
  }
  
  # 7. 组间两两比较（手动方法，不依赖pairwiseAdonis）
  cat("\n=== 组间两两比较 ===\n")
  
  pairwise_pvalues <- data.frame()
  
  for(i in 1:(length(group_levels)-1)) {
    for(j in (i+1):length(group_levels)) {
      group1 <- group_levels[i]
      group2 <- group_levels[j]
      
      # 筛选当前比较的两组样本
      subset_data <- nmds_coord %>% filter(group %in% c(group1, group2))
      
      if(nrow(subset_data) >= 3) {  # 确保每组至少有2个样本
        subset_dist <- as.dist(as.matrix(bc_dist_subset)[subset_data$sample, subset_data$sample])
        
        # 进行PERMANOVA检验
        tryCatch({
          pairwise_test <- adonis2(subset_dist ~ group, data = subset_data, permutations = 999)
          p_val <- pairwise_test$Pr[1]
          f_val <- pairwise_test$F[1]
          r2_val <- pairwise_test$R2[1]
          
          pairwise_pvalues <- rbind(pairwise_pvalues, 
                                    data.frame(Group1 = group1, Group2 = group2, 
                                               p_value = p_val, F_value = f_val, R2 = r2_val))
          cat(sprintf("  %s vs %s: p = %.4f, F = %.3f, R² = %.3f\n", 
                      group1, group2, p_val, f_val, r2_val))
        }, error = function(e) {
          cat(sprintf("  %s vs %s: 计算失败 - %s\n", group1, group2, e$message))
        })
      } else {
        cat(sprintf("  %s vs %s: 样本数不足，跳过\n", group1, group2))
      }
    }
  }
  
  # 8. 可视化 - 带显著性标注（严格分组顺序）
  cat("\n=== 生成带显著性标注的图形 ===\n")
  
  group_col <- c(IS = "#00FF7F", AS = "#FF1493", ES = "#FFA500", NS = "#00FFFF")
  
  # 统计标签
  if(!is.na(p_manova)) {
    stat_lab <- sprintf("stress = %.3f\nOverall PERMANOVA p = %.3g\nR² = %.3f", 
                        stress, p_manova, r2_manova)
  } else {
    stat_lab <- sprintf("stress = %.3f", stress)
  }
  
  # 图形1: NMDS1箱线图 - 带显著性标注（严格分组顺序）
  p_nmds1 <- ggplot(nmds_coord, aes(x = group, y = NMDS1, fill = group)) +
    geom_boxplot(width = 0.5, alpha = 0.9, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2, shape = 21, color = "black") +
    scale_fill_manual(values = group_col) +
    scale_x_discrete(limits = group_levels) +  # 强制使用指定顺序
    labs(title = "NMDS1 Coordinate with Pairwise Significance (Bray-Curtis)",
         subtitle = "Group order: IS (Initial Seep Age) → AS → ES → NS",
         x = "Group", y = "NMDS1") +
    theme_bw(base_size = 15) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 添加显著性标注
  if(nrow(pairwise_pvalues) > 0) {
    # 计算y轴位置
    y_range <- range(nmds_coord$NMDS1)
    y_max <- max(nmds_coord$NMDS1)
    step <- (y_range[2] - y_range[1]) * 0.15
    
    # 为每对比较添加标注
    for(i in 1:nrow(pairwise_pvalues)) {
      group1 <- pairwise_pvalues$Group1[i]
      group2 <- pairwise_pvalues$Group2[i]
      p_val <- pairwise_pvalues$p_value[i]
      
      # 确定显著性水平
      sig_label <- ifelse(p_val < 0.001, "***",
                          ifelse(p_val < 0.01, "**",
                                 ifelse(p_val < 0.05, "*", "ns")))
      
      if(p_val < 0.05) {  # 只标注显著的结果
        y_position <- y_max + step * i
        
        p_nmds1 <- p_nmds1 + 
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
  }
  
  # 添加统计信息标注
  p_nmds1 <- p_nmds1 + 
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = stat_lab, size = 4, colour = "black")
  
  # 图形2: NMDS散点图 - 在标题中显示显著性信息
  scatter_title <- if(nrow(pairwise_pvalues) > 0) {
    sig_pairs <- pairwise_pvalues %>% filter(p_value < 0.05)
    if(nrow(sig_pairs) > 0) {
      sig_text <- paste(apply(sig_pairs, 1, function(x) {
        sprintf("%s-%s(p=%.3f)", x[1], x[2], as.numeric(x[3]))
      }), collapse = "; ")
      paste("NMDS Plot - Significant pairs:", sig_text)
    } else {
      "NMDS Plot (Bray-Curtis) - No significant pairs"
    }
  } else {
    "NMDS Plot (Bray-Curtis)"
  }
  
  p_nmds_scatter <- ggplot(nmds_coord, aes(x = NMDS1, y = NMDS2, color = group)) +
    geom_point(size = 3) +
    scale_color_manual(values = group_col, 
                       breaks = group_levels,  # 确保图例顺序正确
                       labels = c("IS (Initial Seep Age)", "AS", "ES", "NS")) +
    labs(title = scatter_title,
         subtitle = "Group order: IS → AS → ES → NS",
         x = paste("NMDS1 (stress =", stress, ")"),
         y = "NMDS2",
         color = "Group") +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  if(length(unique(nmds_coord$group)) >= 2) {
    p_nmds_scatter <- p_nmds_scatter + stat_ellipse(level = 0.68)
  }
  
  # 9. 保存图形（PDF和SVG矢量图）
  ggsave("beta-NMDS1-box-with-significance.pdf", p_nmds1, width = 10, height = 8)
  ggsave("beta-NMDS1-box-with-significance.svg", p_nmds1, width = 10, height = 8)
  ggsave("beta-NMDS-scatter-with-significance.pdf", p_nmds_scatter, width = 12, height = 8)
  ggsave("beta-NMDS-scatter-with-significance.svg", p_nmds_scatter, width = 12, height = 8)
  
  # 10. 保存显著性结果表格
  if(nrow(pairwise_pvalues) > 0) {
    # 添加显著性标注
    pairwise_pvalues <- pairwise_pvalues %>%
      mutate(Significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**", 
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ))
    
    write.csv(pairwise_pvalues, "pairwise_permanova_results.csv", row.names = FALSE)
    cat("两两比较结果已保存至: pairwise_permanova_results.csv\n")
  }
  
  # 11. 输出总结
  cat("\n=== 分析完成 ===\n")
  cat("📊 NMDS stress值:", stress, "\n")
  if(!is.na(p_manova)) {
    cat("📈 整体PERMANOVA - p值:", p_manova, "R²:", r2_manova, "\n")
  }
  cat("👥 有效样本数:", nrow(nmds_coord), "\n")
  cat("🔄 分组显示顺序: IS (Initial Seep Age) → AS → ES → NS\n")
  cat("🎨 生成文件:\n")
  cat("   - beta-NMDS1-box-with-significance.pdf/.svg (带显著性标注的箱线图)\n")
  cat("   - beta-NMDS-scatter-with-significance.pdf/.svg (散点图)\n")
  if(nrow(pairwise_pvalues) > 0) {
    cat("   - pairwise_permanova_results.csv (两两比较结果表格)\n")
  }
  
  # 返回完整结果
  return(list(
    nmds_coord = nmds_coord,
    bc_dist = bc_dist,
    stress = stress,
    p_manova = p_manova,
    r2_manova = r2_manova,
    pairwise_results = pairwise_pvalues,
    plot_box = p_nmds1,
    plot_scatter = p_nmds_scatter
  ))
}

# 使用方法：
# 1. 安装依赖包: install_if_missing()
# 2. 运行分析: results <- beta_diversity_analysis_with_significance()
# 3. 查看结果: print(results$plot_box)
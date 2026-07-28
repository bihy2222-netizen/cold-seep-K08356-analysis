# 0. 工作目录 ---------------------------------------------------
setwd("/Users/catherine/Downloads/metaphlan 作图/class 级别作图-柱状图")

# 检查当前目录文件
cat("当前目录中的文件：\n")
list.files()

# 1. 装包 & 加载 ---------------------------------------------------------
pkgs <- c("tidyverse", "vegan", "ggplot2")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
}))
lapply(pkgs, library, character.only = TRUE)

# 2. 直接读取CSV文件 -----------------------------------------------------
# 读取物种表
tax_tab <- read.csv("class34.csv", row.names = 1, check.names = FALSE) %>%
  as.matrix()

# 检查分组文件是否存在，如果不存在则创建示例分组
if(file.exists("class_group.csv")) {
  group_tab <- read.csv("class_group.csv", check.names = FALSE)
} else {
  cat("分组文件不存在，基于样本名称创建示例分组...\n")
  # 基于样本名称创建分组
  samples <- colnames(tax_tab)
  # 假设分组为IS, AS, ES, NS，根据样本数量均匀分配
  groups <- rep(c("IS", "AS", "ES", "NS"), length.out = length(samples))
  group_tab <- data.frame(
    "sample-id" = samples,
    group = groups,
    check.names = FALSE
  )
  # 保存示例分组文件供参考
  write.csv(group_tab, "class_group_example.csv", row.names = FALSE)
  cat("已创建示例分组文件：class_group_example.csv\n")
}

# 3. Alpha 多样性计算 ----------------------------------------------------
Shannon  <- vegan::diversity(tax_tab, index = "shannon", MARGIN = 2)
Simpson  <- vegan::diversity(tax_tab, index = "simpson", MARGIN = 2)
Richness <- vegan::specnumber(tax_tab, MARGIN = 2)

# 计算Chao1指数
obs_chao_ace <- lapply(seq_len(ncol(tax_tab)), function(i) {
  tmp <- vegan::estimateR(ceiling(tax_tab[, i]))
  data.frame(S.obs = tmp["S.obs"], S.chao1 = tmp["S.chao1"], S.ace = tmp["S.ace"])
}) %>% bind_rows()

alpha_df <- data.frame(
  sample   = colnames(tax_tab),
  Shannon  = Shannon,
  Simpson  = Simpson,
  Richness = Richness,
  Chao     = obs_chao_ace$S.chao1
) %>%
  left_join(group_tab, by = c("sample" = "sample-id")) %>%
  filter(!is.na(group))  # 去掉NA分组

# 显示分组信息
cat("分组信息：\n")
table(alpha_df$group)

# 4. 统计检验函数 --------------------------------------------------------
calculate_pvalues <- function(data, value_var, group_var) {
  groups <- unique(data[[group_var]])
  pvalues <- list()
  
  for(i in 1:(length(groups)-1)) {
    for(j in (i+1):length(groups)) {
      group1 <- groups[i]
      group2 <- groups[j]
      
      data1 <- data[data[[group_var]] == group1, value_var]
      data2 <- data[data[[group_var]] == group2, value_var]
      
      # 使用Wilcoxon检验（非参数检验）
      pval <- wilcox.test(data1, data2)$p.value
      
      comparison <- paste(group1, "vs", group2)
      pvalues[[comparison]] <- pval
    }
  }
  return(pvalues)
}

# 5. 计算各指数的p值 -----------------------------------------------------
indices <- c("Shannon", "Simpson", "Chao")
pvalue_results <- list()

for(index in indices) {
  pvalue_results[[index]] <- calculate_pvalues(alpha_df, index, "group")
}

# 6. 导出p值结果 ---------------------------------------------------------
pvalue_df <- do.call(rbind, lapply(names(pvalue_results), function(index) {
  comparisons <- names(pvalue_results[[index]])
  data.frame(
    Index = index,
    Comparison = comparisons,
    Pvalue = unlist(pvalue_results[[index]][comparisons]),
    Significance = ifelse(unlist(pvalue_results[[index]][comparisons]) < 0.001, "***",
                          ifelse(unlist(pvalue_results[[index]][comparisons]) < 0.01, "**",
                                 ifelse(unlist(pvalue_results[[index]][comparisons]) < 0.05, "*", "ns")))
  )
}))

write.csv(pvalue_df, "alpha_diversity_pvalues.csv", row.names = FALSE)

# 7. 生成公式和说明文档 ---------------------------------------------------
formula_doc <- data.frame(
  Index = c("Shannon", "Simpson", "Chao1"),
  Formula = c(
    "H = -∑(p_i × ln(p_i))，其中p_i是第i个物种的相对丰度",
    "D = 1 - ∑(p_i^2)，其中p_i是第i个物种的相对丰度", 
    "S_chao1 = S_obs + (n1^2)/(2×n2)，其中S_obs是观测物种数，n1是单例物种数，n2是双例物种数"
  ),
  Description = c(
    "Shannon指数：衡量群落多样性，考虑物种丰富度和均匀度",
    "Simpson指数：衡量群落优势度，值越小多样性越高",
    "Chao1指数：估计群落中实际物种总数，考虑稀有物种"
  )
)

write.csv(formula_doc, "alpha_diversity_formulas.csv", row.names = FALSE)

# 8. 设置颜色 -----------------------------------------------------------
group_col <- c(IS = "#00FF7F", AS = "#FF1493", ES = "#FFA500", NS = "#00FFFF")

# 9. 手动添加p值标注的函数 ----------------------------------------------
add_pvalue_annotation <- function(plot, pvalues, y_position_ratio = 0.95) {
  # 获取绘图数据
  plot_data <- ggplot_build(plot)
  y_range <- plot_data$layout$panel_params[[1]]$y.range
  y_span <- y_range[2] - y_range[1]
  
  # 计算标注位置
  max_y <- max(alpha_df[[pvalues$Index[1]]], na.rm = TRUE)
  base_y <- max_y * 1.1
  
  # 为每个比较添加标注
  comparisons <- unique(pvalues$Comparison)
  
  for(comp in comparisons) {
    comp_pvals <- pvalues[pvalues$Comparison == comp, ]
    
    # 获取组名和位置
    groups <- strsplit(comp, " vs ")[[1]]
    group_levels <- levels(factor(alpha_df$group))
    x1 <- which(group_levels == groups[1])
    x2 <- which(group_levels == groups[2])
    
    for(i in 1:nrow(comp_pvals)) {
      index <- comp_pvals$Index[i]
      sig <- comp_pvals$Significance[i]
      
      # 计算该指数在数据中的位置
      index_data <- alpha_df[[index]]
      index_max <- max(index_data, na.rm = TRUE)
      y_pos <- base_y + (which(indices == index) - 1) * index_max * 0.15
      
      # 添加线段和标注
      plot <- plot +
        geom_segment(x = x1, xend = x2, y = y_pos, yend = y_pos, 
                     color = "black", linewidth = 0.5) +
        geom_text(x = (x1 + x2) / 2, y = y_pos * 1.02, label = sig, 
                  size = 4, fontface = "bold")
    }
  }
  
  return(plot)
}

# 10. 创建单个指数的箱线图 -----------------------------------------------
create_alpha_plot <- function(data, index, ylab, pvalues_df) {
  # 筛选该指数的p值
  index_pvals <- pvalues_df[pvalues_df$Index == index, ]
  
  # 创建基础图
  p <- ggplot(data, aes(x = group, y = .data[[index]], fill = group)) +
    geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2, shape = 21, color = "black") +
    scale_fill_manual(values = group_col) +
    labs(title = paste(ylab, "Diversity Index"),
         x = "Group", y = ylab) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # 添加p值标注
  if(nrow(index_pvals) > 0) {
    p <- add_pvalue_annotation(p, index_pvals)
  }
  
  return(p)
}

# 11. 生成三个指数的图（保存为PDF和SVG矢量图）---------------------------
p_shannon <- create_alpha_plot(alpha_df, "Shannon", "Shannon", pvalue_df)
p_simpson <- create_alpha_plot(alpha_df, "Simpson", "Simpson", pvalue_df)
p_chao <- create_alpha_plot(alpha_df, "Chao", "Chao1", pvalue_df)

# 保存为PDF（矢量图）
ggsave("Shannon_diversity_with_pvalues.pdf", p_shannon, width = 10, height = 8)
ggsave("Simpson_diversity_with_pvalues.pdf", p_simpson, width = 10, height = 8)
ggsave("Chao1_diversity_with_pvalues.pdf", p_chao, width = 10, height = 8)

# 同时保存为SVG（另一种矢量图格式）
ggsave("Shannon_diversity_with_pvalues.svg", p_shannon, width = 10, height = 8)
ggsave("Simpson_diversity_with_pvalues.svg", p_simpson, width = 10, height = 8)
ggsave("Chao1_diversity_with_pvalues.svg", p_chao, width = 10, height = 8)

# 可选：同时保存高分辨率PNG用于预览
ggsave("Shannon_diversity_with_pvalues.png", p_shannon, width = 10, height = 8, dpi = 300)
ggsave("Simpson_diversity_with_pvalues.png", p_simpson, width = 10, height = 8, dpi = 300)
ggsave("Chao1_diversity_with_pvalues.png", p_chao, width = 10, height = 8, dpi = 300)

# 12. 三个指数的分面图（保存为矢量图）-----------------------------------
alpha_long <- alpha_df %>%
  pivot_longer(c(Shannon, Simpson, Chao), 
               names_to = "index", values_to = "value") %>%
  filter(!is.na(value))

# 设置因子水平以便正确排序
alpha_long$index <- factor(alpha_long$index, 
                           levels = c("Shannon", "Simpson", "Chao"),
                           labels = c("Shannon", "Simpson", "Chao1"))

p_total <- ggplot(alpha_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, shape = 21, color = "black") +
  scale_fill_manual(values = group_col) +
  facet_wrap(~index, scales = "free_y", nrow = 1) +
  labs(title = "Alpha Diversity Indices Comparison",
       x = "Group", y = "Diversity Index Value") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 保存分面图为矢量图
ggsave("combined_alpha_diversity.pdf", p_total, width = 14, height = 8)
ggsave("combined_alpha_diversity.svg", p_total, width = 14, height = 8)
# 可选：保存PNG用于预览
ggsave("combined_alpha_diversity.png", p_total, width = 14, height = 8, dpi = 300)

# 13. 导出最终的alpha多样性数据 -------------------------------------------
write.csv(alpha_df, "final_alpha_diversity_data.csv", row.names = FALSE)

# 14. 创建汇总报告PDF ----------------------------------------------------
# 安装并加载gridExtra包来组合多个表格到PDF
if(!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library(gridExtra)

# 创建汇总表格的PDF
pdf("alpha_diversity_summary_tables.pdf", width = 11, height = 8.5)

# 添加标题页
grid::grid.newpage()
grid::grid.text("Alpha Diversity Analysis Summary", 
                gp = grid::gpar(fontsize = 20, fontface = "bold"),
                y = 0.8)

# 添加分组信息
group_summary <- table(alpha_df$group)
group_df <- data.frame(Group = names(group_summary), Count = as.numeric(group_summary))
grid::grid.text("Sample Group Distribution", gp = grid::gpar(fontsize = 16), y = 0.6)
gridExtra::grid.table(group_df, rows = NULL)

# 添加多样性指数描述
grid::grid.newpage()
grid::grid.text("Alpha Diversity Index Formulas", 
                gp = grid::gpar(fontsize = 16, fontface = "bold"),
                y = 0.9)
gridExtra::grid.table(formula_doc, rows = NULL)

# 添加p值结果
grid::grid.newpage()
grid::grid.text("Statistical Test Results (P-values)", 
                gp = grid::gpar(fontsize = 16, fontface = "bold"),
                y = 0.9)
gridExtra::grid.table(pvalue_df, rows = NULL)

dev.off()

cat("\n分析完成！\n")
cat("生成的文件：\n")
cat("矢量图文件（可无限放大不失真）：\n")
cat("1. Shannon_diversity_with_pvalues.pdf/.svg - Shannon指数矢量图\n")
cat("2. Simpson_diversity_with_pvalues.pdf/.svg - Simpson指数矢量图\n")
cat("3. Chao1_diversity_with_pvalues.pdf/.svg - Chao1指数矢量图\n")
cat("4. combined_alpha_diversity.pdf/.svg - 三个指数组合矢量图\n")
cat("5. alpha_diversity_summary_tables.pdf - 汇总表格PDF报告\n\n")

cat("数据文件：\n")
cat("6. final_alpha_diversity_data.csv - 最终的alpha多样性数据\n")
cat("7. alpha_diversity_pvalues.csv - 组间比较的p值结果\n") 
cat("8. alpha_diversity_formulas.csv - 多样性指数公式说明\n")

cat("预览文件（PNG格式）：\n")
cat("9. Shannon_diversity_with_pvalues.png - Shannon指数预览图\n")
cat("10. Simpson_diversity_with_pvalues.png - Simpson指数预览图\n")
cat("11. Chao1_diversity_with_pvalues.png - Chao1指数预览图\n")
cat("12. combined_alpha_diversity.png - 三个指数组合预览图\n")

if(!file.exists("class_group.csv")) {
  cat("13. class_group_example.csv - 示例分组文件（请根据实际情况修改此文件）\n")
}

cat("\n矢量图格式说明：\n")
cat("- PDF: 适合学术出版和打印，在所有平台上兼容性好\n")
cat("- SVG: 适合网页展示和进一步编辑，可在Illustrator等软件中编辑\n")
cat("- PNG: 适合快速预览和日常使用\n")


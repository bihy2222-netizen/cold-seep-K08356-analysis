# 加载包
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(stringr)
library(cowplot)
library(gridExtra)
library(RColorBrewer)

# 设置工作目录
setwd("/Users/catherine/Downloads/砷循环作图热图和直方图/")

cat("==================================================\n")
cat("开始处理碳氮硫砷功能基因数据\n")
cat("==================================================\n")

# 读取CSV文件
cat("\n读取CNSAS.csv文件...\n")

if (file.exists("CNSAS.csv")) {
  cat("找到CNSAS.csv文件，开始读取...\n")
  
  # 读取CSV文件，保留原始列名
  tax_tab <- read.csv("CNSAS.csv", header = TRUE, stringsAsFactors = FALSE, 
                      check.names = FALSE, encoding = "UTF-8")
  
  if (!is.null(tax_tab)) {
    cat("成功读取CNSAS.csv文件！\n")
  }
} else {
  cat("错误：找不到CNSAS.csv文件！\n")
  tax_tab <- NULL
}

# 检查数据
if (!is.null(tax_tab) && nrow(tax_tab) > 0) {
  
  cat("\n==================================================\n")
  cat("成功读取数据！\n")
  cat("==================================================\n")
  cat("原始数据维度:", dim(tax_tab), "\n")
  
  # 标准化列名
  colnames(tax_tab)[1] <- "KO"
  colnames(tax_tab)[2] <- "Taxon"
  
  # 移除包含#DIV/0!的行
  tax_tab <- tax_tab[!apply(tax_tab, 1, function(x) any(grepl("#DIV/0!", x))), ]
  
  # 获取样本列
  sample_cols <- c("IS", "AS", "ES", "NS")
  sample_cols <- sample_cols[sample_cols %in% colnames(tax_tab)]
  
  cat("\n样本列:", paste(sample_cols, collapse=", "), "\n")
  cat("清理后数据维度:", dim(tax_tab), "\n")
  
  # 转换数值列
  cat("\n转换数值列...\n")
  for (col in sample_cols) {
    tax_tab[[col]] <- as.numeric(as.character(tax_tab[[col]]))
    tax_tab[[col]][is.na(tax_tab[[col]])] <- 0
  }
  
  # ====================================================
  # 定义功能基因
  # ====================================================
  
  # 创建原始基因名到标准化基因名的映射
  tax_tab$KO_clean <- trimws(tax_tab$KO)
  
  # 碳功能基因
  c_genes <- c("pmoA", "pmoB", "pmoC", "mcrA", "mcrB")
  
  # 氮功能基因
  n_genes <- c("nirB", "narG", "nifH", "amoA")
  
  # 硫功能基因
  s_genes <- c("dsrA", "dsrB", "soxB")
  
  # 砷功能基因
  as_genes <- c("acr3", "arsC1", "arsM")
  
  # 所有功能基因分类 - 按C、N、S、AS顺序
  gene_categories <- list(
    "Carbon" = c_genes,
    "Nitrogen" = n_genes,
    "Sulfur" = s_genes,
    "Arsenic" = as_genes
  )
  
  # 找到每个类别中实际存在的基因
  available_genes <- list()
  for (cat in names(gene_categories)) {
    available_genes[[cat]] <- intersect(gene_categories[[cat]], unique(tax_tab$KO_clean))
  }
  
  cat("\n==================================================\n")
  cat("各循环可用功能基因 (按C、N、S、AS顺序):\n")
  for (cat in names(available_genes)) {
    cat(sprintf("%s: %d个基因\n", cat, length(available_genes[[cat]])))
    if (length(available_genes[[cat]]) > 0) {
      cat("  ", paste(available_genes[[cat]], collapse=", "), "\n")
    }
  }
  cat("==================================================\n")
  
  # ====================================================
  # 生成鲜艳颜色的函数（每个循环独立生成，避免重复）
  # ====================================================
  generate_bright_colors <- function(n, seed = NULL) {
    # 使用鲜艳的调色板
    bright_palettes <- list(
      # RColorBrewer鲜艳调色板
      brewer.pal(9, "Set1"),
      brewer.pal(8, "Set2"),
      brewer.pal(12, "Set3"),
      brewer.pal(8, "Dark2"),
      brewer.pal(8, "Accent"),
      brewer.pal(12, "Paired"),
      # 自定义鲜艳颜色
      c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF",
        "#FF4500", "#32CD32", "#1E90FF", "#FF1493", "#FFD700", "#00FA9A",
        "#FF6346", "#7CFC00", "#4169E1", "#FF69B4", "#F0E68C", "#40E0D0",
        "#FF8C00", "#98FB98", "#87CEFA", "#DA70D6", "#BDB76B", "#48D1CC",
        "#DC143C", "#32CD32", "#00BFFF", "#FF00FF", "#FFD700", "#7FFFD4",
        "#FF3030", "#7CFC00", "#4876FF", "#FF6EB4", "#FFC125", "#00E5EE",
        "#FF8247", "#00CD66", "#436EEE", "#EE6AA7", "#EEDD82", "#7FFFD4",
        "#CD5C5C", "#F4A460", "#D2B48C", "#9ACD32", "#40E0D0", "#EE82EE")
    )
    
    all_colors <- unique(unlist(bright_palettes))
    
    if (length(all_colors) < n) {
      # 使用HSV生成额外颜色
      golden_angle <- 0.618033988749895
      hues <- (seq(0, n - 1) * golden_angle) %% 1
      hsv_colors <- hsv(hues, s = runif(n, 0.8, 1), v = runif(n, 0.8, 1))
      all_colors <- c(all_colors, hsv_colors)
    }
    
    # 如果不指定seed，则完全随机
    if (!is.null(seed)) {
      set.seed(seed)
    }
    all_colors <- sample(all_colors)
    
    return(all_colors[1:n])
  }
  
  # ====================================================
  # 设置样本顺序
  # ====================================================
  sample_cols_ordered <- sample_cols  # IS, AS, ES, NS
  
  # ====================================================
  # 为每个循环独立分配颜色
  # ====================================================
  cycle_colors <- list()
  cycle_taxa <- list()
  
  for (cat in names(available_genes)) {
    genes <- available_genes[[cat]]
    
    # 获取该类循环中所有出现的分类单元
    cat_taxa <- unique(tax_tab[tax_tab$KO_clean %in% genes, "Taxon"])
    cycle_taxa[[cat]] <- cat_taxa
    
    # 为每个循环独立生成颜色（使用不同的seed确保颜色不同）
    seed_map <- list("Carbon" = 123, "Nitrogen" = 456, "Sulfur" = 789, "Arsenic" = 101)
    cycle_colors[[cat]] <- generate_bright_colors(length(cat_taxa), seed = seed_map[[cat]])
    names(cycle_colors[[cat]]) <- cat_taxa
    
    cat(sprintf("\n%s循环: %d个分类单元，颜色分配完成\n", cat, length(cat_taxa)))
  }
  
  # ====================================================
  # 创建每个循环的图形（去掉柱子黑线，调整灰色框位置）
  # ====================================================
  
  # 存储每个循环的组合图
  cycle_combined_plots <- list()
  
  for (cat in names(available_genes)) {
    genes <- available_genes[[cat]]
    if (length(genes) == 0) next
    
    cat(sprintf("\n==================================================\n"))
    cat(sprintf("处理%s循环 (%d个基因)\n", cat, length(genes)))
    cat(sprintf("==================================================\n"))
    
    cat_plots <- list()
    
    for (i in 1:length(genes)) {
      gene <- genes[i]
      cat(sprintf("  处理基因 %d/%d: %s\n", i, length(genes), gene))
      
      gene_data <- tax_tab[tax_tab$KO_clean == gene, ]
      if (nrow(gene_data) == 0) {
        cat(sprintf("    警告: %s 没有数据\n", gene))
        next
      }
      
      # 转换为长格式
      plot_data <- gene_data %>%
        pivot_longer(
          cols = all_of(sample_cols_ordered),
          names_to = "Sample",
          values_to = "Abundance"
        )
      
      # 计算百分比
      plot_data <- plot_data %>%
        group_by(Sample) %>%
        mutate(Percentage = Abundance / sum(Abundance) * 100) %>%
        ungroup()
      
      # 将Sample转换为因子
      plot_data$Sample <- factor(plot_data$Sample, levels = sample_cols_ordered)
      
      # 按分类单元的总丰度排序
      taxon_order <- plot_data %>%
        group_by(Taxon) %>%
        summarise(Total = sum(Abundance)) %>%
        arrange(desc(Total)) %>%
        pull(Taxon)
      
      plot_data$Taxon <- factor(plot_data$Taxon, levels = taxon_order)
      
      # 创建图形 - 去掉柱子黑线，调整灰色框位置
      p <- ggplot(plot_data, aes(x = Sample, y = Percentage, fill = Taxon)) +
        # 去掉柱子黑线 (去掉color和linewidth参数)
        geom_col(width = 0.7) +
        # 灰色框放在柱子顶部
        annotate("rect",
                 xmin = 0.5, xmax = length(sample_cols_ordered) + 0.5,
                 ymin = 102, ymax = 108,  # 上移灰色框
                 fill = "#E0E0E0", color = "black", linewidth = 0.3) +
        annotate("text",
                 x = length(sample_cols_ordered)/2 + 0.5,
                 y = 105,  # 调整文字位置到灰色框中间
                 label = gene,
                 size = 5,
                 fontface = "bold") +
        scale_fill_manual(values = cycle_colors[[cat]]) +
        scale_y_continuous(breaks = seq(0, 100, by = 20),
                           labels = paste0(seq(0, 100, by = 20), "%"),
                           limits = c(0, 108),  # 扩大Y轴范围容纳灰色框
                           expand = c(0, 0)) +
        scale_x_discrete(expand = c(0, 0.5)) +
        labs(y = NULL, x = NULL) +
        theme_minimal() +
        theme(
          axis.line = element_line(color = "black", linewidth = 0.2),
          axis.ticks = element_line(color = "black", linewidth = 0.2),
          axis.ticks.length = unit(0.15, "cm"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
          axis.text.y = element_text(size = 9, color = "black"),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.2),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          legend.position = "none",
          plot.margin = margin(2, 2, 2, 2)
        )
      
      cat_plots[[i]] <- p
    }
    
    if (length(cat_plots) == 0) {
      cat(sprintf("\n%s循环没有生成任何图形\n", cat))
      next
    }
    
    # ====================================================
    # 组合该类循环的图形（水平排列，间隔0.5）
    # ====================================================
    
    # 在每个图形之间添加间隔
    if (length(cat_plots) > 1) {
      # 创建间隔图
      spacer <- ggplot() + theme_void() + theme(plot.margin = margin(0, 2.5, 0, 2.5))
      
      # 在图形之间插入间隔
      plots_with_spacers <- list()
      for (j in 1:length(cat_plots)) {
        plots_with_spacers[[length(plots_with_spacers) + 1]] <- cat_plots[[j]]
        if (j < length(cat_plots)) {
          plots_with_spacers[[length(plots_with_spacers) + 1]] <- spacer
        }
      }
      
      # 计算相对宽度：每个基因图3.5cm，间隔0.5cm
      rel_widths <- c()
      for (j in 1:length(cat_plots)) {
        rel_widths <- c(rel_widths, 3.5)
        if (j < length(cat_plots)) {
          rel_widths <- c(rel_widths, 0.5)
        }
      }
      
      cat_grid <- plot_grid(plotlist = plots_with_spacers, 
                            nrow = 1, 
                            rel_widths = rel_widths,
                            align = "h")
    } else {
      cat_grid <- cat_plots[[1]]
    }
    
    # ====================================================
    # 创建该类循环的图例
    # ====================================================
    
    # 获取该类循环的分类单元和颜色
    cat_taxa <- cycle_taxa[[cat]]
    cat_colors <- cycle_colors[[cat]]
    
    # 按字母顺序排序图例
    cat_taxa_sorted <- sort(cat_taxa)
    
    legend_data <- data.frame(
      Sample = rep(sample_cols_ordered[1], length(cat_taxa_sorted)),
      Percentage = rep(100/length(cat_taxa_sorted), length(cat_taxa_sorted)),
      Taxon = cat_taxa_sorted
    )
    
    legend_plot <- ggplot(legend_data, aes(x = Sample, y = Percentage, fill = Taxon)) +
      geom_col() +
      scale_fill_manual(values = cat_colors[cat_taxa_sorted], 
                        name = paste(cat, "Cycle Taxa"),
                        guide = guide_legend(ncol = 1,
                                             keyheight = 0.5,
                                             keywidth = 0.8,
                                             byrow = TRUE)) +
      theme_void() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 7, color = "black"),
        legend.key.size = unit(0.4, "cm"),
        legend.spacing.y = unit(0.05, "cm"),
        legend.box.margin = margin(0, 0, 0, 5)
      )
    
    legend_grob <- get_legend(legend_plot)
    
    # ====================================================
    # 组合图形和图例
    # ====================================================
    
    # 计算宽度：基因图总宽度 + 图例宽度
    genes_width <- length(cat_plots) * 3.5 + (length(cat_plots) - 1) * 0.5
    legend_width <- 5
    total_width <- genes_width + legend_width
    
    # 创建左侧标签（循环名称）
    label_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = paste(cat, "Genes"), 
               angle = 90, size = 5, fontface = "bold", color = "black") +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))
    
    row_with_label <- plot_grid(label_plot, cat_grid, 
                                nrow = 1, 
                                rel_widths = c(0.5, genes_width))
    
    combined_plot <- plot_grid(row_with_label, legend_grob, 
                               nrow = 1, 
                               rel_widths = c(genes_width + 0.5, legend_width))
    
    # 存储
    cycle_combined_plots[[cat]] <- list(
      plot = combined_plot,
      width = total_width,
      height = 12,
      n_genes = length(cat_plots)
    )
    
    # ====================================================
    # 保存单个循环的图
    # ====================================================
    
    filename_pdf <- sprintf("%s_Cycle_Plot.pdf", cat)
    filename_png <- sprintf("%s_Cycle_Plot.png", cat)
    
    ggsave(filename_pdf, combined_plot, 
           width = total_width, height = 12, units = "cm", dpi = 600, limitsize = FALSE)
    ggsave(filename_png, combined_plot, 
           width = total_width, height = 12, units = "cm", dpi = 600, bg = "white", limitsize = FALSE)
    
    cat(sprintf("\n✅ %s循环图形保存完成！宽度: %.1f cm, %d个基因\n", 
                cat, total_width, length(cat_plots)))
  }
  
  # ====================================================
  # 组合所有循环 - 一列组合图（C、N、S、AS顺序）
  # ====================================================
  
  cat("\n==================================================\n")
  cat("组合所有循环图形 - 一列（C、N、S、AS顺序）\n")
  cat("==================================================\n")
  
  # 按C、N、S、AS顺序提取图形
  final_order <- c("Carbon", "Nitrogen", "Sulfur", "Arsenic")
  all_plots <- list()
  all_heights <- c()
  
  for (cat in final_order) {
    if (cat %in% names(cycle_combined_plots)) {
      all_plots[[cat]] <- cycle_combined_plots[[cat]]$plot
      all_heights <- c(all_heights, cycle_combined_plots[[cat]]$height)
    }
  }
  
  if (length(all_plots) > 0) {
    # 垂直排列所有循环，添加少量间隔
    combined_all_onecol <- plot_grid(plotlist = all_plots, 
                                     ncol = 1, 
                                     align = "v",
                                     rel_heights = rep(1, length(all_plots)))
    
    # 计算总高度
    total_height_onecol <- sum(all_heights)
    
    # 取最大宽度作为组合图宽度
    max_width_onecol <- max(sapply(cycle_combined_plots, function(x) x$width))
    
    # 保存一列组合图
    ggsave("All_Cycles_Combined_OneColumn.pdf", combined_all_onecol, 
           width = max_width_onecol, height = total_height_onecol, 
           units = "cm", dpi = 600, limitsize = FALSE)
    ggsave("All_Cycles_Combined_OneColumn.png", combined_all_onecol, 
           width = max_width_onecol, height = total_height_onecol, 
           units = "cm", dpi = 600, bg = "white", limitsize = FALSE)
    
    cat(sprintf("\n✅ 一列组合图保存完成！宽度: %.1f cm, 高度: %.1f cm\n", 
                max_width_onecol, total_height_onecol))
  }
  
  # ====================================================
  # 组合所有循环 - 两列组合图
  # ====================================================
  
  cat("\n==================================================\n")
  cat("组合所有循环图形 - 两列（C+N 左列，S+AS 右列）\n")
  cat("==================================================\n")
  
  # 分成两列：左列(C+N)，右列(S+AS)
  left_col_cats <- c("Carbon", "Nitrogen")
  right_col_cats <- c("Sulfur", "Arsenic")
  
  left_plots <- list()
  right_plots <- list()
  left_heights <- c()
  right_heights <- c()
  
  # 收集左列图形
  for (cat in left_col_cats) {
    if (cat %in% names(cycle_combined_plots)) {
      left_plots[[cat]] <- cycle_combined_plots[[cat]]$plot
      left_heights <- c(left_heights, cycle_combined_plots[[cat]]$height)
    }
  }
  
  # 收集右列图形
  for (cat in right_col_cats) {
    if (cat %in% names(cycle_combined_plots)) {
      right_plots[[cat]] <- cycle_combined_plots[[cat]]$plot
      right_heights <- c(right_heights, cycle_combined_plots[[cat]]$height)
    }
  }
  
  if (length(left_plots) > 0 && length(right_plots) > 0) {
    # 垂直排列左列
    left_column <- plot_grid(plotlist = left_plots, 
                             ncol = 1, 
                             align = "v",
                             rel_heights = rep(1, length(left_plots)))
    
    # 垂直排列右列
    right_column <- plot_grid(plotlist = right_plots, 
                              ncol = 1, 
                              align = "v",
                              rel_heights = rep(1, length(right_plots)))
    
    # 水平排列两列
    combined_all_twocol <- plot_grid(left_column, right_column, 
                                     nrow = 1,
                                     rel_widths = c(1, 1))
    
    # 计算两列组合图的尺寸
    max_width_left <- if(length(left_plots) > 0) 
      max(sapply(cycle_combined_plots[left_col_cats[left_col_cats %in% names(cycle_combined_plots)]], function(x) x$width)) else 0
    max_width_right <- if(length(right_plots) > 0) 
      max(sapply(cycle_combined_plots[right_col_cats[right_col_cats %in% names(cycle_combined_plots)]], function(x) x$width)) else 0
    
    twocol_width <- max_width_left + max_width_right
    twocol_height <- max(sum(left_heights), sum(right_heights))
    
    # 保存两列组合图
    ggsave("All_Cycles_Combined_TwoColumn.pdf", combined_all_twocol, 
           width = twocol_width, height = twocol_height, 
           units = "cm", dpi = 600, limitsize = FALSE)
    ggsave("All_Cycles_Combined_TwoColumn.png", combined_all_twocol, 
           width = twocol_width, height = twocol_height, 
           units = "cm", dpi = 600, bg = "white", limitsize = FALSE)
    
    cat(sprintf("\n✅ 两列组合图保存完成！宽度: %.1f cm, 高度: %.1f cm\n", 
                twocol_width, twocol_height))
  }
  
  # ====================================================
  # 总结
  # ====================================================
  cat("\n==================================================\n")
  cat("所有图形生成完毕！\n")
  cat("==================================================\n")
  cat("\n生成的文件：\n")
  cat("\n【单个循环图】\n")
  cat("  - Carbon_Cycle_Plot.pdf/png\n")
  cat("  - Nitrogen_Cycle_Plot.pdf/png\n")
  cat("  - Sulfur_Cycle_Plot.pdf/png\n")
  cat("  - Arsenic_Cycle_Plot.pdf/png\n")
  cat("\n【组合图】\n")
  cat("  - All_Cycles_Combined_OneColumn.pdf/png （一列，C、N、S、AS顺序）\n")
  cat("  - All_Cycles_Combined_TwoColumn.pdf/png （两列，C+N 左列，S+AS 右列）\n")
  cat("\n图形特点：\n")
  cat("  - ✅ 柱子上的黑线已移除\n")
  cat("  - ✅ 灰色框和基因名称上移，参考原代码位置\n")
  cat("  - 每个循环有自己的图例，颜色独立分配不重复\n")
  cat("  - 图形顺序：Carbon > Nitrogen > Sulfur > Arsenic\n")
  cat("  - 每个功能基因柱子图之间间隔0.5cm\n")
  cat("  - 省略中间坐标轴名称，缩紧空间\n")
  cat("  - 每个基因方块宽度: 3.5cm\n")
  cat("  - 颜色: 每个循环独立使用鲜艳配色方案\n")
  cat("==================================================\n")
  
} else {
  cat("\n❌ 数据读取失败。请检查文件格式。\n")
  cat("请确保CNSAS.csv文件存在于当前目录\n")
  cat("\n当前工作目录:", getwd(), "\n")
  cat("目录中的文件:\n")
  print(list.files())
}


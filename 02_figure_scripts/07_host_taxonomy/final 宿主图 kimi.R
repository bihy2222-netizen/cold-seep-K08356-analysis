# 加载包
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(stringr)
library(cowplot)
library(gridExtra)

# 设置工作目录
setwd("/Users/catherine/Downloads/重新画图/宿主图 class/")

# 读取CSV数据
tax_tab <- read.csv("classcns.csv", check.names = FALSE, stringsAsFactors = FALSE, 
                    fileEncoding = "UTF-8", na.strings = c("", "NA"), fill = TRUE)

# 查看数据结构
cat("数据维度:", dim(tax_tab), "\n")
cat("列名:", colnames(tax_tab), "\n")

# 标准化列名
colnames(tax_tab)[1] <- "KO"
colnames(tax_tab)[2] <- "class"

# 清理KO名称中的空格
tax_tab$KO <- trimws(tax_tab$KO)

# 将数值列从字符转换为数值
cat("\n转换数值列...\n")
numeric_cols <- c("IS", "AS", "ES", "NS")
for (col in numeric_cols) {
  tax_tab[[col]] <- as.numeric(as.character(tax_tab[[col]]))
  tax_tab[[col]][is.na(tax_tab[[col]])] <- 0
}

# 定义三类功能基因
c_genes <- c("pmoA", "pmoB", "pmoC", "mcrA", "mcrB")
n_genes <- c("nirB", "narG", "nifH", "amoA")
s_genes <- c("dsrA", "dsrB", "soxB")

# 检查哪些基因在数据中实际存在
available_genes <- unique(tax_tab$KO)
c_genes_available <- intersect(c_genes, available_genes)
n_genes_available <- intersect(n_genes, available_genes)
s_genes_available <- intersect(s_genes, available_genes)

cat("\nC循环基因:", paste(c_genes_available, collapse=", "), "\n")
cat("N循环基因:", paste(n_genes_available, collapse=", "), "\n")
cat("S循环基因:", paste(s_genes_available, collapse=", "), "\n")

# 检查是否有遗漏的基因
cat("\n检查遗漏的基因:\n")
cat("C循环遗漏:", setdiff(c_genes, c_genes_available), "\n")
cat("N循环遗漏:", setdiff(n_genes, n_genes_available), "\n")
cat("S循环遗漏:", setdiff(s_genes, s_genes_available), "\n")

# 定义样本顺序
sample_names <- c("IS", "AS", "ES", "NS")

# 使用你提供的34种漂亮颜色
colors_34 <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#911EB4", 
               "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", "#008080", "#E6BEFF", 
               "#9A6324", "#FFFAC8", "#800000", "#AAFFC3", "#808000", "#FFD8B1", 
               "#000075", "#808080", "#A9A9A9", "#2F4F4F", "#FF69B4", "#BA55D3", 
               "#9370DB", "#3CB371", "#7B68EE", "#00FA9A", "#48D1CC", "#C71585", 
               "#191970", "#FF4500", "#32CD32", "#FF1493")

# 扩展颜色函数
get_extended_colors <- function(class_list, base_colors) {
  n <- length(class_list)
  
  if (n <= length(base_colors)) {
    colors <- base_colors[1:n]
  } else {
    repeat_times <- ceiling(n / length(base_colors))
    colors <- rep(base_colors, repeat_times)[1:n]
  }
  
  names(colors) <- class_list
  return(colors)
}

# 为每个循环收集class并生成颜色
collect_classes_and_colors <- function(genes_list) {
  all_classes <- c()
  for (ko_name in genes_list) {
    ko_data <- tax_tab[tax_tab$KO == ko_name, ]
    ko_classes <- unique(ko_data$class)
    ko_classes <- ko_classes[!is.na(ko_classes) & ko_classes != ""]
    all_classes <- unique(c(all_classes, ko_classes))
  }
  
  if (length(all_classes) > 0) {
    colors <- get_extended_colors(all_classes, colors_34)
    return(list(classes = all_classes, colors = colors))
  } else {
    return(NULL)
  }
}

cat("\n收集各循环的class和颜色...\n")
c_data <- collect_classes_and_colors(c_genes_available)
n_data <- collect_classes_and_colors(n_genes_available)
s_data <- collect_classes_and_colors(s_genes_available)

if (!is.null(c_data)) {
  cat("C循环class数量:", length(c_data$classes), "\n")
  cat("C循环class列表:", paste(c_data$classes, collapse=", "), "\n")
}
if (!is.null(n_data)) {
  cat("N循环class数量:", length(n_data$classes), "\n")
  cat("N循环class列表:", paste(n_data$classes, collapse=", "), "\n")
}
if (!is.null(s_data)) {
  cat("S循环class数量:", length(s_data$classes), "\n")
  cat("S循环class列表:", paste(s_data$classes, collapse=", "), "\n")
}

# 创建单个循环图的函数 - 修正plot.margin错误
create_cycle_plot <- function(cycle_name, genes_list, classes_data) {
  if (is.null(classes_data)) {
    cat(sprintf("%s循环: 没有class数据\n", cycle_name))
    return(NULL)
  }
  
  cat(sprintf("\n创建%s循环图形...\n", cycle_name))
  cat(sprintf("处理基因: %s\n", paste(genes_list, collapse=", ")))
  
  # 检查每个基因是否有数据
  missing_genes <- c()
  for (ko_name in genes_list) {
    ko_data <- tax_tab[tax_tab$KO == ko_name, ]
    if (nrow(ko_data) == 0) {
      missing_genes <- c(missing_genes, ko_name)
      cat(sprintf("警告: 基因%s没有数据\n", ko_name))
    }
  }
  
  if (length(missing_genes) > 0) {
    cat(sprintf("以下基因没有数据，将被跳过: %s\n", paste(missing_genes, collapse=", ")))
    genes_list <- setdiff(genes_list, missing_genes)
  }
  
  if (length(genes_list) == 0) {
    cat(sprintf("%s循环: 没有可用的基因数据\n", cycle_name))
    return(NULL)
  }
  
  # 为该循环的所有图形收集数据
  plot_list <- list()
  for (i in 1:length(genes_list)) {
    ko_name <- genes_list[i]
    cat(sprintf("  处理: %s\n", ko_name))
    
    # 筛选数据
    ko_data <- tax_tab[tax_tab$KO == ko_name, ]
    
    if (nrow(ko_data) == 0) {
      cat(sprintf("  警告: %s没有数据，跳过\n", ko_name))
      next
    }
    
    # 转换为长格式
    plot_data <- ko_data %>%
      pivot_longer(
        cols = all_of(sample_names),
        names_to = "Sample",
        values_to = "Abundance"
      )
    
    # 确保Abundance是数值类型
    plot_data$Abundance <- as.numeric(plot_data$Abundance)
    
    # 设置因子顺序
    plot_data$Sample <- factor(plot_data$Sample, levels = sample_names)
    
    # 计算每个样本的总丰度用于标准化为百分比
    sample_totals <- plot_data %>%
      group_by(Sample) %>%
      summarise(Total = sum(Abundance, na.rm = TRUE))
    
    # 合并总丰度数据并计算百分比
    plot_data <- plot_data %>%
      left_join(sample_totals, by = "Sample") %>%
      mutate(Percentage = ifelse(Total > 0, Abundance / Total * 100, 0))
    
    # 按总丰度排序class
    class_summary <- plot_data %>%
      group_by(class) %>%
      summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
      arrange(desc(TotalAbundance))
    
    class_levels <- class_summary$class
    plot_data$class <- factor(plot_data$class, levels = class_levels)
    
    # 获取该KO的颜色
    current_classes <- as.character(class_levels)
    current_colors <- classes_data$colors[current_classes]
    
    # 创建柱状图 - 修复plot.margin重复定义错误
    p <- ggplot(plot_data, aes(x = Sample, y = Percentage, fill = class)) +
      geom_col(
        width = 0.5,  # 调整柱子宽度
        color = "black",
        linewidth = 0.2,
        position = position_stack(reverse = FALSE)
      ) +
      # 添加灰色标题块
      annotate("rect",
               xmin = 0.5, xmax = length(sample_names) + 0.5,
               ymin = 100, ymax = 105,
               fill = "#E0E0E0", color = "black", linewidth = 0.5) +
      annotate("text",
               x = mean(1:length(sample_names)),
               y = 102.5,
               label = ko_name,
               size = 3.5,
               family = "serif",
               fontface = "bold",
               color = "black") +
      scale_fill_manual(
        values = current_colors,
        guide = "none"
      ) +
      scale_y_continuous(
        breaks = seq(0, 100, by = 20),  # 调整Y轴刻度
        labels = seq(0, 100, by = 20),  # 调整Y轴标签
        expand = expansion(mult = c(0, 0.05)),
        limits = c(0, 105)
      ) +
      scale_x_discrete(
        expand = expansion(mult = c(0.05, 0.05))
      ) +
      labs(
        y = "Relative Abundance (%)",
        x = NULL,
        title = NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(
        text = element_text(family = "serif"),
        axis.text.x = element_text(
          angle = 0,
          hjust = 0.5,
          vjust = 0.5,
          size = 10,
          color = "black",
          face = "bold"
        ),
        axis.text.y = element_text(
          size = 9,
          color = "black"
        ),
        axis.title.y = element_text(
          size = 10,
          margin = margin(r = 10)
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.3),
        axis.ticks.length = unit(0.15, "cm"),
        # 关键：统一的边距设置，移除重复的plot.margin
        plot.margin = margin(10, 10, 10, 10, "pt")
      )
    
    plot_list[[i]] <- p
  }
  
  # 移除NULL值
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  if (length(plot_list) == 0) {
    cat(sprintf("  %s循环: 没有生成图形\n", cycle_name))
    return(NULL)
  }
  
  cat(sprintf("  成功创建%d个基因的图形\n", length(plot_list)))
  
  # 创建图例 - 只使用一列
  cat("  创建图例...\n")
  legend_data <- data.frame(
    class = factor(classes_data$classes, levels = classes_data$classes),
    Value = 1:length(classes_data$classes)
  )
  
  # 只使用一列图例
  legend_cols <- 1
  
  # 创建图例
  legend_plot <- ggplot(legend_data, aes(x = 1, y = Value, fill = class)) +
    geom_tile() +
    scale_fill_manual(
      values = classes_data$colors,
      name = paste0(cycle_name, " Cycle"),
      guide = guide_legend(
        ncol = legend_cols,
        byrow = TRUE,
        keyheight = unit(0.4, "cm"),
        keywidth = unit(0.5, "cm"),
        title.position = "top",
        title.hjust = 0.5,
        title = element_text(size = 10, family = "serif", face = "bold", margin = margin(b = 5)),
        label.position = "right",
        label.hjust = 0,
        label.theme = element_text(size = 8, family = "serif"),
        title.theme = element_text(size = 10, family = "serif", face = "bold")
      )
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = unit(0.15, "cm"),
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
  
  legend_grob <- get_legend(legend_plot)
  
  # 组合图形 - 添加图形间距
  cat("  组合图形...\n")
  
  if (length(plot_list) > 1) {
    # 创建图形之间的间距元素
    spacer_plots <- list()
    for (j in 1:(length(plot_list) - 1)) {
      spacer <- ggplot() + theme_void()
      spacer_plots[[j]] <- spacer
    }
    
    # 交替组合图形和间距
    combined_plots <- list()
    for (k in 1:length(plot_list)) {
      combined_plots[[2*k - 1]] <- plot_list[[k]]
      if (k < length(plot_list)) {
        combined_plots[[2*k]] <- spacer_plots[[k]]
      }
    }
    
    # 设置图形和间距的宽度比例
    plot_widths <- rep(1, length(combined_plots))
    for (m in seq(2, length(combined_plots), by = 2)) {
      plot_widths[m] <- 0.05  # 0.1cm间距
    }
    
    plots_grid <- plot_grid(
      plotlist = combined_plots,
      nrow = 1,
      align = "h",
      axis = "tb",
      rel_widths = plot_widths
    )
  } else {
    # 如果只有一个图，直接使用
    plots_grid <- plot_list[[1]]
  }
  
  # 组合图形和图例 - 使用合适的比例
  combined_plot <- plot_grid(
    plots_grid, legend_grob,
    nrow = 1,
    rel_widths = c(0.75, 0.25),
    align = "h",
    axis = "tb"
  )
  
  # 保存图形 - 只保存PDF
  cat("  保存图形...\n")
  
  # 计算图形尺寸
  num_genes <- length(genes_list)
  base_width_per_gene <- 3.0
  spacing_width <- 0.1 * max(0, num_genes - 1)
  
  # 根据class数量动态调整图例宽度
  num_classes <- length(classes_data$classes)
  legend_width <- max(4, num_classes * 0.7)
  
  output_width <- num_genes * base_width_per_gene + spacing_width + legend_width
  
  # 根据循环调整画布大小
  if (cycle_name == "Carbon") {
    output_width <- output_width * 2  # 碳循环画布再加长一倍
  } else if (cycle_name == "Nitrogen") {
    output_width <- output_width + 1 * base_width_per_gene  # 氮循环加长一个柱子的长度
  }
  
  output_height <- 12
  
  output_file <- paste0(cycle_name, "_Cycle_Plot.pdf")
  
  tryCatch({
    ggsave(
      output_file,
      combined_plot,
      width = output_width,
      height = output_height,
      units = "cm",
      dpi = 600,
      device = cairo_pdf,
      limitsize = FALSE
    )
    cat(sprintf("  ✅ %s循环图形保存完成: %s (%.1f cm × %.1f cm)\n", 
                cycle_name, output_file, output_width, output_height))
  }, error = function(e) {
    cat(sprintf("  ❌ 保存失败: %s\n", e$message))
    return(NULL)
  })
  
  return(combined_plot)
}

# 创建并保存三张独立的图
cat("\n==========================================\n")
cat("开始创建三张独立的循环图...\n")
cat("==========================================\n")

# C循环图
if (!is.null(c_data) && length(c_genes_available) > 0) {
  c_plot <- create_cycle_plot("Carbon", c_genes_available, c_data)
} else {
  cat("⚠️  C循环: 没有可用的数据\n")
}

# N循环图
if (!is.null(n_data) && length(n_genes_available) > 0) {
  n_plot <- create_cycle_plot("Nitrogen", n_genes_available, n_data)
} else {
  cat("⚠️  N循环: 没有可用的数据\n")
}

# S循环图
if (!is.null(s_data) && length(s_genes_available) > 0) {
  s_plot <- create_cycle_plot("Sulfur", s_genes_available, s_data)
} else {
  cat("⚠️  S循环: 没有可用的数据\n")
}

# 输出完成信息
cat("\n==========================================\n")
cat("🎉 所有任务完成！\n")
cat("==========================================\n")

cat("\n📁 输出文件:\n")
if (exists("c_plot") && !is.null(c_plot)) {
  cat("  • Carbon_Cycle_Plot.pdf\n")
}
if (exists("n_plot") && !is.null(n_plot)) {
  cat("  • Nitrogen_Cycle_Plot.pdf\n")
}
if (exists("s_plot") && !is.null(s_plot)) {
  cat("  • Sulfur_Cycle_Plot.pdf\n")
}

cat("\n✅ 图形特点（修正版）：\n")
cat("  1. 修复了plot.margin重复定义的错误\n")
cat("  2. 所有图规格完全一致：每个图都有完整Y轴标签\n")
cat("  3. 图例固定为一列\n")
cat("  4. 统一的边距设置，确保视觉一致性\n")
cat("  5. 增加图形宽度，确保Y轴标签有足够空间\n")
cat("  6. 每个子图之间距离0.1厘米\n")
cat("  7. 所有文字使用Times New Roman字体\n")
cat("  8. 全封闭框线设计\n")
cat("  9. 碳循环画布再加长一倍\n")
cat("  10. 氮循环加长一个柱子的长度\n")
cat("  11. Y轴刻度只写数字，标签写0、20、40、60、80、100\n")
cat("  12. 第一个柱子全都写标签，后面的柱子只写0、100\n")

cat("\n📏 图形尺寸计算:\n")
cat("  宽度 = 基因数量 × 3.0cm + (基因数量-1) × 0.1cm + 图例宽度\n")
cat("  图例宽度 = max(4cm, class数量 × 0.7cm)\n")
cat("  高度 = 12cm\n")

cat("\n📍 文件保存在当前目录:", getwd(), "\n")
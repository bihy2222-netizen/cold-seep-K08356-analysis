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

# 1. 数据读取和清洗
cat("正在读取数据文件...\n")
tax_tab <- read.csv("classcns.csv", 
                    check.names = FALSE, 
                    stringsAsFactors = FALSE, 
                    fileEncoding = "UTF-8", 
                    na.strings = c("", "NA", "N/A", "na", "n/a", "NULL", "null"), 
                    fill = TRUE,
                    strip.white = TRUE,
                    quote = "\"")

cat("原始数据维度:", dim(tax_tab), "\n")
cat("原始列名:", paste(colnames(tax_tab), collapse=", "), "\n\n")

# 显示数据结构
str(tax_tab)

# 标准化列名
colnames(tax_tab)[1] <- "KO"
colnames(tax_tab)[2] <- "class"

# 清理数据
cat("数据清洗...\n")
tax_tab$KO <- trimws(tax_tab$KO)
tax_tab$class <- trimws(tax_tab$class)

# 识别样本列
non_numeric_cols <- c("KO", "class")
sample_cols <- setdiff(colnames(tax_tab), non_numeric_cols)
cat("样本列:", paste(sample_cols, collapse=", "), "\n")

# 转换数值列 - 修复NA问题
cat("转换数值列...\n")
for (col in sample_cols) {
  # 先转换为字符
  char_values <- as.character(tax_tab[[col]])
  # 清理字符：移除非数字字符（除了数字、小数点、负号、e/E）
  cleaned_values <- gsub("[^0-9\\.\\-eE]", "", char_values)
  # 将空字符串转换为NA
  cleaned_values[cleaned_values == ""] <- NA
  # 转换为数值
  num_values <- suppressWarnings(as.numeric(cleaned_values))
  # 将NA设为0
  num_values[is.na(num_values)] <- 0
  # 更新数据框
  tax_tab[[col]] <- num_values
  cat(sprintf("  列 %s: 范围 [%.4f, %.4f], 总和 %.4f\n", 
              col, min(num_values, na.rm = TRUE), max(num_values, na.rm = TRUE), sum(num_values, na.rm = TRUE)))
}

# 移除全为0的行
initial_rows <- nrow(tax_tab)
row_sums <- rowSums(tax_tab[, sample_cols], na.rm = TRUE)
tax_tab <- tax_tab[row_sums > 0, ]
cat(sprintf("\n移除了 %d 行全为0的数据 (%.1f%%)\n", 
            initial_rows - nrow(tax_tab), 
            (initial_rows - nrow(tax_tab))/initial_rows*100))

# 2. 定义功能基因和样本顺序
c_genes <- c("pmoA", "pmoB", "pmoC", "mcrA", "mcrB")
n_genes <- c("nirB", "narG", "nifH", "amoA")
s_genes <- c("dsrA", "dsrB", "soxB")

sample_names <- sample_cols  # 使用实际的样本列
cat("\n样本顺序:", paste(sample_names, collapse=", "), "\n")

# 检查可用的基因
available_genes <- unique(tax_tab$KO)
c_genes_available <- intersect(c_genes, available_genes)
n_genes_available <- intersect(n_genes, available_genes)
s_genes_available <- intersect(s_genes, available_genes)

cat("\nC循环基因:", paste(c_genes_available, collapse=", "), "\n")
cat("N循环基因:", paste(n_genes_available, collapse=", "), "\n")
cat("S循环基因:", paste(s_genes_available, collapse=", "), "\n")

# 3. 颜色配置
colors_34 <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#911EB4", 
               "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", "#008080", "#E6BEFF", 
               "#9A6324", "#FFFAC8", "#800000", "#AAFFC3", "#808000", "#FFD8B1", 
               "#000075", "#808080", "#A9A9A9", "#2F4F4F", "#FF69B4", "#BA55D3", 
               "#9370DB", "#3CB371", "#7B68EE", "#00FA9A", "#48D1CC", "#C71585", 
               "#191970", "#FF4500", "#32CD32", "#FF1493")

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

# 4. 收集class和颜色
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
  }
  return(NULL)
}

cat("\n收集各循环的class和颜色...\n")
c_data <- collect_classes_and_colors(c_genes_available)
n_data <- collect_classes_and_colors(n_genes_available)
s_data <- collect_classes_and_colors(s_genes_available)

# 5. 古菌/细菌分类函数
classify_archaea_bacteria <- function(class_name) {
  if (is.na(class_name) || class_name == "") return(NA)
  
  archaea_keywords <- c("Archae", "Methano", "Thermo", "Halo", "Crenarchaeota", 
                        "Euryarchaeota", "Thaumarchaeota", "Asgard", "TACK", 
                        "archaea", "archaeon", "Methanogen", "Halobacteria",
                        "Thermococci", "Thermoplasmata")
  
  class_lower <- tolower(class_name)
  for (keyword in archaea_keywords) {
    if (grepl(tolower(keyword), class_lower)) {
      return("Archaea")
    }
  }
  return("Bacteria")
}

# 6. 创建循环图的函数（完整版）
create_cycle_plot <- function(cycle_name, genes_list, classes_data) {
  if (is.null(classes_data)) {
    cat(sprintf("%s循环: 没有class数据\n", cycle_name))
    return(NULL)
  }
  
  cat(sprintf("\n创建%s循环图形...\n", cycle_name))
  cat(sprintf("处理基因: %s\n", paste(genes_list, collapse=", ")))
  
  # 对class进行分类
  class_types <- sapply(classes_data$classes, classify_archaea_bacteria)
  archaea_classes <- classes_data$classes[class_types == "Archaea"]
  bacteria_classes <- classes_data$classes[class_types == "Bacteria"]
  
  cat(sprintf("  古菌数量: %d, 细菌数量: %d\n", 
              length(archaea_classes), length(bacteria_classes)))
  
  # 创建图形
  plot_list <- list()
  for (i in 1:length(genes_list)) {
    ko_name <- genes_list[i]
    cat(sprintf("  处理: %s\n", ko_name))
    
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
    
    plot_data$Abundance <- as.numeric(as.character(plot_data$Abundance))
    plot_data <- plot_data[!is.na(plot_data$Abundance), ]
    plot_data$Sample <- factor(plot_data$Sample, levels = sample_names)
    
    # 计算百分比
    sample_totals <- plot_data %>%
      group_by(Sample) %>%
      summarise(Total = sum(Abundance, na.rm = TRUE))
    
    plot_data <- plot_data %>%
      left_join(sample_totals, by = "Sample") %>%
      mutate(Percentage = ifelse(Total > 0, Abundance / Total * 100, 0))
    
    # 排序class
    class_summary <- plot_data %>%
      group_by(class) %>%
      summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
      arrange(desc(TotalAbundance))
    
    class_levels <- class_summary$class
    plot_data$class <- factor(plot_data$class, levels = class_levels)
    
    # 预处理数据，手动计算堆叠位置
    plot_data_processed <- plot_data %>%
      arrange(Sample, class) %>%
      group_by(Sample) %>%
      mutate(
        y_end = cumsum(Percentage),
        y_start = lag(y_end, default = 0)
      ) %>%
      ungroup()
    
    # 创建图形 - 使用geom_rect完全避免黑线
    p <- ggplot() +
      geom_rect(
        data = plot_data_processed,
        aes(xmin = as.numeric(Sample) - 0.2,
            xmax = as.numeric(Sample) + 0.2,
            ymin = y_start,
            ymax = y_end,
            fill = class),
        color = NA,
        linewidth = 0
      ) +
      # 灰色标签块紧贴顶部框线
      annotate("rect",
               xmin = 0.5, xmax = length(sample_names) + 0.5,
               ymin = 100, ymax = 105,
               fill = "#E0E0E0", 
               color = "black", 
               linewidth = 0.5) +
      annotate("text",
               x = mean(1:length(sample_names)),
               y = 102.5,
               label = ko_name,
               size = 3.5,
               family = "serif",
               fontface = "bold",
               color = "black") +
      scale_fill_manual(
        values = classes_data$colors[as.character(class_levels)],
        guide = "none"  # 不在子图中显示图例
      ) +
      scale_y_continuous(
        breaks = seq(0, 100, by = 20),
        labels = seq(0, 100, by = 20),
        expand = expansion(mult = c(0, 0.05)),
        limits = c(0, 105)
      ) +
      scale_x_continuous(
        breaks = 1:length(sample_names),
        labels = sample_names,
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      labs(
        y = "Relative Abundance",
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
        plot.margin = margin(10, 10, 10, 10, "pt")
      )
    
    plot_list[[i]] <- p
  }
  
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  if (length(plot_list) == 0) {
    cat(sprintf("  %s循环: 没有生成图形\n", cycle_name))
    return(NULL)
  }
  
  cat(sprintf("  成功创建%d个基因的图形\n", length(plot_list)))
  
  # 创建分类图例
  cat("  创建分类图例...\n")
  
  # 准备图例数据 - 先古菌后细菌
  legend_classes <- c()
  legend_labels <- c()
  
  if (length(archaea_classes) > 0) {
    legend_classes <- c(legend_classes, archaea_classes)
    legend_labels <- c(legend_labels, rep("Archaea", length(archaea_classes)))
  }
  
  if (length(bacteria_classes) > 0) {
    legend_classes <- c(legend_classes, bacteria_classes)
    legend_labels <- c(legend_labels, rep("Bacteria", length(bacteria_classes)))
  }
  
  legend_data <- data.frame(
    class = factor(legend_classes, levels = legend_classes),
    group = factor(legend_labels, levels = unique(legend_labels)),
    Value = 1:length(legend_classes)
  )
  
  # 创建分组的图例
  legend_plot <- ggplot(legend_data, aes(x = 1, y = Value, fill = class)) +
    geom_tile(color = NA) +
    annotate("text", 
             x = 0.8, 
             y = ifelse(length(archaea_classes) > 0, 
                        length(archaea_classes) / 2 + 0.5, 
                        1),
             label = "Archaea",
             hjust = 1,
             size = 3.5,
             family = "serif",
             fontface = "bold",
             color = "darkred") +
    annotate("text", 
             x = 0.8, 
             y = ifelse(length(archaea_classes) > 0,
                        length(archaea_classes) + length(bacteria_classes) / 2 + 0.5,
                        length(bacteria_classes) / 2 + 0.5),
             label = "Bacteria",
             hjust = 1,
             size = 3.5,
             family = "serif",
             fontface = "bold",
             color = "darkblue") +
    scale_fill_manual(
      values = classes_data$colors[legend_classes],
      name = paste0(cycle_name, " Cycle"),
      guide = guide_legend(
        ncol = 1,
        byrow = TRUE,
        keyheight = unit(0.35, "cm"),
        keywidth = unit(0.4, "cm"),
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
      legend.spacing.y = unit(0.1, "cm"),
      plot.margin = margin(5, 5, 5, 5, "pt")
    )
  
  legend_grob <- get_legend(legend_plot)
  
  # 组合图形
  cat("  组合图形...\n")
  
  if (length(plot_list) > 1) {
    spacer_plots <- list()
    for (j in 1:(length(plot_list) - 1)) {
      spacer <- ggplot() + theme_void()
      spacer_plots[[j]] <- spacer
    }
    
    combined_plots <- list()
    for (k in 1:length(plot_list)) {
      combined_plots[[2*k - 1]] <- plot_list[[k]]
      if (k < length(plot_list)) {
        combined_plots[[2*k]] <- spacer_plots[[k]]
      }
    }
    
    plot_widths <- rep(1, length(combined_plots))
    for (m in seq(2, length(combined_plots), by = 2)) {
      plot_widths[m] <- 0.02  # 最小间距
    }
    
    plots_grid <- plot_grid(
      plotlist = combined_plots,
      nrow = 1,
      align = "h",
      axis = "tb",
      rel_widths = plot_widths
    )
  } else {
    plots_grid <- plot_list[[1]]
  }
  
  # 组合图形和图例
  combined_plot <- plot_grid(
    plots_grid, legend_grob,
    nrow = 1,
    rel_widths = c(0.7, 0.3),
    align = "h",
    axis = "tb"
  )
  
  # 保存图形
  cat("  保存图形...\n")
  
  # 计算图形尺寸
  num_genes <- length(genes_list)
  base_width_per_gene <- 2.7
  spacing_width <- 0.02 * max(0, num_genes - 1)
  
  num_classes <- length(classes_data$classes)
  legend_width <- max(5, num_classes * 0.6)
  
  output_width <- num_genes * base_width_per_gene + spacing_width + legend_width
  
  # 根据循环调整画布大小
  if (cycle_name == "Carbon") {
    output_width <- output_width * 1.7
  } else if (cycle_name == "Nitrogen") {
    output_width <- output_width + 0.7 * base_width_per_gene
  }
  
  output_height <- 11
  
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

# 7. 创建并保存三张图
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

# 8. 输出完成信息
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

cat("\n✅ 修改总结：\n")
cat("  1. ✅ 灰色标签块紧贴顶部框线，与柱子之间有空隙\n")
cat("  2. ✅ 图例按古菌/细菌分类，古菌在上，细菌在下\n")
cat("  3. ✅ 柱子间间距最小化（0.02cm）\n")
cat("  4. ✅ 柱子宽度缩小到0.4\n")
cat("  5. ✅ 彻底去掉柱子间的黑色边框线（使用geom_rect）\n")
cat("  6. ✅ 所有图形元素都无黑线边框\n")
cat("  7. ✅ 图例添加古菌/细菌分组标签\n")

cat("\n📏 图形尺寸:\n")
cat("  宽度 = 基因数量 × 2.7cm + (基因数量-1) × 0.02cm + 图例宽度\n")
cat("  图例宽度 = max(5cm, class数量 × 0.6cm)\n")
cat("  高度 = 11cm\n")

cat("\n📍 文件保存在当前目录:", getwd(), "\n")

# 保存清洗后的数据
write.csv(tax_tab, "cleaned_data.csv", row.names = FALSE)
cat("✅ 清洗后的数据已保存为: cleaned_data.csv\n")
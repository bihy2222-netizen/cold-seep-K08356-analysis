# 加载包
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(stringr)
library(cowplot)
library(gridExtra)

# 设置工作目录
setwd("/Users/catherine/Downloads/重新画图/宿主图/")

# 读取CSV数据
tax_tab <- read.csv("huizong.csv", check.names = FALSE, stringsAsFactors = FALSE, 
                    fileEncoding = "UTF-8", na.strings = c("", "NA"), fill = TRUE)

# 查看数据结构
cat("数据维度:", dim(tax_tab), "\n")
cat("列名:", colnames(tax_tab), "\n")

# 标准化列名
colnames(tax_tab)[1] <- "KO"
colnames(tax_tab)[2] <- "Order"

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

# 定义样本顺序
sample_names <- c("IS", "AS", "ES", "NS")

# 使用Set3色板，适合分类数据
nature_colors <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
  "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
  "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
  "#66A61E", "#E6AB02", "#A6761D", "#666666"
)

# 扩展配色方案
get_extended_colors <- function(order_list, base_colors) {
  n <- length(order_list)
  
  if (n <= length(base_colors)) {
    colors <- base_colors[1:n]
  } else {
    # 使用viridis色板生成更多颜色
    library(viridis)
    colors <- viridis(n, option = "D", begin = 0.1, end = 0.9)
  }
  
  names(colors) <- order_list
  return(colors)
}

# 收集所有Order用于统一配色
cat("\n收集所有Order...\n")
all_orders <- c()
for (ko_name in available_genes) {
  ko_data <- tax_tab[tax_tab$KO == ko_name, ]
  ko_orders <- unique(ko_data$Order)
  ko_orders <- ko_orders[!is.na(ko_orders) & ko_orders != ""]
  all_orders <- unique(c(all_orders, ko_orders))
}

cat("总共发现", length(all_orders), "个不同的Order\n")

# 生成统一的颜色
order_colors <- get_extended_colors(all_orders, nature_colors)

# 修改绘图函数
create_ko_plot <- function(ko_name, tax_data, order_colors, sample_names, show_x_labels = TRUE) {
  # 筛选特定功能基因的数据
  ko_data <- tax_data[tax_tab$KO == ko_name, ]
  
  if (nrow(ko_data) == 0) {
    return(NULL)
  }
  
  # 获取该KO的所有Order
  ko_orders <- unique(ko_data$Order)
  ko_orders <- ko_orders[!is.na(ko_orders) & ko_orders != ""]
  
  if (length(ko_orders) == 0) {
    return(NULL)
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
  
  # 按总丰度排序Order（按所有样本的总和）
  order_summary <- plot_data %>%
    group_by(Order) %>%
    summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
    arrange(desc(TotalAbundance))
  
  order_levels <- order_summary$Order
  plot_data$Order <- factor(plot_data$Order, levels = order_levels)
  
  # 只保留当前KO有的Order颜色
  current_colors <- order_colors[as.character(order_levels)]
  
  # 创建简洁的柱状图
  p <- ggplot(plot_data, aes(x = Sample, y = Percentage, fill = Order)) +
    geom_col(
      width = 0.35,  # 柱子更细（原来的一半）
      color = "black",  # 柱子边框
      linewidth = 0.1,
      position = position_stack(reverse = FALSE)
    ) +
    # 添加灰色标题方块（覆盖整个X轴范围）
    annotate("rect",
             xmin = 0.5, xmax = length(sample_names) + 0.5,
             ymin = 100, ymax = 105,
             fill = "#E0E0E0", color = "black", linewidth = 0.5) +
    annotate("text",
             x = mean(1:length(sample_names)),
             y = 102.5,
             label = ko_name,
             size = 3.2,
             fontface = "bold",
             color = "black") +
    scale_fill_manual(
      values = current_colors,
      guide = "none"  # 不显示图例，使用统一的图例
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, by = 25),
      labels = paste0(seq(0, 100, by = 25), "%"),
      expand = expansion(mult = c(0, 0.05)),  # 扩展Y轴上限以容纳标题块
      limits = c(0, 105)  # 增加上限以容纳标题块
    ) +
    scale_x_discrete(
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    labs(
      y = NULL,
      x = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      # 坐标轴
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5,
        size = 9,
        color = "black",
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 8,
        color = "black"
      ),
      axis.title = element_blank(),
      
      # 网格线
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # 边框
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # 坐标轴线
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      
      # 边距
      plot.margin = margin(2, 2, 2, 2, "pt")
    )
  
  return(p)
}

# 为每个功能基因创建图形
cat("\n为每个功能基因创建图形...\n")

# C循环图形
c_plot_list <- list()
for (i in 1:length(c_genes_available)) {
  ko_name <- c_genes_available[i]
  cat(sprintf("  处理C循环: %s\n", ko_name))
  plot <- create_ko_plot(ko_name, tax_tab, order_colors, sample_names, TRUE)
  if (!is.null(plot)) {
    c_plot_list[[i]] <- plot
  }
}

# N循环图形
n_plot_list <- list()
for (i in 1:length(n_genes_available)) {
  ko_name <- n_genes_available[i]
  cat(sprintf("  处理N循环: %s\n", ko_name))
  plot <- create_ko_plot(ko_name, tax_tab, order_colors, sample_names, TRUE)
  if (!is.null(plot)) {
    n_plot_list[[i]] <- plot
  }
}

# S循环图形
s_plot_list <- list()
for (i in 1:length(s_genes_available)) {
  ko_name <- s_genes_available[i]
  cat(sprintf("  处理S循环: %s\n", ko_name))
  plot <- create_ko_plot(ko_name, tax_tab, order_colors, sample_names, TRUE)
  if (!is.null(plot)) {
    s_plot_list[[i]] <- plot
  }
}

cat(sprintf("\n成功生成 %d 个C循环图形\n", length(c_plot_list)))
cat(sprintf("成功生成 %d 个N循环图形\n", length(n_plot_list)))
cat(sprintf("成功生成 %d 个S循环图形\n", length(s_plot_list)))

# 创建图例函数
create_cycle_legend <- function(plot_list, genes_available, cycle_name, order_colors) {
  if (length(plot_list) == 0) return(NULL)
  
  # 收集该循环中出现的所有Order
  cycle_orders <- c()
  for (ko_name in genes_available) {
    ko_data <- tax_tab[tax_tab$KO == ko_name, ]
    ko_orders <- unique(ko_data$Order)
    ko_orders <- ko_orders[!is.na(ko_orders) & ko_orders != ""]
    cycle_orders <- unique(c(cycle_orders, ko_orders))
  }
  
  if (length(cycle_orders) == 0) return(NULL)
  
  # 创建图例数据
  legend_data <- data.frame(
    Order = factor(cycle_orders, levels = cycle_orders),
    Value = 1:length(cycle_orders)
  )
  
  # 计算图例列数
  legend_cols <- if(length(cycle_orders) > 25) {
    4
  } else if(length(cycle_orders) > 15) {
    3
  } else if(length(cycle_orders) > 8) {
    2
  } else {
    1
  }
  
  # 创建图例
  p <- ggplot(legend_data, aes(x = 1, y = Value, fill = Order)) +
    geom_tile() +
    scale_fill_manual(
      values = order_colors[cycle_orders],
      name = paste0(cycle_name, " Cycle"),
      guide = guide_legend(
        ncol = legend_cols,
        byrow = TRUE,
        keyheight = unit(0.3, "cm"),
        keywidth = unit(0.4, "cm"),
        title.position = "top",
        title.hjust = 0.5,
        title = element_text(size = 9, face = "bold", margin = margin(b = 3)),
        label.position = "right",
        label.hjust = 0,
        label.theme = element_text(size = 6.5),
        title.theme = element_text(size = 9, face = "bold")
      )
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "left",
      plot.margin = margin(0, 5, 0, 5, "pt")
    )
  
  return(p)
}

# 创建三种循环的图例
cat("\n创建图例...\n")
c_legend_plot <- create_cycle_legend(c_plot_list, c_genes_available, "C", order_colors)
n_legend_plot <- create_cycle_legend(n_plot_list, n_genes_available, "N", order_colors)
s_legend_plot <- create_cycle_legend(s_plot_list, s_genes_available, "S", order_colors)

# 组合每行的图形和图例
cat("\n组合图形和图例...\n")

# 组合C循环行
c_plots_grid <- if(length(c_plot_list) > 0) {
  plot_grid(
    plotlist = c_plot_list,
    nrow = 1,
    align = "h",
    axis = "tb",
    rel_widths = rep(1, length(c_plot_list))
  )
} else NULL

# 组合N循环行
n_plots_grid <- if(length(n_plot_list) > 0) {
  plot_grid(
    plotlist = n_plot_list,
    nrow = 1,
    align = "h",
    axis = "tb",
    rel_widths = rep(1, length(n_plot_list))
  )
} else NULL

# 组合S循环行
s_plots_grid <- if(length(s_plot_list) > 0) {
  plot_grid(
    plotlist = s_plot_list,
    nrow = 1,
    align = "h",
    axis = "tb",
    rel_widths = rep(1, length(s_plot_list))
  )
} else NULL

# 为每行添加Y轴标签
add_y_label <- function(plot_grid, label) {
  if (is.null(plot_grid)) return(NULL)
  
  label_plot <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = label, angle = 90, 
             size = 4, fontface = "bold") +
    theme_void() +
    theme(
      plot.margin = margin(0, 5, 0, 0, "pt")
    )
  
  combined <- plot_grid(
    label_plot, plot_grid,
    nrow = 1,
    rel_widths = c(0.03, 0.97),  # Y轴标签占3%
    align = "h",
    axis = "tb"
  )
  
  return(combined)
}

c_row_with_label <- add_y_label(c_plots_grid, "C cycle")
n_row_with_label <- add_y_label(n_plots_grid, "N cycle")
s_row_with_label <- add_y_label(s_plots_grid, "S cycle")

# 组合每行和对应的图例
combine_row_with_legend <- function(row_plot, legend_plot, row_height = 1) {
  if (is.null(row_plot) || is.null(legend_plot)) return(row_plot)
  
  # 提取图例
  legend_grob <- get_legend(legend_plot)
  
  # 组合图形和图例
  combined <- plot_grid(
    row_plot, legend_grob,
    nrow = 1,
    rel_widths = c(0.7, 0.3),  # 图形占70%，图例占30%
    align = "h",
    axis = "tb"
  )
  
  return(combined)
}

# 组合所有行
c_row_combined <- combine_row_with_legend(c_row_with_label, c_legend_plot)
n_row_combined <- combine_row_with_legend(n_row_with_label, n_legend_plot)
s_row_combined <- combine_row_with_legend(s_row_with_label, s_legend_plot)

# 收集所有行
all_rows <- list()
if (!is.null(c_row_combined)) all_rows[[length(all_rows) + 1]] <- c_row_combined
if (!is.null(n_row_combined)) all_rows[[length(all_rows) + 1]] <- n_row_combined
if (!is.null(s_row_combined)) all_rows[[length(all_rows) + 1]] <- s_row_combined

# 垂直排列所有行
if (length(all_rows) > 0) {
  # 计算相对高度
  rel_heights <- c()
  if (!is.null(c_row_combined)) rel_heights <- c(rel_heights, 1.2)  # C循环行，高度增加
  if (!is.null(n_row_combined)) rel_heights <- c(rel_heights, 1.2)  # N循环行，高度增加
  if (!is.null(s_row_combined)) rel_heights <- c(rel_heights, 1.2)  # S循环行，高度增加
  
  final_plot <- plot_grid(
    plotlist = all_rows,
    ncol = 1,
    rel_heights = rel_heights,
    align = "v",
    axis = "lr"
  )
  
  # 保存最终图形
  cat("\n保存最终图形...\n")
  
  # 计算输出尺寸（增加高度和宽度）
  total_genes <- length(c_genes_available) + length(n_genes_available) + length(s_genes_available)
  output_width <- 15 + total_genes * 3  # 拉宽画布
  output_height <- 15 * 3  # 高度增加一倍（原来是20）
  
  # 高质量PDF输出
  ggsave(
    "Multi_Cycle_Right_Legends.pdf",
    final_plot,
    width = output_width,
    height = output_height,
    units = "cm",
    dpi = 600,
    device = cairo_pdf,
    limitsize = FALSE
  )
  
  # PNG输出
  ggsave(
    "Multi_Cycle_Right_Legends.png",
    final_plot,
    width = output_width,
    height = output_height,
    units = "cm",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  
  # 输出完成信息
  cat("\n🎉 任务完成！\n")
  cat("📁 输出文件:\n")
  cat("  - Multi_Cycle_Right_Legends.pdf\n")
  cat("  - Multi_Cycle_Right_Legends.png\n")
  cat("\n✅ 图形特点：\n")
  cat("  1. 三行布局：C循环、N循环、S循环各一行\n")
  cat("  2. 每个子图都有外边框\n")
  cat("  3. 灰色标题块在柱子顶部\n")
  cat("  4. 所有X轴都标注IS/AS/ES/NS\n")
  cat("  5. 柱子宽度减少一半\n")
  cat("  6. 图形高度增加一倍\n")
  cat("  7. 图例在右侧，与对应行对齐\n")
  cat("  8. 使用统一的Set3色板配色\n")
  cat("\n📊 统计信息：\n")
  cat("  C循环基因: ", length(c_genes_available), "个\n")
  cat("  N循环基因: ", length(n_genes_available), "个\n")
  cat("  S循环基因: ", length(s_genes_available), "个\n")
  cat("  总Order数量: ", length(all_orders), "\n")
  cat("  输出尺寸: ", output_width, "cm × ", output_height, "cm\n")
  
} else {
  cat("\n❌ 错误：没有成功生成任何图形！\n")
}
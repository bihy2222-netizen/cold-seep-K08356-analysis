# 加载包
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(stringr)
library(cowplot)
library(gridExtra)

# 设置工作目录
setwd("/Users/catherine/Downloads/画功能基因所在宿主的图/")

# 读取CSV数据
tax_tab <- read.csv("C-order.csv", check.names = FALSE, stringsAsFactors = FALSE, 
                    fileEncoding = "UTF-8", na.strings = c("", "NA"), fill = TRUE)

# 找出有效列（非空列名）
valid_cols <- colnames(tax_tab)
valid_cols <- valid_cols[!is.na(valid_cols) & valid_cols != ""]

# 只保留有效列
tax_tab <- tax_tab[, valid_cols]

# 根据数据结构，第一列是KO，第二列是Order
colnames(tax_tab)[1] <- "KO"
colnames(tax_tab)[2] <- "Order"

# 查看唯一的KO和Order
ko_list <- unique(tax_tab$KO)
order_list <- unique(tax_tab$Order)
order_list <- order_list[!is.na(order_list) & order_list != ""]

# 提取样本列（从第3列开始）
sample_cols <- colnames(tax_tab)[3:ncol(tax_tab)]

# 定义样本顺序
sample_order <- c(
  "SY365BB-0-4", "SY365BB-4-8", "SY365BB-8-12",
  "SY366YB-0-4", "SY366YB-4-8", "SY366YB-8-12",
  "SY366YW-0-4", "SY366YW-4-8", "SY366YW-8-12",
  "SY368YW-0-4", "SY368YW-4-8", "SY368YW-8-12",
  "SY456YB-0-4", "SY456YB-4-8", "SY456YB-8-12",
  "SY457BB-0-4", "SY457BB-4-8", "SY457BB-8-12",
  "SY459WG-0-4", "SY459WG-4-8", "SY459WG-8-12",
  "S4_12-15", "S4_9-12",
  "SQ_58_0-4", "SQ_58_4-8", "SQ_58_8-12",
  "SQ_81_0-4", "SQ_81_4-8", "SQ_81_8-12",
  "S2_0-3", "S2_3-6", "S2_12-15",
  "S1_0-3", "S1_6-9", "S1_9-12",
  "C2_0-6", "C2_6-12", "C2_12-18",
  "C1_0-6", "C1_6-12", "C1_12-18",
  "C3_0-6", "C3_6-12", "C3_12-18",
  "S13_0-2", "S14_4-6", "S15_8-10",
  "ES_2_0-6",
  "S3_0-3", "S3_6-9", "S3_9-12",
  "R2111_N300_0-10", "R2111_N500_0-10", 
  "R2111_S300_0-10", "R2111_S500_0-10",
  "NS_0-6"
)

# 只保留在样本顺序中且在数据中的样本
sample_names <- intersect(sample_order, sample_cols)

# 定义分组
group_info <- data.frame(
  Sample = sample_names,
  Group = NA_character_,
  stringsAsFactors = FALSE
)

# 分配分组
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
                "SQ_58_0-4", "SQ_58_4-8", "SQ_58_8-12",
                "SQ_81_0-4", "SQ_81_4-8", "SQ_81_8-12")

ES_samples <- c("C1_0-6", "C1_6-12", "C1_12-18",
                "C2_0-6", "C2_6-12", "C2_12-18",
                "C3_0-6", "C3_6-12", "C3_12-18",
                "ES_2_0-6",
                "S13_0-2", "S14_4-6", "S15_8-10")

NS_samples <- c("S3_0-3", "S3_6-9", "S3_9-12",
                "NS_0-6",
                "R2111_N300_0-10", "R2111_N500_0-10", 
                "R2111_S300_0-10", "R2111_S500_0-10")

group_info$Group[group_info$Sample %in% IS_samples] <- "IS"
group_info$Group[group_info$Sample %in% AS_samples] <- "AS"
group_info$Group[group_info$Sample %in% ES_samples] <- "ES"
group_info$Group[group_info$Sample %in% NS_samples] <- "NS"

# 设置因子顺序
group_info$Group <- factor(group_info$Group, levels = c("IS", "AS", "ES", "NS"))
group_info$Sample <- factor(group_info$Sample, levels = sample_names)

# 使用您指定的分组颜色
group_colors <- c(
  "IS" = "#00FF7F",   # 春绿色
  "AS" = "#9E39E6",   # 紫色
  "ES" = "#FFA500",   # 橙色
  "NS" = "#00FFFF"    # 青色
)

# 使用您提供的配色方案
custom_colors <- c(
  "#8CD0C3",  # 青色
  "#FAF5B5",  # 淡黄色
  "#BCB9D8",  # 淡紫色
  "#F18072",  # 珊瑚红
  "#80B1D2",  # 淡蓝色
  "#F9B063",  # 橙色
  "#B3D46B",  # 黄绿色
  "#F7CBDF",  # 淡粉色
  "#D7D7D5",  # 灰色
  "#BA7FB5"   # 紫色
)

# 扩展配色方案
get_extended_colors <- function(order_list, base_colors) {
  n <- length(order_list)
  
  if (n <= length(base_colors)) {
    colors <- base_colors[1:n]
  } else {
    # 使用色轮生成更多颜色
    hues <- seq(15, 375, length = n + 1)[1:n]
    colors <- hcl(h = hues, c = 100, l = 70)
  }
  
  names(colors) <- order_list
  return(colors)
}

# 获取所有Order并生成颜色
all_orders <- unique(tax_tab$Order)
all_orders <- all_orders[!is.na(all_orders) & all_orders != ""]
order_colors <- get_extended_colors(all_orders, custom_colors)

# 修改绘图函数
create_ko_plot <- function(ko_name, tax_data, group_info, order_colors) {
  cat(sprintf("  处理功能基因: %s\n", ko_name))
  
  # 筛选特定功能基因的数据
  ko_data <- tax_data[tax_data$KO == ko_name, ]
  
  if (nrow(ko_data) == 0) {
    warning(sprintf("没有找到功能基因 %s 的数据", ko_name))
    return(NULL)
  }
  
  # 获取该KO的所有Order
  ko_orders <- unique(ko_data$Order)
  ko_orders <- ko_orders[!is.na(ko_orders) & ko_orders != ""]
  
  if (length(ko_orders) == 0) {
    warning(sprintf("功能基因 %s 没有有效的Order数据", ko_name))
    return(NULL)
  }
  
  # 手动转换为长格式
  plot_data_list <- list()
  for (i in 1:length(sample_names)) {
    sample_name <- sample_names[i]
    if (sample_name %in% colnames(ko_data)) {
      temp_data <- data.frame(
        KO = ko_data$KO,
        Order = ko_data$Order,
        Sample = sample_name,
        Abundance = as.numeric(ko_data[[sample_name]]),
        stringsAsFactors = FALSE
      )
      plot_data_list[[i]] <- temp_data
    }
  }
  
  plot_data <- do.call(rbind, plot_data_list)
  
  # 将NA值替换为0
  plot_data$Abundance[is.na(plot_data$Abundance)] <- 0
  
  # 添加分组信息
  plot_data$Sample <- factor(plot_data$Sample, levels = sample_names)
  plot_data <- merge(plot_data, group_info, by = "Sample")
  
  # 按总丰度排序Order
  order_summary <- aggregate(Abundance ~ Order, data = plot_data, sum)
  order_summary <- order_summary[order(-order_summary$Abundance), ]
  order_levels <- order_summary$Order
  
  plot_data$Order <- factor(plot_data$Order, levels = order_levels)
  
  # 计算每个样本的总丰度用于标准化
  sample_totals <- aggregate(Abundance ~ Sample, data = plot_data, sum)
  colnames(sample_totals)[2] <- "Total"
  
  plot_data <- merge(plot_data, sample_totals, by = "Sample")
  plot_data$Relative <- ifelse(plot_data$Total > 0, 
                               plot_data$Abundance / plot_data$Total * 100, 0)
  
  # 获取该KO的Order颜色
  plot_order_colors <- order_colors[as.character(order_levels)]
  
  # 计算分组位置
  unique_samples <- unique(plot_data[, c("Sample", "Group")])
  unique_samples <- unique_samples[order(unique_samples$Sample), ]
  
  group_positions <- data.frame()
  for (grp in levels(unique_samples$Group)) {
    group_samples <- unique_samples$Sample[unique_samples$Group == grp]
    if (length(group_samples) > 0) {
      start_idx <- which(sample_names == group_samples[1])[1] - 0.5
      end_idx <- which(sample_names == group_samples[length(group_samples)])[1] + 0.5
      label_x <- (which(sample_names == group_samples[1])[1] + 
                    which(sample_names == group_samples[length(group_samples)])[1]) / 2
      
      group_positions <- rbind(group_positions, data.frame(
        Group = grp,
        start = start_idx,
        end = end_idx,
        label_x = label_x,
        n_samples = length(group_samples)
      ))
    }
  }
  
  # 创建柱状图
  p <- ggplot() +
    # 柱状图
    geom_col(
      data = plot_data,
      aes(x = Sample, y = Relative, fill = Order),
      width = 0.8,
      color = "white",
      linewidth = 0.1,
      na.rm = TRUE
    ) +
    # 顶部分组色块 - 拉宽一些
    geom_rect(
      data = group_positions,
      aes(xmin = start + 0.1, xmax = end - 0.1, ymin = 100, ymax = 103),
      fill = group_colors[as.character(group_positions$Group)],
      alpha = 0.8
    ) +
    # 分组标签
    geom_text(
      data = group_positions,
      aes(x = label_x, y = 101.5, label = Group),
      size = 3.5,
      fontface = "bold",
      color = "black",
      vjust = 0.5
    ) +
    # 使用自定义配色
    scale_fill_manual(
      values = plot_order_colors,
      name = "Host Order",
      breaks = order_levels,
      labels = order_levels,
      guide = guide_legend(
        ncol = 1,
        keyheight = unit(0.5, "cm"),
        keywidth = unit(0.8, "cm"),
        title.position = "top",
        title.hjust = 0.5,
        title = element_text(face = "bold", size = 10),
        label.position = "right",
        label.hjust = 0,
        label.vjust = 0.5
      )
    ) +
    # Y轴设置 - 限制在0-103%
    scale_y_continuous(
      breaks = c(0, 20, 40, 60, 80, 100),
      labels = c("0%", "20%", "40%", "60%", "80%", "100%"),
      expand = expansion(mult = c(0, 0.03)),  # 稍微扩展以显示分组标签
      limits = c(0, 103)  # 限制在103%
    ) +
    # 主题设置
    theme_minimal(base_size = 12) +
    theme(
      # 坐标轴文本
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 8,
        color = "black",
        margin = margin(t = 3)
      ),
      axis.text.y = element_text(
        size = 10,
        color = "black"
      ),
      
      # 坐标轴标题
      axis.title.x = element_blank(),
      axis.title.y = element_text(
        size = 12, 
        margin = margin(r = 10),
        face = "bold",
        color = "black"
      ),
      
      # 标题 - 移除默认标题
      plot.title = element_blank(),
      
      # 图例
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.text = element_text(
        size = 9,
        color = "black",
        margin = margin(r = 5)
      ),
      legend.title = element_text(
        size = 11, 
        face = "bold",
        color = "black",
        margin = margin(b = 5)
      ),
      legend.key = element_rect(
        color = NA,
        fill = NA
      ),
      legend.key.size = unit(0.4, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.8, "cm"),
      legend.spacing.y = unit(0.2, "cm"),
      legend.box.margin = margin(l = 5),
      legend.background = element_rect(
        fill = "white",
        color = "gray80",
        linewidth = 0.5
      ),
      
      # 网格线
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(
        color = "gray90", 
        linewidth = 0.3,
        linetype = "solid"
      ),
      panel.grid.minor.y = element_blank(),
      
      # 边距和背景
      plot.margin = margin(0.5, 1.5, 0.5, 1, "cm"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      
      # 坐标轴线
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm")
    ) +
    labs(
      y = "Relative Abundance (%)"
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}

# 为每个功能基因创建图形
cat("\n开始为每个功能基因绘图:\n")
plot_list <- list()

for (i in 1:length(ko_list)) {
  cat(sprintf("\n正在绘制第 %d/%d 个功能基因: %s\n", i, length(ko_list), ko_list[i]))
  plot <- create_ko_plot(ko_list[i], tax_tab, group_info, order_colors)
  if (!is.null(plot)) {
    plot_list[[i]] <- plot
  }
}

# 移除NULL值
plot_list <- plot_list[!sapply(plot_list, is.null)]

if (length(plot_list) == 0) {
  stop("没有生成任何图形!")
}

cat(sprintf("\n成功生成 %d 个图形\n", length(plot_list)))

# 创建灰色标题方块（覆盖在分组标签上方）
create_title_box <- function(label, ko_name, group_positions) {
  # 计算灰色方块的宽度（与分组色块一致）
  # 找到第一个和最后一个分组的位置
  if (nrow(group_positions) > 0) {
    min_start <- min(group_positions$start)
    max_end <- max(group_positions$end)
    
    # 创建灰色背景矩形，覆盖整个柱状图区域（包括分组标签）
    bg_grob <- rectGrob(
      x = (min_start + 0.1 + max_end - 0.1) / 2 / length(sample_names),  # 居中
      y = 0.5,
      width = (max_end - min_start - 0.2) / length(sample_names),  # 与分组色块宽度一致
      height = 0.25,  # 高度
      gp = gpar(fill = "#D7D7D5", col = NA, alpha = 0.7)
    )
  } else {
    # 如果没有分组数据，使用默认宽度
    bg_grob <- rectGrob(
      x = 0.5, y = 0.5,
      width = 0.9, height = 0.25,
      gp = gpar(fill = "#D7D7D5", col = NA, alpha = 0.7)
    )
  }
  
  # 创建标签文本
  label_grob <- textGrob(
    label = paste0(label, " ", ko_name),
    x = 0.5, y = 0.5,
    just = "center",
    gp = gpar(fontface = "bold", fontsize = 13, col = "black")
  )
  
  # 组合所有元素
  title_grob <- gTree(children = gList(bg_grob, label_grob))
  return(title_grob)
}

# 获取分组位置信息（用于计算灰色方块宽度）
# 这里我们使用第一个功能基因的分组位置作为参考
ref_plot <- plot_list[[1]]
ref_data <- tax_tab[tax_tab$KO == ko_list[1], ]

# 计算参考分组位置
ref_unique_samples <- unique(group_info[, c("Sample", "Group")])
ref_unique_samples <- ref_unique_samples[order(ref_unique_samples$Sample), ]

ref_group_positions <- data.frame()
for (grp in levels(ref_unique_samples$Group)) {
  group_samples <- ref_unique_samples$Sample[ref_unique_samples$Group == grp]
  if (length(group_samples) > 0) {
    start_idx <- which(sample_names == group_samples[1])[1] - 0.5
    end_idx <- which(sample_names == group_samples[length(group_samples)])[1] + 0.5
    label_x <- (which(sample_names == group_samples[1])[1] + 
                  which(sample_names == group_samples[length(group_samples)])[1]) / 2
    
    ref_group_positions <- rbind(ref_group_positions, data.frame(
      Group = grp,
      start = start_idx,
      end = end_idx,
      label_x = label_x,
      n_samples = length(group_samples)
    ))
  }
}

# 为每个图创建标题方块
title_grobs <- list()
for (i in 1:length(plot_list)) {
  title_grobs[[i]] <- create_title_box(LETTERS[i], ko_list[i], ref_group_positions)
}

# 创建组合图
n_plots <- length(plot_list)
plot_with_titles <- list()

for (i in 1:n_plots) {
  # 创建每个图的组合：标题 + 图表
  title_plot <- ggplot() + 
    theme_void() +
    annotation_custom(title_grobs[[i]], xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    coord_cartesian(clip = "off")
  
  # 组合标题和图表
  combined <- plot_grid(
    title_plot,
    plot_list[[i]],
    ncol = 1,
    rel_heights = c(0.08, 0.92),  # 标题占8%
    align = "v",
    axis = "lr"
  )
  
  plot_with_titles[[i]] <- combined
}

# 垂直排列所有带标题的图
final_combined_plot <- plot_grid(
  plotlist = plot_with_titles,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = rep(1, n_plots)
)

# 保存组合图形
output_width <- 50
output_height <- 11 * n_plots

cat(sprintf("\n保存组合图形 (宽: %dcm, 高: %dcm)\n", output_width, output_height))

# 高质量PDF输出
ggsave(
  "KO_Order_barplots_Final_v2.pdf",
  final_combined_plot,
  width = output_width,
  height = output_height,
  units = "cm",
  dpi = 600,
  device = cairo_pdf
)

# PNG输出
ggsave(
  "KO_Order_barplots_Final_v2.png",
  final_combined_plot,
  width = output_width,
  height = output_height,
  units = "cm",
  dpi = 600,
  bg = "white"
)

# 也单独保存每个图（带标题）
for (i in 1:n_plots) {
  # 保存PDF
  ggsave(
    sprintf("KO_%s_Order_barplot_v2.pdf", ko_list[i]),
    plot_with_titles[[i]],
    width = 40,
    height = 14,
    units = "cm",
    dpi = 600,
    device = cairo_pdf
  )
  
  # 保存PNG
  ggsave(
    sprintf("KO_%s_Order_barplot_v2.png", ko_list[i]),
    plot_with_titles[[i]],
    width = 40,
    height = 14,
    units = "cm",
    dpi = 600,
    bg = "white"
  )
}

# 输出完成信息
cat("\n🎉 所有任务完成！\n")
cat("📁 输出文件:\n")
cat("  - KO_Order_barplots_Final_v2.pdf/png (组合图)\n")
cat("  - KO_<功能基因>_Order_barplot_v2.pdf/png (单图)\n")
cat("\n✅ 修改完成：\n")
cat("  1. 灰色方块覆盖在IS/AS/ES/NS分组标签上方\n")
cat("  2. Y轴高度限制在0-103%\n")
cat("  3. 分组色块拉宽，与灰色方块宽度一致\n")
cat("  4. 分组色块范围: 100-103%\n")
cat("  5. 分组标签位置: 101.5%\n")
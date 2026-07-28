# 加载包
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)

# 设置工作目录
setwd("/Users/catherine/Downloads/metaphlan 作图/class 级别作图-柱状图/")

# 读取CSV数据
tax_tab <- read.csv("class34.csv", check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# 定义样本顺序（严格按照要求的顺序）
sample_order <- c(
  "SY365BB-0-4", "SY365BB-4-8", "SY365BB-8-12",
  "SY366YB-0-4", "SY366YB-4-8", "SY366YB-8-12",
  "SY366YW-0-4", "SY366YW-4-8", "SY366YW-8-12",
  "SY368YW-0-4", "SY368YW-4-8", "SY368YW-8-12",
  "SY456YB-0-4", "SY456YB-4-8", "SY456YB-8-12",
  "SY457BB-0-4", "SY457BB-4-8", "SY457BB-8-12",
  "SY459WG-0-4", "SY459WG-4-8", "SY459WG-8-12",
  "S4_12-15", "S4_9-12",
  "SQ_58_0-4", "SQ_58_4-8", "SQ_58_-8-12",
  "SQ_81_0-4", "SQ_81_4-8", "SQ_81_-8-12",
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

# 定义分组（根据您提供的信息修正）
group_info <- data.frame(
  Sample = sample_order,
  Group = NA_character_,
  stringsAsFactors = FALSE
)

# 分配分组（修正版）
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

# 应用分组
group_info$Group[group_info$Sample %in% IS_samples] <- "IS"
group_info$Group[group_info$Sample %in% AS_samples] <- "AS"
group_info$Group[group_info$Sample %in% ES_samples] <- "ES"
group_info$Group[group_info$Sample %in% NS_samples] <- "NS"

# 检查所有样本都有分组
if (any(is.na(group_info$Group))) {
  missing_samples <- group_info$Sample[is.na(group_info$Group)]
  cat("以下样本没有被分配到任何组:\n")
  cat(paste(missing_samples, collapse = ", "), "\n")
} else {
  cat("所有样本都已成功分配到组\n")
}

cat("\n分组统计:\n")
print(table(group_info$Group))

# 基本数据信息
taxon_names <- tax_tab[, 1]

# 检查sample_order中的样本是否都在数据中
missing_in_data <- setdiff(sample_order, colnames(tax_tab)[-1])
if (length(missing_in_data) > 0) {
  cat("\n警告: 以下样本在sample_order中但不在数据中:\n")
  cat(paste(missing_in_data, collapse = ", "), "\n")
}

# 只保留在sample_order中且在数据中的样本
sample_names <- intersect(sample_order, colnames(tax_tab)[-1])

cat("\n分类群数量:", length(taxon_names), "\n")
cat("样本数量:", length(sample_names), "\n")

# 创建绘图数据
plot_list <- list()
for(i in 1:length(sample_names)) {
  sample_data <- data.frame(
    Sample = sample_names[i],
    Taxon = taxon_names,
    Abundance = as.numeric(tax_tab[[sample_names[i]]])
  )
  plot_list[[i]] <- sample_data
}

plot_data <- do.call(rbind, plot_list)

# 转换为百分比
plot_data <- plot_data %>%
  group_by(Sample) %>%
  mutate(Abundance = Abundance / sum(Abundance) * 100) %>%
  ungroup()

# 添加分组信息
plot_data <- plot_data %>%
  left_join(group_info, by = "Sample")

# 设置因子顺序
plot_data$Taxon <- factor(plot_data$Taxon, levels = taxon_names)
plot_data$Sample <- factor(plot_data$Sample, levels = sample_names)
plot_data$Group <- factor(plot_data$Group, levels = c("IS", "AS", "ES", "NS"))

# 颜色设置 - 34种颜色
colors_34 <- c(
  "#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#911EB4", "#46F0F0", "#F032E6",
  "#BCF60C", "#FABEBE", "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000", "#AAFFC3",
  "#808000", "#FFD8B1", "#000075", "#808080", "#A9A9A9", "#2F4F4F", "#FF69B4", "#BA55D3",
  "#9370DB", "#3CB371", "#7B68EE", "#00FA9A", "#48D1CC", "#C71585", "#191970", "#FF4500",
  "#32CD32", "#8A2BE2"
)

# 计算分组色块的位置（窄条，在顶部）
group_positions <- plot_data %>%
  distinct(Sample, Group) %>%
  arrange(Sample) %>%
  group_by(Group) %>%
  summarize(
    start = which(sample_names == first(Sample))[1] - 0.5,
    end = which(sample_names == last(Sample))[1] + 0.5,
    label_x = (which(sample_names == first(Sample))[1] + which(sample_names == last(Sample))[1]) / 2,
    n_samples = n()
  ) %>%
  ungroup()

# 马卡龙色系（更窄的色块）
group_colors <- c("IS" = "#A8E6CF",    # 马卡龙淡绿色
                  "AS" = "#FFAAA5",    # 马卡龙红色
                  "ES" = "#FFD3B6",    # 马卡龙橙色
                  "NS" = "#AEC6CF")    # 马卡龙蓝色

# 为顶部窄色块准备数据（高度与字号一致）
# 假设字号为6pt，转换为数据单位（100%对应的高度）
label_height <- 2  # 窄条高度，约等于字号高度

# 顶部窄色块数据
top_band_data <- data.frame(
  xmin = group_positions$start,
  xmax = group_positions$end,
  ymin = 100,  # 从100%开始
  ymax = 100 + label_height,  # 窄条高度
  fill_color = group_colors[group_positions$Group]
)

# 创建图形
p <- ggplot() +
  # 先添加柱状图
  geom_col(
    data = plot_data,
    aes(x = Sample, y = Abundance, fill = Taxon),
    width = 0.8,
    color = "white",
    linewidth = 0.1
  ) +
  # 添加顶部窄色块
  geom_rect(
    data = top_band_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = top_band_data$fill_color,
    alpha = 0.8
  ) +
  # 在色块上添加组名标签
  geom_text(
    data = group_positions,
    aes(x = label_x, y = 100 + label_height/2, label = Group),
    size = 5,
    fontface = "bold",
    color = "black",
    vjust = 0.5
  ) +
  # 分类群颜色
  scale_fill_manual(
    values = colors_34,
    name = "Taxonomic Groups",
    breaks = taxon_names,
    labels = taxon_names,
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keyheight = unit(0.5, "cm"),
      keywidth = unit(0.8, "cm"),
      title.position = "top",
      title.hjust = 0.5,
      label.position = "right"
    )
  ) +
  # Y轴设置 - 只显示指定的刻度
  scale_y_continuous(
    breaks = c(0, 20, 40, 60, 80, 100),
    labels = c("0%", "20%", "40%", "60%", "80%", "100%"),
    expand = expansion(mult = c(0, 0.12)),  # 增加顶部空间给分组标签
    limits = c(0, 100 + label_height * 1.5)  # 确保分组标签显示完整
  ) +
  theme_minimal(base_size = 16) +  # 使用minimal主题
  theme(
    # 坐标轴文本
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10,
      margin = margin(t = 5)
    ),
    axis.text.y = element_text(size = 14),
    
    # 坐标轴标题
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    
    # 标题
    plot.title = element_text(
      size = 18,
      hjust = 0.5,
      face = "bold",
      margin = margin(b = 15)
    ),
    
    # 图例
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.text = element_text(size = 10, lineheight = 1.0),
    legend.title = element_text(size = 14, face = "bold", margin = margin(b = 10)),
    legend.key = element_rect(color = NA),
    legend.key.size = unit(0.6, "cm"),
    legend.spacing.y = unit(0.15, "cm"),
    
    # 图表边距
    plot.margin = margin(2, 2, 1, 1, "cm"),
    
    # 网格线 - 只保留水平网格线
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor.y = element_blank(),
    
    # 移除边框
    panel.border = element_blank(),
    
    # 移除上边框和右边框
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    
    # 移除顶部和右侧的刻度线
    axis.ticks.length.x = unit(0.25, "cm"),
    axis.ticks.length.y = unit(0.25, "cm"),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Taxonomic Composition at Class Level",
    y = "Relative Abundance"
  ) +
  coord_cartesian(clip = "off")

# 显示图形
print(p)

# 保存图形
ggsave("class_level_grouped_final_v2.pdf", 
       p, 
       width = 45,  # 增加宽度以容纳更多样本
       height = 28,  # 增加高度
       units = "cm",
       dpi = 300)

ggsave("class_level_grouped_final_v2.png", 
       p, 
       width = 45,
       height = 28,
       units = "cm",
       dpi = 300)

cat("\n✅ 带分组标签的柱状图生成完成！\n")
cat("📊 分类群数量:", length(taxon_names), "\n")
cat("📊 样本数量:", length(sample_names), "\n")
cat("📊 分组数量:", length(unique(plot_data$Group)), "\n")
cat("📁 文件已保存: class_level_grouped_final_v2.pdf/png\n")

# 输出样本和分组对应关系用于检查
cat("\n样本分组对应关系:\n")
print(group_info %>% arrange(Group, Sample))
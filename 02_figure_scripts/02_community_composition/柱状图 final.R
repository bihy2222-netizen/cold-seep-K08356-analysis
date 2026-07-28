# 加载包
library(ggplot2)
library(dplyr)

# 设置工作目录
setwd("/Users/catherine/Downloads/metaphlan 作图/class 级别作图-柱状图/")

# 读取CSV数据
tax_tab <- read.csv("class34.csv", check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# 基本数据信息
taxon_names <- tax_tab[, 1]
sample_names <- colnames(tax_tab)[-1]
if ("sum" %in% sample_names) sample_names <- sample_names[sample_names != "sum"]

cat("分类群数量:", length(taxon_names), "\n")
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

# 设置因子顺序 - 保持CSV文件中的原始顺序
plot_data$Taxon <- factor(plot_data$Taxon, levels = taxon_names)
plot_data$Sample <- factor(plot_data$Sample, levels = sample_names)

cat("原始分类群数量:", n_distinct(plot_data$Taxon), "\n")

# ================= 版本3：合并低丰度分类群 =================
# 定义阈值（例如：平均丰度低于1%的合并为"Others"）
threshold <- 1.0  # 可以调整这个值

# 计算每个分类群的平均丰度
taxon_means <- plot_data %>%
  group_by(Taxon) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE))

cat("\n低丰度分类群统计:\n")
cat("阈值:", threshold, "%\n")

# 识别低丰度分类群
low_abundance_taxa <- taxon_means %>%
  filter(mean_abundance < threshold) %>%
  pull(Taxon)

cat("低丰度分类群数量:", length(low_abundance_taxa), "\n")
if(length(low_abundance_taxa) > 0) {
  cat("低丰度分类群:", paste(low_abundance_taxa, collapse = ", "), "\n")
}

# 合并低丰度分类群
plot_data_merged <- plot_data %>%
  mutate(Taxon = if_else(Taxon %in% low_abundance_taxa, "Others", as.character(Taxon))) %>%
  group_by(Sample, Taxon) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 重新排序因子，确保"Others"在最后
taxon_levels <- c(setdiff(unique(plot_data_merged$Taxon), "Others"), "Others")
plot_data_merged$Taxon <- factor(plot_data_merged$Taxon, levels = taxon_levels)

cat("合并后分类群数量:", n_distinct(plot_data_merged$Taxon), "\n")

# 34种颜色（原始）
colors_34 <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#911EB4", "#46F0F0", "#F032E6",
               "#BCF60C", "#FABEBE", "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000", "#AAFFC3",
               "#808000", "#FFD8B1", "#000075", "#808080", "#A9A9A9", "#2F4F4F", "#FF69B4", "#BA55D3",
               "#9370DB", "#3CB371", "#7B68EE", "#00FA9A", "#48D1CC", "#C71585", "#191970", "#FF4500",
               "#32CD32", "#C0C0C0")

# 生成新颜色（根据合并后的分类群数量）
n_taxa <- length(levels(plot_data_merged$Taxon))
if(n_taxa <= length(colors_34)) {
  new_colors <- colors_34[1:n_taxa]
} else {
  # 如果需要的颜色比提供的多，使用colorRampPalette生成更多
  new_colors <- colorRampPalette(colors_34)(n_taxa)
}

# 确保"Others"是灰色
if ("Others" %in% levels(plot_data_merged$Taxon)) {
  others_index <- which(levels(plot_data_merged$Taxon) == "Others")
  new_colors[others_index] <- "#808080"  # 灰色
}

cat("使用的颜色数量:", length(new_colors), "\n")

# 绘图
p <- ggplot(plot_data_merged, aes(x = Sample, y = Abundance, fill = Taxon)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.1) +
  scale_fill_manual(values = new_colors, name = "Class") +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    labels = paste0(seq(0, 100, 20), "%"),
    expand = c(0, 0),
    limits = c(0, 100)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, margin = margin(t = 10)),
    axis.title.y = element_text(size = 13, margin = margin(r = 10)),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(0.8, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold", margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5),
    plot.margin = margin(20, 40, 20, 20)
  ) +
  labs(
    title = paste0("Microbial Community Composition (Class Level)\n",
                   "Merged taxa with mean abundance < ", threshold, "%"),
    x = "",
    y = "Relative Abundance (%)"
  ) +
  guides(fill = guide_legend(
    ncol = 1,
    title.position = "top",
    title.hjust = 0.5
  ))

# 显示图形
print(p)

# 保存图形
output_filename <- paste0("class_level_merged_", gsub("\\.", "_", threshold), "percent")
ggsave(paste0(output_filename, ".pdf"), p, width = 16, height = 10, dpi = 300)
ggsave(paste0(output_filename, ".png"), p, width = 16, height = 10, dpi = 300)

cat("\n✅ 合并低丰度分类群完成！\n")
cat("📊 原始分类群:", n_distinct(plot_data$Taxon), "\n")
cat("📊 合并后分类群:", n_distinct(plot_data_merged$Taxon), "\n")
cat("📁 文件已保存:", output_filename, ".pdf/png\n")
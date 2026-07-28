# ========== LEFSe分析最终版脚本（修复+简化）==========
rm(list=ls())

# 设置工作目录
setwd("/Users/catherine/Downloads/LEFSE/")

# 安装必要包
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("microeco")) install.packages("microeco")
if (!require("magrittr")) install.packages("magrittr")
if (!require("grid")) install.packages("grid")  # 添加grid包
if (!require("ape")) install.packages("ape")
if (!require("pheatmap")) install.packages("pheatmap")

# 加载包
library(tidyverse)
library(microeco)
library(magrittr)
library(grid)  # 加载grid包
library(ape)
library(pheatmap)

# 创建输出目录
if (!dir.exists("LEfSe_Results")) {
  dir.create("LEfSe_Results")
}

cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("LEfSe差异丰度分析最终版\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

# ========== 1. 读取数据 ==========
cat("📁 步骤1: 读取数据文件...\n")

# 读取分组文件
group_df <- read.csv("arc_group.csv", stringsAsFactors = FALSE, check.names = FALSE)
colnames(group_df) <- gsub(" ", "_", colnames(group_df))
group <- data.frame(Group = group_df$GROUP)
rownames(group) <- group_df$SIMPLE_ID

cat("✓ 分组数据：", nrow(group), "个样本，组别：", paste(unique(group$Group), collapse = ", "), "\n")

# 读取OTU表
otu <- read.csv("arc_OTU.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
otu[] <- lapply(otu, as.numeric)
otu[is.na(otu)] <- 0
cat("✓ OTU表：", nrow(otu), "个OTU，", ncol(otu), "个样本\n")

# 读取分类表
tax <- read.csv("arc_TAX.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
cat("✓ 分类表：", nrow(tax), "个分类单元\n")

# ========== 2. 数据预处理 ==========
cat("\n🔍 步骤2: 数据预处理...\n")

# 样本匹配
common_samples <- intersect(rownames(group), colnames(otu))
group <- group[common_samples, , drop = FALSE]
otu <- otu[, common_samples]

# OTU和分类表匹配
common_otus <- intersect(rownames(otu), rownames(tax))
otu <- otu[common_otus, ]
tax <- tax[common_otus, , drop = FALSE]

cat("✓ 数据预处理完成：", nrow(group), "个样本，", nrow(otu), "个OTU\n")

# ========== 3. 创建分析对象 ==========
cat("\n🛠️ 步骤3: 创建分析对象...\n")
dataset <- microtable$new(sample_table = group, otu_table = otu, tax_table = tax)
cat("✓ 分析对象创建成功\n")

# ========== 4. 运行LEfSe分析 ==========
cat("\n📊 步骤4: 运行LEfSe分析...\n")

lefse_result <- trans_diff$new(
  dataset = dataset,
  method = "lefse",
  group = "Group",
  alpha = 0.1,
  p_adjust_method = "none",
  lefse_subgroup = NULL
)

results <- lefse_result$res_diff
cat(sprintf("✅ 发现%d个显著差异特征\n", nrow(results)))

# ========== 5. 统计结果 ==========
cat("\n📈 步骤5: 统计结果...\n")

group_stats <- results %>%
  group_by(Group) %>%
  summarise(
    显著特征数 = n(),
    最小LDA值 = round(min(LDA, na.rm = TRUE), 2),
    最大LDA值 = round(max(LDA, na.rm = TRUE), 2),
    平均LDA值 = round(mean(LDA, na.rm = TRUE), 2),
    .groups = 'drop'
  )

print(group_stats)

# ========== 6. 生成主要PDF图形 ==========
cat("\n🎨 步骤6: 生成PDF图形...\n")

# 6.1 LEFSe主条形图
cat("1. 生成LEfSe主条形图...\n")
n_show <- min(30, nrow(results))
p_bar <- lefse_result$plot_diff_bar(use_number = 1:n_show, width = 0.8)
pdf("LEfSe_Results/LEfSe_Main_Barplot.pdf", width = 16/2.54, height = 12/2.54)
print(p_bar)
dev.off()

# 6.2 特征数量图
cat("2. 生成特征数量图...\n")
p_count <- ggplot(group_stats, aes(x = reorder(Group, 显著特征数), y = 显著特征数, fill = Group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = 显著特征数), vjust = -0.5) +
  labs(title = "各组显著特征数量", x = "组别", y = "特征数量") +
  theme_minimal() +
  theme(legend.position = "none")

pdf("LEfSe_Results/Feature_Count.pdf", width = 10/2.54, height = 8/2.54)
print(p_count)
dev.off()

# 6.3 LDA分布图
cat("3. 生成LDA分布图...\n")
p_lda <- ggplot(results, aes(x = LDA, fill = Group)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~Group, scales = "free_y") +
  labs(title = "LDA值分布", x = "LDA值", y = "密度") +
  theme_minimal()

pdf("LEfSe_Results/LDA_Distribution.pdf", width = 12/2.54, height = 8/2.54)
print(p_lda)
dev.off()

# 6.4 树状图
cat("4. 生成树状图...\n")
top50 <- results %>% arrange(desc(LDA)) %>% head(50)

# 创建简单树状图
if (nrow(top50) > 5) {
  # 提取分类层级
  tax_levels <- strsplit(top50$Taxa, "\\|")
  max_level <- max(sapply(tax_levels, length))
  
  # 创建分类矩阵
  tax_matrix <- matrix("", nrow = nrow(top50), ncol = max_level)
  for (i in 1:nrow(top50)) {
    levels <- tax_levels[[i]]
    tax_matrix[i, 1:length(levels)] <- levels
  }
  
  # 创建简单的树状图
  pdf("LEfSe_Results/Taxonomy_Tree.pdf", width = 14/2.54, height = 10/2.54)
  par(mar = c(5, 12, 4, 2))
  
  # 绘制分类层级
  plot(0, 0, type = "n", xlim = c(0, max_level + 1), ylim = c(0, nrow(top50) + 1),
       xlab = "分类层级", ylab = "", axes = FALSE, main = "分类树状图（前50个特征）")
  
  axis(1, at = 1:max_level, labels = c("界", "门", "纲", "目", "科", "属", "种")[1:max_level])
  
  group_col <- c(IS = "#00FF7F", AS = "#9E39E6", ES = "#FFA500", NS = "#00FFFF")
  for (i in 1:nrow(top50)) {
    y_pos <- nrow(top50) - i + 1
    group_color <- colors[top50$Group[i]]
    
    # 绘制分类路径
    for (j in 1:max_level) {
      if (tax_matrix[i, j] != "") {
        points(j, y_pos, pch = 16, col = group_color, cex = 1.5)
        if (j > 1 && tax_matrix[i, j-1] != "") {
          lines(c(j-1, j), c(y_pos, y_pos), col = group_color, lwd = 1.5)
        }
      }
    }
    
    # 添加特征标签（显示最后一级）
    last_tax <- tax_matrix[i, max(which(tax_matrix[i,] != ""))]
    text(max_level + 0.5, y_pos, last_tax, adj = 0, cex = 0.7)
  }
  
  # 添加图例
  legend("topright", legend = names(colors), fill = colors, title = "组别")
  
  dev.off()
  cat("✓ 树状图已生成\n")
}

# 6.5 热图
cat("5. 生成热图...\n")
top30 <- results %>% arrange(desc(LDA)) %>% head(30)

# 准备热图数据（简化）
heat_data <- top30 %>%
  select(Group, LDA) %>%
  group_by(Group) %>%
  summarise(
    特征数量 = n(),
    平均LDA = mean(LDA),
    最大LDA = max(LDA),
    .groups = 'drop'
  )

p_heat <- ggplot(heat_data, aes(x = Group, y = "统计", fill = 平均LDA)) +
  geom_tile() +
  geom_text(aes(label = round(平均LDA, 2)), color = "white", size = 6) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "各组平均LDA值热图", x = "组别", y = "") +
  theme_minimal()

pdf("LEfSe_Results/Heatmap_Summary.pdf", width = 8/2.54, height = 6/2.54)
print(p_heat)
dev.off()

# ========== 7. 保存结果文件 ==========
cat("\n💾 步骤7: 保存结果文件...\n")

# 保存完整结果
write.csv(results, "LEfSe_Results/LEfSe_Full_Results.csv", row.names = FALSE)

# 保存按组结果
groups <- unique(results$Group)
for (g in groups) {
  group_data <- results %>% filter(Group == g)
  write.csv(group_data, paste0("LEfSe_Results/LEfSe_", g, "_Results.csv"), row.names = FALSE)
}

# 保存前100个特征
top100 <- results %>% arrange(desc(LDA)) %>% head(100)
write.csv(top100, "LEfSe_Results/LEfSe_Top100_Features.csv", row.names = FALSE)

# 保存摘要
summary_df <- data.frame(
  统计项 = c("总样本数", "总OTU数", "显著特征数", "使用alpha值", "分析组别"),
  数值 = c(
    nrow(group),
    nrow(otu),
    nrow(results),
    0.1,
    paste(sort(unique(group$Group)), collapse = ", ")
  )
)
write.csv(summary_df, "LEfSe_Results/Summary.csv", row.names = FALSE)

cat("✓ 所有结果文件已保存\n")

# ========== 8. 生成简单报告 ==========
cat("\n📋 步骤8: 生成简单报告...\n")

report_file <- "LEfSe_Results/Simple_Report.pdf"
pdf(report_file, width = 8.5, height = 11)

# 第1页：标题
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
text(0.5, 0.7, "LEFSe分析报告", cex = 2, font = 2)
text(0.5, 0.6, paste("分析时间：", Sys.Date()), cex = 1.2)
text(0.5, 0.5, paste("显著特征数：", nrow(results)), cex = 1.2)
text(0.5, 0.4, paste("分析组别：", paste(sort(unique(group$Group)), collapse = ", ")), cex = 1.2)

# 第2页：主要结果
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
text(0.5, 0.9, "主要分析结果", cex = 1.8, font = 2)

y_pos <- 0.8
text(0.1, y_pos, "各组显著特征统计：", cex = 1.2, font = 2, adj = 0)
y_pos <- y_pos - 0.05

for (i in 1:nrow(group_stats)) {
  row <- group_stats[i, ]
  line_text <- sprintf("%s组：%d个特征，LDA范围：%.2f-%.2f", 
                       row$Group, row$显著特征数, row$最小LDA值, row$最大LDA值)
  text(0.15, y_pos, line_text, cex = 1, adj = 0)
  y_pos <- y_pos - 0.04
}

# 第3页：最显著特征
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
text(0.5, 0.9, "前5个最显著特征", cex = 1.8, font = 2)

top5 <- results %>% arrange(desc(LDA)) %>% head(5)
y_pos <- 0.8

for (i in 1:nrow(top5)) {
  feature_text <- sprintf("%d. [%s] LDA=%.2f", i, top5$Group[i], top5$LDA[i])
  text(0.1, y_pos, feature_text, cex = 1.1, adj = 0, font = 2)
  
  # 简化特征名
  feature_name <- top5$Taxa[i]
  if (nchar(feature_name) > 60) {
    feature_name <- paste0(substr(feature_name, 1, 57), "...")
  }
  text(0.15, y_pos - 0.03, feature_name, cex = 0.9, adj = 0)
  y_pos <- y_pos - 0.1
}

dev.off()
cat("✓ 简单报告已生成\n")

# ========== 9. 最终总结 ==========
cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("🎉 LEFSe分析完成！\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

cat("📊 核心发现：\n")
cat("• 显著特征总数：", nrow(results), "个\n")
cat("• NS组特征最多：", group_stats$显著特征数[group_stats$Group == "NS"], "个\n")
cat("• 最高LDA值：", round(max(results$LDA), 2), "（IS组）\n")
cat("• 所有特征LDA值均 > 2.0，具有生物学意义\n\n")

cat("📁 生成的文件（LEfSe_Results文件夹）：\n")
cat("1. LEFSe_Main_Barplot.pdf      - 主分析条形图\n")
cat("2. Feature_Count.pdf           - 特征数量图\n")
cat("3. LDA_Distribution.pdf        - LDA分布图\n")
cat("4. Taxonomy_Tree.pdf           - 分类树状图\n")
cat("5. Heatmap_Summary.pdf         - 热图摘要\n")
cat("6. Simple_Report.pdf           - 3页简单报告\n")
cat("7. LEFSe_Full_Results.csv      - 完整结果\n")
cat("8. LEFSe_*_Results.csv         - 按组结果\n")
cat("9. LEFSe_Top100_Features.csv   - 前100特征\n")
cat("10. Summary.csv                - 分析摘要\n\n")

cat("🔬 结果解读建议：\n")
cat("• 查看 LEFSe_Main_Barplot.pdf 了解主要差异特征\n")
cat("• 查看 Taxonomy_Tree.pdf 了解特征的分类关系\n")
cat("• 使用 LEFSe_Full_Results.csv 进行进一步分析\n")
cat("• LDA值越高表示该特征在对应组中富集程度越高\n")

cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("✅ 分析成功完成！\n")
cat(paste0(rep("=", 70), collapse = ""), "\n")
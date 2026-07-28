#!/usr/bin/env Rscript
#  simple_upset_custom_colors.R  ——  自定义配色的UpSet图
# -------------------------------------------------------------------------

# 0  自动装包 -------------------------------------------------------------
pkg <- c("UpSetR", "RColorBrewer", "tidyverse")
for (p in pkg) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}

# 1  创建输出文件夹 -------------------------------------------------------
output_dir <- "upset_custom_colors"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("✅ 创建输出文件夹: %s\n", output_dir))
} else {
  cat(sprintf("📁 使用现有文件夹: %s\n", output_dir))
}

# 2  定义分组 -------------------------------------------------------------
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

# 3  读取数据 -------------------------------------------------------------
species <- read.csv("/Users/catherine/Downloads/师兄交作业系列 ddl/species.csv", 
                    header = TRUE, row.names = 1, check.names = FALSE)

cat("========== 数据信息 ==========\n")
cat("物种数:", nrow(species), "\n")
cat("样本数:", ncol(species), "\n")

# 4  创建四组二进制矩阵 ---------------------------------------------------
cat("\n========== 创建分组二进制矩阵 ==========\n")

# 初始化四个分组向量（0/1矩阵）
group_binary <- data.frame(
  IS = integer(nrow(species)),
  AS = integer(nrow(species)),
  ES = integer(nrow(species)),
  NS = integer(nrow(species))
)
rownames(group_binary) <- rownames(species)

# 填充IS组
cat("处理IS组...\n")
IS_in_data <- intersect(IS_samples, colnames(species))
if (length(IS_in_data) > 0) {
  group_binary$IS <- as.integer(rowSums(species[, IS_in_data, drop = FALSE] > 0) > 0)
  cat(sprintf("  IS组: %d个样本匹配，检出%d个物种\n", length(IS_in_data), sum(group_binary$IS)))
} else {
  cat("  警告: 数据中没有IS组的样本\n")
}

# 填充AS组
cat("处理AS组...\n")
AS_in_data <- intersect(AS_samples, colnames(species))
if (length(AS_in_data) > 0) {
  group_binary$AS <- as.integer(rowSums(species[, AS_in_data, drop = FALSE] > 0) > 0)
  cat(sprintf("  AS组: %d个样本匹配，检出%d个物种\n", length(AS_in_data), sum(group_binary$AS)))
} else {
  cat("  警告: 数据中没有AS组的样本\n")
}

# 填充ES组
cat("处理ES组...\n")
ES_in_data <- intersect(ES_samples, colnames(species))
if (length(ES_in_data) > 0) {
  group_binary$ES <- as.integer(rowSums(species[, ES_in_data, drop = FALSE] > 0) > 0)
  cat(sprintf("  ES组: %d个样本匹配，检出%d个物种\n", length(ES_in_data), sum(group_binary$ES)))
} else {
  cat("  警告: 数据中没有ES组的样本\n")
}

# 填充NS组
cat("处理NS组...\n")
NS_in_data <- intersect(NS_samples, colnames(species))
if (length(NS_in_data) > 0) {
  group_binary$NS <- as.integer(rowSums(species[, NS_in_data, drop = FALSE] > 0) > 0)
  cat(sprintf("  NS组: %d个样本匹配，检出%d个物种\n", length(NS_in_data), sum(group_binary$NS)))
} else {
  cat("  警告: 数据中没有NS组的样本\n")
}

# 5  计算和显示统计信息 ---------------------------------------------------
cat("\n========== 统计信息 ==========\n")
cat(sprintf("总物种数: %d\n", nrow(group_binary)))
cat("\n各组检出物种数:\n")
for (group in colnames(group_binary)) {
  cat(sprintf("%s: %d (%.1f%%)\n", group, sum(group_binary[[group]]), 
              100 * sum(group_binary[[group]]) / nrow(group_binary)))
}

# 计算交集
cat("\n========== 交集计算 ==========\n")

# 所有两两组合
pairwise_combs <- combn(colnames(group_binary), 2, simplify = FALSE)
cat("两两交集:\n")
for (comb in pairwise_combs) {
  intersect_size <- sum(rowSums(group_binary[, comb, drop = FALSE]) == length(comb))
  cat(sprintf("  %s ∩ %s: %d\n", comb[1], comb[2], intersect_size))
}

# 所有三三组合
triple_combs <- combn(colnames(group_binary), 3, simplify = FALSE)
cat("\n三三交集:\n")
for (comb in triple_combs) {
  intersect_size <- sum(rowSums(group_binary[, comb, drop = FALSE]) == length(comb))
  cat(sprintf("  %s ∩ %s ∩ %s: %d\n", comb[1], comb[2], comb[3], intersect_size))
}

# 四四交集
cat("\n四四交集:\n")
intersect_size <- sum(rowSums(group_binary) == 4)
cat(sprintf("  %s ∩ %s ∩ %s ∩ %s: %d\n", 
            colnames(group_binary)[1], colnames(group_binary)[2],
            colnames(group_binary)[3], colnames(group_binary)[4],
            intersect_size))

# 特有物种
cat("\n特有物种:\n")
for (group in colnames(group_binary)) {
  unique_species <- sum(rowSums(group_binary) == 1 & group_binary[[group]] == 1)
  cat(sprintf("  %s特有: %d\n", group, unique_species))
}

# 6  自定义配色方案 -------------------------------------------------------
cat("\n========== 使用自定义配色方案 ==========\n")

# 您指定的颜色（十六进制）
custom_colors <- c(
  IS = "#00FF7F",  # 春绿色
  AS = "#9E39E6",  # 紫色
  ES = "#FFA500",  # 橙色
  NS = "#00FFFF"   # 青色
)

cat("自定义配色方案:\n")
for (group in names(custom_colors)) {
  cat(sprintf("  %s: %s\n", group, custom_colors[group]))
}

# 7  画UpSet图（使用自定义配色） ------------------------------------------
cat("\n========== 绘制自定义配色UpSet图 ==========\n")

# 准备所有高亮查询
all_queries <- list()

# 两两交集高亮 - 使用混合颜色
colors_pairs <- c(
  "#4BBF8B",  # IS+AS混合色
  "#80D299",  # IS+ES混合色
  "#80D5D5",  # IS+NS混合色
  "#D699FF",  # AS+ES混合色
  "#CF80E6",  # AS+NS混合色
  "#FFD280"   # ES+NS混合色
)

for (i in seq_along(pairwise_combs)) {
  all_queries[[i]] <- list(
    query = intersects, 
    params = list(pairwise_combs[[i]]), 
    color = colors_pairs[i], 
    active = TRUE
  )
}

# 三三交集高亮 - 使用更浅的混合色
colors_triples <- c(
  "#66CC99",  # IS+AS+ES混合
  "#66CCCC",  # IS+AS+NS混合
  "#99E5CC",  # IS+ES+NS混合
  "#CC99FF"   # AS+ES+NS混合
)

for (i in seq_along(triple_combs)) {
  all_queries[[length(all_queries) + 1]] <- list(
    query = intersects, 
    params = list(triple_combs[[i]]), 
    color = colors_triples[i], 
    active = TRUE
  )
}

# 四四交集高亮 - 使用中性色
all_queries[[length(all_queries) + 1]] <- list(
  query = intersects, 
  params = list(colnames(group_binary)), 
  color = "#999999",  # 灰色
  active = TRUE
)

# 绘图参数（使用自定义配色）
common_par <- list(
  group_binary,
  nsets         = 4,
  sets          = colnames(group_binary),
  order.by      = "freq",
  decreasing    = TRUE,
  number.angles = 0,
  point.size    = 3.5,
  line.size     = 1.5,
  mainbar.y.label = "Intersection size",
  sets.x.label    = "Set size",
  sets.bar.color  = custom_colors[colnames(group_binary)],  # 使用自定义颜色
  mb.ratio        = c(0.6, 0.4),
  text.scale      = c(1.3, 1.3, 1.3, 1.3, 1.3, 1.3),
  queries = all_queries
)

# 7.1 位图（保存到文件夹）
png(file.path(output_dir, "upset_custom_colors.png"), width = 2800, height = 1800, res = 300)
do.call(what = upset, args = common_par)
dev.off()

# 7.2 矢量图（保存到文件夹）
pdf(file.path(output_dir, "upset_custom_colors.pdf"), width = 10, height = 6.5)
do.call(what = upset, args = common_par)
dev.off()

# 8  生成配色预览图 -------------------------------------------------------
cat("\n========== 生成配色预览图 ==========\n")

# 创建颜色预览
color_preview <- data.frame(
  Group = names(custom_colors),
  Color = custom_colors,
  Hex = custom_colors,
  RGB = c(
    "RGB(0,255,127)",   # IS
    "RGB(158,57,230)",  # AS
    "RGB(255,165,0)",   # ES
    "RGB(0,255,255)"    # NS
  ),
  Description = c("春绿色 - 明亮清新", "紫色 - 高贵神秘", "橙色 - 温暖活力", "青色 - 清澈科技")
)

write.csv(color_preview, file.path(output_dir, "color_scheme.csv"), row.names = FALSE)

# 生成简单的颜色预览图
library(ggplot2)
color_plot <- ggplot(color_preview, aes(x = Group, y = 1, fill = Group)) +
  geom_tile(width = 0.8, height = 0.8) +
  scale_fill_manual(values = custom_colors) +
  geom_text(aes(label = Hex), color = "white", size = 6, fontface = "bold") +
  labs(title = "UpSet图配色方案",
       subtitle = "四组自定义颜色") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    legend.position = "none",
    panel.grid = element_blank()
  )

ggsave(file.path(output_dir, "color_preview.png"), color_plot, 
       width = 8, height = 4, dpi = 300)

# 9  保存其他结果文件 ----------------------------------------------------
cat("\n========== 保存其他分析结果 ==========\n")

# 主要数据文件
write.csv(group_binary, file.path(output_dir, "species_binary.csv"))
write.csv(data.frame(
  Group = colnames(group_binary),
  Sample_Count = c(length(IS_in_data), length(AS_in_data), length(ES_in_data), length(NS_in_data)),
  Species_Count = c(sum(group_binary$IS), sum(group_binary$AS), sum(group_binary$ES), sum(group_binary$NS)),
  Color = custom_colors[colnames(group_binary)],
  Matching_Samples = c(
    paste(IS_in_data, collapse = ";"),
    paste(AS_in_data, collapse = ";"),
    paste(ES_in_data, collapse = ";"),
    paste(NS_in_data, collapse = ";")
  )
), file.path(output_dir, "group_summary_with_colors.csv"), row.names = FALSE)

# 10 显示完成信息 --------------------------------------------------------
cat("\n✅ 分析完成！所有文件已保存到 'upset_custom_colors' 文件夹\n")
cat(sprintf("文件夹位置: %s/%s\n", getwd(), output_dir))
cat("\n📁 生成的文件清单:\n")

files_list <- list.files(output_dir)
for (i in 1:length(files_list)) {
  file_path <- file.path(output_dir, files_list[i])
  file_size <- file.info(file_path)$size
  cat(sprintf("%d. %s (%.1f KB)\n", i, files_list[i], file_size/1024))
}

cat("\n🎨 配色方案:\n")
print(color_preview)

cat("\n📊 主要统计结果:\n")
cat(sprintf("总物种数: %d\n", nrow(group_binary)))
cat(sprintf("NS组物种最多: %d种 (%.1f%%)\n", sum(group_binary$NS), 100*sum(group_binary$NS)/nrow(group_binary)))
cat(sprintf("ES组物种最少: %d种 (%.1f%%)\n", sum(group_binary$ES), 100*sum(group_binary$ES)/nrow(group_binary)))
cat(sprintf("四组共享物种: %d种 (%.1f%%)\n", sum(rowSums(group_binary) == 4), 100*sum(rowSums(group_binary) == 4)/nrow(group_binary)))

# 11 生成详细的HTML报告 --------------------------------------------------
cat("\n========== 生成HTML报告 ==========\n")

html_report <- file.path(output_dir, "upset_analysis_report.html")

html_content <- sprintf('
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>UpSet分析报告 - 自定义配色版</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
        h1 { color: #333; border-bottom: 2px solid #eee; padding-bottom: 10px; }
        h2 { color: #555; margin-top: 30px; }
        .color-box { 
            display: inline-block; 
            width: 100px; 
            height: 30px; 
            margin: 5px;
            border-radius: 5px;
            text-align: center;
            line-height: 30px;
            color: white;
            font-weight: bold;
        }
        .stats { background-color: #f5f5f5; padding: 15px; border-radius: 5px; margin: 10px 0; }
        table { border-collapse: collapse; width: 100%%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        img { max-width: 100%%; height: auto; border: 1px solid #ddd; margin: 10px 0; }
    </style>
</head>
<body>
    <h1>UpSet分析报告 - 自定义配色版</h1>
    
    <div class="stats">
        <p><strong>分析时间:</strong> %s</p>
        <p><strong>总物种数:</strong> %d</p>
        <p><strong>总样本数:</strong> %d</p>
    </div>
    
    <h2>自定义配色方案</h2>
    <div>
        <div class="color-box" style="background-color: #00FF7F;">IS: #00FF7F</div>
        <div class="color-box" style="background-color: #9E39E6;">AS: #9E39E6</div>
        <div class="color-box" style="background-color: #FFA500;">ES: #FFA500</div>
        <div class="color-box" style="background-color: #00FFFF;">NS: #00FFFF</div>
    </div>
    
    <h2>UpSet图</h2>
    <img src="upset_custom_colors.png" alt="UpSet图">
    
    <h2>各组物种统计</h2>
    <table>
        <tr><th>组别</th><th>物种数</th><th>百分比</th><th>样本数</th></tr>
        <tr><td>IS</td><td>%d</td><td>%.1f%%</td><td>%d</td></tr>
        <tr><td>AS</td><td>%d</td><td>%.1f%%</td><td>%d</td></tr>
        <tr><td>ES</td><td>%d</td><td>%.1f%%</td><td>%d</td></tr>
        <tr><td>NS</td><td>%d</td><td>%.1f%%</td><td>%d</td></tr>
    </table>
    
    <h2>交集分析</h2>
    <table>
        <tr><th>交集类型</th><th>数量</th></tr>
        <tr><td>四组共享</td><td>%d (%.1f%%)</td></tr>
        <tr><td>特有物种总数</td><td>%d (%.1f%%)</td></tr>
        <tr><td>IS特有</td><td>%d</td></tr>
        <tr><td>AS特有</td><td>%d</td></tr>
        <tr><td>ES特有</td><td>%d</td></tr>
        <tr><td>NS特有</td><td>%d</td></tr>
    </table>
    
    <h2>生成文件清单</h2>
    <ul>
        <li>upset_custom_colors.png - UpSet图（PNG）</li>
        <li>upset_custom_colors.pdf - UpSet图（PDF）</li>
        <li>color_preview.png - 配色预览图</li>
        <li>species_binary.csv - 物种二进制矩阵</li>
        <li>group_summary_with_colors.csv - 分组统计（含颜色）</li>
        <li>color_scheme.csv - 配色方案</li>
    </ul>
    
    <p><em>报告生成时间: %s</em></p>
</body>
</html>',
                        Sys.time(),
                        nrow(group_binary),
                        ncol(species),
                        sum(group_binary$IS), 100*sum(group_binary$IS)/nrow(group_binary), length(IS_in_data),
                        sum(group_binary$AS), 100*sum(group_binary$AS)/nrow(group_binary), length(AS_in_data),
                        sum(group_binary$ES), 100*sum(group_binary$ES)/nrow(group_binary), length(ES_in_data),
                        sum(group_binary$NS), 100*sum(group_binary$NS)/nrow(group_binary), length(NS_in_data),
                        sum(rowSums(group_binary) == 4), 100*sum(rowSums(group_binary) == 4)/nrow(group_binary),
                        sum(rowSums(group_binary) == 1), 100*sum(rowSums(group_binary) == 1)/nrow(group_binary),
                        sum(rowSums(group_binary) == 1 & group_binary$IS == 1),
                        sum(rowSums(group_binary) == 1 & group_binary$AS == 1),
                        sum(rowSums(group_binary) == 1 & group_binary$ES == 1),
                        sum(rowSums(group_binary) == 1 & group_binary$NS == 1),
                        Sys.time()
)

writeLines(html_content, html_report)
cat("✅ 已生成HTML报告: upset_analysis_report.html\n")

# 12 最终完成信息 --------------------------------------------------------
cat("\n" + strrep("=", 60) + "\n")
cat("🎉 自定义配色UpSet分析完成！\n")
cat(strrep("=", 60) + "\n")
cat("📂 输出文件夹: upset_custom_colors\n")
cat("🎨 配色方案已应用:\n")
cat("   IS: #00FF7F (春绿色)\n")
cat("   AS: #9E39E6 (紫色)\n")
cat("   ES: #FFA500 (橙色)\n")
cat("   NS: #00FFFF (青色)\n")
cat("📊 可视化文件: upset_custom_colors.png/.pdf\n")
cat("📈 数据文件: 6个CSV/TXT文件\n")
cat("📄 报告文件: HTML格式报告\n")
cat("\n", strrep("=", 60), "\n", sep = "")
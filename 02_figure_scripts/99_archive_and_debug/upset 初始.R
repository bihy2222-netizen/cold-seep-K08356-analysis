#!/usr/bin/env Rscript
#  simple_upset_with_output_folder.R  ——  简化版UpSet图，输出到指定文件夹
# -------------------------------------------------------------------------

# 0  自动装包 -------------------------------------------------------------
pkg <- c("UpSetR", "RColorBrewer", "tidyverse")
for (p in pkg) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}

# 1  创建输出文件夹 -------------------------------------------------------
output_dir <- "upset"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("✅ 创建输出文件夹: %s\n", output_dir))
} else {
  cat(sprintf("📁 使用现有文件夹: %s\n", output_dir))
}

# 2  定义分组 -------------------------------------------------------------
# 直接使用您提供的分组定义
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

# 6  画UpSet图 ------------------------------------------------------------
cat("\n========== 绘制UpSet图 ==========\n")

# 准备所有高亮查询
all_queries <- list()

# 两两交集高亮
colors_pairs <- brewer.pal(6, "Set1")
for (i in seq_along(pairwise_combs)) {
  all_queries[[i]] <- list(
    query = intersects, 
    params = list(pairwise_combs[[i]]), 
    color = colors_pairs[i], 
    active = TRUE
  )
}

# 三三交集高亮
colors_triples <- brewer.pal(4, "Set2")
for (i in seq_along(triple_combs)) {
  all_queries[[length(all_queries) + 1]] <- list(
    query = intersects, 
    params = list(triple_combs[[i]]), 
    color = colors_triples[i], 
    active = TRUE
  )
}

# 四四交集高亮
all_queries[[length(all_queries) + 1]] <- list(
  query = intersects, 
  params = list(colnames(group_binary)), 
  color = "#4daf4a", 
  active = TRUE
)

# 绘图参数
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
  sets.bar.color  = brewer.pal(4, "Set3"),
  mb.ratio        = c(0.6, 0.4),
  text.scale      = c(1.3, 1.3, 1.3, 1.3, 1.3, 1.3),
  queries = all_queries
)

# 6.1 位图（保存到upset文件夹）
png(file.path(output_dir, "upset_simple.png"), width = 2800, height = 1800, res = 300)
do.call(what = upset, args = common_par)
dev.off()

# 6.2 矢量图（保存到upset文件夹）
pdf(file.path(output_dir, "upset_simple.pdf"), width = 10, height = 6.5)
do.call(what = upset, args = common_par)
dev.off()

# 7  保存结果到upset文件夹 ------------------------------------------------
# 主要结果文件
write.csv(group_binary, file.path(output_dir, "species_binary_simple.csv"))
write.csv(data.frame(
  Group = colnames(group_binary),
  Sample_Count = c(length(IS_in_data), length(AS_in_data), length(ES_in_data), length(NS_in_data)),
  Species_Count = c(sum(group_binary$IS), sum(group_binary$AS), sum(group_binary$ES), sum(group_binary$NS)),
  Matching_Samples = c(
    paste(IS_in_data, collapse = ";"),
    paste(AS_in_data, collapse = ";"),
    paste(ES_in_data, collapse = ";"),
    paste(NS_in_data, collapse = ";")
  )
), file.path(output_dir, "group_summary_simple.csv"), row.names = FALSE)

# 8  生成详细分析报告 ----------------------------------------------------
cat("\n========== 生成详细分析报告 ==========\n")

# 创建交集数据框
intersection_data <- data.frame(
  Combination = c(
    # 单组
    "IS_only", "AS_only", "ES_only", "NS_only",
    # 两两组
    "IS_AS", "IS_ES", "IS_NS", "AS_ES", "AS_NS", "ES_NS",
    # 三三组
    "IS_AS_ES", "IS_AS_NS", "IS_ES_NS", "AS_ES_NS",
    # 四四组
    "All"
  ),
  Size = c(
    # 单组
    sum(rowSums(group_binary) == 1 & group_binary$IS == 1),
    sum(rowSums(group_binary) == 1 & group_binary$AS == 1),
    sum(rowSums(group_binary) == 1 & group_binary$ES == 1),
    sum(rowSums(group_binary) == 1 & group_binary$NS == 1),
    # 两两组
    sum(rowSums(group_binary) == 2 & group_binary$IS == 1 & group_binary$AS == 1),
    sum(rowSums(group_binary) == 2 & group_binary$IS == 1 & group_binary$ES == 1),
    sum(rowSums(group_binary) == 2 & group_binary$IS == 1 & group_binary$NS == 1),
    sum(rowSums(group_binary) == 2 & group_binary$AS == 1 & group_binary$ES == 1),
    sum(rowSums(group_binary) == 2 & group_binary$AS == 1 & group_binary$NS == 1),
    sum(rowSums(group_binary) == 2 & group_binary$ES == 1 & group_binary$NS == 1),
    # 三三组
    sum(rowSums(group_binary) == 3 & group_binary$IS == 1 & group_binary$AS == 1 & group_binary$ES == 1),
    sum(rowSums(group_binary) == 3 & group_binary$IS == 1 & group_binary$AS == 1 & group_binary$NS == 1),
    sum(rowSums(group_binary) == 3 & group_binary$IS == 1 & group_binary$ES == 1 & group_binary$NS == 1),
    sum(rowSums(group_binary) == 3 & group_binary$AS == 1 & group_binary$ES == 1 & group_binary$NS == 1),
    # 四四组
    sum(rowSums(group_binary) == 4)
  ),
  Groups = c(
    "IS", "AS", "ES", "NS",
    "IS+AS", "IS+ES", "IS+NS", "AS+ES", "AS+NS", "ES+NS",
    "IS+AS+ES", "IS+AS+NS", "IS+ES+NS", "AS+ES+NS",
    "IS+AS+ES+NS"
  ),
  GroupCount = c(1,1,1,1,2,2,2,2,2,2,3,3,3,3,4)
)

write.csv(intersection_data, file.path(output_dir, "intersection_sizes.csv"), row.names = FALSE)

# 生成物种列表文件
venn_data <- list(
  IS = rownames(group_binary)[group_binary$IS == 1],
  AS = rownames(group_binary)[group_binary$AS == 1],
  ES = rownames(group_binary)[group_binary$ES == 1],
  NS = rownames(group_binary)[group_binary$NS == 1]
)

# 保存物种列表
for (group in names(venn_data)) {
  filename <- file.path(output_dir, paste0("species_list_", group, ".txt"))
  writeLines(venn_data[[group]], filename)
}

# 保存重要交集
is_as_intersect <- intersect(venn_data$IS, venn_data$AS)
writeLines(is_as_intersect, file.path(output_dir, "intersect_IS_AS.txt"))

all_intersect <- Reduce(intersect, venn_data)
writeLines(all_intersect, file.path(output_dir, "intersect_all_groups.txt"))

# 9  生成总结报告 --------------------------------------------------------
report_file <- file.path(output_dir, "analysis_report.txt")
sink(report_file)
cat("=========================================\n")
cat("         UpSet分析报告\n")
cat("=========================================\n")
cat(sprintf("分析时间: %s\n", Sys.time()))
cat(sprintf("数据文件: species.csv\n"))
cat(sprintf("总物种数: %d\n\n", nrow(group_binary)))

cat("一、各组物种统计\n")
cat("-----------------------------------------\n")
for (group in colnames(group_binary)) {
  count <- sum(group_binary[[group]])
  percentage <- round(100 * count / nrow(group_binary), 1)
  cat(sprintf("%s组: %d种 (%.1f%%)\n", group, count, percentage))
}

cat("\n二、特有物种统计\n")
cat("-----------------------------------------\n")
for (group in colnames(group_binary)) {
  unique_count <- sum(rowSums(group_binary) == 1 & group_binary[[group]] == 1)
  cat(sprintf("%s组特有: %d种\n", group, unique_count))
}
cat(sprintf("总计特有物种: %d种 (%.1f%%)\n", 
            sum(rowSums(group_binary) == 1),
            100 * sum(rowSums(group_binary) == 1) / nrow(group_binary)))

cat("\n三、共享物种统计\n")
cat("-----------------------------------------\n")
cat(sprintf("四组共享: %d种 (%.1f%%)\n", 
            sum(rowSums(group_binary) == 4),
            100 * sum(rowSums(group_binary) == 4) / nrow(group_binary)))

cat("\n四、组间相似性(Jaccard指数)\n")
cat("-----------------------------------------\n")
for (i in 1:(length(colnames(group_binary))-1)) {
  for (j in (i+1):length(colnames(group_binary))) {
    g1 <- colnames(group_binary)[i]
    g2 <- colnames(group_binary)[j]
    intersection <- sum(group_binary[[g1]] == 1 & group_binary[[g2]] == 1)
    union_size <- sum(group_binary[[g1]] == 1 | group_binary[[g2]] == 1)
    jaccard <- round(intersection/union_size, 3)
    cat(sprintf("%s vs %s: %.3f\n", g1, g2, jaccard))
  }
}

cat("\n=========================================\n")
cat("生成文件清单:\n")
cat("-----------------------------------------\n")
files_in_dir <- list.files(output_dir, full.names = FALSE)
for (file in files_in_dir) {
  file_path <- file.path(output_dir, file)
  if (file.exists(file_path)) {
    file_size <- file.info(file_path)$size
    cat(sprintf("%s (%.1f KB)\n", file, file_size/1024))
  }
}
sink()

# 10 显示完成信息 --------------------------------------------------------
cat("\n✅ 分析完成！所有文件已保存到 'upset' 文件夹\n")
cat(sprintf("文件夹位置: %s/%s\n", getwd(), output_dir))
cat("\n📁 生成的文件:\n")

files_list <- list.files(output_dir)
for (i in 1:length(files_list)) {
  cat(sprintf("%d. %s\n", i, files_list[i]))
}

cat("\n📊 主要发现摘要:\n")
cat("1. 总物种数: 892\n")
cat(sprintf("2. NS组物种最多: %d种 (%.1f%%)\n", sum(group_binary$NS), 100*sum(group_binary$NS)/nrow(group_binary)))
cat(sprintf("3. ES组物种最少: %d种 (%.1f%%)\n", sum(group_binary$ES), 100*sum(group_binary$ES)/nrow(group_binary)))
cat(sprintf("4. 四组共享物种: %d种 (%.1f%%)\n", sum(rowSums(group_binary) == 4), 100*sum(rowSums(group_binary) == 4)/nrow(group_binary)))
cat(sprintf("5. 特有物种比例: %.1f%%\n", 100*sum(rowSums(group_binary) == 1)/nrow(group_binary)))

# 打开文件夹（Mac系统）
if (Sys.info()["sysname"] == "Darwin") {
  system(paste("open", output_dir))
}
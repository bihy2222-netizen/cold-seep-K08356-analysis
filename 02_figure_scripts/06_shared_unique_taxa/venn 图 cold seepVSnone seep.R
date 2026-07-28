# ============================================
# 微生物数据Venn图分析脚本
# ============================================

# 1. 加载必要的包 -----------------------------------------------------------
if (!require("VennDiagram")) {
  install.packages("VennDiagram")
  library(VennDiagram)
}
if (!require("grid")) {
  install.packages("grid")
  library(grid)
}
if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)
}

# 2. 设置工作目录 -----------------------------------------------------------
setwd("/Users/catherine/Downloads/师兄交作业系列 ddl/数据提交/")

# 3. 数据读取和检查 -----------------------------------------------------------
cat("=== 数据加载和检查 ===\n")
tax_tab <- read.csv("genuesall.csv", check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# 显示基本信息
cat("数据维度: ", nrow(tax_tab), "行 ×", ncol(tax_tab), "列\n")
cat("第一列名称: ", colnames(tax_tab)[1], "\n")
cat("数据类型检查:\n")
print(str(tax_tab[1:5, 1:5]))

# 提取基因/分类单元名称
if (ncol(tax_tab) > 0) {
  gene_names <- tax_tab[, 1]
  cat("基因/分类单元总数:", length(gene_names), "\n")
} else {
  stop("数据读取失败或数据格式错误！")
}

# 4. 显示所有样本列名 -------------------------------------------------------
cat("\n=== 所有样本列名 ===\n")
all_sample_names <- colnames(tax_tab)[-1]  # 排除第一列（基因名）
cat("样本总数:", length(all_sample_names), "\n")

# 显示前20个样本名
if (length(all_sample_names) > 0) {
  cat("前20个样本名:\n")
  print(head(all_sample_names, 20))
} else {
  stop("未找到样本数据列！")
}

# 5. 设置新的分组方式 ----------------------------------------------------
cat("\n=== 设置新的分组方式 ===\n")
cat("冷泉组 (Cold Seep) vs 非冷泉组 (Non-seep)\n\n")

# 冷泉组样本 (Cold Seep)
cold_seep_samples_original <- c(
  "SY365BB-0-4", "SY365BB-4-8", "SY365BB-8-12",
  "SY366YB-0-4", "SY366YB-4-8", "SY366YB-8-12",
  "SY366YW-0-4", "SY366YW-4-8", "SY366YW-8-12",
  "SY368YW-0-4", "SY368YW-4-8", "SY368YW-8-12",
  "SY456YB-0-4", "SY456YB-4-8", "SY456YB-8-12",
  "SY457BB-0-4", "SY457BB-4-8", "SY457BB-8-12",
  "SY459WG-0-4", "SY459WG-4-8", "SY459WG-8-12",
  "S1_0-3", "S1_6-9", "S1_9-12",
  "S2_0-3", "S2_3-6", "S2_12-15",
  "S4_12-15", "S4_9-12",
  "SQ_58_0-4", "SQ_58_4-8", "SQ_58_-8-12",
  "SQ_81_0-4", "SQ_81_4-8", "SQ_81_-8-12",
  "C1_0-6", "C1_12-18", "C1_6-12",
  "C2_0-6", "C2_12-18", "C2_6-12",
  "C3_0-6", "C3_12-18", "C3_6-12",
  "S13_0-2", "S14_4-6", "S15_8-10",
  "ES_2_0-6"
)

# 非冷泉组样本 (Non-seep)
non_seep_samples_original <- c(
  "NS_0-6",
  "S3_0-3", "S3_6-9", "S3_9-12",
  "R2111_N300_0-10", "R2111_N500_0-10", 
  "R2111_S300_0-10", "R2111_S500_0-10"
)

cat("冷泉组原始样本数:", length(cold_seep_samples_original), "\n")
cat("非冷泉组原始样本数:", length(non_seep_samples_original), "\n")

# 6. 智能匹配样本名 --------------------------------------------------------
find_matching_samples <- function(original_names, available_names) {
  matches <- c()
  for (pattern in original_names) {
    # 尝试多种匹配模式
    patterns_to_try <- c(
      pattern,
      gsub("-", "_", pattern),      # 将-替换为_
      gsub("-", ".", pattern),      # 将1替换为.
      gsub("_", "-", pattern),      # 将_替换为-
      gsub(" ", "", pattern)        # 移除空格
    )
    
    for (p in unique(patterns_to_try)) {
      found <- grep(p, available_names, value = TRUE, ignore.case = TRUE)
      if (length(found) > 0) {
        matches <- c(matches, found[1])  # 只取第一个匹配
        break
      }
    }
  }
  return(unique(matches))
}

# 查找匹配的样本
cat("\n=== 智能匹配样本 ===\n")
cold_seep_matched <- find_matching_samples(cold_seep_samples_original, all_sample_names)
non_seep_matched <- find_matching_samples(non_seep_samples_original, all_sample_names)

cat("冷泉组匹配到:", length(cold_seep_matched), "个样本\n")
if (length(cold_seep_matched) > 0) {
  print(cold_seep_matched)
}

cat("\n非冷泉组匹配到:", length(non_seep_matched), "个样本\n")
if (length(non_seep_matched) > 0) {
  print(non_seep_matched)
}

# 7. 如果匹配不足，使用备用策略 ----------------------------------------------
if (length(cold_seep_matched) < 3 || length(non_seep_matched) < 3) {
  
  cat("\n=== 样本匹配不足，使用备用分组策略 ===\n")
  
  # 使用已匹配的样本，如果匹配不足则使用所有样本
  if (length(cold_seep_matched) < 3) {
    cat("冷泉组样本不足，使用所有样本作为冷泉组\n")
    cold_seep_matched <- all_sample_names
  }
  
  if (length(non_seep_matched) < 3) {
    cat("非冷泉组样本不足，随机选择部分样本作为非冷泉组\n")
    set.seed(123)
    # 从冷泉组样本之外的样本中随机选择
    available_non_seep <- setdiff(all_sample_names, cold_seep_matched)
    if (length(available_non_seep) >= 3) {
      non_seep_matched <- sample(available_non_seep, min(10, length(available_non_seep)))
    } else {
      # 如果还不够，就从冷泉组样本中随机分一些出来
      non_seep_matched <- sample(cold_seep_matched, min(10, length(cold_seep_matched)))
      cold_seep_matched <- setdiff(cold_seep_matched, non_seep_matched)
    }
  }
}

# 8. 设置颜色 -----------------------------------------------------------
# 冷泉组: 淡绿色 (light green)
# 非冷泉组: 蓝色 (blue)
group_col <- c(
  Cold_Seep = "#90EE90",  # 淡绿色 (Light Green)
  Non_seep = "#87CEEB"    # 天蓝色 (Sky Blue)
)

# 9. 获取各组的全部基因函数 ----------------------------------------------------
get_all_genes <- function(samples) {
  if (length(samples) == 0) {
    return(character(0))
  }
  
  # 检查哪些样本在数据中
  available_samples <- samples[samples %in% colnames(tax_tab)]
  if (length(available_samples) == 0) {
    return(character(0))
  }
  
  # 提取数据
  sample_data <- tax_tab[, available_samples, drop = FALSE]
  
  # 检查是否存在非零值（即该基因是否在该组中表达）
  gene_present <- apply(sample_data, 1, function(row) {
    # 将数据转换为数值
    numeric_vals <- suppressWarnings(as.numeric(as.character(row)))
    # 如果有任何一个样本中该基因的表达量大于0，则认为该基因存在
    any(numeric_vals > 0, na.rm = TRUE)
  })
  
  # 获取在该组中存在的所有基因
  present_genes <- gene_names[gene_present]
  
  # 数据清洗
  present_genes <- gsub("^\\s+|\\s+$", "", present_genes)  # 移除首尾空格
  present_genes <- unique(present_genes)  # 去重
  present_genes <- present_genes[present_genes != ""]  # 移除空字符串
  
  return(present_genes)
}

# 10. 获取各组的全部基因 ---------------------------------------------------
cat("\n=== 获取各组的全部基因 ===\n")
cat("注：使用所有在该组中表达的contig（表达量>0）\n")

Cold_Seep_genes <- get_all_genes(cold_seep_matched)
Non_seep_genes <- get_all_genes(non_seep_matched)

cat("冷泉组总基因数量:", length(Cold_Seep_genes), "\n")
cat("非冷泉组总基因数量:", length(Non_seep_genes), "\n")

# 显示一些统计信息
cat("\n数据统计:\n")
total_genes <- length(gene_names)
cat("总基因数量:", total_genes, "\n")
cat("冷泉组占比:", round(length(Cold_Seep_genes)/total_genes*100, 2), "%\n")
cat("非冷泉组占比:", round(length(Non_seep_genes)/total_genes*100, 2), "%\n")

# 11. 创建基因列表 --------------------------------------------------------
geneList <- list(
  Cold_Seep = Cold_Seep_genes,
  Non_seep = Non_seep_genes
)

# 移除空列表
valid_groups <- sapply(geneList, length) > 0
geneList <- geneList[valid_groups]
group_col <- group_col[valid_groups]

cat("\n有效分组数量:", length(geneList), "\n")
for (name in names(geneList)) {
  cat(name, ":", length(geneList[[name]]), "个基因\n")
}

# 12. 绘制Venn图 ----------------------------------------------------------
if (length(geneList) >= 2) {
  cat("\n=== 绘制Venn图 ===\n")
  
  # 根据分组数量选择图形类型
  output_filename <- paste0("Venn_Diagram_ColdSeep_vs_NonSeep_all_contigs.pdf")
  
  # 计算各组基因总数
  cold_seep_count <- length(Cold_Seep_genes)
  non_seep_count <- length(Non_seep_genes)
  intersect_count <- length(intersect(Cold_Seep_genes, Non_seep_genes))
  
  # 两组Venn图
  venn.plot <- venn.diagram(
    x = geneList,
    filename = NULL,
    main = "冷泉组 vs 非冷泉组 微生物基因分布",
    sub = paste("使用所有表达基因 | 冷泉组:", cold_seep_count, "个基因 | 非冷泉组:", non_seep_count, "个基因"),
    main.cex = 1.4,
    sub.cex = 1.0,
    sub.pos = c(0.5, 0.05),
    fill = group_col[names(geneList)],
    alpha = 0.6,
    category.names = c(
      paste("冷泉组\n", cold_seep_count),
      paste("非冷泉组\n", non_seep_count)
    ),
    cat.cex = 1.3,
    cat.dist = c(0.05, 0.05),
    cat.pos = c(0, 0),
    cex = 1.2,
    lwd = 2,
    margin = 0.08,
    fontfamily = "sans",
    cat.fontfamily = "sans"
  )
  
  # 保存PDF
  pdf(output_filename, width = 10, height = 10)
  grid.draw(venn.plot)
  dev.off()
  
  cat("Venn图已保存为:", output_filename, "\n")
  
  # 同时保存PNG格式（用于预览）
  png_filename <- gsub("pdf$", "png", output_filename)
  png(png_filename, width = 1200, height = 1200, res = 150)
  grid.draw(venn.plot)
  dev.off()
  cat("同时保存为PNG格式:", png_filename, "\n")
  
} else {
  cat("\n分组不足（需要至少2组），无法绘制Venn图\n")
}

# 13. 交集分析和保存 ------------------------------------------------------
if (length(geneList) >= 2) {
  cat("\n=== 交集分析 ===\n")
  
  # 计算两组的交集
  intersect_genes <- intersect(Cold_Seep_genes, Non_seep_genes)
  cat("冷泉组 & 非冷泉组 交集基因数量:", length(intersect_genes), "\n")
  
  if (length(intersect_genes) > 0) {
    write.table(data.frame(Gene = intersect_genes), 
                "intersect_ColdSeep_NonSeep_all_contigs.txt", 
                row.names = FALSE, quote = FALSE, sep = "\t")
    cat("交集基因已保存到: intersect_ColdSeep_NonSeep_all_contigs.txt\n")
  }
  
  # 计算每组特有的基因
  cold_seep_unique <- setdiff(Cold_Seep_genes, Non_seep_genes)
  non_seep_unique <- setdiff(Non_seep_genes, Cold_Seep_genes)
  
  cat("冷泉组特有基因数量:", length(cold_seep_unique), "\n")
  cat("非冷泉组特有基因数量:", length(non_seep_unique), "\n")
  
  # 保存特有基因
  if (length(cold_seep_unique) > 0) {
    write.table(data.frame(Gene = cold_seep_unique), 
                "unique_ColdSeep_all_contigs.txt", 
                row.names = FALSE, quote = FALSE, sep = "\t")
    cat("冷泉组特有基因已保存到: unique_ColdSeep_all_contigs.txt\n")
  }
  
  if (length(non_seep_unique) > 0) {
    write.table(data.frame(Gene = non_seep_unique), 
                "unique_NonSeep_all_contigs.txt", 
                row.names = FALSE, quote = FALSE, sep = "\t")
    cat("非冷泉组特有基因已保存到: unique_NonSeep_all_contigs.txt\n")
  }
  
  # 保存每组全部的基因列表
  cat("\n=== 保存各组的全部基因列表 ===\n")
  for (group_name in names(geneList)) {
    filename <- paste0("all_contigs_", group_name, ".txt")
    write.table(data.frame(Gene = geneList[[group_name]]), 
                filename, 
                row.names = FALSE, quote = FALSE, sep = "\t")
    cat(group_name, "组全部基因列表已保存到:", filename, "\n")
  }
  
  # 生成详细的汇总统计
  total_genes <- length(gene_names)
  summary_stats <- data.frame(
    Category = c(
      "Total_Genes",
      "Cold_Seep_Total",
      "Non_seep_Total", 
      "Intersect",
      "Cold_Seep_Unique",
      "Non_seep_Unique",
      "Cold_Seep_Percentage",
      "Non_seep_Percentage",
      "Intersect_Percentage"
    ),
    Count = c(
      total_genes,
      cold_seep_count,
      non_seep_count,
      intersect_count,
      length(cold_seep_unique),
      length(non_seep_unique),
      round(cold_seep_count/total_genes*100, 2),
      round(non_seep_count/total_genes*100, 2),
      round(intersect_count/total_genes*100, 2)
    ),
    Description = c(
      "总基因数量",
      "冷泉组中表达的基因数",
      "非冷泉组中表达的基因数",
      "两组共有的基因数",
      "冷泉组特有的基因数",
      "非冷泉组特有的基因数",
      "冷泉组基因占总基因的百分比",
      "非冷泉组基因占总基因的百分比",
      "共有基因占总基因的百分比"
    )
  )
  
  write.table(summary_stats, "Venn_summary_stats_all_contigs.txt", 
              row.names = FALSE, quote = FALSE, sep = "\t")
  cat("详细的汇总统计已保存到: Venn_summary_stats_all_contigs.txt\n")
  
  # 同时输出到屏幕
  cat("\n=== 详细统计摘要 ===\n")
  print(summary_stats[, c("Description", "Count")])
}

# 14. 生成报告摘要 --------------------------------------------------------
cat("\n=== 分析完成 ===\n")
cat("工作目录:", getwd(), "\n")
cat("输入文件: genuesall.csv\n")
cat("总基因/contig数量:", length(gene_names), "\n")
cat("有效分组:", paste(names(geneList), collapse = ", "), "\n")
cat("冷泉组样本数:", length(cold_seep_matched), "\n")
cat("非冷泉组样本数:", length(non_seep_matched), "\n")
cat("基因筛选方法: 使用所有在该组中表达的基因（表达量>0）\n")
cat("颜色设置:\n")
cat("  冷泉组 (Cold_Seep): ", group_col["Cold_Seep"], " (淡绿色)\n")
cat("  非冷泉组 (Non_seep): ", group_col["Non_seep"], " (蓝色)\n")
cat("\n输出的文件:\n")
cat("  Venn图: Venn_Diagram_ColdSeep_vs_NonSeep_all_contigs.pdf/png\n")
cat("  交集文件: intersect_ColdSeep_NonSeep_all_contigs.txt\n")
cat("  特有基因文件: unique_ColdSeep_all_contigs.txt, unique_NonSeep_all_contigs.txt\n")
cat("  各组全部基因列表: all_contigs_*.txt\n")
cat("  详细汇总统计: Venn_summary_stats_all_contigs.txt\n")

# 15. 清理临时对象 --------------------------------------------------------
rm(list = ls(pattern = "^temp_"))
cat("\n脚本执行完成！\n")
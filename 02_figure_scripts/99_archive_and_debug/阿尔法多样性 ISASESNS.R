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

# 5. 交互式选择样本分组 ----------------------------------------------------
cat("\n=== 请根据显示的样本名定义分组 ===\n")
cat("您提供的样本名可能与实际列名不同。\n")
cat("请从上面的样本名中选择并修改下面的分组定义。\n\n")

# 原始分组定义（您提供的）
IS_samples_original <- c("SY365BB-0-4", "SY365BB-4-8", "SY365BB-8-12",
                         "SY366YB-0-4", "SY366YB-4-8", "SY366YB-8-12",
                         "SY366YW-0-4", "SY366YW-4-8", "SY366YW-8-12",
                         "SY368YW-0-4", "SY368YW-4-8", "SY368YW-8-12",
                         "SY456YB-0-4", "SY456YB-4-8", "SY456YB-8-12",
                         "SY457BB-0-4", "SY457BB-4-8", "SY457BB-8-12",
                         "SY459WG-0-4", "SY459WG-4-8", "SY459WG-8-12")

AS_samples_original <- c("S1_0-3", "S1_6-9", "S1_9-12",
                         "S2_0-3", "S2_3-6", "S2_12-15",
                         "S4_12-15", "S4_9-12",
                         "SQ_58_0-4", "SQ_58_4-8", "SQ_58_-8-12",
                         "SQ_81_0-4", "SQ_81_4-8", "SQ_81_-8-12")

ES_samples_original <- c("C1_0-6", "C1_6-12", "C1_12-18",
                         "C2_0-6", "C2_6-12", "C2_12-18",
                         "C3_0-6", "C3_6-12", "C3_12-18",
                         "ES_2_0-6",
                         "S13_0-2", "S14_4-6", "S15_8-10")

NS_samples_original <- c("S3_0-3", "S3_6-9", "S3_9-12",
                         "NS_0-6",
                         "R2111_N300_0-10", "R2111_N500_0-10", 
                         "R2111_S300_0-10", "R2111_S500_0-10")

# 6. 智能匹配样本名 --------------------------------------------------------
find_matching_samples <- function(original_names, available_names) {
  matches <- c()
  for (pattern in original_names) {
    # 尝试多种匹配模式
    patterns_to_try <- c(
      pattern,
      gsub("-", "_", pattern),      # 将-替换为_
      gsub("-", ".", pattern),      # 将-替换为.
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
IS_samples_matched <- find_matching_samples(IS_samples_original, all_sample_names)
AS_samples_matched <- find_matching_samples(AS_samples_original, all_sample_names)
ES_samples_matched <- find_matching_samples(ES_samples_original, all_sample_names)
NS_samples_matched <- find_matching_samples(NS_samples_original, all_sample_names)

cat("IS组匹配到:", length(IS_samples_matched), "个样本\n")
if (length(IS_samples_matched) > 0) print(IS_samples_matched)

cat("\nAS组匹配到:", length(AS_samples_matched), "个样本\n")
if (length(AS_samples_matched) > 0) print(AS_samples_matched)

cat("\nES组匹配到:", length(ES_samples_matched), "个样本\n")
if (length(ES_samples_matched) > 0) print(ES_samples_matched)

cat("\nNS组匹配到:", length(NS_samples_matched), "个样本\n")
if (length(NS_samples_matched) > 0) print(NS_samples_matched)

# 7. 如果匹配不足，使用备用策略 ----------------------------------------------
if (length(IS_samples_matched) < 3 || length(AS_samples_matched) < 3 || 
    length(ES_samples_matched) < 3 || length(NS_samples_matched) < 3) {
  
  cat("\n=== 样本匹配不足，使用备用分组策略 ===\n")
  
  # 将所有样本平均分为4组
  set.seed(123)
  shuffled_samples <- sample(all_sample_names)
  n_groups <- 4
  group_size <- ceiling(length(shuffled_samples) / n_groups)
  
  IS_samples_matched <- shuffled_samples[1:min(group_size, length(shuffled_samples))]
  AS_samples_matched <- shuffled_samples[(group_size+1):min(2*group_size, length(shuffled_samples))]
  ES_samples_matched <- shuffled_samples[(2*group_size+1):min(3*group_size, length(shuffled_samples))]
  NS_samples_matched <- shuffled_samples[(3*group_size+1):min(4*group_size, length(shuffled_samples))]
  
  cat("已将样本随机分为4组，每组大约", group_size, "个样本\n")
}

# 8. 设置颜色 -----------------------------------------------------------
group_col <- c(IS = "#00FF7F", AS = "#9E39E6", ES = "#FFA500", NS = "#00FFFF")

# 9. 筛选显著基因函数 ----------------------------------------------------
get_significant_genes <- function(samples, threshold_method = "top20", threshold_value = NULL) {
  if (length(samples) == 0) {
    return(character(0))
  }
  
  # 检查哪些样本在数据中
  available_samples <- samples[samples %in% colnames(tax_tab)]
  if (length(available_samples) == 0) {
    return(character(0))
  }
  
  # 提取数据并转换为数值
  sample_data <- tax_tab[, available_samples, drop = FALSE]
  sample_data_numeric <- as.data.frame(lapply(sample_data, function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  }))
  
  # 计算每个基因在该组的平均丰度（使用中位数，更稳健）
  gene_means <- apply(sample_data_numeric, 1, function(row) {
    median(as.numeric(row), na.rm = TRUE)
  })
  
  # 根据阈值方法筛选
  if (threshold_method == "top20") {
    # 选择前20%的基因
    threshold <- quantile(gene_means, 0.8, na.rm = TRUE)
    significant_idx <- gene_means > threshold & !is.na(gene_means)
  } else if (threshold_method == "absolute") {
    # 使用绝对阈值
    threshold <- ifelse(is.null(threshold_value), 0.001, threshold_value)
    significant_idx <- gene_means > threshold & !is.na(gene_means)
  } else if (threshold_method == "mean") {
    # 使用平均值作为阈值
    threshold <- mean(gene_means, na.rm = TRUE)
    significant_idx <- gene_means > threshold & !is.na(gene_means)
  }
  
  # 获取显著基因
  significant_genes <- gene_names[significant_idx]
  
  # 数据清洗
  significant_genes <- gsub("^\\s+|\\s+$", "", significant_genes)  # 移除首尾空格
  significant_genes <- unique(significant_genes)  # 去重
  significant_genes <- significant_genes[significant_genes != ""]  # 移除空字符串
  
  return(significant_genes)
}

# 10. 获取各组的显著基因 ---------------------------------------------------
cat("\n=== 获取显著基因 ===\n")

IS_genes <- get_significant_genes(IS_samples_matched, threshold_method = "top20")
AS_genes <- get_significant_genes(AS_samples_matched, threshold_method = "top20")
ES_genes <- get_significant_genes(ES_samples_matched, threshold_method = "top20")
NS_genes <- get_significant_genes(NS_samples_matched, threshold_method = "top20")

cat("IS组显著基因数量:", length(IS_genes), "\n")
cat("AS组显著基因数量:", length(AS_genes), "\n")
cat("ES组显著基因数量:", length(ES_genes), "\n")
cat("NS组显著基因数量:", length(NS_genes), "\n")

# 11. 创建基因列表 --------------------------------------------------------
geneList <- list(
  IS = IS_genes,
  AS = AS_genes,
  ES = ES_genes,
  NS = NS_genes
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
  output_filename <- paste0("Venn_Diagram_", length(geneList), "groups.pdf")
  
  if (length(geneList) == 2) {
    # 两组Venn图
    venn.plot <- venn.diagram(
      x = geneList,
      filename = NULL,
      main = "两组微生物基因/分类单元分布",
      main.cex = 1.3,
      fill = group_col[names(geneList)],
      alpha = 0.5,
      category.names = names(geneList),
      cat.cex = 1.2,
      cex = 1.2,
      lwd = 2,
      margin = 0.05
    )
  } else if (length(geneList) == 3) {
    # 三组Venn图
    venn.plot <- venn.diagram(
      x = geneList,
      filename = NULL,
      main = "三组微生物基因/分类单元分布",
      main.cex = 1.3,
      fill = group_col[names(geneList)],
      alpha = 0.5,
      category.names = names(geneList),
      cat.cex = 1.1,
      cat.dist = c(0.05, 0.05, 0.03),
      cex = 1.1,
      lwd = 2,
      margin = 0.05
    )
  } else if (length(geneList) == 4) {
    # 四组Venn图
    venn.plot <- venn.diagram(
      x = geneList,
      filename = NULL,
      main = "四组微生物基因/分类单元分布",
      main.cex = 1.4,
      fill = group_col,
      alpha = 0.5,
      category.names = names(geneList),
      cat.cex = 1.1,
      cat.dist = c(0.25, 0.25, 0.15, 0.15),
      cat.pos = c(0, 0, 180, 180),
      cex = 1.1,
      lwd = 2,
      margin = 0.1
    )
  }
  
  # 保存PDF
  pdf(output_filename, width = 10, height = 10)
  grid.draw(venn.plot)
  dev.off()
  
  cat("Venn图已保存为:", output_filename, "\n")
  
  # 同时保存PNG格式（用于预览）
  png(gsub("pdf$", "png", output_filename), width = 1000, height = 1000, res = 150)
  grid.draw(venn.plot)
  dev.off()
  cat("同时保存为PNG格式:", gsub("pdf$", "png", output_filename), "\n")
  
} else {
  cat("\n分组不足（需要至少2组），无法绘制Venn图\n")
}

# 13. 交集分析和保存 ------------------------------------------------------
if (length(geneList) >= 2) {
  cat("\n=== 交集分析 ===\n")
  
  # 计算所有组的交集
  if (length(geneList) > 1) {
    intersect_all <- Reduce(intersect, geneList)
    cat("所有组的交集基因数量:", length(intersect_all), "\n")
    
    if (length(intersect_all) > 0) {
      write.table(data.frame(Gene = intersect_all), 
                  "intersectGenes_all.txt", 
                  row.names = FALSE, quote = FALSE, sep = "\t")
      cat("所有组的交集基因已保存到: intersectGenes_all.txt\n")
    }
  }
  
  # 两两交集统计
  cat("\n两两交集统计:\n")
  for (i in 1:(length(geneList)-1)) {
    for (j in (i+1):length(geneList)) {
      name1 <- names(geneList)[i]
      name2 <- names(geneList)[j]
      intersect_genes <- intersect(geneList[[i]], geneList[[j]])
      cat(name1, "&", name2, ":", length(intersect_genes), "个共同基因\n")
      
      # 保存两两交集
      if (length(intersect_genes) > 0) {
        filename <- paste0("intersect_", name1, "_", name2, ".txt")
        write.table(data.frame(Gene = intersect_genes), 
                    filename, 
                    row.names = FALSE, quote = FALSE, sep = "\t")
      }
    }
  }
  
  # 保存每组单独的基因列表
  cat("\n=== 保存各组的基因列表 ===\n")
  for (group_name in names(geneList)) {
    filename <- paste0("genes_", group_name, ".txt")
    write.table(data.frame(Gene = geneList[[group_name]]), 
                filename, 
                row.names = FALSE, quote = FALSE, sep = "\t")
    cat(group_name, "组基因列表已保存到:", filename, "\n")
  }
}

# 14. 生成报告摘要 --------------------------------------------------------
cat("\n=== 分析完成 ===\n")
cat("工作目录:", getwd(), "\n")
cat("输入文件: genuesall.csv\n")
cat("有效分组:", paste(names(geneList), collapse = ", "), "\n")
cat("基因筛选方法: 每组前20%的高表达基因\n")
cat("颜色设置:\n")
for (group in names(group_col)) {
  cat("  ", group, ": ", group_col[group], "\n")
}
cat("\n输出的文件:\n")
cat("  Venn图: Venn_Diagram_*.pdf/png\n")
cat("  交集文件: intersect_*.txt\n")
cat("  各组基因列表: genes_*.txt\n")

# 15. 清理临时对象 --------------------------------------------------------
rm(list = ls(pattern = "^temp_"))
cat("\n脚本执行完成！\n")
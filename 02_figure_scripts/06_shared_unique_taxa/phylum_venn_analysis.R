# ============================================
# phylum.txt UTF-16LE 编码Venn图分析脚本
# ============================================

# 1. 加载包
if (!require("VennDiagram")) install.packages("VennDiagram"); library(VennDiagram)
if (!require("grid")) install.packages("grid"); library(grid)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)

# 2. 设置工作目录
setwd("/Users/catherine/Downloads/师兄交作业系列 ddl/数据提交/")

# 3. 读取数据（使用UTF-16LE编码）
cat("正在读取phylum.txt（UTF-16LE编码）...\n")
try({
  # 方法1: 先转换编码再读取
  raw_data <- readBin("phylum.txt", what = "raw", n = file.size("phylum.txt"))
  utf8_data <- iconv(list(raw_data), from = "UTF-16LE", to = "UTF-8")
  writeLines(utf8_data, "phylum_temp_utf8.txt")
  
  # 读取转换后的文件
  full_data <- read.delim("phylum_temp_utf8.txt", 
                         check.names = FALSE,
                         stringsAsFactors = FALSE)
  
  cat("✓ 数据读取成功！\n")
  cat("数据维度:", dim(full_data), "\n")
}, silent = TRUE)

if (!exists("full_data")) {
  # 方法2: 直接读取
  try({
    con <- file("phylum.txt", encoding = "UTF-16LE")
    lines <- readLines(con)
    close(con)
    
    # 写入临时文件
    writeLines(lines, "phylum_temp_utf8.txt")
    
    # 读取
    full_data <- read.delim("phylum_temp_utf8.txt", 
                           check.names = FALSE,
                           stringsAsFactors = FALSE)
    
    cat("✓ 数据读取成功！\n")
    cat("数据维度:", dim(full_data), "\n")
  }, silent = TRUE)
}

if (!exists("full_data")) {
  stop("无法读取数据，请检查文件编码")
}

# 4. 显示数据结构
cat("\n=== 数据结构 ===\n")
cat("门水平分类单元数:", nrow(full_data), "\n")
cat("样本数:", ncol(full_data) - 1, "\n")

# 显示门名称
cat("\n前20个门水平分类单元:\n")
phylum_names <- full_data[, 1]
for (i in 1:min(20, length(phylum_names))) {
  cat(i, ": ", phylum_names[i], "\n", sep = "")
}

# 显示样本名称
cat("\n样本名称:\n")
sample_names <- colnames(full_data)[-1]
print(sample_names)

# 5. 选择要分析的菌门
cat("\n=== 选择分析目标 ===\n")
cat("请从上面的列表中选择要分析的菌门\n")

# 默认选择第一个门，您可以根据需要修改
target_phylum <- phylum_names[1]
cat("默认选择:", target_phylum, "\n")
cat("如果要分析其他菌门，请修改 target_phylum 变量\n")

# 6. 查找目标菌门
target_indices <- which(phylum_names == target_phylum)
if (length(target_indices) == 0) {
  cat("未找到", target_phylum, "，使用第一个门\n")
  target_indices <- 1
  target_phylum <- phylum_names[1]
}

cat("分析目标:", target_phylum, "\n")

# 7. 定义样本分组
cat("\n=== 样本分组 ===\n")

# 您的样本分组定义
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

# 8. 匹配实际样本
cat("\n=== 匹配样本 ===\n")

match_samples <- function(target_samples, all_samples) {
  matched <- c()
  for (sample in target_samples) {
    # 创建正则表达式模式
    pattern <- gsub("([().])", "\\\\\\1", sample)  # 转义特殊字符
    pattern <- gsub("[-_]", ".", pattern)  # 将-和_替换为通配符
    
    matches <- grep(pattern, all_samples, value = TRUE, ignore.case = TRUE)
    if (length(matches) > 0) {
      matched <- c(matched, matches[1])
    }
  }
  return(unique(matched))
}

IS_actual <- match_samples(IS_samples, sample_names)
AS_actual <- match_samples(AS_samples, sample_names)
ES_actual <- match_samples(ES_samples, sample_names)
NS_actual <- match_samples(NS_samples, sample_names)

cat("IS组匹配:", length(IS_actual), "个样本\n")
cat("AS组匹配:", length(AS_actual), "个样本\n")
cat("ES组匹配:", length(ES_actual), "个样本\n")
cat("NS组匹配:", length(NS_actual), "个样本\n")

# 9. 分析每组中的目标菌门
cat("\n=== 分析各组中的", target_phylum, " ===\n")

analyze_group <- function(samples, target_idx) {
  if (length(samples) == 0) return(FALSE)
  
  # 提取数据
  group_data <- full_data[target_idx, samples, drop = FALSE]
  
  # 转换为数值
  numeric_values <- sapply(group_data, function(x) as.numeric(as.character(x)))
  
  # 判断是否存在（> 0）
  any(numeric_values > 0, na.rm = TRUE)
}

# 检查目标菌门在各组中的存在情况
is_present <- list(
  IS = if (length(IS_actual) > 0) analyze_group(IS_actual, target_indices) else FALSE,
  AS = if (length(AS_actual) > 0) analyze_group(AS_actual, target_indices) else FALSE,
  ES = if (length(ES_actual) > 0) analyze_group(ES_actual, target_indices) else FALSE,
  NS = if (length(NS_actual) > 0) analyze_group(NS_actual, target_indices) else FALSE
)

cat(target_phylum, "在各组中的存在情况:\n")
for (group in names(is_present)) {
  cat(group, ":", ifelse(is_present[[group]], "存在", "不存在"), "\n")
}

# 10. 创建Venn图数据
cat("\n=== 创建Venn图数据 ===\n")

# 创建各组列表
group_list <- list()
if (is_present$IS) group_list$IS <- target_phylum
if (is_present$AS) group_list$AS <- target_phylum
if (is_present$ES) group_list$ES <- target_phylum
if (is_present$NS) group_list$NS <- target_phylum

cat("有效分组数:", length(group_list), "\n")

# 11. 绘制Venn图
if (length(group_list) >= 2) {
  cat("\n=== 绘制Venn图 ===\n")
  
  # 设置颜色
  group_col <- c(IS = "#00FF7F", AS = "#FF1493", ES = "#FFA500", NS = "#00FFFF")
  group_col <- group_col[names(group_list)]
  
  # 清理文件名
  clean_name <- gsub("[^[:alnum:]]", "_", target_phylum)
  output_pdf <- paste0("Phylum_", clean_name, "_Venn.pdf")
  
  # 根据分组数量绘制
  if (length(group_list) == 2) {
    venn.plot <- venn.diagram(
      x = group_list,
      filename = NULL,
      main = paste("Phylum:", target_phylum),
      fill = group_col,
      alpha = 0.6,
      category.names = names(group_list),
      cat.cex = 1.2,
      cex = 1.5,
      lwd = 2,
      margin = 0.05
    )
  } else if (length(group_list) == 3) {
    venn.plot <- venn.diagram(
      x = group_list,
      filename = NULL,
      main = paste("Phylum:", target_phylum),
      fill = group_col,
      alpha = 0.6,
      category.names = names(group_list),
      cat.cex = 1.1,
      cat.dist = c(0.06, 0.06, 0.04),
      cex = 1.3,
      lwd = 2,
      margin = 0.08
    )
  } else if (length(group_list) == 4) {
    venn.plot <- venn.diagram(
      x = group_list,
      filename = NULL,
      main = paste("Phylum:", target_phylum),
      fill = group_col,
      alpha = 0.6,
      category.names = names(group_list),
      cat.cex = 1.0,
      cat.dist = c(0.28, 0.28, 0.18, 0.18),
      cat.pos = c(0, 0, 180, 180),
      cex = 1.2,
      lwd = 2,
      margin = 0.1
    )
  }
  
  # 保存图形
  pdf(output_pdf, width = 10, height = 10)
  grid.draw(venn.plot)
  dev.off()
  
  cat("Venn图已保存为:", output_pdf, "\n")
  
  # 保存PNG
  png(gsub("pdf$", "png", output_pdf), width = 1000, height = 1000, res = 150)
  grid.draw(venn.plot)
  dev.off()
  cat("PNG格式已保存:", gsub("pdf$", "png", output_pdf), "\n")
  
} else {
  cat("\n分组不足（需要至少2组），无法绘制Venn图\n")
  cat("建议分析多个菌门\n")
}

# 12. 清理临时文件
if (file.exists("phylum_temp_utf8.txt")) {
  file.remove("phylum_temp_utf8.txt")
}

cat("\n=== 分析完成 ===\n")
cat("请查看生成的文件:\n")
cat("1. 分析报告: phylum_analysis_results.txt\n")
cat("2. Venn图: Phylum_*_Venn.pdf/png\n")
cat("3. 转换后的文件: phylum_utf8_fixed.txt\n")


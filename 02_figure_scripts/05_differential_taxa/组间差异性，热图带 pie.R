# 清理工作空间
rm(list = ls())
gc()

# 设置工作目录到桌面
setwd("~/Desktop")
cat("工作目录设置为:", getwd(), "\n")

# 设置系统语言为英文避免编码问题
Sys.setlocale("LC_ALL", "en_US.UTF-8")

# 加载必要的包
if (!require("corrgram")) install.packages("corrgram")
if (!require("viridis")) install.packages("viridis")

library(corrgram)
library(viridis)

cat("Packages loaded successfully!\n")

# 读取CSV数据
file_path <- "/Users/catherine/Downloads/董西洋铁-砷-汞等金属循环文献-冷泉/角鲨烯基于 contig 寻找的功能基因/squalene.csv"

cat("Reading file:", file_path, "\n")

# 读取数据
data <- read.csv(file_path, 
                 header = TRUE, 
                 stringsAsFactors = FALSE,
                 fileEncoding = "UTF-8",
                 check.names = FALSE)

cat("Data read successfully!\n")
cat("数据维度 - 行:", nrow(data), "列:", ncol(data), "\n")

# 提取基因名（第一列）
gene_names <- data[, 1]

# 提取P值数据（排除第一列）
pvalue_data <- data[, -1, drop = FALSE]

# 只需要前8列（因为这是一个8x8的对称矩阵）
# 从你之前的数据片段看，这是一个8x8的P值矩阵
pvalue_data <- pvalue_data[, 1:8]

# 设置行名和列名
rownames(pvalue_data) <- gene_names
colnames(pvalue_data) <- gene_names

# 转换为矩阵
pvalue_matrix <- as.matrix(pvalue_data)

cat("\n=== 处理后的数据 ===\n")
cat("P值矩阵维度:", dim(pvalue_matrix), "\n")
cat("是否为方阵:", nrow(pvalue_matrix) == ncol(pvalue_matrix), "\n")
cat("基因列表:\n")
print(gene_names)

# 显示P值矩阵
cat("\nP值矩阵:\n")
print(round(pvalue_matrix, 4))

# 创建相关性矩阵（基于P值）
create_correlation_matrix <- function(p_matrix) {
  n <- nrow(p_matrix)
  
  # 初始化相关性矩阵（对角线为1）
  cor_mat <- matrix(1, nrow = n, ncol = n)
  
  # 设置行名和列名
  gene_names <- rownames(p_matrix)
  rownames(cor_mat) <- gene_names
  colnames(cor_mat) <- gene_names
  
  # 根据P值计算相关性
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        p_val <- p_matrix[i, j]
        
        if (!is.na(p_val) && is.numeric(p_val)) {
          # 转换P值为相关性强度
          # P值越小，相关性越强
          if (p_val < 0.001) {
            cor_mat[i, j] <- 0.95  # 强相关
          } else if (p_val < 0.01) {
            cor_mat[i, j] <- 0.80  # 中等相关
          } else if (p_val < 0.05) {
            cor_mat[i, j] <- 0.65  # 弱相关
          } else {
            cor_mat[i, j] <- 0.10  # 不相关
          }
        }
      }
    }
  }
  
  return(cor_mat)
}

# 创建相关性矩阵
cor_matrix <- create_correlation_matrix(pvalue_matrix)

cat("\n=== 相关性矩阵 ===\n")
cat("维度:", dim(cor_matrix), "\n")
cat("值范围:", round(range(cor_matrix), 3), "\n")

# 显示相关性矩阵
cat("\n相关性矩阵:\n")
print(round(cor_matrix, 3))

# 创建简短的基因标签（用于图形显示）
short_labels <- gene_names
for (i in seq_along(gene_names)) {
  if (nchar(gene_names[i]) > 15) {
    short_labels[i] <- substr(gene_names[i], 1, 12)
  }
}

cat("\n使用的基因标签:\n")
for (i in 1:length(gene_names)) {
  cat(i, ".", gene_names[i], "->", short_labels[i], "\n")
}

# 方法1：使用自定义对角面板函数
cat("\n=== 方法1：使用自定义对角面板 ===\n")

# 自定义对角面板函数
custom_diag_panel <- function(x = 0, varname = "", ...) {
  # 获取当前面板的索引
  panel_index <- sys.nframe() - 4
  
  if (panel_index > 0 && panel_index <= length(short_labels)) {
    current_gene <- short_labels[panel_index]
  } else {
    current_gene <- varname
  }
  
  # 清空面板
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = "white", border = NA)
  
  # 显示基因名
  text(0.5, 0.5, 
       labels = current_gene, 
       cex = 0.8,
       font = 2,
       col = "black")
}

# 绘制PDF
tryCatch({
  pdf("Gene_Correlation_Corrgram.pdf", width = 10, height = 9)
  
  par(mar = c(3, 3, 5, 3), family = "sans", bg = "white")
  
  # 绘制corrgram
  corrgram(cor_matrix,
           order = FALSE,
           lower.panel = panel.shade,
           upper.panel = panel.pie,
           diag.panel = custom_diag_panel,
           main = "Gene Correlation Analysis\n(Based on P-values)",
           col.regions = colorRampPalette(viridis(10)),
           cex.main = 1.4,
           gap = 0.2)
  
  dev.off()
  cat("✓ Gene_Correlation_Corrgram.pdf 创建成功\n")
}, error = function(e) {
  cat("✗ 方法1失败:", e$message, "\n")
})

# 方法2：使用最简单的corrgram（不使用对角面板）
cat("\n=== 方法2：简单corrgram（无对角面板）===\n")

tryCatch({
  pdf("Gene_Correlation_Simple.pdf", width = 9, height = 8)
  
  par(mar = c(3, 3, 5, 3), family = "sans", bg = "white")
  
  # 使用NULL作为对角面板
  corrgram(cor_matrix,
           order = FALSE,
           lower.panel = panel.shade,
           upper.panel = panel.shade,
           diag.panel = NULL,
           main = "Gene Correlation Heatmap",
           col.regions = colorRampPalette(c("blue", "white", "red")),
           cex.main = 1.3,
           gap = 0.15)
  
  dev.off()
  cat("✓ Gene_Correlation_Simple.pdf 创建成功\n")
}, error = function(e) {
  cat("✗ 方法2失败:", e$message, "\n")
})

# 方法3：使用panel.minmax作为对角面板
cat("\n=== 方法3：使用panel.minmax ===\n")

tryCatch({
  pdf("Gene_Correlation_MinMax.pdf", width = 10, height = 9)
  
  par(mar = c(3, 3, 5, 3), family = "sans", bg = "white")
  
  corrgram(cor_matrix,
           order = FALSE,
           lower.panel = panel.shade,
           upper.panel = panel.cor,
           diag.panel = panel.minmax,
           main = "Gene Correlation with Min/Max Values",
           col.regions = colorRampPalette(viridis(10)),
           cex.main = 1.4,
           gap = 0.2)
  
  dev.off()
  cat("✓ Gene_Correlation_MinMax.pdf 创建成功\n")
}, error = function(e) {
  cat("✗ 方法3失败:", e$message, "\n")
})

# 方法4：使用热图作为替代
cat("\n=== 方法4：使用热图替代 ===\n")

if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)

tryCatch({
  pdf("Gene_Correlation_Heatmap.pdf", width = 10, height = 9)
  
  # 创建更美观的颜色方案
  heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # 绘制热图
  pheatmap(cor_matrix,
           color = heatmap_colors,
           main = "Gene Correlation Heatmap\n(Based on P-values)",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_row = 10,
           fontsize_col = 10,
           fontsize_number = 8)
  
  dev.off()
  cat("✓ Gene_Correlation_Heatmap.pdf 创建成功\n")
}, error = function(e) {
  cat("✗ 方法4失败:", e$message, "\n")
})

# 创建数据表格（与之前相同）
cat("\n=== 创建数据表格 ===\n")

# 创建基因对表格
n_genes <- length(gene_names)
gene_pairs <- data.frame()

for (i in 1:(n_genes-1)) {
  for (j in (i+1):n_genes) {
    gene1 <- gene_names[i]
    gene2 <- gene_names[j]
    p_val <- pvalue_matrix[i, j]
    
    if (!is.na(p_val) && is.numeric(p_val)) {
      significance <- ""
      
      if (p_val < 0.001) {
        significance <- "***"
      } else if (p_val < 0.01) {
        significance <- "**"
      } else if (p_val < 0.05) {
        significance <- "*"
      }
      
      gene_pairs <- rbind(gene_pairs, data.frame(
        Gene1 = gene1,
        Gene2 = gene2,
        P_value = p_val,
        Significance = significance,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 保存基因对表格
if (nrow(gene_pairs) > 0) {
  gene_pairs <- gene_pairs[order(gene_pairs$P_value), ]
  write.csv(gene_pairs, "Gene_Correlation_Pairs.csv", row.names = FALSE)
  cat("✓ Gene_Correlation_Pairs.csv 保存成功 (", nrow(gene_pairs), " 对基因)\n", sep = "")
}

# 创建汇总表格
summary_df <- data.frame(
  Gene = gene_names,
  Significant_Pairs = 0,
  Min_P_value = NA,
  stringsAsFactors = FALSE
)

for (i in 1:n_genes) {
  p_vals <- pvalue_matrix[i, -i]
  p_vals <- p_vals[!is.na(p_vals)]
  
  if (length(p_vals) > 0) {
    sig_count <- sum(p_vals < 0.05)
    min_p <- min(p_vals)
    
    summary_df$Significant_Pairs[i] <- sig_count
    summary_df$Min_P_value[i] <- min_p
  }
}

write.csv(summary_df, "Gene_Correlation_Summary.csv", row.names = FALSE)
cat("✓ Gene_Correlation_Summary.csv 保存成功\n")

# 创建图例说明
cat("\n=== 创建图例说明 ===\n")

legend_content <- paste(
  "=================================================",
  "         基因相关性分析结果说明",
  "=================================================",
  "",
  paste("生成时间: ", Sys.time()),
  paste("工作目录: ", getwd()),
  "",
  "生成的文件:",
  "1. Gene_Correlation_Corrgram.pdf - 标准corrgram图",
  "2. Gene_Correlation_Simple.pdf   - 简化版热图",
  "3. Gene_Correlation_MinMax.pdf   - 带最小最大值显示",
  "4. Gene_Correlation_Heatmap.pdf  - 传统热图",
  "5. Gene_Correlation_Pairs.csv    - 基因对P值数据",
  "6. Gene_Correlation_Summary.csv  - 基因汇总统计",
  "",
  "图形解读:",
  "- 对角线: 基因名称",
  "- 下三角区域: 颜色阴影表示相关性强度",
  "  * 颜色越深表示相关性越强",
  "  * 红色表示正相关，蓝色表示负相关",
  "- 上三角区域: 饼图表示相关性强度",
  "  * 填充比例越高表示相关性越强",
  "",
  "P值显著性:",
  "- *** : P < 0.001 (极显著)",
  "- **  : P < 0.01  (非常显著)",
  "- *   : P < 0.05  (显著)",
  "- 无星号: P >= 0.05 (不显著)",
  "",
  "分析的基因 (共8个):",
  paste(paste0(1:length(gene_names), ". ", gene_names), collapse = "\n"),
  "",
  "注意:",
  "1. 相关性值是基于P值计算的",
  "2. P值越小，计算出的相关性值越大",
  "3. 对角线上的值始终为1（基因与自身的相关性）",
  "",
  "=================================================",
  sep = "\n"
)

writeLines(legend_content, "README_说明.txt", useBytes = TRUE)
cat("✓ README_说明.txt 保存成功\n")

# 显示结果
cat("\n=================================================\n")
cat("           分析完成！\n")
cat("=================================================\n")
cat("生成的文件:\n")
cat("-------------------------------------------------\n")

files_to_check <- c("Gene_Correlation_Corrgram.pdf",
                    "Gene_Correlation_Simple.pdf",
                    "Gene_Correlation_MinMax.pdf",
                    "Gene_Correlation_Heatmap.pdf",
                    "Gene_Correlation_Pairs.csv",
                    "Gene_Correlation_Summary.csv",
                    "README_说明.txt")

for (file in files_to_check) {
  if (file.exists(file)) {
    size_kb <- round(file.info(file)$size / 1024, 1)
    cat(sprintf("✓ %-35s %5.1f KB\n", file, size_kb))
  } else {
    cat(sprintf("✗ %-35s (未找到)\n", file))
  }
}

cat("=================================================\n")
cat("\n重要提示:\n")
cat("1. 请查看桌面上的PDF文件\n")
cat("2. Gene_Correlation_Corrgram.pdf 是主要的corrgram图\n")
cat("3. 如果corrgram有问题，Gene_Correlation_Heatmap.pdf 是很好的替代\n")
cat("4. 详细数据在CSV文件中\n")
cat("\n文件位置: ", getwd(), "\n")

# 尝试打开第一个PDF
if (file.exists("Gene_Correlation_Corrgram.pdf")) {
  cat("\n正在打开PDF文件...\n")
  if (Sys.info()["sysname"] == "Darwin") {
    system("open 'Gene_Correlation_Corrgram.pdf'")
  }
}

cat("\n分析完成！\n")
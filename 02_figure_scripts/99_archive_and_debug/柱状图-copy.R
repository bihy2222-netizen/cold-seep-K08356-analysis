# 简单直接的版本
input_file <- "/Users/catherine/Downloads/师兄交作业系列 ddl/species.csv"
output_file <- "/Users/catherine/Downloads/师兄交作业系列 ddl/species_upset_binary.csv"

# 读取数据
cat("正在读取数据...\n")
data <- read.csv(input_file, row.names = 1, check.names = FALSE)

# 显示数据信息
cat("数据维度:", dim(data), "\n")
cat("数据类型示例:\n")
print(str(data[1:3, 25:28]))

# 直接转换：大于0的值变为1
cat("\n正在转换数据...\n")
binary_data <- data
binary_data[binary_data > 0] <- 1

# 验证转换
cat("\n转换后数据示例:\n")
print(binary_data[1:5, 25:28])

# 统计
total <- nrow(binary_data) * ncol(binary_data)
ones <- sum(binary_data == 1, na.rm = TRUE)
zeros <- sum(binary_data == 0, na.rm = TRUE)

cat("\n统计信息:\n")
cat("总数据点:", total, "\n")
cat("1的数量:", ones, sprintf("(%.2f%%)", ones/total*100), "\n")
cat("0的数量:", zeros, sprintf("(%.2f%%)", zeros/total*100), "\n")

# 保存
write.csv(binary_data, output_file, row.names = TRUE)
cat("\n文件已保存:", output_file, "\n")
# 最小版本
setwd("/Users/catherine/Downloads/metaphlan 作图/class 级别作图-柱状图/")

# 清理文件
raw <- readBin("class34.csv", raw(), file.info("class34.csv")$size)
clean <- gsub("\x00", "", rawToChar(raw))
writeLines(strsplit(clean, "\n")[[1]], "clean.csv")

# 读取数据
d <- read.csv("clean.csv", check.names = FALSE, stringsAsFactors = FALSE)
cat("数据：", nrow(d), "行，", ncol(d), "列\n")

# 创建二进制矩阵
binary <- data.frame(Taxon = d[, 1])
for(i in 2:ncol(d)) binary[[colnames(d)[i]]] <- ifelse(d[, i] > 0, 1, 0)

# 绘图
install.packages("UpSetR")
library(UpSetR)

lst <- list()
for(i in 2:ncol(binary)) {
  taxa <- binary$Taxon[binary[, i] == 1]
  if(length(taxa) > 0) lst[[colnames(binary)[i]]] <- as.character(taxa)
}

pdf("Upset.pdf", 12, 8)
upset(fromList(lst), order.by = "freq")
dev.off()
cat("完成！图在Upset.pdf\n")
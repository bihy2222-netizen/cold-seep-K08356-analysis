# phylum.txt 编码修复脚本
setwd("/Users/catherine/Downloads/师兄交作业系列 ddl/数据提交/")

# 方法1: 尝试自动转换
try({
  # 读取原始数据
  raw_data <- readBin("phylum.txt", what = "raw", n = file.size("phylum.txt"))
  
  # 尝试多种编码
  encodings <- c("UTF-8", "GBK", "GB18030", "latin1", "CP936")
  
  for (enc in encodings) {
    try({
      text_data <- iconv(raw_data, from = enc, to = "UTF-8")
      writeLines(text_data, "phylum_utf8.txt")
      cat("成功使用", enc, "编码转换，保存为 phylum_utf8.txt\n")
      break
    }, silent = TRUE)
  }
}, silent = TRUE)

# 方法2: 使用 system 命令（在终端中执行）
cat("\n也可以在终端中执行以下命令:\n")
cat("iconv -f GBK -t UTF-8 phylum.txt > phylum_utf8.txt\n")
cat("或:\n")
cat("nkf -w --overwrite phylum.txt  # 如果安装了nkf工具\n")


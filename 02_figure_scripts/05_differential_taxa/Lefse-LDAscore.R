rm(list=ls())
#调用R包
library(tidyverse)
library(microeco)
library(magrittr)
#读取数据
otu <- read.table("arc_OTU.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
group <- read.table("arc_group.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
tax <- read.table("arc_TAX.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
#创建microeco包可识别的整合对象
dataset <- microtable$new(sample_table = group,
                          otu_table = otu,
                          tax_table = tax)


#开始LEfse分析
lefse <- trans_diff$new(dataset = dataset,
                        method = "lefse",
                        group = "Group",
                        alpha = 0.1,
                        p_adjust_method = "none",  # 同时关闭 p 值校正
                        lefse_subgroup = NULL)


# 查看分析结果
head(lefse$res_diff)


# 生成图并赋值
p <- lefse$plot_diff_bar(use_number = 1:30,
                         width = 0.8,
                         group_order = c("C12", "C13"))

# 保存为 SVG 文件
ggsave("arc_LEfSe_0.1_top30_diff_bar.svg", plot = p, width = 8, height = 6)

print(p)
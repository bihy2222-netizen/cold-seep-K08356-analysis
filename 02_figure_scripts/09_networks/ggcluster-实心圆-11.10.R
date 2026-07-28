# 清理工作空间
rm(list = ls())

# 载入必要R包
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(ggClusterNet)
library(Biostrings)

# ---读取数据---
# OTU丰度表（行为OTU，列为样品）
otutab = read.csv("./pathway-tax-OTUS.csv", row.names = 1)

# 分类注释表（行为OTU，列为门、纲、目等）
taxonomy = read.csv("./pathway-tax-taxoS.csv", row.names = 1)

# 分组信息（每个OTU的模块）
netClu = read.csv("pathway-tax-groupS.csv", header = TRUE, stringsAsFactors = TRUE)
netClu$ID <- as.character(netClu$ID)
netClu$group <- as.factor(netClu$group)

# 构建phyloseq对象
ps = phyloseq(
  otu_table(as.matrix(otutab), taxa_are_rows = TRUE),
  tax_table(as.matrix(taxonomy))
)

# ---计算相关矩阵---
result = corMicro(
  ps = ps,
  N = 0,                                # 使用所有OTU
  method.scale = "TMM",                 # 标准化方法
  r.threshold = 0.6,
  p.threshold = 0.05,
  method = "spearman"
)

# ---提取相关矩阵---
cor = result[[1]]
ps_net = result[[3]]

# ---OTU表提取并转置---
otu_table = ps_net %>%
  vegan_otu() %>%
  t() %>%
  as.data.frame()

# ---构建节点坐标---
result2 = model_filled_circle(
  cor = cor,
  culxy = TRUE,
  da = NULL,
  nodeGroup = netClu,
  mi.size = 2,
  zoom = 0.6
)

node = result2[[1]]
rownames(node) <- node$elements  # 设置行名为元素名，便于匹配

# ---注释节点信息---
tax_mat = as.matrix(taxonomy)  # 避免与函数 tax_table 冲突
nodes = nodeadd(
  plotcord = node,
  otu_table = otu_table,
  tax_table = tax_mat
)

# ---构建边---
edge = edgeBuild(
  cor = cor,
  node = node
)

pnet <- ggplot() +
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = as.factor(cor)),
               data = edge, linewidth = 0.3, alpha = 0.6) +
  geom_point(aes(X1, X2, fill = Phylum, size = mean), pch = 21, data = nodes) +
  scale_colour_brewer(palette = "Set1") +
  scale_size_continuous(
    name = "Relative Abundance",
    range = c(2, 20),
    breaks = c(1, 10, 100, 1000, 10000),
    labels = c("1", "10", "100", "1000", "10000")
    #limits = c(1, 10000)  # 也可以用 range(nodes$mean)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5), title = "Tax"),
    size = guide_legend(title = "Relative Abundance")
  ) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme_void() +
  theme(
    legend.background = element_rect(colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

# ---保存图形---
ggsave("mags_network_plotS.pdf", pnet, width = 12, height = 9)
ggsave("mags_network_plotS.svg", pnet, width = 12, height = 9, device = "svg")

# ---展示图形---
print(pnet)



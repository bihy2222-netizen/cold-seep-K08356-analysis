# ============================================================
# LEfSe 分析：仅绘制分支图（Cladogram），保留 MetaPhlAn 前缀
# 基于 MetaPhlAn 合并丰度表（threshold = 2.0）
# 请重启 R 会话后运行本脚本
# ============================================================

# 0. 生成样本信息文件 ------------------------------------------
abund_file <- "/Users/catherine/Downloads/mgmksz/merged_abundance_table.txt"

con <- file(abund_file, "r")
tmp_lines <- readLines(con, n = 2)
close(con)

group_line <- tmp_lines[1]
header_line <- tmp_lines[2]
header <- strsplit(header_line, "\t")[[1]]
header <- trimws(header)
samples <- header[-1]
samples_clean <- gsub("_metaphlan$", "", samples)

group_header <- trimws(strsplit(group_line, "\t")[[1]])
file_groups <- group_header[-1]
file_groups <- trimws(file_groups)

sample_info <- data.frame(
  Sample = samples_clean,
  stringsAsFactors = FALSE
)
sample_info$Site <- gsub("_.*", "", samples_clean)
sample_info$Type <- ifelse(grepl("_[0-9]+_[0-9]+", samples_clean), "Sediment", "Seawater")
if (length(file_groups) == length(samples_clean) &&
    tolower(trimws(group_header[1])) %in% c("group", "class", "type")) {
  sample_info$Group <- file_groups
} else {
  sample_info$Group <- paste(sample_info$Site, sample_info$Type, sep = "_")
}
rownames(sample_info) <- sample_info$Sample

write.csv(sample_info, "sample_info.csv", row.names = TRUE)
print("样本信息文件已生成，前几行：")
print(head(sample_info))

# 1. 加载包 ----------------------------------------------------
local_lib <- file.path(getwd(), "R_libs")
if (!dir.exists(local_lib)) dir.create(local_lib, recursive = TRUE)
.libPaths(c(local_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("R6", quietly = TRUE)) install.packages("R6")
if (!requireNamespace("microeco", quietly = TRUE)) install.packages("microeco")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("magrittr", quietly = TRUE)) install.packages("magrittr")
if (!requireNamespace("tidytree", quietly = TRUE)) install.packages("tidytree")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ggtree", quietly = TRUE)) BiocManager::install("ggtree", ask = FALSE, update = FALSE)

library(R6)
library(microeco)
library(ggplot2)
library(magrittr)
library(tidytree)
library(ggtree)

# 2. 数据准备 -------------------------------------------------
group_file <- "sample_info.csv"

abund_raw <- read.table(abund_file, sep = "\t", header = TRUE, 
                        row.names = 1, check.names = FALSE, 
                        skip = 1, comment.char = "", quote = "")

colnames(abund_raw) <- gsub("_metaphlan$", "", colnames(abund_raw))

# 确保 OTU 表为数值数据框
otu_mat <- as.data.frame(abund_raw)
otu_mat[] <- lapply(otu_mat, as.numeric)

sample_info <- read.csv(group_file, row.names = 1, check.names = FALSE)

common_samples <- intersect(colnames(otu_mat), rownames(sample_info))
if (length(common_samples) == 0) {
  stop("没有匹配的样本！请检查丰度表列名与样本信息文件行名。")
}
otu_mat <- otu_mat[, common_samples, drop = FALSE]
sample_info <- sample_info[common_samples, , drop = FALSE]

# 3. 解析分类表（保留 MetaPhlAn 前缀，形成标准7列）-------------
parse_metaphlan_with_prefix <- function(tax_strings) {
  res_list <- lapply(tax_strings, function(tax_str) {
    levels <- strsplit(tax_str, "\\|")[[1]]
    out <- setNames(rep(NA_character_, 7), 
                    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
    for (lv in levels) {
      if (grepl("^k__", lv)) out["Kingdom"] <- lv
      if (grepl("^p__", lv)) out["Phylum"]  <- lv
      if (grepl("^c__", lv)) out["Class"]   <- lv
      if (grepl("^o__", lv)) out["Order"]   <- lv
      if (grepl("^f__", lv)) out["Family"]  <- lv
      if (grepl("^g__", lv)) out["Genus"]   <- lv
      if (grepl("^s__", lv)) out["Species"] <- lv
    }
    return(out)
  })
  tax_df <- as.data.frame(do.call(rbind, res_list), stringsAsFactors = FALSE)
  rownames(tax_df) <- tax_strings
  
  # 清理可能残留的 "<NA>" 字符串
  tax_df[] <- lapply(tax_df, function(col) {
    col[col == "<NA>"] <- NA_character_
    col
  })
  return(tax_df)
}

tax_table <- parse_metaphlan_with_prefix(rownames(otu_mat))
print("分类表解析完成，前几行：")
print(head(tax_table))

# 4. 构造 microtable 对象 -------------------------------------
dataset <- tryCatch({
  microtable$new(
    otu_table = otu_mat,
    sample_table = sample_info,
    tax_table = tax_table
  )
}, error = function(e) {
  stop("创建 microtable 对象失败: ", e$message)
})

# 注意：由于分类表已包含标准7列并带有前缀，无需再执行 tidy_taxonomy()
# 若仍想检查，可取消下一行注释，但可能引发错误
# dataset$tidy_taxonomy()

print("microtable 对象创建成功：")
print(dataset)

# 5. LEfSe 分析 -----------------------------------------------
available_ranks <- colnames(dataset$tax_table)
target_rank <- if ("Genus" %in% available_ranks) "Genus" else available_ranks[length(available_ranks)]
message("使用分类等级: ", target_rank)

lefse <- tryCatch({
  trans_diff$new(
    dataset = dataset,
    method = "lefse",
    group = "Group",
    alpha = 0.05,
    lefse_subgroup = NULL,
    taxa_level = "all",
    filter_thres = 0.0001,
    p_adjust_method = "none"
  )
}, error = function(e) {
  stop("LEfSe 分析失败: ", e$message)
})

write.csv(lefse$res_diff, "lefse_res_diff.csv", row.names = TRUE)
print("差异分析完成，结果已保存。")

# 6. 绘制 Nature 风格分支图 ------------------------------------
group_order <- unique(sample_info$Group)
n_groups <- length(group_order)

# 色盲友好、高对比、低饱和度：适合期刊印刷和屏幕阅读
nature_palette <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9", "#F0E442", "#000000"
)
if (n_groups <= length(nature_palette)) {
  my_colors <- nature_palette[seq_len(n_groups)]
} else {
  my_colors <- grDevices::hcl.colors(n_groups, palette = "Dark 3")
}
names(my_colors) <- group_order

nature_theme <- theme_void(base_family = "Helvetica") +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(8, 10, 8, 8, unit = "mm"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold", colour = "grey10"),
    legend.text = element_text(size = 9, colour = "grey15"),
    legend.key.size = unit(4.2, "mm"),
    legend.spacing.y = unit(1.2, "mm"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, colour = "grey10"),
    plot.subtitle = element_text(size = 9.5, hjust = 0.5, colour = "grey35",
                                 margin = margin(b = 4, unit = "mm"))
  )

p_clado <- lefse$plot_diff_cladogram(
  use_taxa_num = 200,
  use_feature_num = 40,
  clade_label_level = 5,
  group_order = group_order,
  color = unname(my_colors),
  branch_size = 0.38,
  alpha = 0.12,
  clade_label_size = 2.7
) +
  nature_theme +
  labs(
    title = "LEfSe cladogram",
    subtitle = "Discriminative microbial lineages among sampling groups",
    colour = "Enriched group",
    fill = "Enriched group"
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 4, alpha = 1)),
    fill = guide_legend(override.aes = list(alpha = 1))
  )

nature_pdf <- function(filename, ...) {
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(filename, family = "Helvetica", ...)
  } else {
    grDevices::pdf(filename, family = "Helvetica", useDingbats = FALSE, ...)
  }
}

ggsave("LDA_cladogram_Nature_200.pdf", p_clado,
       width = 180, height = 150, units = "mm", device = nature_pdf)
ggsave("LDA_cladogram_Nature_200.png", p_clado,
       width = 180, height = 150, units = "mm", dpi = 600, bg = "white")

print("Nature风格分支图绘制完成：LDA_cladogram_Nature_200.pdf / .png")



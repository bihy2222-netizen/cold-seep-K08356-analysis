# ==================== 10. 额外的统计分析和可视化 ====================

# 10.1 查看每个基因在不同样本中的分布
cat("\n")
cat(rep("-", 60), "\n", sep = "")
cat("各基因在不同样本中的平均丰度:\n")
cat(rep("-", 60), "\n", sep = "")

gene_avg_abundance <- data.frame(
  Gene_ID = rownames(abundance_matrix),
  Average_Abundance = rowMeans(abundance_matrix, na.rm = TRUE)
) %>%
  arrange(desc(Average_Abundance))

print(gene_avg_abundance)

# ==================== 重要发现：基因ID匹配问题 ====================

cat("\n")
cat(rep("=", 60), "\n", sep = "")
cat("重要发现：基因ID匹配问题\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("问题分析：\n")
cat("1. 所有基因的门分类都是'Unknown'\n")
cat("2. 这意味着基因ID与KO-MAG数据中的KO_ID不匹配\n")
cat("3. 可能原因：基因丰度数据中的基因ID是基因名称，而不是KO编号\n\n")

cat("基因丰度数据中的基因ID:\n")
print(analyzed_genes)

cat("\nKO-MAG数据中的KO_ID示例:\n")
print(head(unique(ko_mag_data$KO_ID), 10))

# ==================== 修复：尝试不同的匹配方式 ====================

cat("\n")
cat(rep("=", 60), "\n", sep = "")
cat("尝试修复：使用基因名称匹配\n")
cat(rep("=", 60), "\n\n", sep = "")

# 10.2 重新尝试为基因找到门分类信息
# 假设基因丰度数据中的"Gene_ID"实际上是基因名称而不是KO编号

# 创建一个基因名称到KO编号的映射（如果可能）
gene_name_mapping <- list(
  "Squalene_synthase" = "K00801",  # 这是示例，需要根据实际情况调整
  "Squalene_monooxygenase" = "K00498",
  "FPP_synthase_euk" = "K00791",
  "FPP_synthase_prok" = "K00787",
  "OSC_lanosterol_cyclase" = "K01852",
  "HMG-CoA_synthase" = "K01641",
  "HMG-CoA_reductase" = "K00021",
  "IspG_MEP_pathway" = "K03526"
)

# 尝试使用基因名称在描述中查找
cat("尝试在描述信息中查找基因名称...\n")
for (gene_name in analyzed_genes) {
  # 在描述中搜索基因名称
  matching_descriptions <- ko_mag_data %>%
    filter(grepl(gene_name, Description, ignore.case = TRUE)) %>%
    distinct(KO_ID, Description)
  
  if (nrow(matching_descriptions) > 0) {
    cat("\n基因:", gene_name, "\n")
    cat("找到的匹配项:\n")
    print(matching_descriptions)
    
    # 使用第一个匹配的KO_ID
    suggested_ko <- matching_descriptions$KO_ID[1]
    gene_name_mapping[[gene_name]] <- suggested_ko
    cat("建议的KO编号:", suggested_ko, "\n")
  } else {
    cat("基因", gene_name, ": 未在描述中找到匹配项\n")
  }
}

# 10.3 重新创建基因注释信息
cat("\n重新创建基因注释信息...\n")

# 创建新的基因注释数据框
gene_annotation_new <- data.frame(
  Gene_Name = analyzed_genes,
  Suggested_KO_ID = sapply(analyzed_genes, function(x) {
    if (!is.null(gene_name_mapping[[x]]) && gene_name_mapping[[x]] != "") {
      return(gene_name_mapping[[x]])
    } else {
      return("Unknown")
    }
  }),
  stringsAsFactors = FALSE
)

# 使用建议的KO_ID查找门分类
gene_phylum_info_new <- data.frame(
  Gene_Name = analyzed_genes,
  Phylum = "Unknown",
  stringsAsFactors = FALSE
)

for (i in 1:length(analyzed_genes)) {
  gene_name <- analyzed_genes[i]
  suggested_ko <- gene_annotation_new$Suggested_KO_ID[i]
  
  if (suggested_ko != "Unknown") {
    # 使用建议的KO_ID查找门分类
    gene_records <- ko_mag_tax %>%
      filter(KO_ID == suggested_ko) %>%
      filter(!is.na(Phylum)) %>%
      filter(Phylum != "Unknown")
    
    if (nrow(gene_records) > 0) {
      # 找到最频繁的门分类
      phylum_counts <- table(gene_records$Phylum)
      most_common_phylum <- names(phylum_counts)[which.max(phylum_counts)]
      gene_phylum_info_new$Phylum[i] <- most_common_phylum
      gene_phylum_info_new$KO_ID[i] <- suggested_ko
      
      cat(gene_name, " (KO:", suggested_ko, ") 的主要门分类:", most_common_phylum, 
          "(出现在", max(phylum_counts), "个MAG中)\n")
    } else {
      # 如果找不到，尝试在描述中搜索
      desc_matches <- ko_mag_data %>%
        filter(grepl(gene_name, Description, ignore.case = TRUE))
      
      if (nrow(desc_matches) > 0) {
        # 获取所有匹配记录的门分类
        all_phyla <- ko_mag_tax %>%
          filter(KO_ID %in% desc_matches$KO_ID) %>%
          filter(!is.na(Phylum)) %>%
          filter(Phylum != "Unknown")
        
        if (nrow(all_phyla) > 0) {
          phylum_counts <- table(all_phyla$Phylum)
          most_common_phylum <- names(phylum_counts)[which.max(phylum_counts)]
          gene_phylum_info_new$Phylum[i] <- most_common_phylum
          gene_phylum_info_new$KO_ID[i] <- desc_matches$KO_ID[1]
          
          cat(gene_name, " (通过描述匹配) 的主要门分类:", most_common_phylum, 
              "(出现在", max(phylum_counts), "个MAG中)\n")
        }
      }
    }
  }
}

# 10.4 如果还是找不到，使用更宽松的匹配
cat("\n进行更宽松的匹配...\n")

# 创建基因名称的关键词
gene_keywords <- list(
  "Squalene_synthase" = c("squalene", "synthase"),
  "Squalene_monooxygenase" = c("squalene", "monooxygenase"),
  "FPP_synthase" = c("farnesyl", "synthase"),
  "OSC" = c("oxidosqualene", "cyclase"),
  "HMG-CoA" = c("HMG-CoA", "hydroxymethylglutaryl"),
  "IspG" = c("IspG", "MEP")
)

for (i in 1:length(analyzed_genes)) {
  if (gene_phylum_info_new$Phylum[i] == "Unknown") {
    gene_name <- analyzed_genes[i]
    
    # 尝试关键词匹配
    keywords <- if (!is.null(gene_keywords[[gene_name]])) {
      gene_keywords[[gene_name]]
    } else {
      # 将基因名称拆分为单词
      strsplit(gsub("_", " ", gene_name), " ")[[1]]
    }
    
    # 在描述中搜索关键词
    desc_matches <- ko_mag_data %>%
      filter(grepl(paste(keywords, collapse = "|"), Description, ignore.case = TRUE))
    
    if (nrow(desc_matches) > 0) {
      # 获取所有匹配记录的门分类
      all_phyla <- ko_mag_tax %>%
        filter(KO_ID %in% desc_matches$KO_ID) %>%
        filter(!is.na(Phylum)) %>%
        filter(Phylum != "Unknown")
      
      if (nrow(all_phyla) > 0) {
        phylum_counts <- table(all_phyla$Phylum)
        most_common_phylum <- names(phylum_counts)[which.max(phylum_counts)]
        gene_phylum_info_new$Phylum[i] <- most_common_phylum
        gene_phylum_info_new$KO_ID[i] <- desc_matches$KO_ID[1]
        gene_phylum_info_new$Match_Method[i] <- "关键词匹配"
        
        cat(gene_name, " (通过关键词匹配) 的主要门分类:", most_common_phylum, 
            "(出现在", max(phylum_counts), "个MAG中)\n")
      }
    }
  }
}

# 10.5 更新注释数据框
annotation_df_new <- data.frame(
  Phylum = gene_phylum_info_new$Phylum,
  row.names = analyzed_genes
)

cat("\n=== 更新后的基因门分类分布 ===\n")
print(table(annotation_df_new$Phylum))

# 10.6 重新创建热图（如果找到了门分类）- 改为PDF输出
if (any(annotation_df_new$Phylum != "Unknown")) {
  cat("\n重新创建带门分类注释的热图（PDF格式）...\n")
  
  # 准备颜色方案
  unique_phyla_new <- unique(annotation_df_new$Phylum)
  n_phyla_new <- length(unique_phyla_new)
  
  if (n_phyla_new <= 8) {
    phylum_colors_new <- brewer.pal(max(3, n_phyla_new), "Set2")[1:n_phyla_new]
  } else {
    phylum_colors_new <- colorRampPalette(brewer.pal(8, "Set2"))(n_phyla_new)
  }
  
  names(phylum_colors_new) <- unique_phyla_new
  annotation_colors_new <- list(Phylum = phylum_colors_new)
  
  # 创建热图 - 保存为PDF
  pdf("角鲨烯功能基因相关性热图_修复带分类.pdf", width = 8, height = 6)
  pheatmap(gene_corr,
           color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
           clustering_method = "complete",
           annotation_row = annotation_df_new,
           annotation_colors = annotation_colors_new,
           main = "角鲨烯功能基因相关性热图 (带菌门分类注释)",
           fontsize = 10,
           fontsize_row = 10,
           fontsize_col = 10,
           display_numbers = TRUE,
           number_format = "%.2f")
  dev.off()
  
  cat("新的热图已保存到: 角鲨烯功能基因相关性热图_修复带分类.pdf\n")
}

# 10.7 创建基因丰度的箱线图 - 改为PDF输出
library(ggplot2)

# 准备数据
abundance_long <- as.data.frame(abundance_matrix) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(cols = -Gene_ID, names_to = "Sample", values_to = "Abundance") %>%
  left_join(annotation_df_new %>% rownames_to_column("Gene_ID"), by = "Gene_ID")

# 创建箱线图
p_boxplot <- ggplot(abundance_long, aes(x = reorder(Gene_ID, Abundance, median), 
                                        y = Abundance, 
                                        fill = Phylum)) +
  geom_boxplot() +
  labs(title = "角鲨烯功能基因在各样本中的丰度分布",
       x = "基因",
       y = "丰度",
       fill = "菌门") +  # 修改图例标题
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

# 如果有门分类信息，添加颜色
if (any(annotation_df_new$Phylum != "Unknown")) {
  p_boxplot <- p_boxplot + scale_fill_manual(values = phylum_colors_new)
}

# 保存为PDF
ggsave("基因丰度分布箱线图.pdf", p_boxplot, width = 10, height = 6)
cat("\n基因丰度箱线图已保存到: 基因丰度分布箱线图.pdf\n")

# 10.8 查看最丰富的样本
cat("\n")
cat(rep("-", 60), "\n", sep = "")
cat("各样本中角鲨烯基因总丰度排名:\n")
cat(rep("-", 60), "\n", sep = "")

sample_totals <- colSums(abundance_matrix, na.rm = TRUE)
sample_ranking <- data.frame(
  Sample = names(sample_totals),
  Total_Abundance = sample_totals
) %>%
  arrange(desc(Total_Abundance))

print(head(sample_ranking, 10))

# 10.9 创建菌门-基因映射关系可视化
if (any(annotation_df_new$Phylum != "Unknown")) {
  cat("\n创建菌门-基因映射关系图...\n")
  
  # 准备菌门-基因映射数据
  phylum_gene_mapping <- annotation_df_new %>%
    rownames_to_column("Gene") %>%
    group_by(Phylum) %>%
    summarise(
      Genes = paste(Gene, collapse = ", "),
      Gene_Count = n(),
      .groups = "drop"
    ) %>%
    filter(Phylum != "Unknown")
  
  if (nrow(phylum_gene_mapping) > 0) {
    # 创建菌门-基因映射条形图
    p_mapping <- ggplot(phylum_gene_mapping, 
                        aes(x = reorder(Phylum, Gene_Count), y = Gene_Count, fill = Phylum)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = Gene_Count), hjust = -0.3, size = 3) +
      labs(title = "角鲨烯功能基因在各菌门中的分布",
           x = "菌门",
           y = "基因数量",
           fill = "菌门") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5)) +
      scale_fill_manual(values = phylum_colors_new[names(phylum_colors_new) %in% phylum_gene_mapping$Phylum]) +
      coord_flip()
    
    # 保存为PDF
    ggsave("菌门-基因映射分布.pdf", p_mapping, width = 8, height = 6)
    cat("菌门-基因映射分布图已保存到: 菌门-基因映射分布.pdf\n")
    
    # 打印详细映射关系
    cat("\n=== 菌门-基因详细映射关系 ===\n")
    for (i in 1:nrow(phylum_gene_mapping)) {
      cat(sprintf("%d. %s (包含 %d 个基因):\n", 
                  i, phylum_gene_mapping$Phylum[i], phylum_gene_mapping$Gene_Count[i]))
      cat("   基因列表: ", phylum_gene_mapping$Genes[i], "\n\n")
    }
  }
}

# 10.10 保存修复后的基因-门对应关系
gene_phylum_final <- data.frame(
  Gene_Name = analyzed_genes,
  Suggested_KO_ID = ifelse(is.na(gene_phylum_info_new$KO_ID), "", gene_phylum_info_new$KO_ID),
  Phylum = gene_phylum_info_new$Phylum,
  Match_Method = ifelse("Match_Method" %in% colnames(gene_phylum_info_new), 
                        gene_phylum_info_new$Match_Method, 
                        "直接匹配")
)

write.csv(gene_phylum_final, "基因-菌门对应关系_修复.csv", row.names = FALSE)
cat("\n修复后的基因-菌门对应关系已保存到: 基因-菌门对应关系_修复.csv\n")

# 10.11 创建菌门相关性网络图（可选）
if (any(annotation_df_new$Phylum != "Unknown")) {
  cat("\n创建菌门相关性网络图...\n")
  
  # 按菌门分组计算平均丰度
  phylum_abundance <- abundance_long %>%
    group_by(Phylum, Sample) %>%
    summarise(Phylum_Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Sample, values_from = Phylum_Abundance)
  
  # 转换为矩阵
  phylum_matrix <- as.matrix(phylum_abundance[, -1])
  rownames(phylum_matrix) <- phylum_abundance$Phylum
  
  # 计算菌门间的相关性
  phylum_corr <- cor(t(phylum_matrix), use = "complete.obs")
  
  # 创建菌门相关性热图
  pdf("菌门间相关性热图.pdf", width = 7, height = 6)
  pheatmap(phylum_corr,
           color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
           clustering_method = "complete",
           main = "各菌门间角鲨烯基因平均丰度相关性",
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize = 12,
           fontsize_row = 10,
           fontsize_col = 10)
  dev.off()
  
  cat("菌门相关性热图已保存到: 菌门间相关性热图.pdf\n")
}

# 10.12 保存最终的报告 - 简化版本，避免复杂的gridExtra问题
cat("\n生成详细分析报告...\n")

# 创建简化的文本报告
sink("角鲨烯基因分析报告_详细.txt")

cat(rep("=", 70), "\n", sep = "")
cat("角鲨烯功能基因详细分析报告\n")
cat(rep("=", 70), "\n\n", sep = "")
cat("分析日期:", format(Sys.Date(), "%Y年%m月%d日"), "\n")
cat("分析时间:", format(Sys.time(), "%H:%M:%S"), "\n\n")
cat("基于菌门分类的基因功能映射分析\n\n")

cat("1. 数据概况\n")
cat(rep("-", 70), "\n", sep = "")
cat("分析的基因数量:", nrow(gene_corr), "\n")
cat("样本数量:", ncol(abundance_matrix), "\n")
cat("MAG数量:", nrow(mag_taxonomy), "\n")
cat("KO-MAG对应记录:", nrow(ko_mag_data), "\n")
cat("鉴定到的菌门数量:", length(unique(annotation_df_new$Phylum[annotation_df_new$Phylum != "Unknown"])), "\n\n")

cat("2. 基因基本信息及菌门映射\n")
cat(rep("-", 70), "\n", sep = "")
for (i in 1:length(analyzed_genes)) {
  gene <- analyzed_genes[i]
  phylum <- annotation_df_new[gene, "Phylum"]
  avg_abundance <- gene_avg_abundance$Average_Abundance[gene_avg_abundance$Gene_ID == gene]
  match_method <- ifelse(i <= nrow(gene_phylum_final), gene_phylum_final$Match_Method[i], "未知")
  
  cat(i, ". ", gene, "\n", sep = "")
  cat("   菌门分类: ", phylum, "\n", sep = "")
  cat("   平均丰度: ", round(avg_abundance, 3), "\n", sep = "")
  cat("   匹配方法: ", match_method, "\n", sep = "")
  
  # 尝试查找描述
  if (i <= nrow(gene_phylum_final) && 
      gene_phylum_final$Suggested_KO_ID[i] != "" && 
      !is.na(gene_phylum_final$Suggested_KO_ID[i])) {
    ko_id <- gene_phylum_final$Suggested_KO_ID[i]
    gene_desc <- ko_mag_data %>%
      filter(KO_ID == ko_id) %>%
      slice(1)
    
    if (nrow(gene_desc) > 0) {
      cat("   建议的KO编号: ", ko_id, "\n", sep = "")
      cat("   描述: ", gene_desc$Description[1], "\n", sep = "")
    }
  }
  cat("\n")
}

cat("3. 样本角鲨烯基因总丰度排名（前10）\n")
cat(rep("-", 70), "\n", sep = "")
for (i in 1:min(10, nrow(sample_ranking))) {
  cat(i, ". ", sample_ranking$Sample[i], ": ", 
      round(sample_ranking$Total_Abundance[i], 2), "\n", sep = "")
}

cat("\n4. 菌门-基因映射统计\n")
cat(rep("-", 70), "\n", sep = "")
if (any(annotation_df_new$Phylum != "Unknown")) {
  phylum_gene_mapping <- annotation_df_new %>%
    rownames_to_column("Gene") %>%
    group_by(Phylum) %>%
    summarise(
      Genes = paste(Gene, collapse = ", "),
      Gene_Count = n(),
      .groups = "drop"
    ) %>%
    filter(Phylum != "Unknown") %>%
    arrange(desc(Gene_Count))
  
  for (i in 1:nrow(phylum_gene_mapping)) {
    cat(sprintf("%d. %s (包含 %d 个基因):\n", 
                i, phylum_gene_mapping$Phylum[i], phylum_gene_mapping$Gene_Count[i]))
    cat("   基因列表: ", phylum_gene_mapping$Genes[i], "\n\n")
  }
} else {
  cat("未找到有效的菌门分类信息\n")
}

cat("5. 基因间相关性分析\n")
cat(rep("-", 70), "\n", sep = "")
cat("基因之间的相关性矩阵:\n\n")

# 格式化显示相关性矩阵
corr_formatted <- round(gene_corr, 3)
print(corr_formatted)

cat("\n最强正相关性（前5）:\n")
corr_df <- as.data.frame(as.table(gene_corr))
colnames(corr_df) <- c("Gene1", "Gene2", "Correlation")
corr_df <- corr_df %>%
  filter(Gene1 != Gene2) %>%
  arrange(desc(Correlation))

for (i in 1:min(5, nrow(corr_df))) {
  cat("  ", corr_df$Gene1[i], " - ", corr_df$Gene2[i], ": ", 
      round(corr_df$Correlation[i], 3), "\n", sep = "")
}

cat("\n最强负相关性（前5）:\n")
corr_df_neg <- corr_df %>%
  arrange(Correlation)

for (i in 1:min(5, nrow(corr_df_neg))) {
  cat("  ", corr_df_neg$Gene1[i], " - ", corr_df_neg$Gene2[i], ": ", 
      round(corr_df_neg$Correlation[i], 3), "\n", sep = "")
}

cat("\n6. 生成的文件列表\n")
cat(rep("-", 70), "\n", sep = "")
cat("1. 角鲨烯功能基因相关性热图_修复带分类.pdf\n")
cat("2. 基因丰度分布箱线图.pdf\n")
cat("3. 菌门-基因映射分布.pdf\n")
cat("4. 菌门间相关性热图.pdf\n")
cat("5. 基因-菌门对应关系_修复.csv\n")
cat("6. 角鲨烯基因分析报告_详细.txt (本文件)\n")

sink()

# 另外创建一个简单的PDF总结报告
if (requireNamespace("ggplot2", quietly = TRUE) && any(annotation_df_new$Phylum != "Unknown")) {
  cat("\n创建PDF总结报告...\n")
  
  pdf("角鲨烯基因分析总结.pdf", width = 11, height = 8.5)
  
  # 页面1: 标题
  plot.new()
  text(0.5, 0.7, "角鲨烯功能基因分析总结报告", 
       cex = 2, font = 2)
  text(0.5, 0.6, paste("分析日期:", format(Sys.Date(), "%Y年%m月%d日")), 
       cex = 1.2)
  text(0.5, 0.55, paste("分析时间:", format(Sys.time(), "%H:%M:%S")), 
       cex = 1.2)
  text(0.5, 0.45, "基于菌门分类的基因功能映射分析", 
       cex = 1.5, font = 2)
  
  # 页面2: 数据概况
  plot.new()
  text(0.1, 0.9, "数据概况", cex = 1.8, font = 2, adj = 0)
  
  summary_text <- c(
    paste("分析的基因数量:", nrow(gene_corr)),
    paste("样本数量:", ncol(abundance_matrix)),
    paste("MAG数量:", nrow(mag_taxonomy)),
    paste("KO-MAG对应记录:", nrow(ko_mag_data)),
    paste("鉴定到的菌门数量:", length(unique(annotation_df_new$Phylum[annotation_df_new$Phylum != "Unknown"]))),
    "",
    "生成的文件:",
    "1. 角鲨烯功能基因相关性热图_修复带分类.pdf",
    "2. 基因丰度分布箱线图.pdf",
    "3. 菌门-基因映射分布.pdf",
    "4. 菌门间相关性热图.pdf",
    "5. 基因-菌门对应关系_修复.csv",
    "6. 角鲨烯基因分析报告_详细.txt"
  )
  
  for (i in 1:length(summary_text)) {
    text(0.1, 0.8 - i*0.05, summary_text[i], cex = 1, adj = 0)
  }
  
  # 页面3: 菌门分布
  plot.new()
  text(0.1, 0.9, "菌门-基因映射统计", cex = 1.8, font = 2, adj = 0)
  
  if (any(annotation_df_new$Phylum != "Unknown")) {
    phylum_stats <- table(annotation_df_new$Phylum)
    phylum_stats <- phylum_stats[names(phylum_stats) != "Unknown"]
    
    y_pos <- 0.8
    for (i in 1:length(phylum_stats)) {
      phylum_name <- names(phylum_stats)[i]
      gene_count <- phylum_stats[i]
      
      # 获取该菌门的基因列表
      genes_in_phylum <- rownames(annotation_df_new)[annotation_df_new$Phylum == phylum_name]
      
      text(0.1, y_pos, paste(phylum_name, ": ", gene_count, "个基因", sep = ""), 
           cex = 1, adj = 0)
      y_pos <- y_pos - 0.04
      
      if (length(genes_in_phylum) > 0) {
        genes_text <- paste(genes_in_phylum, collapse = ", ")
        text(0.15, y_pos, paste("基因: ", genes_text), cex = 0.8, adj = 0)
        y_pos <- y_pos - 0.03
      }
      
      if (y_pos < 0.1) {
        break  # 防止超出页面
      }
    }
  } else {
    text(0.1, 0.5, "未找到有效的菌门分类信息", cex = 1.2, adj = 0)
  }
  
  dev.off()
  cat("PDF总结报告已保存到: 角鲨烯基因分析总结.pdf\n")
}

cat("\n")
cat(rep("*", 60), "\n", sep = "")
cat("*             所有分析已完成！请检查生成的文件。              *\n")
cat(rep("*", 60), "\n", sep = "")
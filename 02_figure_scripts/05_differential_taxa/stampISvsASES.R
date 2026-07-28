# 0 环境 & 路径 -------------------------------------------------------------
# 安装必要的包
p_list = c("tidyverse", "ggsci", "magrittr", "ggh4x", "rstatix", 
           "ggsignif", "ggpubr", "ggnewscale", "patchwork", 
           "reshape2", "ggplot2", "broom", "ggrepel")

for(p in p_list){
  if (!requireNamespace(p, quietly = TRUE)){
    install.packages(p)
  }
  library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
}

setwd("/Users/catherine/Downloads/stamp 作图")

# 1 读 OTU 表和分组表 --------------------------------------------------------
# 读取物种数据
cat("读取数据...\n")
otu <- read_csv("species-stamp.csv") %>% 
  column_to_rownames("...1")

# 读取分组信息
grp <- read_csv("group.csv") %>% 
  mutate(sample = str_trim(simple)) %>% 
  select(sample, GROUP)

# 2 数据预处理 --------------------------------------------------------------
# 转置数据并合并分组信息
data <- otu %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  left_join(grp, by = "sample") %>% 
  filter(!is.na(GROUP)) %>%          
  mutate(GROUP = factor(GROUP))

# 检查所有分组
cat("\n数据中的分组类别:\n")
print(table(data$GROUP))

# 3 定义所有要比较的组对 -------------------------------------------------------
# 定义所有可能的IS与其他组的比较
comparisons <- list(
  IS_vs_NS = c("IS", "NS"),
  IS_vs_AS = c("IS", "AS"),
  IS_vs_ES = c("IS", "ES")
)

# 只保留实际存在的比较
available_groups <- unique(data$GROUP)
valid_comparisons <- list()

cat("\n检查可用的比较组:\n")
for (comp_name in names(comparisons)) {
  grp1 <- comparisons[[comp_name]][1]
  grp2 <- comparisons[[comp_name]][2]
  
  if (grp1 %in% available_groups & grp2 %in% available_groups) {
    # 检查每组是否有足够样本
    n1 <- sum(data$GROUP == grp1)
    n2 <- sum(data$GROUP == grp2)
    
    if (n1 >= 3 & n2 >= 3) {  # 每组至少3个样本
      valid_comparisons[[comp_name]] <- c(grp1, grp2)
      cat(sprintf("✓ %s: %s (n=%d) vs %s (n=%d)\n", 
                  comp_name, grp1, n1, grp2, n2))
    } else {
      cat(sprintf("✗ %s: 跳过 - 样本数不足 (%s: %d, %s: %d)\n", 
                  comp_name, grp1, n1, grp2, n2))
    }
  } else {
    missing_groups <- c()
    if (!(grp1 %in% available_groups)) missing_groups <- c(missing_groups, grp1)
    if (!(grp2 %in% available_groups)) missing_groups <- c(missing_groups, grp2)
    cat(sprintf("✗ %s: 跳过 - 缺少分组: %s\n", 
                comp_name, paste(missing_groups, collapse = ", ")))
  }
}

if (length(valid_comparisons) == 0) {
  stop("没有有效的组对可以比较！")
}

# 4 数据分析函数 -----------------------------------------------------------
perform_stamp_analysis <- function(grp1, grp2, comparison_name, top_n = 10) {
  cat(sprintf("\n%s: 开始分析 %s vs %s\n", 
              comparison_name, grp1, grp2))
  
  # 4.1 筛选数据
  sub_data <- data %>% 
    filter(GROUP %in% c(grp1, grp2)) %>% 
    mutate(GROUP = factor(GROUP, levels = c(grp1, grp2)))
  
  # 4.2 数据清理
  numeric_cols <- sub_data %>% 
    select(-sample, -GROUP) %>% 
    names()
  
  clean_data <- sub_data %>% 
    mutate(across(all_of(numeric_cols), ~ {
      x <- ifelse(is.na(.), 0, .)
      x <- ifelse(is.infinite(x), max(x[is.finite(x)], na.rm = TRUE), x)
      return(x)
    }))
  
  # 检查是否需要转换为百分比
  max_val <- max(as.matrix(clean_data %>% select(all_of(numeric_cols))), na.rm = TRUE)
  if (max_val < 1) {
    clean_data <- clean_data %>% 
      mutate(across(all_of(numeric_cols), ~ . * 100))
  }
  
  # 4.3 过滤低丰度物种 (平均丰度 > 0.5%)
  mean_abundance <- clean_data %>% 
    select(all_of(numeric_cols)) %>% 
    summarise(across(everything(), mean)) %>% 
    unlist()
  
  keep_species <- names(mean_abundance[mean_abundance > 0.5])
  
  if (length(keep_species) == 0) {
    cat("警告：没有物种满足平均丰度>0.5%的条件，使用平均丰度>0.1%的物种\n")
    keep_species <- names(mean_abundance[mean_abundance > 0.1])
  }
  
  if (length(keep_species) == 0) {
    cat("警告：仍然没有物种满足条件，使用所有物种\n")
    keep_species <- numeric_cols
  }
  
  sub_data_filtered <- clean_data %>% 
    select(sample, GROUP, all_of(keep_species))
  
  cat(sprintf("  分析物种数: %d\n", length(keep_species)))
  
  # 4.4 安全Wilcoxon检验函数
  safe_wilcox_test <- function(x, group) {
    # 检查数据有效性
    valid_x <- x[!is.na(x) & !is.infinite(x)]
    if (length(valid_x) < 6) {  # 每组至少3个样本
      return(tibble(statistic = NA_real_, p.value = NA_real_))
    }
    
    if (length(unique(valid_x)) <= 1) {
      return(tibble(statistic = NA_real_, p.value = NA_real_))
    }
    
    tryCatch({
      grp_factor <- factor(group)
      group_counts <- table(grp_factor)
      if (any(group_counts < 2)) return(tibble(statistic = NA_real_, p.value = NA_real_))
      
      result <- wilcox.test(x ~ grp_factor, conf.int = TRUE, exact = FALSE, na.rm = TRUE)
      medians <- tapply(x, grp_factor, median, na.rm = TRUE)
      
      broom::tidy(result) %>% 
        mutate(
          estimate1 = medians[1],
          estimate2 = medians[2]
        )
    }, error = function(e) {
      return(tibble(statistic = NA_real_, p.value = NA_real_))
    })
  }
  
  # 4.5 进行Wilcoxon检验
  group_var <- sub_data_filtered$GROUP
  
  # 分批处理避免内存问题
  cat("  进行Wilcoxon检验...\n")
  diff_results_list <- list()
  batch_size <- 30
  n_species <- length(keep_species)
  n_batches <- ceiling(n_species / batch_size)
  
  if (interactive()) {
    pb <- txtProgressBar(min = 0, max = n_batches, style = 3)
  }
  
  for (i in 1:n_batches) {
    start_idx <- (i-1) * batch_size + 1
    end_idx <- min(i * batch_size, n_species)
    current_species <- keep_species[start_idx:end_idx]
    
    batch_results <- sub_data_filtered %>% 
      select(all_of(current_species)) %>% 
      map_dfr(~ safe_wilcox_test(.x, group_var), .id = "species")
    
    diff_results_list[[i]] <- batch_results
    if (interactive()) {
      setTxtProgressBar(pb, i)
    }
  }
  if (interactive()) {
    close(pb)
  }
  
  # 合并结果
  diff_results <- bind_rows(diff_results_list) %>% 
    filter(!is.na(p.value)) %>% 
    mutate(
      adj.p = p.adjust(p.value, method = "BH"),
      median_diff = estimate2 - estimate1,
      log2FC = log2((estimate2 + 0.001) / (estimate1 + 0.001)),  # 避免除零
      direction = ifelse(median_diff > 0, paste0(grp2, " > ", grp1), 
                         paste0(grp1, " > ", grp2))
    ) %>% 
    arrange(p.value)
  
  cat(sprintf("  成功检验物种数: %d\n", nrow(diff_results)))
  
  # 4.6 筛选Top显著物种
  if (nrow(diff_results) > 0) {
    if (nrow(diff_results) < top_n) {
      actual_top_n <- nrow(diff_results)
      cat(sprintf("  注意：只有 %d 个物种有检验结果，显示全部\n", actual_top_n))
      significant_species <- diff_results$species
    } else {
      significant_species <- diff_results %>% 
        arrange(adj.p) %>% 
        head(top_n) %>% 
        pull(species)
      cat(sprintf("  显示Top %d显著物种\n", top_n))
    }
  } else {
    significant_species <- character(0)
    cat("  警告：没有获得有效的检验结果\n")
  }
  
  # 4.7 准备绘图数据
  if (length(significant_species) > 0) {
    bar_data <- sub_data_filtered %>% 
      select(all_of(significant_species), GROUP) %>% 
      pivot_longer(-GROUP, names_to = "Species", values_to = "Value") %>% 
      group_by(Species, GROUP) %>% 
      summarise(
        Mean = mean(Value, na.rm = TRUE),
        SE = sd(Value, na.rm = TRUE) / sqrt(n()),
        Median = median(Value, na.rm = TRUE),
        .groups = "drop"
      ) %>% 
      mutate(
        Species = factor(Species, levels = rev(significant_species))
      )
    
    scatter_data <- diff_results %>% 
      filter(species %in% significant_species) %>% 
      mutate(
        Species = factor(species, levels = rev(significant_species)),
        p.value.display = case_when(
          adj.p < 0.001 ~ "< 0.001",
          adj.p < 0.01 ~ sprintf("%.4f", adj.p),
          TRUE ~ sprintf("%.3f", adj.p)
        ),
        Significance = case_when(
          adj.p < 0.001 ~ "***",
          adj.p < 0.01 ~ "**",
          adj.p < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      ) %>% 
      arrange(desc(Species))
    
    cat("  绘图数据准备完成\n")
  } else {
    bar_data <- NULL
    scatter_data <- NULL
  }
  
  # 4.8 绘图函数 - 修复版本
  create_stamp_plot <- function(bar_data, scatter_data, grp1, grp2, comparison_name, top_n) {
    if (is.null(bar_data) || nrow(bar_data) == 0) {
      cat("  跳过绘图：没有有效数据\n")
      return(NULL)
    }
    
    # 颜色方案
    colors <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")  # IS, NS, ES, AS
    names(colors) <- c("IS", "NS", "ES", "AS")
    plot_colors <- colors[c(grp1, grp2)]
    
    # 获取物种显示名（简化长名称）
    species_display_names <- sapply(levels(bar_data$Species), function(name) {
      # 如果是完整的分类路径，取最后的物种名
      if (grepl("s__", name)) {
        parts <- strsplit(name, " ")[[1]]
        species_part <- parts[length(parts)]
        # 提取物种名
        if (grepl("s__", species_part)) {
          return(gsub("s__", "", species_part))
        }
      }
      # 如果名称太长，截断
      if (nchar(name) > 50) {
        return(paste0(substr(name, 1, 47), "..."))
      }
      return(name)
    })
    
    # 创建交替背景色函数
    create_striped_background <- function(n_items) {
      # 计算条纹的数量
      n_stripes <- floor((n_items + 1) / 2)
      xmin_seq <- seq(0.5, n_items + 0.5, by = 2)
      xmax_seq <- seq(1.5, n_items + 0.5, by = 2)
      
      # 确保序列长度一致
      min_length <- min(length(xmin_seq), length(xmax_seq))
      xmin_seq <- xmin_seq[1:min_length]
      xmax_seq <- xmax_seq[1:min_length]
      
      list(xmin = xmin_seq, xmax = xmax_seq)
    }
    
    # 获取物种数量
    n_species <- length(unique(bar_data$Species))
    bg_stripes <- create_striped_background(n_species)
    
    # 4.8.1 左侧条形图 - 修复版
    p_left <- ggplot(bar_data, aes(x = Species, y = Mean, fill = GROUP)) +
      # 交替背景色 - 修复版本
      annotate("rect", 
               xmin = bg_stripes$xmin,
               xmax = bg_stripes$xmax,
               ymin = -Inf, ymax = Inf,
               fill = "grey95", alpha = 0.3) +
      geom_bar(stat = "identity", 
               position = position_dodge(0.85), 
               width = 0.7, 
               color = "black", 
               linewidth = 0.15) +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                    position = position_dodge(0.85),
                    width = 0.22,
                    size = 0.22) +
      coord_flip() +
      scale_x_discrete(labels = species_display_names) +
      scale_fill_manual(values = plot_colors,
                        labels = c(grp1, grp2)) +
      labs(x = "", 
           y = "Mean proportion (%)",
           fill = "Group") +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.ticks = element_line(color = 'black', size = 0.3),
        axis.line = element_line(colour = "black", size = 0.3),
        axis.title.x = element_text(colour = 'black', 
                                    size = 10, 
                                    face = "bold",
                                    margin = margin(t = 8)),
        axis.text.x = element_text(colour = 'black', size = 9),
        axis.text.y = element_text(colour = 'black', size = 8, face = "italic"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(t = 0, b = 5),
        plot.margin = margin(15, 10, 15, 15)
      )
    
    # 4.8.2 中间散点图 - 修复版
    bg_stripes_scatter <- create_striped_background(nrow(scatter_data))
    
    p_center <- ggplot(scatter_data, 
                       aes(x = Species, y = median_diff, 
                           ymin = conf.low, ymax = conf.high)) +
      # 交替背景色
      annotate("rect", 
               xmin = bg_stripes_scatter$xmin,
               xmax = bg_stripes_scatter$xmax,
               ymin = -Inf, ymax = Inf,
               fill = "grey95", alpha = 0.3) +
      geom_hline(yintercept = 0, 
                 linetype = 'dashed', 
                 color = 'grey60', 
                 linewidth = 0.4) +
      geom_errorbar(width = 0.15, 
                    size = 0.3,
                    color = "black") +
      geom_point(aes(fill = median_diff > 0), 
                 shape = 21, 
                 size = 2.5, 
                 stroke = 0.4) +
      # 添加显著性标记
      geom_text(aes(label = Significance), 
                size = 3.5, 
                y = ifelse(scatter_data$median_diff >= 0, 
                           scatter_data$conf.high + 0.5,
                           scatter_data$conf.low - 0.5),
                vjust = 0.5, 
                fontface = "bold") +
      coord_flip() +
      scale_x_discrete(labels = species_display_names) +
      scale_fill_manual(values = c("FALSE" = plot_colors[1], "TRUE" = plot_colors[2]),
                        guide = "none") +
      labs(x = "", 
           y = paste0("Difference in median (%)\n", grp2, " - ", grp1)) +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.ticks = element_line(color = 'black', size = 0.3),
        axis.line = element_line(colour = "black", size = 0.3),
        axis.title.x = element_text(colour = 'black', 
                                    size = 10, 
                                    face = "bold",
                                    margin = margin(t = 8)),
        axis.text.x = element_text(colour = 'black', size = 9),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(15, 10, 15, 10)
      )
    
    # 4.8.3 右侧p值标签 - 修复版
    bg_stripes_pvalue <- create_striped_background(nrow(scatter_data))
    
    p_right <- ggplot(scatter_data, aes(x = Species)) +
      # 交替背景色
      annotate("rect", 
               xmin = bg_stripes_pvalue$xmin,
               xmax = bg_stripes_pvalue$xmax,
               ymin = -Inf, ymax = Inf,
               fill = "grey95", alpha = 0.3) +
      geom_text(aes(y = 0.35, label = p.value.display), 
                size = 3.2, 
                hjust = 0,
                fontface = "bold") +
      geom_text(aes(x = median(1:nrow(scatter_data)), 
                    y = 0.85), 
                label = "Adj. p-value", 
                srt = 90,
                fontface = "bold", 
                size = 3.5) +
      coord_flip() +
      scale_x_discrete(labels = species_display_names) +
      ylim(c(0, 1)) +
      theme_void() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(15, 15, 15, 10)
      )
    
    # 4.8.5 组合图形
    combined_plot <- (p_left | p_center | p_right) + 
      plot_layout(widths = c(3.8, 2.2, 1.5),  # 调整宽度比例
                  guides = 'collect') &
      theme(legend.position = "top",
            legend.box.margin = margin(0, 0, 0, 0))
    
    # 添加标题
    title_text <- paste0("STAMP Plot: ", comparison_name, " (Top ", top_n, " Species)")
    combined_plot <- combined_plot + 
      plot_annotation(
        title = title_text,
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5,
                                                margin = margin(b = 10)))
      )
    
    return(combined_plot)
  }
  
  # 4.9 创建图形
  if (length(significant_species) > 0) {
    plot <- create_stamp_plot(bar_data, scatter_data, grp1, grp2, comparison_name, top_n)
  } else {
    plot <- NULL
    cat("  跳过绘图：没有显著物种\n")
  }
  
  # 4.10 准备返回结果
  results <- list(
    comparison_name = comparison_name,
    grp1 = grp1,
    grp2 = grp2,
    n_samples = nrow(sub_data_filtered),
    n_species_analyzed = length(keep_species),
    n_species_tested = nrow(diff_results),
    diff_results = diff_results,
    top_species = significant_species,
    bar_data = bar_data,
    scatter_data = scatter_data,
    plot = plot
  )
  
  cat(sprintf("  %s分析完成！\n", comparison_name))
  return(results)
}

# 5 执行所有比较分析 -------------------------------------------------------
cat("\n")
cat(paste0(rep("=", 70), collapse = ""))
cat("\n开始STAMP分析\n")
cat(paste0(rep("=", 70), collapse = ""))
cat("\n")

all_results <- list()

for (comp_name in names(valid_comparisons)) {
  grp_pair <- valid_comparisons[[comp_name]]
  results <- perform_stamp_analysis(grp_pair[1], grp_pair[2], comp_name, top_n = 10)
  
  if (!is.null(results$plot)) {
    all_results[[comp_name]] <- results
    
    # 单独保存每个比较的图形
    output_file <- paste0("STAMP_", comp_name, "_Top10.pdf")
    
    # 计算合适的图形尺寸
    n_species <- length(results$top_species)
    plot_height <- max(6, n_species * 0.35 + 3)  # 调整高度
    plot_width <- 18  # 增加宽度
    
    tryCatch({
      ggsave(output_file, 
             results$plot, 
             width = plot_width, 
             height = plot_height,
             dpi = 300)
      cat(sprintf("✓ 图形已保存: %s (%.1f x %.1f 英寸)\n", 
                  output_file, plot_width, plot_height))
    }, error = function(e) {
      cat(sprintf("✗ 保存图形失败: %s\n", e$message))
    })
    
    # 保存统计结果
    stats_file <- paste0("STAMP_stats_", comp_name, "_Top10.csv")
    
    if (!is.null(results$diff_results) && nrow(results$diff_results) > 0) {
      stats_output <- results$diff_results %>% 
        filter(species %in% results$top_species) %>% 
        select(
          Species = species,
          Median_grp1 = estimate1,
          Median_grp2 = estimate2,
          Median_Difference = median_diff,
          Log2FC = log2FC,
          CI_lower = conf.low,
          CI_upper = conf.high,
          W_statistic = statistic,
          P_value = p.value,
          FDR_adjusted_p = adj.p,
          Direction = direction
        ) %>% 
        mutate(
          Significance = case_when(
            FDR_adjusted_p < 0.001 ~ "***",
            FDR_adjusted_p < 0.01 ~ "**",
            FDR_adjusted_p < 0.05 ~ "*",
            TRUE ~ "ns"
          )
        ) %>% 
        arrange(FDR_adjusted_p)
      
      write_csv(stats_output, stats_file)
      cat(sprintf("✓ 统计结果已保存: %s\n", stats_file))
    }
    
    cat("\n")
  }
}

# 6 创建汇总报告 -----------------------------------------------------------
if (length(all_results) > 0) {
  cat("\n")
  cat(paste0(rep("=", 70), collapse = ""))
  cat("\nSTAMP分析汇总报告\n")
  cat(paste0(rep("=", 70), collapse = ""))
  cat("\n")
  
  # 6.1 创建汇总表格
  summary_table <- tibble(
    Comparison = character(),
    Group1 = character(),
    Group2 = character(),
    Samples_Group1 = integer(),
    Samples_Group2 = integer(),
    Total_Samples = integer(),
    Species_Analyzed = integer(),
    Species_Tested = integer(),
    Top10_Shown = integer(),
    Min_adj_p = character(),
    Max_abs_diff = character(),
    Most_Sig_Species = character()
  )
  
  for (comp_name in names(all_results)) {
    results <- all_results[[comp_name]]
    
    # 计算各组样本数
    samples_grp1 <- sum(data$GROUP == results$grp1)
    samples_grp2 <- sum(data$GROUP == results$grp2)
    
    # 获取统计信息
    if (nrow(results$diff_results) > 0) {
      min_p <- min(results$diff_results$adj.p, na.rm = TRUE)
      max_diff <- max(abs(results$diff_results$median_diff), na.rm = TRUE)
      
      # 获取最显著物种
      if (length(results$top_species) > 0) {
        most_sig <- results$diff_results %>% 
          filter(species == results$top_species[1])
        if (nrow(most_sig) > 0) {
          most_sig_species <- paste0(
            most_sig$species[1], 
            " (adj.p=", 
            ifelse(most_sig$adj.p[1] < 0.001, "<0.001", 
                   sprintf("%.4f", most_sig$adj.p[1])), 
            ")"
          )
        } else {
          most_sig_species <- "None"
        }
      } else {
        most_sig_species <- "None"
      }
    } else {
      min_p <- NA
      max_diff <- NA
      most_sig_species <- "None"
    }
    
    summary_table <- summary_table %>% 
      add_row(
        Comparison = comp_name,
        Group1 = results$grp1,
        Group2 = results$grp2,
        Samples_Group1 = samples_grp1,
        Samples_Group2 = samples_grp2,
        Total_Samples = results$n_samples,
        Species_Analyzed = results$n_species_analyzed,
        Species_Tested = results$n_species_tested,
        Top10_Shown = length(results$top_species),
        Min_adj_p = ifelse(is.na(min_p), "NA", 
                           ifelse(min_p < 0.001, "<0.001", sprintf("%.4f", min_p))),
        Max_abs_diff = ifelse(is.na(max_diff), "NA", sprintf("%.2f%%", max_diff)),
        Most_Sig_Species = most_sig_species
      )
  }
  
  # 6.2 打印汇总表格
  cat("\n分析汇总:\n")
  print(summary_table, n = Inf)
  
  # 保存汇总表格
  tryCatch({
    write_csv(summary_table, "STAMP_Analysis_Summary.csv")
    cat("\n✓ 汇总表格已保存: STAMP_Analysis_Summary.csv\n")
  }, error = function(e) {
    cat(sprintf("\n✗ 保存汇总表格失败: %s\n", e$message))
  })
  
  # 6.3 显示简要总结
  cat("\n")
  cat(paste0(rep("-", 70), collapse = ""))
  cat("\n分析完成总结:\n")
  cat(paste0(rep("-", 70), collapse = ""))
  cat("\n")
  
  for (comp_name in names(all_results)) {
    results <- all_results[[comp_name]]
    cat(sprintf("\n%s:\n", comp_name))
    cat(sprintf("  • 比较: %s (n=%d) vs %s (n=%d)\n", 
                results$grp1, sum(data$GROUP == results$grp1),
                results$grp2, sum(data$GROUP == results$grp2)))
    cat(sprintf("  • 分析物种数: %d\n", results$n_species_analyzed))
    cat(sprintf("  • 显著物种数 (adj.p < 0.05): %d\n", 
                sum(results$diff_results$adj.p < 0.05, na.rm = TRUE)))
    cat(sprintf("  • 输出文件: STAMP_%s_Top10.pdf\n", comp_name))
  }
  
  cat("\n")
  cat(paste0(rep("=", 70), collapse = ""))
  cat("\n🎉 所有STAMP分析完成！\n")
  cat(paste0(rep("=", 70), collapse = ""))
  cat("\n✓ PDF图形文件已保存\n")
  cat("✓ CSV统计文件已保存\n")
  cat("✓ 汇总报告已生成\n")
  cat("✓ 所有输出文件保存在当前目录\n")
  
} else {
  cat("\n⚠️ 警告：没有产生任何有效结果！\n")
  cat("请检查数据格式和分组信息。\n")
}

# 7 清理工作空间 ---------------------------------------------------------
# 询问用户是否清理工作空间
cat("\n是否清理工作空间？(y/n): ")
response <- readline()

if (tolower(response) == "y") {
  keep_vars <- c("data", "all_results", "summary_table")
  rm(list = setdiff(ls(), keep_vars))
  cat("工作空间已清理，保留主要数据对象。\n")
} else {
  cat("保留所有工作空间对象。\n")
}
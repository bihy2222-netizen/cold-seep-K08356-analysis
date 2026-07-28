# 0. 工作目录（自己改）--------------------------------------------------
setwd("/Users/catherine/Downloads/metaphlan 作图/class 级别作图-柱状图")

# 1. 装包 & 加载 ---------------------------------------------------------
pkgs <- c("tidyverse", "vegan", "ggplot2")
invisible(lapply(pkgs, \(p) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)))
lapply(pkgs, library, character.only = TRUE)

# 2. 读物种表（UTF-16 BOM + 第一列变行名）--------------------------------
tax_tab <- readr::read_tsv("class34.txt",
                           locale = readr::locale(encoding = "UTF-16LE"),
                           show_col_types = FALSE) %>%
  tibble::column_to_rownames(var = names(.)[1]) %>%
  as.matrix()                                         # vegan 需要矩阵

# 3. 读分组表 ------------------------------------------------------------
group_tab <- readr::read_tsv("class_group.txt",
                             locale = readr::locale(encoding = "UTF-16LE"),
                             show_col_types = FALSE)

# 4. Alpha 多样性 ---------------------------------------------------------
Shannon  <- vegan::diversity(tax_tab, index = "shannon", MARGIN = 2)
Simpson  <- vegan::diversity(tax_tab, index = "simpson",  MARGIN = 2)
Richness <- vegan::specnumber(tax_tab, MARGIN = 2)

obs_chao_ace <- lapply(seq_len(ncol(tax_tab)), \(i) {
  tmp <- vegan::estimateR(ceiling(tax_tab[, i]))
  tibble::tibble(S.obs = tmp["S.obs"], S.chao1 = tmp["S.chao1"], S.ace = tmp["S.ace"])
}) %>% bind_rows()

alpha_df <- tibble::tibble(
  sample   = colnames(tax_tab),
  Shannon  = Shannon,
  Simpson  = Simpson,
  Richness = Richness,
  Obs      = obs_chao_ace$S.obs,
  Chao     = obs_chao_ace$S.chao1,
  Ace      = obs_chao_ace$S.ace,
  Pielou   = ifelse(Richness > 0 & Shannon >= 0, Shannon / log(Richness), NA_real_)
) %>%
  left_join(group_tab, by = c("sample" = "sample-id"))

# 5. 导出结果 -------------------------------------------------------------
write.csv(alpha_df, "group1_alpha_data_class数据.csv", row.names = FALSE)

# 6. Shannon 箱线图（去掉 NA/Inf 行）-------------------------------------
p_shannon <- ggplot(alpha_df %>% filter(!is.na(Shannon)),
                    aes(x = group, y = Shannon, fill = group)) +
  geom_boxplot(width = 0.5, alpha = 0.9) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c('#ca0020','#A39AF3','#F3C051',
                               '#92c5de','#0571b0','#0D99D4')) +
  labs(title = "Shannon diversity index") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("group2-Shannon-class.png", p_shannon, width = 16, height = 9, dpi = 300)

# 7. 所有指数分面图（去掉 NA/Inf 点）-------------------------------------
alpha_long <- alpha_df %>%
  pivot_longer(c(Shannon:Pielou), names_to = "index", values_to = "value") %>%
  filter(!is.na(value), !is.infinite(value))   # 只丢这一个缺失值

p_total <- ggplot(alpha_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.5, alpha = 0.9) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c('#ca0020','#A39AF3','#F3C051',
                               '#92c5de','#0571b0','#0D99D4')) +
  facet_wrap(~index, ncol = 1, scales = "free_y") +
  labs(title = "Alpha diversity") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("group4-total_alpha_class.png", p_total, width = 9, height = 16, dpi = 300)


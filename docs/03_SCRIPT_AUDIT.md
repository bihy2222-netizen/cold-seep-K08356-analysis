# 脚本清点与版本审计

## 清点结果

- 原始 ZIP 中共有 61 个 R/Shell 脚本文件。
- 按内容哈希识别到 57 份不同脚本内容。
- 内容完全相同但出现在多个目录的副本已合并保留一份；不同内容但同名的脚本使用 `root_` 或 `final_` 前缀区分。
- v4 主流程为当前推荐版；v1–v3、`.bak` 和个人命令笔记放入 `99_legacy_and_versions/`。

## 推荐版本

| 类别 | 推荐文件 |
|---|---|
| 宏基因组主流程 | `01_metagenomic_workflow/00_master_pipeline/metagenome_pipeline_v4.sh` |
| K08356 流程说明 | `01_metagenomic_workflow/00_master_pipeline/K08356_pipeline_v4_使用说明.md` |
| MopB/AioA/IdrA 建树 | `01_metagenomic_workflow/02_phylogeny/01_run_MopB_AioA_IdrA_phylogeny.sh` |
| 群落组成 | `02_figure_scripts/02_community_composition/柱状图最终脚本——带分组色块.R` |
| α 多样性 | `02_figure_scripts/03_alpha_diversity/阿尔法多样性_shannon_wills_FDRcorrect.R` |
| β 多样性 | `02_figure_scripts/04_beta_diversity/β多样性升级脚本ISASESNS.R` |
| UpSet | `02_figure_scripts/06_shared_unique_taxa/upsetfinal.R` |
| 宿主分类 | `02_figure_scripts/07_host_taxonomy/宿主图 final 清洗数据加柱状图无黑线.R` |
| C/N/S/As 网络 | `02_figure_scripts/09_networks/make_cnsas_habitat_nested_networks.R` |

## 已发现的问题

1. 大部分早期 R 脚本写死了 `/Users/catherine/Downloads/...`，换电脑或目录后会报错。
2. 部分脚本运行时自动安装 R 包，不适合服务器批处理；建议后续统一用 `renv`。
3. 早期 Mantel 脚本包含模拟环境数据，不能直接当作真实统计结果。
4. 文件名中的“final”“finally”“升级”不能单独证明它是论文最终版本，因此本仓库同时依据修改日期、代码内容和统计处理进行推荐。
5. 原始流程笔记曾包含服务器地址和疑似登录口令，仓库版本已替换为 `[REDACTED_*]`。
6. 本次 ZIP 不包含近期讨论的若干 K08356 专题图最终脚本，缺失项已记录在分析—作图对应表。


# Cold-seep K08356 metagenome analysis

本仓库整理南海冷泉四类生境（IS、AS、ES、NS）的宏基因组分析脚本，并将每一步分析结果与后续画图脚本按顺序对应。

## 仓库定位

- `01_metagenomic_workflow/`：先运行的宏基因组主流程，包括质控、组装、MetaPhlAn、分箱、MAG 质控、基因丰度、KOfam、K08356/AioA/IdrA 证据链和系统发育。
- `02_figure_scripts/`：分析结束后运行的画图脚本，按论文结果顺序整理。
- `03_metadata_examples/`：压缩包中附带的小型分组与样品信息示例。
- `docs/`：完整执行顺序、分析—作图对应关系及脚本审计记录。

## 第一部分：宏基因组分析顺序

当前推荐主脚本：

```text
01_metagenomic_workflow/00_master_pipeline/metagenome_pipeline_v4.sh
```

脚本默认 `DRY_RUN=1`，只打印命令，不直接运行。先初始化并核查：

```bash
export PROJECT_DIR=/path/to/metagenome_project
bash 01_metagenomic_workflow/00_master_pipeline/metagenome_pipeline_v4.sh init
bash 01_metagenomic_workflow/00_master_pipeline/metagenome_pipeline_v4.sh check
bash 01_metagenomic_workflow/00_master_pipeline/metagenome_pipeline_v4.sh dry-run-plan
```

确认样品表、数据库、Conda 环境和输出目录后，才对单个步骤设置 `DRY_RUN=0`。不要用真实运行模式一次性执行 `all`。

完整推荐顺序见 [docs/01_ANALYSIS_ORDER.md](docs/01_ANALYSIS_ORDER.md)。

Contig-level K08356/AioA 追溯树脚本按顺序保存在：

```text
01_metagenomic_workflow/02_phylogeny/02_extract_contig_level_K08356_faa.py
01_metagenomic_workflow/02_phylogeny/03_run_contig_level_K08356_phylogeny.sh
01_metagenomic_workflow/02_phylogeny/04_extract_group_derep_K08356_from_kofamscan.sh
```

## 第二部分：画图顺序

画图从采样地图开始，随后是群落组成、α/β 多样性、差异类群、共有/特有类群、宿主分类、环境关联和网络图：

```text
02_figure_scripts/
├── 01_sampling_map/
├── 02_community_composition/
├── 03_alpha_diversity/
├── 04_beta_diversity/
├── 05_differential_taxa/
├── 06_shared_unique_taxa/
├── 07_host_taxonomy/
├── 08_environment_association/
├── 09_networks/
└── 99_archive_and_debug/
```

分析输出与画图脚本的逐项对应见 [docs/02_ANALYSIS_FIGURE_CROSSWALK.md](docs/02_ANALYSIS_FIGURE_CROSSWALK.md)。

## 重要说明

1. v4 是当前推荐流程；v1–v3 和个人命令笔记只用于方法追溯，不建议直接执行。
2. K08356、AioA 和 IdrA 的功能分类必须联合 HMM、系统发育和基因邻域证据；不能只凭 KOfam 或单一 HMM 命中定义功能。
3. `MAG TPM` 表示 reads 回贴得到的宿主 MAG 丰度，不等同于分箱来源数量。
4. 宏基因组 TPM 表示基因组层面的功能潜力，不代表转录活性。
5. 多数旧 R 脚本保留了原电脑的 `setwd()`。正式运行前需改为项目路径；这些脚本已按用途归档，但没有假设用户本地数据布局并强行改写。
6. 原压缩包中的服务器地址和疑似登录口令已经脱敏。不要把密码、Token 或服务器登录信息写入 Git 仓库。

## 未纳入仓库的内容

为避免仓库臃肿和泄露本机状态，以下内容未上传：

- 本机 `R_libs/` 包目录；
- `.RData`、`.Rhistory`、`.DS_Store` 和 `__MACOSX`；
- 已生成的 PNG、临时结果和 Word 文档；
- 大型原始 reads、组装结果、数据库和 MAG 文件。

这些文件应通过 Conda、renv 或独立数据存储管理，不应提交到 Git。

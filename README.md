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
01_metagenomic_workflow/02_phylogeny/05_extract_MAG_proteins_by_hmmsearch.sh
```

统一参考的 K08356-bearing MAG 丰度脚本保存在：

```text
01_metagenomic_workflow/03_abundance/02_run_unified_MAG48_coverm_tpm.sh
01_metagenomic_workflow/03_abundance/03_audit_unified_MAG48_tpm.R
```

四个生境的 56 个样本必须全部映射到同一套 48-MAG 参考。旧的分生境
derep 目录结果不能直接横向比较；旧脚本名仅保留为统一流程的兼容入口。

后续 K08356/AioA 蛋白树主流程见：

```text
docs/04_STANDARD_PROTEIN_TREE_WORKFLOW.md
01_metagenomic_workflow/02_phylogeny/06_run_standard_K08356_protein_tree_workflow.sh
```

分组 derep K08356 候选 MAG 的 METABOLIC-G 代谢重建流程见：

```text
docs/05_STANDARD_METABOLIC_WORKFLOW.md
01_metagenomic_workflow/04_metabolism/METABOLIC_group_derep49/
```

该数据集包含 49 条 K08356 候选蛋白，但对应 48 个唯一 MAG；METABOLIC-G
使用每个 MAG 的完整 `*.cds.faa` 蛋白组，因此标准输入为 48 个文件。

分组参考TPM与锁定代谢矩阵的丰度加权整合见：

```text
02_figure_scripts/10_k08356_metabolism/abundance_group_reference/
```

该分析严格核验56样本×48目标MAG，但IS/AS/ES/NS使用的背景参考分别含
315/128/252/156个MAG。Clade 1–3混合参考，只作描述性展示；Clade 4的20个
MAG全部来自同一IS参考，可进行同参考的样本生境比较，但仍需站位/深度校正。

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
├── 10_k08356_metabolism/
└── 99_archive_and_debug/
```

分析输出与画图脚本的逐项对应见 [docs/02_ANALYSIS_FIGURE_CROSSWALK.md](docs/02_ANALYSIS_FIGURE_CROSSWALK.md)。

## 重要说明

1. v4 是当前推荐流程；v1–v3 和个人命令笔记只用于方法追溯，不建议直接执行。
2. K08356、AioA 和 IdrA 的功能分类必须联合 HMM、系统发育和基因邻域证据；不能只凭 KOfam 或单一 HMM 命中定义功能。
3. `MAG TPM` 表示 reads 回贴得到的宿主 MAG 丰度，不等同于分箱来源数量。
4. 宏基因组 TPM 表示基因组层面的功能潜力，不代表转录活性。
5. 四生境丰度比较以 56 个样本为统计重复，并使用同一 48-MAG 参考；48 个
   MAG 不是环境重复。
6. 多数旧 R 脚本保留了原电脑的 `setwd()`。正式运行前需改为项目路径；这些脚本已按用途归档，但没有假设用户本地数据布局并强行改写。
7. 原压缩包中的服务器地址和疑似登录口令已经脱敏。不要把密码、Token 或服务器登录信息写入 Git 仓库。

## 未纳入仓库的内容

为避免仓库臃肿和泄露本机状态，以下内容未上传：

- 本机 `R_libs/` 包目录；
- `.RData`、`.Rhistory`、`.DS_Store` 和 `__MACOSX`；
- 已生成的 PNG、临时结果和 Word 文档；
- 大型原始 reads、组装结果、数据库和 MAG 文件。

这些文件应通过 Conda、renv 或独立数据存储管理，不应提交到 Git。

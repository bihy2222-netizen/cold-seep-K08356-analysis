# 宏基因组分析执行顺序

## 0. 项目初始化与元数据

主脚本命令：

```bash
bash 01_metagenomic_workflow/00_master_pipeline/metagenome_pipeline_v4.sh init
bash 01_metagenomic_workflow/00_master_pipeline/metagenome_pipeline_v4.sh check
bash 01_metagenomic_workflow/00_master_pipeline/metagenome_pipeline_v4.sh versions
```

核心输入为 `samples.tsv`，至少核查样品名、生境、站位、岩芯、深度层、双端 reads 路径和 SRA accession。四类生境样本量应与最终统计口径保持一致。

## 1. 原始数据获取与质控

```bash
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh download
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh qc
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh fastp
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh qc-clean
```

输出：原始/过滤后 FastQC、MultiQC 报告及 clean reads。

## 2. 组装与组装质量

```bash
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh assembly
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh quast
```

输出：各样品 contigs、QUAST 统计。正式分析时记录最短 contig 阈值和组装参数。

## 3. 群落组成

```bash
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh metaphlan
```

输出：MetaPhlAn 物种丰度表。该表进入群落组成、α/β 多样性、差异类群、Venn/UpSet 等分析。

## 4. 分箱、MAG 质控与分类

```bash
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh binning
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh refinement
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh bin-stats
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh collect-bins
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh mag-qc
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh drep
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh gtdb
```

推荐证据顺序：MetaWRAP refinement → CheckM2/GUNC → 整体 dRep → GTDB-Tk。宿主分类图应使用去冗余代表 MAG 的统一分类表。

## 5. 基因预测、去冗余基因集与丰度

```bash
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh genes
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh catalog
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh salmon
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh coverm
bash 01_metagenomic_workflow/03_abundance/01_run_group_derep_coverm_tpm.sh
```

- Prodigal：contig/MAG 基因预测；
- catalog：建立统一非冗余基因集；
- Salmon：gene/contig 层 TPM；
- CoverM：MAG reads recruitment 丰度。

gene TPM 与 MAG TPM 是两条互补证据链，不应混为一个统计单位。

## 6. 功能注释

```bash
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh kofam
```

输出：KOfamScan 注释及 K08356 初始候选。C/N/S/As 功能基因丰度需要按“样品 × 基因”汇总后再进行四生境统计。

## 7. K08356/AioA/IdrA 参考整合与候选筛选

```bash
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh reference-integrate
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh k08356-check-refs
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh k08356-search
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh k08356-hmms
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh k08356-dual
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh k08356-calibrate
```

阈值校准完成前不要直接执行分类。论文作者 HMM 的 `-T 640` 不能机械套用到新建 HMM。

## 8. DMSOR 系统发育

```bash
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh k08356-tree
bash 01_metagenomic_workflow/02_phylogeny/01_run_MopB_AioA_IdrA_phylogeny.sh
bash 01_metagenomic_workflow/02_phylogeny/03_run_contig_level_K08356_phylogeny.sh
bash 01_metagenomic_workflow/02_phylogeny/04_extract_group_derep_K08356_from_kofamscan.sh
bash 01_metagenomic_workflow/02_phylogeny/05_extract_MAG_proteins_by_hmmsearch.sh
bash 01_metagenomic_workflow/02_phylogeny/06_run_standard_K08356_protein_tree_workflow.sh
```

树中应包含 AioA、IdrA、ArxA、ArrA 以及其他 DMSOR 近缘参考，避免只在 AioA–IdrA 小范围内错误定名。
后续 K08356/AioA 蛋白树以 `docs/04_STANDARD_PROTEIN_TREE_WORKFLOW.md` 为主流程说明。

其中 contig-level K08356/AioA 追溯树使用 `*.cdhit.kofamscan.best.txt`
中的 K08356 best hit，从对应 `*.cdhit.cds.faa` 抽提蛋白，先生成
`aio-A_contig_level/*.K08356.faa`，再合并为
`aio-A_contig_level/tree/all_contig_level_K08356.faa`。

分组 derep MAG 的严格 K08356 抽提则使用
`04_extract_group_derep_K08356_from_kofamscan.sh`：先合并每组
`*.kofamscan.txt`，保留第一列为 `*` 的 KOfam 可靠注释，再筛
`K08356` 并从对应 `*.cds.faa` 抽蛋白。

如果已有某个功能基因的专用 HMM，也可以用
`05_extract_MAG_proteins_by_hmmsearch.sh` 直接从 MAG 蛋白全集抽提候选。
这一路线适合 dmdA、IdrA 这类需要专用模型追踪的基因，但最终仍建议和
KOfam、系统发育、邻域证据交叉比较。

## 9. 基因邻域与最终证据表

```bash
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh neighborhood
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh neighborhood-markers
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh evidence
```

最终分类需联合：

1. HMM 得分、覆盖度和 score delta；
2. 系统发育位置与支持率；
3. 同一 contig 上的基因顺序、方向、间距和边缘截断状态。

对 IdrA 相关序列，建议区分“with DIRM-like synteny”和“without complete DIRM-like synteny”，避免把 partial neighborhood 写成完整 `idrABP1P2`。

## 10. 统计分析

```bash
DRY_RUN=0 bash .../metagenome_pipeline_v4.sh stats-template
```

推荐统计单位为样品。常用处理：

- `log10(TPM + 1)`；
- Kruskal–Wallis + Wilcoxon/Dunn；
- 多重比较使用 BH-FDR；
- Hedges’ g 与 log2FC 分开解释；
- 有站位/深度嵌套时优先考虑 mixed-effects model。

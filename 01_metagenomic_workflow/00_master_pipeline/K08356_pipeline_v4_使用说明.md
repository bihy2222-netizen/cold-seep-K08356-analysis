# K08356/AioA/IdrA 宏基因组流程 v4 使用说明

## 定位

v4 是在 v3 基础上增加“参考整合”的版本。它适合保存为当前主流程说明：

- v3 说明可作为历史归档；
- v4 脚本是当前推荐运行版本；
- HMM、树和邻域仍需联合判定，不能只凭 HMM 结果直接定义功能。

对应脚本：

```bash
metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh
```

## 核心原则

K08356/MopB 候选蛋白最终分类必须联合三类证据：

1. AioA/IdrA 双 HMM 得分、target coverage、model coverage 和 score delta；
2. 包含 AioA、IdrA、ArxA、ArrA、DmsA/DorA、NarG/NapA、PsrA/PhsA、TtrA、SerA 等 DMSOR 参考蛋白的系统发育树；
3. 候选基因上下游各 10 个 ORF 的基因邻域。

HMM 结果只能用于初筛或支持证据，不能单独证明砷氧化或碘酸盐还原功能。

## 四类推荐命名

| Clade | 推荐名称 |
|---|---|
| Clade 1 | canonical_AioA_associated |
| Clade 2 | AioA_proximal_uncharacterized_DMSOR |
| Clade 3 | IdrA_associated |
| Clade 4 | IdrA_proximal_uncharacterized_DMSOR |

## 运行前准备

设置项目路径：

```bash
export PROJECT_DIR=/home/ps/ps1/data/bihongyu/cold_seep/metagenome_project
```

初始化目录和模板：

```bash
bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh init
```

需要填写或检查：

```text
metadata/samples.tsv
metadata/reference_protein_manifest.tsv
```

v4 参考清单模板：

```text
metadata/reference_protein_manifest.v4_template.tsv
```

`reference_protein_manifest.tsv` 必须包含：

```text
sequence_id
family
role
accession
evidence_level
source
fasta_file
include_in_aioa_seed
include_in_idra_seed
include_in_calibration
include_in_tree
note
```

## 参考整合

先准备原始参考 FASTA，例如：

```text
references/raw/aioA_refs.faa
references/raw/idrA_refs.faa
references/raw/dmsor_outgroups.faa
```

然后在 `metadata/reference_protein_manifest.tsv` 中逐条记录参考序列，并标记它属于：

- AioA seed；
- IdrA seed；
- independent calibration；
- DMSOR tree reference。

运行参考整合：

```bash
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh reference-integrate
```

该步骤会生成：

```text
references/aioA_seed.faa
references/idrA_seed.faa
references/dmsor_family_references.faa
references/integrated/calibration_references.faa
references/integrated/reference_family_role_counts.tsv
references/integrated/integrated_reference_md5.tsv
references/integrated/integrated_reference_seqkit_stats.tsv
```

随后检查：

```bash
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh k08356-check-refs
```

## 推荐执行顺序

先 dry run：

```bash
bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh dry-run-plan
```

真实运行时逐步执行，不建议一键跑到底：

```bash
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh check
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh genes
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh catalog
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh kofam
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh reference-integrate
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh k08356-check-refs
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh k08356-search
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh k08356-hmms
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh k08356-dual
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh k08356-calibrate
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh k08356-tree
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh neighborhood
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh neighborhood-markers
```

检查参考得分、树和邻域后，人工填写并锁定：

```text
08_annotation/hmm_calibration/aioA_idrA_thresholds.tsv
```

只有其中 `status` 不再是 `PLACEHOLDER`，且阈值经过独立阳性/阴性参考验证后，才运行：

```bash
DRY_RUN=0 bash metagenome_reproducible_pipeline_v4_reference_integrated_20260721.sh k08356-classify
```

## HMM 阈值原则

- 论文中的 bit score 640 只适用于论文作者的原始 IdrA HMM，不能直接套用到自建模型。
- 同时检查 domain score、target coverage、model coverage 和 AioA-IdrA score delta。
- 校准集应包含独立 AioA/IdrA 阳性参考和其他 DMSOR 阴性/近缘参考。
- seed 序列和 calibration 序列应尽量分开，避免循环验证。
- 如果 IdrA seed 数量过少，模型只能作为探索性工具，不能直接作为正式分类器。

## 邻域判定

- Rieske/PF00355 加相邻至少两个 PF03150：支持 IdrA-like `idrABP1P2` 邻域。
- Rieske/PF00355 存在但相邻 PF03150 缺失：支持 AioA-like `aioBA` 邻域。
- 候选靠近 contig 边缘时，P1/P2 未检出不能解释为生物学缺失。
- Clade 3 应同时获得 IdrA 参考分支和 `idrABP1P2` 支持；否则降级为 IdrA-proximal uncharacterized。

## 归档建议

`K08356_pipeline_v3_使用说明.md` 建议保留为历史版本归档；正式项目记录使用本 v4 说明。

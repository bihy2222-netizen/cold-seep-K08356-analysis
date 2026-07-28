# K08356/AioA/IdrA 宏基因组流程 v3

## 目的

该流程用于南海冷泉宏基因组的标准处理，以及 K08356/MopB 候选蛋白的分层分类。最终分类必须联合三类证据：

1. AioA/IdrA 双 HMM 得分及覆盖度；
2. 含多类 DMSOR 参考蛋白的系统发育树；
3. 候选基因上下游各 10 个 ORF 的基因邻域。

HMM 结果只能用于初筛，不能单独证明砷氧化或碘酸盐还原功能。

## 最终四类

| Clade | 推荐名称 |
|---|---|
| Clade 1 | canonical_AioA_associated |
| Clade 2 | AioA_proximal_uncharacterized_DMSOR |
| Clade 3 | IdrA_associated |
| Clade 4 | IdrA_proximal_uncharacterized_DMSOR |

## 运行前准备

先设置项目路径：

```bash
export PROJECT_DIR=/home/ps/ps1/data/bihongyu/cold_seep/metagenome_project
```

初始化目录和模板：

```bash
bash metagenome_reproducible_pipeline_v3_20260721.sh init
```

填写：

- `metadata/samples.tsv`
- `metadata/reference_protein_manifest.tsv`

准备参考文件：

```text
references/PF00384.hmm
references/PF00355.hmm
references/PF03150.hmm
references/aioA_seed.faa
references/idrA_seed.faa
references/dmsor_family_references.faa
```

其中 `dmsor_family_references.faa` 应包含可靠的 AioA、IdrA、ArxA、ArrA、DmsA/DorA、NarG/NapA、PsrA/PhsA、TtrA、SerA 等参考序列。

## 推荐执行顺序

先做 dry run：

```bash
bash metagenome_reproducible_pipeline_v3_20260721.sh dry-run-plan
```

真实运行单步时添加 `DRY_RUN=0`：

```bash
DRY_RUN=0 bash metagenome_reproducible_pipeline_v3_20260721.sh check
DRY_RUN=0 bash metagenome_reproducible_pipeline_v3_20260721.sh genes
DRY_RUN=0 bash metagenome_reproducible_pipeline_v3_20260721.sh catalog
DRY_RUN=0 bash metagenome_reproducible_pipeline_v3_20260721.sh kofam
DRY_RUN=0 bash metagenome_reproducible_pipeline_v3_20260721.sh k08356-check-refs
DRY_RUN=0 bash metagenome_reproducible_pipeline_v3_20260721.sh k08356-search
DRY_RUN=0 bash metagenome_reproducible_pipeline_v3_20260721.sh k08356-hmms
DRY_RUN=0 bash metagenome_reproducible_pipeline_v3_20260721.sh k08356-dual
DRY_RUN=0 bash metagenome_reproducible_pipeline_v3_20260721.sh k08356-calibrate
DRY_RUN=0 bash metagenome_reproducible_pipeline_v3_20260721.sh k08356-tree
DRY_RUN=0 bash metagenome_reproducible_pipeline_v3_20260721.sh neighborhood
DRY_RUN=0 bash metagenome_reproducible_pipeline_v3_20260721.sh neighborhood-markers
```

检查参考得分、树和邻域后，人工填写并锁定：

```text
08_annotation/hmm_calibration/aioA_idrA_thresholds.tsv
```

只有其中 `status` 不再是 `PLACEHOLDER`，且阈值经过独立阳性/阴性参考验证后，才运行：

```bash
DRY_RUN=0 bash metagenome_reproducible_pipeline_v3_20260721.sh k08356-classify
```

## HMM阈值原则

- 论文中的 bit score 640 只适用于论文作者的原始 IdrA HMM，不能直接套用到自建模型。
- 同时检查 domain score、target coverage、model coverage 和 AioA-IdrA score delta。
- 校准集应包含独立的 AioA/IdrA 阳性参考和其他 DMSOR 阴性/近缘参考。
- 如果 IdrA seed 只有 IR-12 和 SCT 两条序列，模型只能用于探索，不能直接作为正式分类器。

## 邻域判定

- Rieske/PF00355 加相邻至少两个 PF03150：支持 IdrA-like `idrABP1P2` 邻域。
- Rieske/PF00355 存在但相邻 PF03150 缺失：支持 AioA-like `aioBA` 邻域。
- 候选靠近 contig 边缘时，P1/P2 未检出不能解释为生物学缺失。
- Clade 3 应同时获得 IdrA参考分支和 `idrABP1P2` 支持；否则降级为 IdrA-proximal uncharacterized。

## 重要修正

v3 在Prodigal输出后给所有基因、蛋白和GFF feature添加`样品|原始ID`前缀，避免不同样品重复出现`k141_*`导致串样；同时修复多GFF坐标合并、KOfam显著命中ID解析、多domain覆盖、阈值覆盖和非GNU awk兼容性问题。

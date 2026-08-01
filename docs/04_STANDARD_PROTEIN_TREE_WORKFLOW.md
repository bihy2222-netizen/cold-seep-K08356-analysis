# 标准蛋白树流程：分组 derep K08356/AioA

本文档作为后续跑 K08356/AioA 蛋白树的主流程。除非明确做探索性对照，优先使用这里的严格流程。

## 适用范围

- 输入对象：IS、AS、ES、NS 四个分组 derep MAG。
- 候选定义：KOfamScan 可靠注释行，即第一列为 `*`，且 KO 列为 `K08356`。
- 蛋白来源：每个分组 derep MAG 目录中的 `*.cds.faa`。
- 主要输出：严格 K08356 蛋白集合、大树、小树。

## 服务器默认路径

```text
/home/ps/ps1/data/bihongyu/cold_seep/illu/binning
```

四个分组 derep 目录：

```text
IS/IS_50_10_dRep/dereplicated_genomes
AS/AS_50_10_dRep/dereplicated_genomes
ES/ES_50_10_dRep/dereplicated_genomes
NS/NS_50_10_dRep/dereplicated_genomes
```

每个目录中应有：

```text
*.kofamscan.txt
*.cds.faa
```

## Step 1：严格抽提 K08356 蛋白

使用脚本：

```text
01_metagenomic_workflow/02_phylogeny/04_extract_group_derep_K08356_from_kofamscan.sh
```

核心逻辑：

```bash
cat *.kofamscan.txt > all_bin.kofamscan.txt
awk '$1 == "*"' all_bin.kofamscan.txt > all_bin.kofamscan_filtered.txt
awk '$1 == "*" && $3 == "K08356" {print $2}' all_bin.kofamscan.txt > K08356.ids
seqkit grep -f K08356.ids *.cds.faa > K08356.faa
```

整理版脚本会额外输出分组、`bin_id`、`gene_id`、threshold、score、E-value 和 definition。

注意：从 MAG 的 `*.cds.faa` 抽序列时不要加 `seqkit grep -n`，因为 FASTA header 后面带 Prodigal 坐标；默认按 sequence ID 匹配才能抽到。

当前严格抽提结果：

```text
IS  38
AS  3
ES  6
NS  2
ALL 49
```

结果目录：

```text
/home/ps/ps1/data/bihongyu/cold_seep/illu/binning/K08356_four_group_annotation/strict_kofam_best
```

关键文件：

```text
all_groups.K08356_hits.tsv
all_groups.K08356.faa
```

## Step 2：构建大树

大树用于 DMSOR 范围定位，参考序列来自旧整体树：

```text
/home/ps/ps1/data/bihongyu/cold_seep/ASproteintree/20260323/all_proteins.faa
```

该文件包含旧整体 derep 的 48 条 MAG 和 261 条 DMSOR 参考。标准流程会先去掉旧 48 条 MAG，只保留 261 条参考，再合并新的 49 条分组 derep K08356：

```text
49 条分组 derep K08356 + 261 条 DMSOR 参考 = 310 条
```

输出目录：

```text
strict_kofam_best/tree_big_DMSOR_group_derep49_refs261
```

## Step 3：构建小树

小树用于重点比较 AioA、IdrA 和 Unknown clade，参考序列来自：

```text
/home/ps/ps1/data/bihongyu/cold_seep/ASproteintree/MAG+idratree4-small-smalltree/all_proteins.faa
```

该文件包含旧 48 条 MAG 和 136 条 AioA/IdrA/Unknown 参考。标准流程会去掉旧 48 条 MAG，只保留 136 条参考，再合并新的 49 条：

```text
49 条分组 derep K08356 + 136 条 AioA/IdrA/Unknown 参考 = 185 条
```

输出目录：

```text
strict_kofam_best/tree_small_AioA_IdrA_Unknown_group_derep49_refs136
```

## Step 4：MAFFT、trimAl、IQ-TREE 参数

两棵树统一使用：

```bash
mafft --auto --quiet --anysymbol input.fasta > aligned.fasta

trimal \
  -in aligned.fasta \
  -out trimmed.fasta \
  -gt 0.8 \
  -st 0.001 \
  -cons 60

nohup iqtree \
  -s trimmed.fasta \
  -m MFP \
  -B 1000 \
  --alrt 1000 \
  -T AUTO \
  > iqtree.log \
  2>&1 &
```

非标准氨基酸统一替换为 `X`。

## 一键标准流程

使用脚本：

```text
01_metagenomic_workflow/02_phylogeny/06_run_standard_K08356_protein_tree_workflow.sh
```

服务器上运行：

```bash
cd /home/ps/ps1/data/bihongyu/cold_seep/illu/binning
bash /path/to/06_run_standard_K08356_protein_tree_workflow.sh
```

如果只准备输入和 trim 后比对，不启动 IQ-TREE：

```bash
RUN_IQTREE=0 bash /path/to/06_run_standard_K08356_protein_tree_workflow.sh
```

如果只跑小树：

```bash
RUN_BIG_TREE=0 bash /path/to/06_run_standard_K08356_protein_tree_workflow.sh
```

如果已经抽好了 49 条，不想重新抽提：

```bash
RUN_EXTRACT=0 bash /path/to/06_run_standard_K08356_protein_tree_workflow.sh
```

## 方法选择

K08356/AioA 主结果使用 KOfam `*` 严格法。之前的宽松 `score >= 100`、`evalue <= 1e-5`、长度筛选版本只作为探索性对照，不作为最终候选定义。

IdrA、dmdA 等 KO 覆盖不足或需要专用模型的基因，可以使用 HMM 直搜脚本：

```text
01_metagenomic_workflow/02_phylogeny/05_extract_MAG_proteins_by_hmmsearch.sh
```

最终功能判断仍需联合系统发育位置、HMM 证据和基因邻域。

## 小树 49 条的 640 阈值和基因岛补充流程

本仓库另整理了小树 `tree_small_AioA_IdrA_Unknown_group_derep49_refs136` 的标准化补充模块：

```text
01_metagenomic_workflow/02_phylogeny/small_tree_derep49_standardized_workflow
```

该模块包含 2026-08-01 在服务器 `/home/ps/ps1/data/bihongyu/cold_seep` 下重跑的 HMM `-T 640` 结果、按 iTOL SVG 顺序整理的基因岛箭头图、以及 IdrA-PP/gene-island 预测表。复跑本地图表：

```bash
bash 01_metagenomic_workflow/02_phylogeny/small_tree_derep49_standardized_workflow/run_workflow.sh
```

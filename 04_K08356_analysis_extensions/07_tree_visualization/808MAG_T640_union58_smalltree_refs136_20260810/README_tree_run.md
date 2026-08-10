# 808 derepMAG 三模型 T640 并集 58 条候选蛋白小树

## 结论与口径

本目录是 **808 个 derepMAG 完整蛋白库的三模型 T640 高分候选并集**系统树结果。候选集合包含 58 条蛋白、来自 56 个 MAG；它们是 IdrA-specific、AioA-specific 和 combined 三模型高分候选的并集，**不能把58条全部称为已确认 IdrA**。

本树与旧“49条人工分类蛋白 + 136条参考”的核心 K08356 小树不同，也与720条T100宽松候选及contig候选无关。本次输入为58条候选 + 原小树同一套136条参考，共194条唯一FASTA ID。

系统树落点只提供亲缘关系证据，不能单独替代功能鉴定。尤其在当前没有NCBI HitTable的情况下，不应凭BLAST名称作最终功能判定；后续仍需检查基因邻域和其他独立证据。

## 树构建与验收

- 比对：MAFFT 7.525，`--auto --thread 16`
- 修剪：trimAl 1.4.rev15，`-automated1`
- 建树：IQ-TREE 2.3.4，`-m MFP -B 1000 --alrt 1000 -T AUTO`
- BIC最佳模型：`LG+R6`
- 支持度：1000次ultrafast bootstrap完成；1000次SH-like aLRT完成
- 最终treefile和contree均含194个tips
- FASTA首字段ID无重复
- 完成时间：2026-08-10T11:39:33+08:00

输入SHA-256：

```text
0141c95e1319e990bac69d1d9ee0bc81d2f6aacc235f215d34df7f254a790434  T640_union58.faa
199df840a52014ecb86d56782ca064a0e2557e150c9cb5b7ab993e6bf3f69125  smalltree_refs136.faa
999c785680fe7e033da2195791bf37918a2a479ed95214d71595790c50549329  MAG_T640_union58_plus_refs136.faa
```

## 主要文件

- `02_tree/MAG_T640_union58_plus_refs136.treefile`：最终ML树，节点标签为SH-aLRT/UFBoot
- `02_tree/MAG_T640_union58_plus_refs136.contree`：1000次UFBoot共识树
- `02_tree/MAG_T640_union58_plus_refs136.iqtree`：完整模型、参数和结果报告
- `04_tables/union58_three_model_score_coverage_matrix.tsv`：58条候选的三模型score/coverage原表
- `04_tables/candidate58_tree_placement_summary.tsv`：候选ID、MAG、三模型score/coverage、最近参考及树上距离
- `05_itol/MAG_T640_union58_plus_refs136_for_iTOL.treefile`：iTOL上传树
- `05_itol/itol_COLORSTRIP_candidate_vs_reference_groups.txt`：候选/参考组颜色条
- `06_figures/MAG_T640_union58_plus_refs136_tree.{pdf,svg,png}`：Times New Roman字体概览图
- `03_audit/`：计数、输入哈希、重复ID及运行摘要

## 落点表的保守解释

最近参考按最终treefile的patristic distance计算。58条候选中，最近参考类别为 IdrA reference 47条、AioA reference 3条、Unknown-clade reference 8条。这只是“最近参考类别”，不是自动功能命名；同一候选所在更深层分支、节点支持度、序列质量、HMM证据和基因邻域需要联合解释。

运行期间IQ-TREE对部分经验模型报告“state frequencies归一化”警告；作业正常完成，最终报告、树和支持度文件均完整。绘图时Fontconfig曾提示缓存目录不可写，但Times New Roman已成功嵌入PDF，不影响树数据。

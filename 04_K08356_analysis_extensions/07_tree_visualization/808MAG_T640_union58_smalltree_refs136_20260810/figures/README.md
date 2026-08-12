# 58-candidate four-clade tree with gene-island neighborhoods

本图基于最新的 `vj6w2EeEGYAiMv6aDe8mIg.svg` 系统树，将58条候选逐条与三模型HMM结果、最近参考落点和基因邻域证据整合，并在右侧绘制以候选A基因为中心的±10 CDS基因岛。

## 四类与颜色

- Canonical AioA clade：3条，绿色 `#009E73`
- Unknown DMSOR clade：8条，橙色 `#D55E00`
- IdrA-associated clade：27条，蓝色 `#0072B2`
- IdrA clade：20条，紫色 `#7B3294`

其中48条与旧严格QC证据表重叠，沿用旧表分类；10条新T640候选中，8条根据IdrA参考落点且缺少完整B+2P邻域证据暂列为IdrA-associated，2条根据Unknown-clade参考落点暂列为Unknown DMSOR。新10条分类均是保守工作分类，不是最终功能确认。

## Protein subfamily 图例

- Molybdopterin oxidoreductase：目标AioA/IdrA-related大亚基
- Rieske [2Fe-2S] small subunit：IdrB/AioB-related小亚基
- IdrP-like accessory protein
- Cytochrome c / peroxidase
- Redox / oxidoreductase
- Transporter
- Arsenic resistance protein
- Sulfur oxidation protein
- Other annotated protein
- Non-conserved protein

命名和配色参考蛋白亚家族式基因岛图，但只使用当前证据能够支持的类别。旧严格QC邻域中的KO/product分类得到保留；缺少可靠注释的ORF统一标为`Non-conserved protein`，不按图形位置推测具体功能。

基因邻域来自每个MAG的Prodigal蛋白文件，提取目标基因上下游各10个CDS；B/P-like标注由既有IdrB/AioB与IdrP-like联合参考库的DIAMOND搜索支持。图中基因位置和方向来自Prodigal坐标。

## 主要文件

- `MAG_T640_union58_four_clades_gene_islands.svg`：可编辑主图
- `MAG_T640_union58_four_clades_gene_islands.png`：高分辨率预览
- `../04_tables/candidate58_four_clade_integrated_evidence.tsv`：58条四类整合证据表
- `../04_neighborhood/union58_gene_island_manifest.tsv`：671个邻域基因绘图清单
- `../04_neighborhood/union58_neighborhood_long.tsv`：原始±10 CDS邻域表

## 解释限制

四类为当前系统树和邻域证据下的工作分类。系统树位置或单次相似性命中均不能单独完成最终功能鉴定；尤其新10条候选仍建议结合更完整的功能注释、基因岛边界、蛋白结构域和必要的人工审阅进一步确认。

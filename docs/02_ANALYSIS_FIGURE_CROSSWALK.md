# 分析步骤与画图脚本对应表

| 顺序 | 上游分析/输入 | 推荐图形 | 脚本目录 | 状态 |
|---|---|---|---|---|
| Fig.01 | 站位经纬度、GEBCO | 采样地图 | `02_figure_scripts/01_sampling_map/` | 有 GMT 与 R 两版 |
| Fig.02 | MetaPhlAn 合并丰度表 | 四生境群落组成柱状图 | `02_figure_scripts/02_community_composition/` | 以“带分组色块”版为主 |
| Fig.03 | MetaPhlAn 物种/门水平丰度 | Shannon、Simpson、Chao1 | `02_figure_scripts/03_alpha_diversity/` | 优先 FDR corrected 版 |
| Fig.04 | MetaPhlAn 丰度表 + 分组表 | NMDS、PCoA、PERMANOVA | `02_figure_scripts/04_beta_diversity/` | 优先 IS/AS/ES/NS 升级版 |
| Fig.05 | 物种丰度与分组 | LEfSe、LDA、STAMP、差异热图 | `02_figure_scripts/05_differential_taxa/` | 多版并存，需按论文面板选择 |
| Fig.06 | 各生境 presence/absence | Venn、UpSet | `02_figure_scripts/06_shared_unique_taxa/` | `upsetfinal.R` 为推荐版 |
| Fig.07 | MAG 分类表、功能基因宿主表 | 宿主门/纲/目组成 | `02_figure_scripts/07_host_taxonomy/` | 优先“清洗数据加柱状图无黑线”版 |
| Fig.08 | 环境因子 + gene/MAG TPM | Spearman/Mantel 倒三角热图 | `02_figure_scripts/08_environment_association/` | 旧模板含示例数据，需换成真实 IS 数据 |
| Fig.09 | C/N/S/As TPM 或 MAG 模块 | 生境嵌套网络、MAG 网络 | `02_figure_scripts/09_networks/` | 可输出 SVG |
| Fig.10 | 48 个分组 derep K08356 MAG 的 METABOLIC-G 结果 | 四生境 × 四 clade 代谢气泡图 | `02_figure_scripts/10_k08356_metabolism/` | 已按 49 条序列/48 个 MAG 更新并验证 |

## K08356 专题图的对应关系

| 分析输出 | 应绘制的图 | 当前压缩包情况 |
|---|---|---|
| 56 样品 K08356 gene TPM + relative abundance | 顶部 56 样品组成条 + TPM 主面板 | 未找到最近讨论的 `make_reference_style_K08356_figure.R` |
| 四分支宿主 MAG TPM | 四生境 × 四分支丰度图 | 未找到 `make_contig_MAG_TPM_combo_figure.R` 的当前版 |
| DMSOR 系统发育树 | AioA/IdrA/unknown/DIRM-like clade 树 | 含建树 Shell，未含最终树注释绘图脚本 |
| 严格 QC 邻域证据表 | 基因箭头图 | 未找到最终箭头图绘图脚本 |
| MAG module/function reconstruction score and prevalence | 四生境 × 四 clade 代谢气泡图 | 已补入 METABOLIC-G v4.0 新结果、脚本、SVG 和绘图数据 |
| IS 营养盐数据 | Spearman/Mantel 图 | 含早期 Mantel 模板，未确认是否为最终 21 样品版 |

上述缺失并不表示分析不存在，只表示这些最近使用的脚本不在本次 `Rstudio.zip` 内。后续应将最新版本补入相应目录，而不是用旧模板覆盖。

## 统计解释约束

- Hedges’ g 图：横轴是标准化效应量，不称为火山图。
- Volcano：横轴为 log2FC，纵轴为 `−log10(FDR)`。
- 气泡大小若表示原始平均 TPM，应在图例中明确，不能称为标准化丰度。
- `mean gene-set coverage` 与 `MAG prevalence` 是两个指标，需要分别说明。
- 未经转录组或酶学验证时，使用“丰度”“潜力”“associated”，不使用“活跃”“主导反应”。

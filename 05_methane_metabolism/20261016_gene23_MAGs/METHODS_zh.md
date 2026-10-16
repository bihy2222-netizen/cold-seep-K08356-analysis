# 冷泉甲烷循环的代谢分析：详细方法与统计口径

## 1. 分析范围与数据层级

本次归档的是已有候选筛查结果的23KO重整与作图，不在本目录重新运行蛋白预测、HMM搜索、GTDB、CheckM2、dRep或参考蛋白树。上游1449个原始MAG用于候选筛查，808个全局去冗余代表用于原代表注释汇总；二者不是同一个统计母集。经上游证据审查选出的16个mcrA阳性代表中，13个属于原808代表，另3个是同簇补充阳性MAG。

原808代表的分类、质量和dRep簇信息见tables/cold_seep_all808_taxonomy_quality.tsv；其他必要上游输入归档于inputs/。重整脚本仅对这些证据表进行连接、筛选和汇总，不把上游分类结果当作本次新增树证据。

## 2. 目标基因及名称

| 功能组 | KO | 图中名称 |
|---|---|---|
| 甲烷生成或活化 | K00399、K00401、K00402 | mcrA、mcrB、mcrG |
| 甲基转移 | K00577–K00584 | mtrA–mtrH |
| 一碳转化 | K00320、K00319、K01499、K00672 | mer、mtd、mch、ftr |
| 甲酰甲烷呋喃相关转化/电子传递 | K00200、K00201、K00202、K00203、K00205 | fwdA/fmdA、fwdB/fmdB、fwdC/fmdC、fwdD/fmdD、fwdF/fmdF |
| 辅酶再生与电子传递 | K03388、K03389、K03390 | hdrA2、hdrB2、hdrC2 |

精确一对一映射、功能说明与KEGG链接见tables/KO_gene_mapping.tsv。本批次不纳入K00204。斜线表示同一KO不能在此分辨Fwd/Fmd，并非检出两条独立蛋白。不得把早期mcrA sequence groups a–f与ANME-2a/b/c名称混用。

## 3. 接受命中的规则

prepare.py依次读取四类证据：原代表的严格KO结果；screen_status等于strict的MCR HMM结果；3个补充代表已完成严格筛选的通路HMM结果；上游已接受、且其他证据表尚未包含的McrA记录。最后一类保留McrA_evidence_tier而不伪造得分。原代表同一MAG/蛋白/KO已有KO证据时，HMM不重复增加记录；补充来源单独标识。

本阶段没有新设统一HMM分值阈值，也没有把所有relaxed候选一并接受。上游表保留score、阈值或覆盖度等可用字段供追溯，无法由本快照完全重演原始HMM搜索。蛋白-KO明细以MAG ID、protein ID和KO三元组去重。不同KO即使指向同一蛋白也保留证据记录，因此“蛋白-KO记录数”与“独立蛋白数”概念不同。

## 4. MAG携带与宿主分类

携带者定义为至少有一个目标KO的接受命中；此定义不要求mcrA或完整23KO组合。每个MAG的gene_names及KO_list按CH4→CO2概念途径顺序列出。bin_filename由MAG ID附加.fasta生成，是对应文件名标签，本归档未包含bin序列文件，也未重新核对文件是否存在。

MAG×基因矩阵为二元检出：一个MAG至少一条接受ORF对应某KO则为1，否则为0。ORF计数为该MAG/KO下独立protein ID的集合大小，不直接等同真实基因拷贝数。MAG不完整或蛋白未注释均可能造成未检出。

原808代表使用其自身GTDB分类。3个补充MAG使用同一dRep簇原代表的分类代理，taxonomy_evidence明确标记，质量不推算而保留NA。门类缺失以Unclassified显示；“门类标签数”不是对全部标签进行新分类学认证。347条携带记录对应344个dRep簇，补充记录不得当作3个新增独立物种簇。

## 5. mcrA阳性代表与分支

重点图沿用inputs/tables/mcrA_positive_dereplicated_MAGs.tsv和selected_MAG_McrA_classification.tsv。16个阳性代表包括ANME-1 2个、ANME-2a 2个、ANME-2c 2个、未解析8个、产甲烷近缘未解析2个。未解析类别不强行拆为2a、2b、2c，不用共享代谢基因把细菌或产甲烷近缘类群自动归为ANME。

本次无树图不显示节点支持度，也不构造不存在的系统发育关系；分支标签是沿用上游结果的分组注释。后续严格的ANME分型仍须结合参考蛋白树、节点支持度、宿主分类及代谢背景。

## 6. 四生境丰度

样品次序以tables/sample_metadata.tsv为准，共NS 8、IS 21、AS 14、ES 13个样品。保留实际样品名，不以旧示例分组替换。三种数量独立保存：MAG或dRep簇数量、蛋白ORF/命中数量、样品丰度。

全群落KO丰度来自inputs/local_audit/local_community_pathway_KO_56samples.tsv，输出community_gene_abundance.tsv，仅更换行标签及排列，保留来源数值与量纲。阳性MAG基因TPM来自inputs/tables/gene_TPM_56samples_mcrA_positive_MAGs.tsv，输出mcrA_MAG_gene_abundance.tsv，仅包含13个直接定量代表的基因贡献。逐MAG丰度见selected_MAG_TPM_56samples.tsv。3个补充MAG的丰度NA不使用同簇其他MAG的丰度代替。

Salmon索引缺失3条已注释ORF会影响mtrB/mtrF汇总，相关值保留NA，详见RESULTS_zh.md。主携带表对原代表显示的丰度可用状态描述其既有MAG定量背景；该表并未输出全部344个携带MAG的逐样品丰度，不应把其余共享基因携带者加入阳性MAG丰度矩阵。

## 7. 两版无树检出图与丰度图

图01参考用户所给图片的复合体排列：Mcr A/G/B、Mtr H/G/F/A/B/C/D/E、Mer、Mtd、Mch、Ftr、Fwd/Fmd C/A/B/D/F、Hdr2。图02按概念CH4→CO2路径：Mcr A/B/G、Mtr A–H、Mer、Mtd、Mch、Ftr、Fwd/Fmd A/B/C/D/F，最后独立列出Hdr2。两图仅改变列顺序，不改变检出值。

彩色格为检出，灰色斜线为未检出；颜色对应上游分支分组，未解析类别合并为灰蓝色。完整度为可用CheckM2值的整数显示，代理MAG为NA。顶部箭头是碳转化示意，不描述所有电子、辅酶或离子化学计量，不把Hdr列为末端碳反应。未推断基因邻接、操纵子或额外旁系同源拷贝。

图03为全群落23基因×56样品，图04为13个阳性代表的基因贡献，两者使用各自色标的log10(x+1)，不作跨图绝对色值比较；灰色为缺失定量，零值保留。图05为全部347条携带记录的分页检出矩阵。输出PNG、PDF、SVG（分页图仅PDF）。

## 8. 复现与核验

在当前目录运行python3 prepare.py重整表格，再运行python3 verify.py核验关键统计和集合一致性；安装requirements.txt中的绘图库后运行python3 plot.py重绘。prepare.py需要已归档的808代表分类质量表；不从远程服务器动态抓取。Excel是导出快照，不由上述命令重建。依赖未锁定版本，不能保证不同版本的字体布局逐像素一致。

关键核验：仅23KO且不含K00204；808个原代表；344个原携带代表与3个补充记录；347条记录、344个独立簇；16个重点代表与各表ID一致；56样品与四生境计数；蛋白-KO记录数1876；所有携带记录至少检出一个目标基因。manifest_sha256.json用于归档文件完整性校验。

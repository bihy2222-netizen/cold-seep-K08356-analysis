# 表格、字段和缺失值说明

所有TSV第一行为字段名；分隔符为制表符。MAG_ID和protein_ID视为字符串。NA表示数据不可用，不能替换为零。0/1仅用于检出，TPM中的0表示来源定量为0。

## KO_gene_mapping.tsv

KO、gene_name、functional_group、function、CH4_CO2_order、reference_order、KEGG_source。

数据行数：23；字段数：7。

`KO`, `gene_name`, `functional_group`, `function`, `CH4_CO2_order`, `reference_order`, `KEGG_source`

## MAG_gene_ORF_counts.tsv

每MAG/基因的独立protein ID数，未经装配邻域复核，不等同真实拷贝数。

数据行数：347；字段数：24。

`MAG_ID`, `mcrA`, `mcrB`, `mcrG`, `mtrA`, `mtrB`, `mtrC`, `mtrD`, `mtrE`, `mtrF`, `mtrG`, `mtrH`, `mer`, `mtd`, `mch`, `ftr`, `fwdA/fmdA`, `fwdB/fmdB`, `fwdC/fmdC`, `fwdD/fmdD`, `fwdF/fmdF`, `hdrA2`, `hdrB2`, `hdrC2`

## MAG_taxonomy.tsv

每MAG的域、门、纲、目、科、属、种，full_GTDB_taxonomy、taxonomy_evidence及可用完整度/污染度。代理分类不是本MAG自身GTDB结果。

数据行数：347；字段数：13。

`MAG_ID`, `domain`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `full_GTDB_taxonomy`, `taxonomy_evidence`, `completeness_percent`, `contamination_percent`, `cohort`

## all_carrier_detection_CH4_CO2.tsv

347条全部携带记录的基因检出矩阵，基因列为0/1；reference_style与CH4_CO2使用同一检出值、不同列顺序。

数据行数：347；字段数：26。

`MAG_ID`, `domain`, `phylum`, `mcrA`, `mcrB`, `mcrG`, `mtrA`, `mtrB`, `mtrC`, `mtrD`, `mtrE`, `mtrF`, `mtrG`, `mtrH`, `mer`, `mtd`, `mch`, `ftr`, `fwdA/fmdA`, `fwdB/fmdB`, `fwdC/fmdC`, `fwdD/fmdD`, `fwdF/fmdF`, `hdrA2`, `hdrB2`, `hdrC2`

## all_carrier_detection_reference_style.tsv

347条全部携带记录的基因检出矩阵，基因列为0/1；reference_style与CH4_CO2使用同一检出值、不同列顺序。

数据行数：347；字段数：26。

`MAG_ID`, `domain`, `phylum`, `mcrA`, `mcrG`, `mcrB`, `mtrH`, `mtrG`, `mtrF`, `mtrA`, `mtrB`, `mtrC`, `mtrD`, `mtrE`, `mer`, `mtd`, `mch`, `ftr`, `fwdC/fmdC`, `fwdA/fmdA`, `fwdB/fmdB`, `fwdD/fmdD`, `fwdF/fmdF`, `hdrA2`, `hdrB2`, `hdrC2`

## carrying_MAGs.tsv

347条携带记录；MAG_ID、bin_filename、cohort、dRep_cluster、mcrA_positive_selected（0/1）、target_gene_types（不同KO种类数）、gene_names、KO_list、protein_records（独立protein ID数）、abundance_status。

数据行数：347；字段数：10。

`MAG_ID`, `bin_filename`, `cohort`, `dRep_cluster`, `mcrA_positive_selected`, `target_gene_types`, `gene_names`, `KO_list`, `protein_records`, `abundance_status`

## cold_seep_all808_taxonomy_quality.tsv

808原代表的MAG_ID、classification、completeness、contamination、dRep_cluster；表格重整的分类质量输入。

数据行数：808；字段数：5。

`MAG_ID`, `classification`, `completeness`, `contamination`, `dRep_cluster`

## community_gene_abundance.tsv

23基因×56样品；gene_name、KO及实际样品名列，原始来源KO丰度量纲。

数据行数：23；字段数：58。

`gene_name`, `KO`, `NS_0-6`, `R2111_N300_0-10`, `R2111_N500_0-10`, `R2111_S300_0-10`, `R2111_S500_0-10`, `S3_0-3`, `S3_6-9`, `S3_9-12`, `SY365BB-0-4`, `SY365BB-4-8`, `SY365BB-8-12`, `SY366YB-0-4`, `SY366YB-4-8`, `SY366YB-8-12`, `SY366YW-0-4`, `SY366YW-4-8`, `SY366YW-8-12`, `SY368YW-0-4`, `SY368YW-4-8`, `SY368YW-8-12`, `SY456YB-0-4`, `SY456YB-4-8`, `SY456YB-8-12`, `SY457BB-0-4`, `SY457BB-4-8`, `SY457BB-8-12`, `SY459WG-0-4`, `SY459WG-4-8`, `SY459WG-8-12`, `S1_0-3`, `S1_6-9`, `S1_9-12`, `S2_0-3`, `S2_12-15`, `S2_3-6`, `S4_12-15`, `S4_9-12`, `SQ_58_0-4`, `SQ_58_4-8`, `SQ_58_8-12`, `SQ_81_0-4`, `SQ_81_4-8`, `SQ_81_8-12`, `C1_0-6`, `C1_12-18`, `C1_6-12`, `C2_0-6`, `C2_12-18`, `C2_6-12`, `C3_0-6`, `C3_12-18`, `C3_6-12`, `ES_2_0-6`, `S13_0-2`, `S14_4-6`, `S15_8-10`

## mcrA_MAG_gene_abundance.tsv

23基因×56样品；仅13个直接定量阳性代表的基因TPM贡献；受缺失索引影响的汇总保留NA。

数据行数：23；字段数：58。

`gene_name`, `KO`, `NS_0-6`, `R2111_N300_0-10`, `R2111_N500_0-10`, `R2111_S300_0-10`, `R2111_S500_0-10`, `S3_0-3`, `S3_6-9`, `S3_9-12`, `SY365BB-0-4`, `SY365BB-4-8`, `SY365BB-8-12`, `SY366YB-0-4`, `SY366YB-4-8`, `SY366YB-8-12`, `SY366YW-0-4`, `SY366YW-4-8`, `SY366YW-8-12`, `SY368YW-0-4`, `SY368YW-4-8`, `SY368YW-8-12`, `SY456YB-0-4`, `SY456YB-4-8`, `SY456YB-8-12`, `SY457BB-0-4`, `SY457BB-4-8`, `SY457BB-8-12`, `SY459WG-0-4`, `SY459WG-4-8`, `SY459WG-8-12`, `S1_0-3`, `S1_6-9`, `S1_9-12`, `S2_0-3`, `S2_12-15`, `S2_3-6`, `S4_12-15`, `S4_9-12`, `SQ_58_0-4`, `SQ_58_4-8`, `SQ_58_8-12`, `SQ_81_0-4`, `SQ_81_4-8`, `SQ_81_8-12`, `C1_0-6`, `C1_12-18`, `C1_6-12`, `C2_0-6`, `C2_12-18`, `C2_6-12`, `C3_0-6`, `C3_12-18`, `C3_6-12`, `ES_2_0-6`, `S13_0-2`, `S14_4-6`, `S15_8-10`

## mcrA_positive_detection_CH4_CO2.tsv

16个mcrA阳性重点代表的基因检出矩阵，基因列为0/1；reference_style与CH4_CO2使用同一检出值、不同列顺序。

数据行数：16；字段数：27。

`MAG_ID`, `clade`, `phylum`, `completeness_percent`, `mcrA`, `mcrB`, `mcrG`, `mtrA`, `mtrB`, `mtrC`, `mtrD`, `mtrE`, `mtrF`, `mtrG`, `mtrH`, `mer`, `mtd`, `mch`, `ftr`, `fwdA/fmdA`, `fwdB/fmdB`, `fwdC/fmdC`, `fwdD/fmdD`, `fwdF/fmdF`, `hdrA2`, `hdrB2`, `hdrC2`

## mcrA_positive_detection_reference_style.tsv

16个mcrA阳性重点代表的基因检出矩阵，基因列为0/1；reference_style与CH4_CO2使用同一检出值、不同列顺序。

数据行数：16；字段数：27。

`MAG_ID`, `clade`, `phylum`, `completeness_percent`, `mcrA`, `mcrG`, `mcrB`, `mtrH`, `mtrG`, `mtrF`, `mtrA`, `mtrB`, `mtrC`, `mtrD`, `mtrE`, `mer`, `mtd`, `mch`, `ftr`, `fwdC/fmdC`, `fwdA/fmdA`, `fwdB/fmdB`, `fwdD/fmdD`, `fwdF/fmdF`, `hdrA2`, `hdrB2`, `hdrC2`

## phylum_summary.tsv

按域/门分别统计原代表、补充记录、独立簇、mcrA阳性代表数。补充记录不能直接加作独立簇。

数据行数：43；字段数：6。

`domain`, `phylum`, `original808_carrier_MAGs`, `supplement_carrier_MAGs`, `unique_dRep_clusters`, `selected_mcrA_positive_MAGs`

## protein_KO_hits.tsv

MAG_ID、protein_ID、KO、gene_name、evidence、score、E_value、cohort；每MAG/蛋白/KO一行。

数据行数：1876；字段数：8。

`MAG_ID`, `protein_ID`, `KO`, `gene_name`, `evidence`, `score`, `E_value`, `cohort`

## sample_metadata.tsv

56样品的sample、habitat、site；作为全部定量表和分组图的样品顺序依据。

数据行数：56；字段数：3。

`sample`, `habitat`, `site`

## selected_MAG_TPM_56samples.tsv

16个阳性MAG的clade、abundance_status及56样品MAG TPM；3个补充MAG无直接TPM。

数据行数：16；字段数：59。

`MAG_ID`, `clade`, `abundance_status`, `NS_0-6`, `R2111_N300_0-10`, `R2111_N500_0-10`, `R2111_S300_0-10`, `R2111_S500_0-10`, `S3_0-3`, `S3_6-9`, `S3_9-12`, `SY365BB-0-4`, `SY365BB-4-8`, `SY365BB-8-12`, `SY366YB-0-4`, `SY366YB-4-8`, `SY366YB-8-12`, `SY366YW-0-4`, `SY366YW-4-8`, `SY366YW-8-12`, `SY368YW-0-4`, `SY368YW-4-8`, `SY368YW-8-12`, `SY456YB-0-4`, `SY456YB-4-8`, `SY456YB-8-12`, `SY457BB-0-4`, `SY457BB-4-8`, `SY457BB-8-12`, `SY459WG-0-4`, `SY459WG-4-8`, `SY459WG-8-12`, `S1_0-3`, `S1_6-9`, `S1_9-12`, `S2_0-3`, `S2_12-15`, `S2_3-6`, `S4_12-15`, `S4_9-12`, `SQ_58_0-4`, `SQ_58_4-8`, `SQ_58_8-12`, `SQ_81_0-4`, `SQ_81_4-8`, `SQ_81_8-12`, `C1_0-6`, `C1_12-18`, `C1_6-12`, `C2_0-6`, `C2_12-18`, `C2_6-12`, `C3_0-6`, `C3_12-18`, `C3_6-12`, `ES_2_0-6`, `S13_0-2`, `S14_4-6`, `S15_8-10`

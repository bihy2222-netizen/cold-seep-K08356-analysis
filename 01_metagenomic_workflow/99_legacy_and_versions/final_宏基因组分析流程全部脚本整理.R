# ============================================================
# 宏基因组分析流程全过程整理版
# ============================================================
#
# 用途：
#   这是一个“流程记录 + 命令模板”文件，不是直接运行的 R 分析脚本。
#   所有服务器/conda/shell 命令均以注释形式保存，避免在 RStudio 中误点 Source 后
#   直接执行下载、删除、移动、建树等耗时或危险操作。
#
# 建议用法：
#   1. 在 RStudio 中按章节查看。
#   2. 将需要的命令复制到 Linux 服务器终端运行。
#   3. 先替换 PROJECT_DIR、THREADS、SAMPLE、输入文件名等占位符。
#
# 本次标准化修正：
#   - conde/condla/booconda/fasta-dump/cal/rawdedta 等拼写错误已改正。
#   - Shell 多行命令统一使用反斜杠 `\` 续行。
#   - 将危险命令 rm/mv 默认保留为注释模板，运行前需手动确认路径。
#   - 统一样本命名示例为 SAMPLE_1.fastq.gz / SAMPLE_2.fastq.gz。
#   - 明确 MetaPhlAn、assembly、binning、gene prediction、annotation、MAG tree 的顺序。
#
# ============================================================
# 0. 推荐项目目录结构
# ============================================================
#
# PROJECT_DIR/
#   00_raw/                 # 原始 fastq/sra
#   01_qc/                  # FastQC/MultiQC
#   02_clean/               # 去接头/质控后的 reads，可选
#   03_assembly/            # MEGAHIT 拼接结果
#   04_metaphlan/           # 物种组成
#   05_binning/             # MetaWRAP binning/refinement
#   06_genes/               # Prodigal/CD-HIT 非冗余基因集
#   07_salmon/              # 基因丰度定量
#   08_annotation/          # KOfamScan/HMM 等功能注释
#   09_MAGs/                # dRep/CheckM2/GTDB-Tk/MAG 树
#   logs/
#
# 常用变量示例，复制到 bash 中执行：
#
#   PROJECT_DIR=/home/ps/ps1/data/bihongyu/cold_seep/illu
#   THREADS=56
#   SAMPLE=SY457YW-8-12
#
# ============================================================
# 1. 原始数据质控：FastQC + MultiQC
# ============================================================
#
#   conda activate fastqc
#   mkdir -p "$PROJECT_DIR/01_qc/fastqc"
#   nohup fastqc "$PROJECT_DIR"/00_raw/*.fastq.gz \
#     -o "$PROJECT_DIR/01_qc/fastqc" \
#     -t "$THREADS" \
#     > "$PROJECT_DIR/logs/fastqc.log" 2>&1 &
#
#   conda activate multiqc
#   multiqc "$PROJECT_DIR/01_qc/fastqc" \
#     -o "$PROJECT_DIR/01_qc/multiqc"
#
# 检查：
#   - 每个样本是否有 *_fastqc.html。
#   - multiqc_report.html 中接头、低质量碱基、GC 分布是否异常。
#
# ============================================================
# 2. NCBI/SRA 数据下载与转换
# ============================================================
#
# 安装环境示例：
#
#   conda create -n sratools -c bioconda -c conda-forge sra-tools
#   conda activate sratools
#
# 下载单个 SRR：
#
#   mkdir -p "$PROJECT_DIR/00_raw/sra"
#   prefetch SRR19605702 -O "$PROJECT_DIR/00_raw/sra"
#
# 下载 SRR 列表：
#
#   prefetch --min-size 0 --max-size 100G \
#     --option-file SRR_Acc_List.txt \
#     -O "$PROJECT_DIR/00_raw/sra"
#
# 转换为 fastq。推荐 fasterq-dump；如需压缩，再用 pigz/gzip：
#
#   mkdir -p "$PROJECT_DIR/00_raw/fastq"
#   fasterq-dump "$PROJECT_DIR/00_raw/sra/SRR19605702" \
#     -O "$PROJECT_DIR/00_raw/fastq" \
#     --split-3 \
#     -e 16
#
#   pigz -p 16 "$PROJECT_DIR"/00_raw/fastq/*.fastq
#
# 说明：
#   - split-3 同时适用于单端和双端数据。
#   - 双端一般得到 *_1.fastq 和 *_2.fastq；单端得到一个 fastq。
#
# ============================================================
# 3. 组装：MEGAHIT
# ============================================================
#
# 单样本示例：
#
#   conda activate megahit
#   SAMPLE=NS_0-6
#   megahit \
#     -1 "$PROJECT_DIR/00_raw/fastq/${SAMPLE}_1.fastq.gz" \
#     -2 "$PROJECT_DIR/00_raw/fastq/${SAMPLE}_2.fastq.gz" \
#     --k-min 21 --k-max 141 --k-step 10 \
#     -t 80 \
#     -o "$PROJECT_DIR/03_assembly/${SAMPLE}_assembly"
#
# 注意：
#   - 原稿中 `-0` 应为 `-o`。
#   - 文件名包含空格会导致命令失败，建议统一改为下划线命名。
#   - 每个样本输出的 contig 通常为 `${SAMPLE}_assembly/final.contigs.fa`。
#
# ============================================================
# 4. 物种组成：MetaPhlAn
# ============================================================
#
# 单样本示例：
#
#   conda activate metaphlan
#   mkdir -p "$PROJECT_DIR/04_metaphlan"
#   SAMPLE=NS_0-6
#   metaphlan \
#     "$PROJECT_DIR/00_raw/fastq/${SAMPLE}_1.fastq.gz","$PROJECT_DIR/00_raw/fastq/${SAMPLE}_2.fastq.gz" \
#     --bowtie2out "$PROJECT_DIR/04_metaphlan/${SAMPLE}.bowtie2.bz2" \
#     --nproc 60 \
#     --input_type fastq \
#     -o "$PROJECT_DIR/04_metaphlan/${SAMPLE}_metaphlan.txt"
#
# 合并所有样本：
#
#   cd "$PROJECT_DIR/04_metaphlan"
#   merge_metaphlan_tables.py *_metaphlan.txt > merged_abundance_table.txt
#
# 后续 R/LEfSe 分析可读取：
#   - merged_abundance_table.txt
#
# ============================================================
# 5. 分箱：MetaWRAP binning
# ============================================================
#
# 单样本示例：
#
#   conda activate metawrap
#   SAMPLE=C1_0-6
#   mkdir -p "$PROJECT_DIR/05_binning"
#   nohup metawrap binning \
#     -o "$PROJECT_DIR/05_binning/${SAMPLE}_bins" \
#     -t 50 \
#     -a "$PROJECT_DIR/03_assembly/${SAMPLE}_assembly/final.contigs.fa" \
#     --metabat2 --maxbin2 --concoct \
#     "$PROJECT_DIR/00_raw/fastq/${SAMPLE}_1.fastq.gz" \
#     "$PROJECT_DIR/00_raw/fastq/${SAMPLE}_2.fastq.gz" \
#     > "$PROJECT_DIR/logs/${SAMPLE}.metawrap_binning.log" 2>&1 &
#
# 注意：
#   - 不要重复运行同一个输出目录，例如原稿中 SY365BB-0-4 出现重复。
#   - final.contigs.fa 应使用对应样本的 assembly 结果。
#
# ============================================================
# 6. Bin refinement：MetaWRAP bin_refinement
# ============================================================
#
#   conda activate metawrap
#   SAMPLE=C1_0-6
#   metawrap bin_refinement \
#     -o "$PROJECT_DIR/05_binning/${SAMPLE}_bins/BIN_REFINEMENT" \
#     -t 40 \
#     -A "$PROJECT_DIR/05_binning/${SAMPLE}_bins/metabat2_bins/" \
#     -B "$PROJECT_DIR/05_binning/${SAMPLE}_bins/maxbin2_bins/" \
#     -C "$PROJECT_DIR/05_binning/${SAMPLE}_bins/concoct_bins/" \
#     -c 70 \
#     -x 10
#
# 常用阈值：
#   - 严格：完整性 >= 70%，污染度 <= 10%
#   - 宽松：完整性 >= 50%，污染度 <= 10%
#
# ============================================================
# 7. 基因预测：Prodigal
# ============================================================
#
# 对 assembly contig 预测 CDS 和蛋白：
#
#   conda activate prodigal
#   SAMPLE=SY457BB-8-12
#   mkdir -p "$PROJECT_DIR/06_genes/CDS" "$PROJECT_DIR/06_genes/AA"
#   prodigal \
#     -i "$PROJECT_DIR/03_assembly/${SAMPLE}_assembly/final.contigs.fa" \
#     -d "$PROJECT_DIR/06_genes/CDS/${SAMPLE}_cds.fa" \
#     -a "$PROJECT_DIR/06_genes/AA/${SAMPLE}_AA.faa" \
#     -p meta \
#     -m
#
# 注意：
#   - 原稿中 `SSY457BB8-12_AA.faa` 多了一个 S，建议统一 `${SAMPLE}_AA.faa`。
#
# ============================================================
# 8. 去冗余：CD-HIT
# ============================================================
#
# CDS 核酸去冗余：
#
#   conda activate cd-hit
#   SAMPLE=SY457BB-8-12
#   cd-hit-est \
#     -i "$PROJECT_DIR/06_genes/CDS/${SAMPLE}_cds.fa" \
#     -o "$PROJECT_DIR/06_genes/CDS/${SAMPLE}_cds.cdhit.fa" \
#     -c 0.95 \
#     -aS 0.9 \
#     -G 0 \
#     -M 0 \
#     -g 1 \
#     -T 100 \
#     > "$PROJECT_DIR/logs/${SAMPLE}.cdhit.log" 2>&1
#
# 蛋白去冗余可用 cd-hit：
#
#   cd-hit \
#     -i "$PROJECT_DIR/06_genes/AA/${SAMPLE}_AA.faa" \
#     -o "$PROJECT_DIR/06_genes/AA/${SAMPLE}_AA.cdhit.faa" \
#     -c 0.95 \
#     -aS 0.9 \
#     -M 0 \
#     -T 100
#
# ============================================================
# 9. 基因丰度定量：Salmon index + quant
# ============================================================
#
# 建索引：
#
#   conda activate salmon
#   SAMPLE=SY457BB-8-12
#   mkdir -p "$PROJECT_DIR/07_salmon/${SAMPLE}"
#   salmon index \
#     -t "$PROJECT_DIR/06_genes/CDS/${SAMPLE}_cds.cdhit.fa" \
#     -i "$PROJECT_DIR/07_salmon/${SAMPLE}/salmon_index" \
#     -p 24
#
# 定量：
#
#   salmon quant \
#     -i "$PROJECT_DIR/07_salmon/${SAMPLE}/salmon_index" \
#     -l A \
#     -1 "$PROJECT_DIR/00_raw/fastq/${SAMPLE}_1.fastq.gz" \
#     -2 "$PROJECT_DIR/00_raw/fastq/${SAMPLE}_2.fastq.gz" \
#     -o "$PROJECT_DIR/07_salmon/${SAMPLE}/salmon_quant" \
#     --validateMappings \
#     --meta \
#     -p 70
#
# ============================================================
# 10. 功能注释：KOfamScan 批处理
# ============================================================
#
# 脚本模板：run_kofamscan.sh
#
#   #!/usr/bin/env bash
#   set -euo pipefail
#
#   MAX_JOBS=4
#   job_count=0
#
#   KOFAM_EXEC="/home/ps/Research/databases/other_database/kofamscan/kofam_scan-1.3.0/exec_annotation"
#   PROFILES="/home/ps/Research/databases/other_database/kofamscan/profiles/"
#   KO_LIST="/home/ps/Research/databases/other_database/kofamscan/ko_list"
#
#   for faa in *.cdhit.faa; do
#     prefix="${faa%.cdhit.faa}"
#     tmpdir="${prefix}_ko_tmp"
#     output="${prefix}.kofamscan.txt"
#     mkdir -p "$tmpdir"
#
#     (
#       echo "开始处理: $faa"
#       "$KOFAM_EXEC" --cpu 48 -E 1e-5 \
#         -f detail-tsv \
#         --profile "$PROFILES" \
#         --ko-list "$KO_LIST" \
#         --tmp-dir "$tmpdir" \
#         -o "$output" \
#         "$faa"
#       echo "完成: $faa -> $output"
#     ) &
#
#     ((job_count++))
#     if [[ "$job_count" -ge "$MAX_JOBS" ]]; then
#       wait
#       job_count=0
#     fi
#   done
#
#   wait
#   echo "全部 KOfamScan 任务完成"
#
# ============================================================
# 11. MAG 汇总、去冗余和质量评估
# ============================================================
#
# 删除 MetaWRAP work_files。危险命令，确认路径后再去掉注释：
#
#   # find "$PROJECT_DIR/05_binning" \
#   #   -type d \
#   #   -path "*/BIN_REFINEMENT/work_files" \
#   #   -exec rm -rf {} +
#
# 汇总 refined bins 到统一目录，建议先 dry-run：
#
#   mkdir -p "$PROJECT_DIR/09_MAGs/all_refined_bins"
#   find "$PROJECT_DIR/05_binning" \
#     -path "*/BIN_REFINEMENT/metawrap_50_10_bins/*.fa" \
#     -print
#
# dRep 去冗余：
#
#   conda activate drep
#   dRep dereplicate "$PROJECT_DIR/09_MAGs/dRep_out" \
#     -g "$PROJECT_DIR/09_MAGs/all_refined_bins"/*.fa \
#     -p 150 \
#     --ignoreGenomeQuality
#
# CheckM2 质量评估：
#
#   conda activate checkm2
#   checkm2 predict \
#     -x fa \
#     -i "$PROJECT_DIR/09_MAGs/dRep_out/dereplicated_genomes" \
#     -o "$PROJECT_DIR/09_MAGs/checkm2_result" \
#     -t 64
#
# ============================================================
# 12. 合并 MetaWRAP bin 统计表
# ============================================================
#
# 脚本模板：combine_bin_stats.sh
#
#   #!/usr/bin/env bash
#   set -euo pipefail
#
#   echo -e "sample\tbin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner" \
#     > all_combined_bins.stats
#
#   for file in */BIN_REFINEMENT/metawrap_50_10_bins.stats; do
#     [[ -f "$file" ]] || continue
#     sample_dir=$(dirname "$(dirname "$file")")
#     sample_name=$(basename "$sample_dir" | sed 's/_bins$//')
#     echo "处理样本: $sample_name"
#
#     awk -v sample="$sample_name" 'NR > 1 {
#       gsub(/^bin\./, "", $1);
#       printf "%s\t%s", sample, $1;
#       for (i = 2; i <= NF; i++) printf "\t%s", $i;
#       printf "\n";
#     }' OFS="\t" "$file" >> all_combined_bins.stats
#   done
#
#   echo "完成: all_combined_bins.stats"
#
# ============================================================
# 13. GTDB-Tk 分类和系统发育树
# ============================================================
#
# GTDB-Tk 分类：
#
#   conda activate gtdbtk
#   gtdb-tk classify_wf \
#     --genome_dir "$PROJECT_DIR/09_MAGs/dRep_out/dereplicated_genomes" \
#     --out_dir "$PROJECT_DIR/09_MAGs/classify_wf_out" \
#     --cpus 48 \
#     --extension fa
#
# 主要分类结果：
#
#   ls -lh "$PROJECT_DIR"/09_MAGs/classify_wf_out/gtdbtk.*.summary.tsv
#
# 基于 GTDB-Tk 输出的 marker MSA 建树：
#
#   cd "$PROJECT_DIR/09_MAGs/classify_wf_out/align"
#   gunzip -f -k gtdbtk.*.user_msa.fasta.gz
#
#   conda activate trimal
#   trimal \
#     -in gtdbtk.bac120.user_msa.fasta \
#     -out gtdbtk.bac120.user_msa.trimal.fasta \
#     -gappyout \
#     -fasta
#
#   trimal \
#     -in gtdbtk.ar53.user_msa.fasta \
#     -out gtdbtk.ar53.user_msa.trimal.fasta \
#     -gappyout \
#     -fasta
#
#   conda activate iqtree
#   iqtree \
#     -s gtdbtk.bac120.user_msa.trimal.fasta \
#     -m MFP \
#     -bb 1000 \
#     -alrt 1000 \
#     -nt AUTO \
#     -mem 80G \
#     -pre gtdbtk.bac120.tree \
#     -seed 12345
#
#   iqtree \
#     -s gtdbtk.ar53.user_msa.trimal.fasta \
#     -m MFP \
#     -bb 1000 \
#     -alrt 1000 \
#     -nt AUTO \
#     -mem 80G \
#     -pre gtdbtk.ar53.tree \
#     -seed 12345
#
# 输出：
#   - gtdbtk.bac120.tree.treefile
#   - gtdbtk.ar53.tree.treefile
#
# ============================================================
# 14. MAG 蛋白树
# ============================================================
#
# 如果蛋白序列中含 U、u 或终止符 *，先替换为 X：
#
#   sed -i '/^>/! s/U/X/g; /^>/! s/u/X/g; /^>/! s/\*/X/g' all_proteins.faa
#
# 多序列比对：
#
#   conda activate mafft
#   nohup mafft \
#     --auto \
#     --quiet \
#     --anysymbol \
#     all_proteins.faa \
#     > all_proteins.aln.fasta 2> mafft.log &
#
# 修剪：
#
#   conda activate trimal
#   trimal \
#     -in all_proteins.aln.fasta \
#     -out all_proteins.trimmed.fasta \
#     -automated1
#
# 建树：
#
#   conda activate iqtree
#   iqtree \
#     -s all_proteins.trimmed.fasta \
#     -m TEST \
#     -bb 1000 \
#     -alrt 1000 \
#     -nt AUTO \
#     -pre all_proteins.tree
#
# ============================================================
# 15. HMM 文件检查示例
# ============================================================
#
# 查看 HMM 文件中的模型名称：
#
#   grep "^NAME" arsenic_gene.hmm
#
# 原稿记录过的砷代谢相关模型名包括：
#   aioA, aioB, aioR, aioS, aioX, ARC3_832,
#   arrA, arrB, arrC1, arrC2, arrD, arrR, arrS,
#   arsA_2491_345, arsB_935, arsC_014, arsC_2689_2691,
#   arsD, arsH_2690, arsI, arsJ, arsM, arsO, arsP, arsR, arsT,
#   arxA, arxB2, arxB, arxC, arxD
#
# ============================================================
# 16. 运行前检查清单
# ============================================================
#
#   - 样本名是否统一，是否没有空格。
#   - fastq 后缀是 .fastq、.fq、.fastq.gz 还是 .fq.gz。
#   - 每一步输入是否来自上一阶段对应样本。
#   - 每个 nohup 是否有独立 log 文件。
#   - rm/mv/find -exec 命令是否先 dry-run 查看过路径。
#   - conda 环境是否安装了对应软件。
#   - 服务器线程数和内存是否足够，避免同时启动过多任务。
#
# ============================================================

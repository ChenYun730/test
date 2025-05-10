# test————定量问题解决方法的探索
### 一、检查多映射reads比例（SAM/BAM文件中的XS标签），过高（>10%）可能提示同源/重复序列问题
```
$ samtools view -F 4 /mnt/alamo01/users/chenyun730/program/test_hisat2/alignment/SRR27961778/SRR27961778_sorted.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ /^NH:i:/) { nh=substr($i,6); if (nh == 1) uniq++; else multi++; } } } END {print "Unique mappings:", uniq; print "Multi-mappings:", multi; print "Percentage multi-mapping:", multi/(uniq+multi)*100}'
Unique mappings: 1415394
Multi-mappings:
Percentage multi-mapping: 0
(R441) chenyun730@mgt01:/mnt/alamo01/users/chenyun730/program/test_hisat2
$ samtools view -F 4 /mnt/alamo01/users/chenyun730/program/test_hisat2/alignment/SRR27961779/SRR27961779_sorted.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ /^NH:i:/) { nh=substr($i,6); if (nh == 1) uniq++; else multi++; } } } END {print "Unique mappings:", uniq; print "Multi-mappings:", multi; print "Percentage multi-mapping:", multi/(uniq+multi)*100}'
Unique mappings: 18106826
Multi-mappings:
Percentage multi-mapping: 0
(R441) chenyun730@mgt01:/mnt/alamo01/users/chenyun730/program/test_hisat2
$ samtools view -F 4 /mnt/alamo01/users/chenyun730/program/test_hisat2/alignment/SRR27961780/SRR27961780_sorted.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ /^NH:i:/) { nh=substr($i,6); if (nh == 1) uniq++; else multi++; } } } END {print "Unique mappings:", uniq; print "Multi-mappings:", multi; print "Percentage multi-mapping:", multi/(uniq+multi)*100}'
Unique mappings: 19704911
Multi-mappings:
Percentage multi-mapping: 0
(R441) chenyun730@mgt01:/mnt/alamo01/users/chenyun730/program/test_hisat2
$ samtools view -F 4 /mnt/alamo01/users/chenyun730/program/test_hisat2/alignment/SRR27961787/SRR27961787_sorted.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ /^NH:i:/) { nh=substr($i,6); if (nh == 1) uniq++; else multi++; } } } END {print "Unique mappings:", uniq; print "Multi-mappings:", multi; print "Percentage multi-mapping:", multi/(uniq+multi)*100}'
Unique mappings: 18680847
Multi-mappings:
Percentage multi-mapping: 0
(R441) chenyun730@mgt01:/mnt/alamo01/users/chenyun730/program/test_hisat2
$ samtools view -F 4 /mnt/alamo01/users/chenyun730/program/test_hisat2/alignment/SRR27961788/SRR27961788_sorted.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ /^NH:i:/) { nh=substr($i,6); if (nh == 1) uniq++; else multi++; } } } END {print "Unique mappings:", uniq; print "Multi-mappings:", multi; print "Percentage multi-mapping:", multi/(uniq+multi)*100}'
Unique mappings: 18875024
Multi-mappings:
Percentage multi-mapping: 0
(R441) chenyun730@mgt01:/mnt/alamo01/users/chenyun730/program/test_hisat2
$ samtools view -F 4 /mnt/alamo01/users/chenyun730/program/test_hisat2/alignment/SRR27961789/SRR27961789_sorted.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ /^NH:i:/) { nh=substr($i,6); if (nh == 1) uniq++; else multi++; } } } END {print "Unique mappings:", uniq; print "Multi-mappings:", multi; print "Percentage multi-mapping:", multi/(uniq+multi)*100}'
Unique mappings: 18252563
Multi-mappings:
Percentage multi-mapping: 0
(R441) chenyun730@mgt01:/mnt/alamo01/users/chenyun730/program/test_hisat2
```
结果显示为0，异常，正常人类和哺乳动物都可能5%左右的比例，认为是比对参数出现问题。此外，在比对时未保留多映射信息是一个原因。

### 二、hisat2+featurecount
```
cd /mnt/alamo01/users/chenyun730/program/test/scripts
vim align.sh
#! /bin/bash
#source /mnt/alamo01/users/chenyun730/bin/micromamba
#micromamba activate R441
SAMPLES=(SRR27961778 SRR27961779 SRR27961780 SRR27961787 SRR27961788 SRR27961789)
for SAMPLE in "${SAMPLES[@]}"; do
hisat2 -p 64 \
  -x /mnt/alamo01/users/chenyun730/program/test_hisat2/homo_sapiens/homo_data/GRCh38 \
  -1 /mnt/alamo01/users/chenyun730/program/test_hisat2/clean_data/${SAMPLE}_cleaned_1.fp.gz \
  -2 /mnt/alamo01/users/chenyun730/program/test_hisat2/clean_data/${SAMPLE}_cleaned_2.fp.gz \
  -S /mnt/alamo01/users/chenyun730/program/test/alignment/${SAMPLE}/${SAMPLE}.sam \
  --no-unal \
  --dta \
  --un-conc-gz /mnt/alamo01/users/chenyun730/program/test/alignment/${SAMPLE}/${SAMPLE}_unmapped.fq.gz \
  -k 10 \
  --no-spliced-alignment \
  --time \
  --phred33 \
   --rg-id ${SAMPLE} \
   --rg "SM:${SAMPLE}" \
  2> /mnt/alamo01/users/chenyun730/program/test/alignment/${SAMPLE}/${SAMPLE}.align.stats
done
```
为保持准确性，定量时先使用gene-id（严格一对一），之后再使用GENCODE/ENSEMBL 官方注释 将 gene_id 映射到 symbol，并处理多对一关系，以此提高准确性
```
vim quantify.sh
#! /bin/bash
#micromamba activate R441
if [[ -z "$(which featureCounts)" ]]; then
    source "/mnt/alamo01/users/chenyun730/bin/micromamba"
    micromamba activate R441
fi
BAM_DIR=/mnt/alamo01/users/chenyun730/program/test/alignment
GTF_FILE=/mnt/alamo01/users/chenyun730/program/test_hisat2/homo_sapiens/homo_data/Homo_sapiens.GRCh38.109.gtf
OUTPUT_DIR=/mnt/alamo01/users/chenyun730/program/test/quantify
SAMPLES=(SRR27961778 SRR27961779 SRR27961780 SRR27961787 SRR27961788 SRR27961789)
bam_files=()
for sample in "${SAMPLES[@]}"; do
    bam_files+=("${BAM_DIR}/${SAMPLES}/${SAMPLES}.sorted.bam")
done
featureCounts -T 8 \
    -a "${GTF_FILE}" \
    -o "${OUTPUT_DIR}/gene_counts_geneid.txt" \
    -g gene_id \
    -p \
    --countReadPairs \
    -s 1 \
    -M \
    -O \
    --fraction \
     "${bam_files[@]}" > "${OUTPUT_DIR}/featurecounts_summary.txt" 2>&1
awk -F'\t' '$3=="gene" {
    match($0, /gene_id "[^"]+"/, gid); gsub(/gene_id "|"/, "", gid[0]);
    match($0, /gene_name "[^"]+"/, gname); gsub(/gene_name "|"/, "", gname[0]);
    if (gid[0] != "" && gname[0] != "") print gid[0]"\t"gname[0]
}' "$GTF_FILE" | sort -u > "${OUTPUT_DIR}/gene_id_to_symbol.tsv"
awk 'NR==FNR {map[$1]=$2; next}
     FNR==1 && NR>FNR {header=$0; next}
     FNR>1 && ($1 in map) {
         print $1"\t"map[$1]"\t"$0
     }' \
    "${OUTPUT_DIR}/gene_id_to_symbol.tsv" \
    <(tail -n +3 "${OUTPUT_DIR}/gene_counts_geneid.txt") > "${OUTPUT_DIR}/gene_counts_with_symbol.tmp"
{
    read -r line
    echo -e "gene_name\tgene_id\tmock_1\tmock_2\tmock_3\tsars2_1\tsars2_2\tsars2_3"
    while IFS=$'\t' read -r gene_id gene_name _ _ _ _ _ _ c1 c2 c3 c4 c5 c6; do
        if [[ "$gene_name" != "" ]]; then
            echo -e "$gene_name\t$gene_id\t$c1\t$c2\t$c3\t$c4\t$c5\t$c6"
        fi
    done
} < "${OUTPUT_DIR}/gene_counts_with_symbol.tmp" > "${OUTPUT_DIR}/gene_count_matrix_final.csv"
cut -f1 "${OUTPUT_DIR}/gene_count_matrix_final.csv" | tail -n +2 | sort | uniq -d > "${OUTPUT_DIR}/symbol_conflict.log"
if [[ -s "${OUTPUT_DIR}/symbol_conflict.log" ]]; then
    echo "Warning: Duplicated gene symbols found. See ${OUTPUT_DIR}/symbol_conflict.log"
fi
echo "Success! Count matrix saved to ${OUTPUT_DIR}/gene_count_matrix_final.csv"

qsub -cwd -V -l cpu=64:mem=64G -q fast -N quantify /mnt/alamo01/users/chenyun730/program/test/scripts/quantify.sh

```
**quantify.sh运行出来的矩阵有问题，另写了一个featurecount.sh重新生成原始counts矩阵**
```
qsub -cwd -V -l cpu=64:mem=64G -q fast -N featurecount /mnt/alamo01/users/chenyun730/program/test/scripts/featurecount.sh

#!/bin/bash

if [[ -z "$(which featureCounts)" ]]; then
    source "/mnt/alamo01/users/chenyun730/bin/micromamba"
    micromamba activate R441
fi
# 定义路径和文件
BAM_DIR=/mnt/alamo01/users/chenyun730/program/test/alignment
GTF_FILE=/mnt/alamo01/users/chenyun730/program/test_hisat2/homo_sapiens/homo_data/Homo_sapiens.GRCh38.109.gtf
OUTPUT_DIR=/mnt/alamo01/users/chenyun730/program/test/featurecount

SAMPLES=(SRR27961778 SRR27961779 SRR27961780 SRR27961787 SRR27961788 SRR27961789)
NEW_NAMES=(mock_1 mock_2 mock_3 sars2_1 sars2_2 sars2_3)

# 准备BAM文件路径
bam_files=()
for sample in "${SAMPLES[@]}"; do
    bam_path="${BAM_DIR}/${sample}/${sample}.sorted.bam"
    if [[ ! -f "$bam_path" ]]; then
        echo "错误：BAM文件不存在 - $bam_path" >&2
        exit 1
    fi
    bam_files+=("$bam_path")
done

echo "正在运行featureCounts..."
# 运行featureCounts并检查是否成功
if ! featureCounts -T 8 \
    -a "$GTF_FILE" \
    -o "${OUTPUT_DIR}/gene_counts_raw.txt" \
    -g gene_id \
    -p \
    --countReadPairs \
    -s 1 \
    -M \
    -O \
    --fraction \
    "${bam_files[@]}" > "${OUTPUT_DIR}/featurecounts_summary.txt" 2>&1
then
        echo "featureCounts运行失败，请检查日志文件：${OUTPUT_DIR}/featurecounts_summary.txt" >&2
    exit 1
fi
echo "处理计数矩阵..."
# 提取并重命名表头
{
    # 打印新表头
    echo -en "gene_id"
    for name in "${NEW_NAMES[@]}"; do
        echo -en "\t$name"
    done
    echo
    
    # 提取计数数据（跳过前两行注释）
    tail -n +3 "${OUTPUT_DIR}/gene_counts_raw.txt" | cut -f1,7-
} > "${OUTPUT_DIR}/gene_count_matrix.txt"

# 生成CSV格式
awk 'BEGIN {OFS=","} {print $0}' "${OUTPUT_DIR}/gene_count_matrix.txt" > "${OUTPUT_DIR}/gene_count_matrix.csv"

echo "分析完成！"
echo "计数矩阵已保存至："
echo "  - ${OUTPUT_DIR}/gene_count_matrix.txt"
echo "  - ${OUTPUT_DIR}/gene_count_matrix.csv"

```

### 核心思路与多映射问题处理方法
1. 整体分析流程
数据准备：加载原始count数据和GTF注释

ID转换：建立gene_id到symbol的精确映射

质量控制：处理缺失注释和低表达基因

差异分析：使用DESeq2进行标准化和统计检验

可视化：通过PCA和热图展示数据结构和差异基因

2. 多映射问题处理策略
优先保留原则：当多个gene_id映射到同一symbol时：

计算每个基因在所有样本中的总表达量

对相同symbol的基因按总表达量降序排列

只保留表达量最高的记录（slice(1)）

科学依据：

高表达基因通常具有更高的生物学重要性

避免人为合并不同基因的表达量（可能来自假基因或旁系同源基因）

保留最可能的蛋白编码基因转录本

3. 关键保障措施
版本号处理：统一去除Ensembl ID的版本后缀（如.10）

缺失值处理：无symbol的基因自动使用gene_id代替

重复检查：每一步都验证是否成功解决重复问题

可视化验证：通过PCA检查批次效应，热图验证差异基因

4. 结果输出
标准化表格：第一列为symbol，后续为样本count的Excel

分析结果：DESeq2统计结果表格

高质量图表：PCA图和热图（PNG格式）

```
library(rtracklayer)    # 用于解析GTF文件
library(DESeq2)        # 差异表达分析
library(ggplot2)       # 基础绘图
library(pheatmap)      # 热图绘制
library(RColorBrewer)  # 颜色配置
library(ggrepel)       # 防止PCA标签重叠
 library(rtracklayer)
library(dplyr)

# 步骤1：加载并预处理count数据 ----------------------------------------
 count_data <- read.table("/mnt/alamo01/users/chenyun730/program/test/featurecount/gene_count_matrix.txt", header = TRUE, sep = "\t")
 head(count_data)
          gene_id mock_1 mock_2 mock_3 sars2_1 sars2_2 sars2_3
1 ENSG00000160072 184.33  190.1  179.5  180.31  184.45  191.58
2 ENSG00000279928   0.00    0.0    0.0    0.20    0.48    0.25

# 步骤2：解析GTF文件构建基因注释映射表 ----------------------------------------
parse_gtf <- function(gtf_path) {
   gtf <- import("/mnt/alamo01/users/chenyun730/program/test_hisat2/homo_sapiens/homo_data/Homo_sapiens.GRCh38.109.gtf") 
   gene_anno <- gtf[gtf$type == "gene"]
   data.frame(
    gene_id = sub("\\..*", "", gene_anno$gene_id), 
    gene_symbol = ifelse(is.na(gene_anno$gene_name), 
                        gene_anno$gene_id,          # 无symbol时使用gene_id
                        gene_anno$gene_name),
    gene_type = gene_anno$gene_biotype,
    stringsAsFactors = FALSE ) %>% distinct()  # 去重
}
gene_map <- parse_gtf("Homo_sapiens.GRCh38.109.gtf")
 head(gene_map)
          gene_id     gene_symbol              gene_type
1 ENSG00000160072          ATAD3B         protein_coding
2 ENSG00000279928        DDX11L17 unprocessed_pseudogene
3 ENSG00000228037 ENSG00000228037                 lncRNA
4 ENSG00000142611          PRDM16         protein_coding
write.csv(gene_map, "gene_annotation.csv", row.names = FALSE)   #保存这个基因注释表到results中

# 步骤3：合并count数据与基因注释 ---------------------------------------- 
merged_data <- merge(count_data, gene_map, by = "gene_id", all.x = TRUE)     # 左连接保留所有 count_data 中的基因

num_missing_symbols <- sum(is.na(merged_data$gene_symbol))
cat("缺失 gene_symbol 的基因数目：", num_missing_symbols, "\n")    # 统计有多少缺失的 gene_symbol

merged_data$gene_symbol <- ifelse(is.na(merged_data$gene_symbol),
                                   merged_data$gene_id,
                                   merged_data$gene_symbol)      # 填补缺失值（用 gene_id 替代）
# 返回结果： 缺失 gene_symbol 的基因数目： 0

# 步骤4：解决多映射问题（关键步骤） ----------------------------------------
# 策略：对于相同symbol的基因，保留表达量总和最高的记录
processed_data <- merged_data %>%
  mutate(total_counts = rowSums(select(., mock_1:sars2_3))) %>%  # 或 where(is.numeric)
  arrange(gene_symbol, desc(total_counts)) %>%
  group_by(gene_symbol) %>%
  slice(1) %>%  # 每组保留表达量最高的一行
  ungroup() %>%
  select(-total_counts)  # 移除临时列
#以上是原本的代码，报错后改为下列代码可运行

 sum(is.na(merged_data$gene_symbol))  # 排除gene_symbol 列有缺失值（NA）
[1] 0
 processed_data <- merged_data %>%
  mutate(total_counts = rowSums(select(., mock_1:sars2_3))) %>%
  group_by(gene_symbol) %>%     #可能是某些基因的 total_counts 相同，导致 slice(1) 不稳定，改用 dplyr::distinct() 直接去重
  arrange(desc(total_counts)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  ungroup() %>%
  select(-total_counts)
 n_distinct(processed_data$gene_symbol) == nrow(processed_data)  # 检查去重后的基因数量，每个基因唯一
[1] TRUE
 any(duplicated(processed_data$gene_symbol))  # 检查是否有重复基因
[1] FALSE
if (any(duplicated(processed_data$gene_symbol))) {
  warning("仍有重复symbol存在，请手动检查！")
}   # 检查重复symbol是否已解决

# 步骤5：输出symbol转换后的Excel文件 ----------------------------------------
# 整理输出格式：symbol为第一列，接着是6个样本count
final_matrix <- processed_data %>%
  select(gene_symbol, mock_1, mock_2, mock_3, sars2_1, sars2_2, sars2_3)
write.csv(final_matrix, "/mnt/alamo01/users/chenyun730/program/test/results/GeneExpression_SymbolCounts.csv", row.names = FALSE, na = "")
write.table(final_matrix, "/mnt/alamo01/users/chenyun730/program/test/results/GeneExpression_SymbolCounts.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")

# 步骤6：进行DEseq2和PCA ----------------------------------------
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
count_data <- read.table("/mnt/alamo01/users/chenyun730/program/test/results/GeneExpression_SymbolCounts.txt", header=TRUE, row.names=1, sep="\t")  #该文件由于多映射等原因不为整数
head(count_data)
counts <- round(read.table("/mnt/alamo01/users/chenyun730/program/test/results/GeneExpression_SymbolCounts.txt", header=TRUE, row.names=1, sep="\t"))  #这里取整，接下来的步骤都按这个矩阵的数值进行
 write.table(counts, file = "gene_counts.txt", sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = NA)

sample_info <- data.frame(
  sample = colnames(counts),
  condition = factor(c(rep("mock", 3), rep("sars2", 3))),
  row.names = colnames(count_data)
)
 dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ condition
)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "sars2", "mock"))
write.csv(as.data.frame(res), "deseq2_results.csv")    #保存了DEseq2的结果
vsd <- vst(dds, blind=FALSE)
norm_counts <- assay(vsd)

# 步骤7：可视化 ----------------------------------------
#PCA分析
pca_data <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color=condition, fill=condition)) +
  geom_point(size=3.5, shape=21, stroke=0.5) +  # 实心圆点，带轮廓
  scale_color_manual(values=c("mock"="blue", "sars2"="red")) +
  scale_fill_manual(values=c("mock"="lightblue", "sars2"="pink")) +
  geom_text_repel(aes(label=name), size=3, box.padding=0.5) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA Plot of Mock vs SARS2") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust=0.5, size=14, face="bold"),
    legend.position = "right"
  )
ggsave("pca_plot.pdf", pca_plot, width=8, height=6)
ggsave("pca_plot.png", pca_plot, width=8, height=6, dpi=300)

#热图
sig_genes <- rownames(res)[which(res$padj < 0.05)]
top_genes <- head(order(res$padj), 30)

heatmap_data <- norm_counts[top_genes, ]
rownames(heatmap_data) <- substr(rownames(heatmap_data), 1, 15) 
annotation_col <- data.frame(
  Condition = factor(c(rep("Mock", 3), rep("SARS2", 3)))
)
rownames(annotation_col) <- colnames(heatmap_data)
ann_colors <- list(
  Condition = c(Mock="blue", SARS2="red")
)
pdf("heatmap.pdf", width=10, height=8)
pheatmap(heatmap_data,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         scale="row",
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         annotation_col=annotation_col,
         annotation_colors=ann_colors,
         show_rownames=TRUE,
         fontsize_row=8,
         main="Top 30 Significant Genes (padj < 0.05)")
dev.off()

png("heatmap.png", width=1000, height=800, res=120)
pheatmap(heatmap_data,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         scale="row",
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         annotation_col=annotation_col,
         annotation_colors=ann_colors,
         show_rownames=TRUE,
         fontsize_row=8,
         main="Top 30 Significant Genes (padj < 0.05)")
dev.off()

```
**进行合并矩阵的分析**
```
geo_counts <- read.csv("/mnt/alamo01/users/chenyun730/program/test_hisat2/geo_matrix/symbol_count_matrix_geo.csv", row.names = NULL, check.names = FALSE)
 counts <- read.table("/mnt/alamo01/users/chenyun730/program/test/results/gene_counts_round.txt", header=TRUE, stringsAsFactors=FALSE, row.names= NULL)
cat("您的矩阵基因数:", nrow(counts), "\n")
cat("GEO矩阵基因数:", nrow(geo_counts), "\n")
cat("共同基因数:", length(intersect(counts$row.names, geo_counts$symbol)), "\n")
您的矩阵基因数: 61255
GEO矩阵基因数: 38879cat
共同基因数: 25583
colnames(counts)[1] <- "symbol"
 merged_counts <- merge(geo_counts, counts, by = "symbol", all = TRUE)
 merged_counts[is.na(merged_counts)] <- 0
write.table(merged_counts, "merged_counts.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.csv(merged_counts, "merged_counts.csv", row.names=FALSE)   # 保存合并矩阵

count_data <- read.table("/mnt/alamo01/users/chenyun730/program/test/results/merged_counts.txt", row.names=1, header=TRUE, sep= "\t")
head(count_data)
sample_info <- data.frame(
  sample = colnames(count_data),
  group = factor(rep(c("geo_mock", "geo_sars2", "local_mock", "local_sars2"), 
                 each = 3)),
  batch = factor(rep(c("geo", "local"), each = 6)),
  row.names = colnames(count_data)
)
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_info,
  design = ~ group
)
 keep <- rowSums(counts(dds) >= 5) >= 3
dds <- dds[keep,]
 dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

pheatmap(
  mat = heatmap_data,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Top 50 Significant Genes (Four Groups Comparison)",
  filename = "heatmap_four_groups.png",
  width = 10,
  height = 8,
  res = 300,
  units = "in"  # 明确指定单位
)
```




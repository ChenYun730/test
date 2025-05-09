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

```

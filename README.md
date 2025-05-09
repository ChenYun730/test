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

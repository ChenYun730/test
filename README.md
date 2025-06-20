***Work Diary***

2025/6/16
```
#From wsy（从下载数据到文章所有样本的表达矩阵）：
/mnt/alamo01/users/chenyun730/download_data/GSE255647/results/run_rnaseq_pipeline_0616.sh
/mnt/alamo01/users/chenyun730/scripts/combine_symbol_counts.R

#输入文件（GSE255647的所有下载的原始数据）
  /mnt/alamo01/users/chenyun730/download_data/GSE255647

#运行结果：
    /mnt/alamo01/users/chenyun730/download_data/GSE255647/results   #生成RNAseq的所有结果及整合矩阵
	 /mnt/alamo01/users/chenyun730/download_data/GSE255647/results/alignment/ gene_count_matrix.csv
	 /mnt/alamo01/users/chenyun730/download_data/GSE255647/results/alignment /gene_count_matrix_symbol_merged.csv
```
2025/6/17-2025/6/18
```
#编写脚本分装矩阵所需的metadata：
/mnt/alamo01/users/chenyun730/scripts/download/ fetch_geo_metadata.R   #脚本
/mnt/alamo01/users/chenyun730/download_data/GSE255647/GSE255647_metadata.json   #结果不理想，只有SRR没有样本信息

#下载GSE147507的数据：
/mnt/alamo01/users/chenyun730/scripts/download/GSE147507_srrdownload_all.sh
```
2025/6/19
```
# 提取GSM元数据：
    /mnt/alamo01/users/chenyun730/scripts/download/ GSE147507_metadata.R  #脚本
		 /mnt/alamo01/users/chenyun730/download_data/GSE147507/ GSE147507.json    #结果符合预期

# 下载时发现有个别样本因为超时等原因下载失败，补充下载了四个样本：
  /mnt/alamo01/users/chenyun730/scripts/download/ add_GSE147507.sh

# 对SRR样本进行合并并改用GSM编号命名：
GSE147507_merged_to_gsm.sh
/mnt/alamo01/users/chenyun730/download_data/GSE147507   #结果

# GSE255647获取了包含GSM和SRR的metadata：
# GSE255647.json和GSE255647_fetch.json合并为GSE255647_metadata.json（包含需要的信息）
		脚本在GSE255647_metadata.R、 GSE255647_ fetch_geo_metadata.R、 GSE255647_ merge_gsm_srr.R
```



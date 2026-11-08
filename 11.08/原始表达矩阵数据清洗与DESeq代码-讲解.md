下面是 run_deseq2.R 的整体结构，我把它拆成 6 个大的条块 来说明每一部分大致做什么。

1. 环境初始化与输出目录设置
对应行： 1–9

声明这是一个 Rscript 脚本。
关闭 R 的默认 stringsAsFactors。
Windows 下把 locale 设成中文 UTF-8，避免中文路径/标签乱码。
静默加载 DESeq2、data.table、ggplot2。
固定随机种子，方便结果可重复。
创建输出目录 deseq2_analysis/，并把当前 R 的 sessionInfo（包版本、R 版本等）写进去，作为可追溯的记录。
2. 读取并审计样本元数据
对应行： 10–22

读取 metadata anno_final_aligned 文件。
检查样本 ID 是否有重复。
整理出后续分析需要的字段：
gender（性别小写）
condition（分组：HD / SLE）
batch（从样本 ID 前 6 位提取的日期批次）
标记三类身份异常：
同一个 sample 对应多个 patient
同一个 patient 跨越了 HD/SLE 两组
patient 字段为空或 #N/A、-
把有问题的样本输出到 identity_exclusions.csv。
3. 读入原始 counts 并做质控
对应行： 23–39

自动查找 gencode*.csv 原始 count 矩阵。
每一行是一个基因，每一列是一个样本的原始整数 counts。
计算每个样本的：
total_counts（总 reads）
detected_genes（检测到表达的基因数）
设定一个宽松的质控门槛：
总 counts ≥ 10,000
检测到基因数 ≥ 1,000
检查 raw count 里的样本和 metadata 里的样本是否对齐，输出没有 metadata 的样本。
把 metadata 和 QC 结果合并，保存为 metadata_audited.csv。
把 count 矩阵转成整数矩阵，保存为 aligned_raw_counts.rds。
画一张 total_counts vs detected_genes 的 QC 散点图保存为 PNG。
4. 核心分析函数 run_one
对应行： 40–96

这是整个脚本的核心，对**每一个细胞类型（cell）**跑一次 DESeq2。主要流程：

建子目录：给当前 analysis + cell 类型创建目录。
样本筛选：
只保留当前 cell、身份无异常、QC 通过的样本。
同一个 patient 有多个样本时，按日期取最早一次；其余输出到 repeated_patient_exclusions.csv。
设置因子：condition（HD 在前、SLE 在后）、gender、batch。
样本量检查：任一组少于 3 个样本就直接跳过。
构建设计矩阵：~ gender + batch + condition（如果 gender 或 batch 只有 1 个水平就自动去掉）。
共线性/残差自由度检查：如果设计矩阵不满秩或剩余自由度不足，也跳过。
基因过滤：保留在较小分组中至少有 min_n 个样本 counts ≥ 10 的基因（默认 DESeq2 过滤）。
跑 DESeq2：
DESeqDataSetFromMatrix 建对象
DESeq(sfType='poscounts') 估计 size factor、离散度、拟合 GLM
支持 DESEQ2_RESUME=1 环境变量来复用已保存的 dds.rds
提取结果：
contrast = condition SLE vs HD
计算 padj_global（先占位，后面再算）
定义 DEG：padj < 0.05 且 |log2FoldChange| ≥ 1
输出结果：
全部基因结果 DESeq2_all_genes.csv
DEG 列表 DEG_FDR05_absLFC1.csv
样本信息（含 size factor）samples.csv
dds.rds、normalized_counts.rds
画图诊断：
plotDispEsts：离散度估计
plotMA：MA 图
PCA：取 top 500 高变基因做 PCA 看样本聚类
Volcano 图
汇总返回值：把当前 cell 的状态、样本数、保留基因数、DEG 数等返回，并把 p 值收集到全局列表 all_p 中，供后面做跨 cell 的 FDR 校正。
5. 外层循环：两次分析策略
对应行： 97–107

对两种比较策略分别跑一遍：

5.1 primary_baseline_SLE_vs_HD
HD 全部保留。
SLE 只保留 Timepoint 为 - 或 baseline 的基线样本。
5.2 sensitivity_naive_SLE_vs_HD
只比较 HD 和初治 SLE（SLE-na_ve）。
对每一种策略，遍历该策略下所有可用的 cell 类型，调用 run_one()，并把结果追加到 summary_rows。每次循环都实时写入 analysis_summary.csv，方便中途查看进度。

如果任何 cell 出现 ERROR，脚本在最后会停下来报错。

6. 全局 FDR 校正与收尾
对应行： 108–117

把所有 cell 的 p 值合并成一张大表。
在每个 analysis 内部，把所有基因 × 所有 cell 的 p 值再用 BH 法校正一次，得到 padj_global（全局 FDR）。
输出：
all_tests_global_FDR.csv.gz：全部检验结果（含全局 padj）
DEG_global_FDR05_absLFC1.csv：基于全局 FDR 和 |LFC|≥1 筛选的 DEG
最后再写一遍 sessionInfo.txt。
打印 COMPLETE。

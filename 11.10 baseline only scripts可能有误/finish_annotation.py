"""Validate organized annotated artifacts and write the directory guide/report."""
from pathlib import Path
import json,hashlib,re
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent.parent
def slug(s):return re.sub(r'[^A-Za-z0-9]+','_',s.replace('-','_minus').replace('+','_plus'))
def mdtable(df):
 return '\n'.join(['| '+' | '.join(df.columns)+' |','| '+' | '.join(['---']*len(df.columns))+' |']+['| '+' | '.join(map(str,row))+' |' for row in df.itertuples(index=False,name=None)])

def main():
 ann=pd.read_csv(ROOT/'01_annotation/gene_annotation.csv')
 prov=json.loads((ROOT/'01_annotation/annotation_provenance.json').read_text(encoding='utf-8'))
 full=pd.read_csv(ROOT/'03_results/all_tests_global_FDR.csv.gz')
 s=pd.read_csv(ROOT/'03_results/analysis_summary_with_global_FDR.csv')
 fi=pd.read_csv(ROOT/'06_reports/figure_inventory.csv')
 labels=pd.read_csv(ROOT/'03_results/volcano_top10_all_comparisons.csv')
 checks=[]
 def check(name,value):checks.append({'check':name,'passed':bool(value)})
 check('all_60649_genes_named',len(ann)==60649 and ann.gene_symbol.notna().all() and not ann.Geneid.duplicated().any())
 check('all_30_primary_comparisons_present',len(fi)==30 and full.groupby(['analysis','cell']).ngroups==30 and set(full.analysis)=={'primary_baseline_SLE_vs_HD'})
 check('unique_gene_test_keys',not full.duplicated(['analysis','cell','Geneid']).any())
 for row in s.itertuples(index=False):
  an,cell,tag=row.analysis,row.cell,slug(row.cell)
  r=pd.read_csv(ROOT/'03_results'/an/tag/'DESeq2_all_genes.csv')
  a=full[(full.analysis==an)&(full.cell==cell)].set_index('Geneid').reindex(r.Geneid)
  check(f'{an}/{tag}/statistics_consistent',len(r)==len(a) and all(np.allclose(r[c],a[d],rtol=1e-13,atol=1e-300,equal_nan=True) for c,d in [('pvalue','pvalue'),('log2FoldChange','log2FoldChange'),('padj_global','padj_global')]))
  n=((r.padj_global<.05)&(r.log2FoldChange.abs()>=1)).sum()
  check(f'{an}/{tag}/DEG_count_unchanged',n==row.global_DEG)
  l=labels[(labels.analysis==an)&(labels.cell==cell)]
  check(f'{an}/{tag}/label_criteria',len(l)<=10 and l.gene_symbol.notna().all() and (l.padj_global<.05).all() and (l.log2FoldChange.abs()>=1).all())
  check(f'{an}/{tag}/figures_exist',all((ROOT/'04_figures'/kind/an/(tag+ext)).exists() and (ROOT/'04_figures'/kind/an/(tag+ext)).stat().st_size>1000 for kind,ext in [('volcano','.pdf'),('volcano','.png'),('diagnostics','.pdf')]))
 validation=pd.DataFrame(checks);validation.to_csv(ROOT/'06_reports/annotation_validation.csv',index=False)
 if not validation.passed.all():raise RuntimeError('Annotation validation failed')
 examples=labels[(labels.analysis=='primary_baseline_SLE_vs_HD')&(labels.cell=='Plasmablasts')][['label_rank','Geneid','gene_symbol','gene_type','log2FoldChange','padj_global']].copy()
 examples.log2FoldChange=examples.log2FoldChange.map(lambda x:f'{x:.3f}')
 examples.padj_global=examples.padj_global.map(lambda x:f'{x:.3g}')
 report=f'''# 基因名称标注与目录整理报告

基于 GENCODE v38（UCSC hg38 官方镜像），为全部原始基因提供名称和基因类型：60,605 个带版本号 ID 精确匹配；44 个 `_PAR_Y` 条目按同版本基础基因关联名称，并保留原 ID 和独立匹配标记。共 60,649 个基因，无未匹配条目。数据库及字段定义原文件、SHA256 和来源说明均保存在 `01_annotation`。

数据库来源：[{prov['source_url']}]({prov['source_url']})。名称为该版本的 GENCODE gene name，不宣称等同于当前最新版 HGNC 名称。原始 Geneid 未被基因名称替换，也没有合并同名基因。

## 表格和模型中的新增信息

所有含 Geneid 的结果表、全局检验表、跨细胞重复基因表、交集表及基因过滤表均增加 `gene_symbol`、`gene_type`、`annotation_match`、`annotation_source`。184 个原结果/筛选表完成标注；Top10 表另行生成。每个 DESeq2 模型的 rowData 及标准化计数的 gene_annotation 属性也补充注释；aligned_raw_counts.rds 新增 gene_annotation 成员。

## 图中的标注规则

- 30 张主分析火山图各有 PDF 和 PNG。横轴为 baseline SLE/HD 的未收缩 log2FC，纵轴为 -log10(全局 BH FDR)。
- 红色为 SLE 较高，蓝色为 SLE 较低，筛选要求全局 FDR<0.05 且 |log2FC|≥1。
- 每图最多10个基因：按全局FDR从小到大，上调、下调各优先5个；一侧不足时由另一侧补足。标签有连线且分开排布。
- 这是一种兼顾上下调的 Top10，未单纯按最大倍数排序。无足够合格基因时不强行补入不显著基因。
- 每个比较的 Top10 表位于对应结果目录，全部标签汇总见 `03_results/volcano_top10_all_comparisons.csv`。
- 30 份诊断 PDF 已重建，包含离散度、MA、PCA、带相同基因名称的火山图；MA 页附该比较的 Top10 名称和效应值。

## 浆母细胞主分析的标注示例

{mdtable(examples)}

## 校验及解释边界

{len(validation)} 项检查全部通过：30项 baseline SLE vs HD 主分析比较完整、基因检验键唯一、表间统计值一致、原汇总DEG数量未变、全部标签满足阈值、图文件完整。文件增加解释信息并调整存放位置。

注释确认原矩阵除 19,955 个 protein_coding 基因外，还含 lncRNA、假基因及其他非编码RNA。因此本轮结果是全部输入基因的差异分析；如果限定蛋白编码基因，应另立分析范围，不能把当前全部行都称为蛋白编码 mRNA。当前保留全部结果，并用 gene_type 明确区分。

身份冲突排除、月份代理批次和治疗混杂等原分析限制仍然适用。图上的基因是筛选后的代表性标签，不意味着其已得到机制验证。
'''
 (ROOT/'06_reports/ANNOTATION_REPORT.md').write_text(report,encoding='utf-8')
 guide='''# SLE 差异表达分析文件入口

| 目录 | 内容 |
| --- | --- |
| 00_scripts | 完整分析、注释、绘图、整理及验证代码 |
| 01_annotation | GENCODE v38 数据库、全基因名称/类型对照表、注释覆盖和来源 |
| 02_quality_control | 文库质控、身份排除、各比较样本/基因筛选及 PCA 分数 |
| 03_results | 差异表达结果、全局 FDR、交集、各比较 Top10 和汇总表 |
| 04_figures/volcano | 带 Top10 基因名称的火山图，PDF 和 PNG |
| 04_figures/diagnostics | 每个比较的离散度、MA、PCA、带名称火山图组成的诊断 PDF |
| 04_figures/overview | 文库 QC 和差异基因数汇总图 |
| 05_models | DESeq2 模型、标准化计数及对齐原始计数，均附基因注释 |
| 06_reports | 原分析报告、注释报告、方法说明、校验和运行日志 |

先看 [本轮注释报告](06_reports/ANNOTATION_REPORT.md)，再看 [原分析结果报告](06_reports/REPORT.md)。

本轮只运行并更新 `primary_baseline_SLE_vs_HD`：全部 HD 对仅 `Timepoint=baseline` 的 SLE；其他 SLE 时点不纳入。`LDG_CD14_plus` 和 `LDG_CD14_minus` 分开保存。

在原始 data 目录运行：

```powershell
# 只更新注释和图，不重新拟合模型
.\\deseq2_baseline_only_analysis\\00_scripts\\refresh_annotations.ps1

# 从原始计数完整重跑：分析、汇总、注释、绘图和分类整理
.\\deseq2_baseline_only_analysis\\00_scripts\\run_pipeline.ps1
```

完整运行脚本先在 `.work` 写入中间结果，再归入上述目录；正式结果无需在 `.work` 查找。同名生成文件会被更新。数据库已缓存在本地，重绘不需要再次下载。原始三个输入文件位于上级 data 目录。
'''
 (ROOT/'README.md').write_text(guide,encoding='utf-8')
 # Keep the statistical report, correcting only its now-outdated annotation statement.
 original=ROOT/'06_reports/REPORT.md'
 if original.exists():
  text=original.read_text(encoding='utf-8')
  text=text.replace('- 没有匹配GENCODE版本的GTF，未做基因符号/基因类型映射，也未作GO/KEGG富集；不能把全部行直接称为蛋白编码mRNA。','- 已补充GENCODE v38基因名称和类型；详见ANNOTATION_REPORT.md。未作GO/KEGG富集，输入包含多种RNA类型。')
  original.write_text(text,encoding='utf-8')
 methods=ROOT/'06_reports/README.md'
 if methods.exists():
  text=methods.read_text(encoding='utf-8')
  if not text.startswith('> 目录已整理'):
   text='> 目录已整理：当前运行入口及文件位置见上级 README.md；以下保留初次分析的方法记录，旧路径不再作为运行入口。基因注释状态已更新，见 ANNOTATION_REPORT.md。\n\n'+text
  methods.write_text(text,encoding='utf-8')
 manifest=[]
 for p in sorted((ROOT/'00_scripts').iterdir()):
  if p.is_file():manifest.append({'file':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
 pd.DataFrame(manifest).to_csv(ROOT/'06_reports/updated_code_manifest_sha256.csv',index=False)
 print(f'PASS: {len(validation)} validation checks; {len(labels)} figure labels; {len(fi)} comparisons')

if __name__=='__main__':main()

"""Check ORA calculations and outputs, then summarize results in Chinese."""
from pathlib import Path
import math,json,hashlib
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
def table(df):
 return '\n'.join(['| '+' | '.join(map(str,df.columns))+' |','| '+' | '.join(['---']*len(df.columns))+' |']+['| '+' | '.join(map(str,r))+' |' for r in df.itertuples(index=False,name=None)])
def main():
 s=pd.read_csv(ROOT/'03_results/enrichment_summary.csv')
 sig=pd.read_csv(ROOT/'03_results/significant_terms_global_FDR05.csv')
 mapping=pd.read_csv(ROOT/'02_gene_sets/ID_mapping_summary.csv')
 inventory=pd.read_csv(ROOT/'05_reports/figure_inventory.csv')
 checks=[]
 def check(n,v):checks.append({'check':n,'passed':bool(v)})
 check('720_families_accounted_for',len(s)==720 and not s.duplicated(['analysis','cell','database','direction']).any())
 check('60_gene_mapping_summaries',len(mapping)==60)
 check('all_foregrounds_within_background',bool((s.annotated_foreground_n<=s.universe_N).all()))
 check('all_significant_terms_match_criteria',bool((sig.padj_analysis_database<.05).all() and (sig.k>=3).all()))
 total=0;valid=True;tail_examples=[]
 for d in pd.read_csv(ROOT/'03_results/all_enrichment_tests.csv.gz',chunksize=100000):
  total+=len(d)
  valid &= bool((d.k>=0).all() and (d.k<=d.K).all() and (d.k<=d.foreground_n).all() and d.pvalue.between(0,1).all() and d.padj_family.between(0,1).all() and d.padj_analysis_database.between(0,1).all())
  if len(tail_examples)<12:
   tail_examples.extend(d[(d.pvalue>1e-30)&(d.pvalue<.8)&(d.k>=3)].head(12-len(tail_examples)).to_dict('records'))
 check('all_contingency_and_probability_bounds',valid)
 check('full_test_count_matches_family_summary',total==s.tested_terms.sum())
 # Exact-integer combinatorics independently verifies selected R phyper tails.
 exact_ok=True
 for row in tail_examples:
  N,K,n,k=(int(row[x]) for x in ['universe_N','K','foreground_n','k'])
  numerator=sum(math.comb(K,j)*math.comb(N-K,n-j) for j in range(k,min(K,n)+1) if 0<=n-j<=N-K)
  ref=numerator/math.comb(N,n)
  exact_ok &= math.isclose(ref,row['pvalue'],rel_tol=1e-10,abs_tol=1e-100)
 check('12_independent_exact_hypergeometric_checks',len(tail_examples)==12 and exact_ok)
 for row in inventory.itertuples(index=False):
  p=ROOT/row.pdf
  check(f'{row.analysis}/{row.cell}/PDF_present',p.exists() and p.stat().st_size>2000)
 check('60_figure_packs_present',len(inventory)==60)
 validation=pd.DataFrame(checks);validation.to_csv(ROOT/'05_reports/validation_checks.csv',index=False)
 if not validation.passed.all():raise RuntimeError('Validation failed; see validation_checks.csv')
 aggregate=s.groupby(['analysis','database','direction'],as_index=False).agg(successful_families=('status',lambda x:int((x=='OK').sum())),global_significant_term_cell_pairs=('significant_global','sum'))
 aggregate.to_csv(ROOT/'03_results/summary_by_analysis_database_direction.csv',index=False,encoding='utf-8-sig')
 primary=sig[(sig.analysis=='primary_baseline_SLE_vs_HD')&sig.direction.isin(['up','down'])]
 recurrent=primary.groupby(['database','direction','term_id','Description'],as_index=False).agg(n_cells=('cell','nunique'),best_global_FDR=('padj_analysis_database','min')).sort_values(['database','n_cells','best_global_FDR'],ascending=[True,False,True])
 recurrent.to_csv(ROOT/'03_results/recurrent_primary_terms.csv',index=False,encoding='utf-8-sig')
 # Overlap is an internal sensitivity check; the cohorts are not independent.
 keys=['cell','database','direction','term_id']
 p=sig[sig.analysis=='primary_baseline_SLE_vs_HD'];n=sig[sig.analysis=='sensitivity_naive_SLE_vs_HD']
 shared=p.merge(n,on=keys,suffixes=('_primary','_naive'))
 shared.to_csv(ROOT/'03_results/shared_primary_naive_terms.csv',index=False,encoding='utf-8-sig')
 examples=recurrent[recurrent.database=='GO_BP'].head(12).copy()
 examples.best_global_FDR=examples.best_global_FDR.map(lambda x:f'{x:.3g}')
 status=s.status.value_counts().to_dict()
 report=f'''# GO / KEGG 富集分析报告

## 运行结果

完成60个原差异表达比较，每个比较分析 up、down、combined 三个基因集及 GO_BP、GO_CC、GO_MF、KEGG 四类数据库，共720个分析组。状态统计：{status}。基因数不足的组明确标记，未伪造富集结果。

共检验 {total:,} 个条目-基因集组合；严格全局FDR筛选得到 {len(sig):,} 个条目-基因集组合。这里不是独立通路数量，相同条目可出现在多个细胞或方向中。完整汇总见 `03_results/enrichment_summary.csv`。

{table(aggregate)}

主分析与未治疗敏感性分析共享 {len(shared):,} 个显著的细胞×数据库×方向×条目组合。两类分析共享受试者，仅作内部稳健性检查。

## 主分析中跨细胞重复的 GO BP 代表条目

按出现细胞数量和最小全局FDR排序；不是互相独立的生物学过程，GO上下级条目高度重叠。

{table(examples)}

## 背景和基因集

输入来自 `deseq2_analysis/03_results/all_tests_global_FDR.csv.gz`，不只使用火山图Top10。差异基因要求原分析全局FDR<0.05、|未收缩log2FC|≥1。

- up：SLE相对HD较高；down：SLE较低；combined：两方向并集。combined不用于判断通路方向。
- 背景从每个比较中 pvalue 和 padj_global 均非缺失的基因开始，再转换为唯一Entrez ID并限制到该数据库/本体中有注释的基因；没有使用全人类基因组作为背景。
- Ensembl移除版本号，PAR_Y关联同一基础基因并按Entrez只计一次。一个Ensembl映射多个Entrez的歧义映射全部排除，不挑选最有利的映射。
- 多个条目映射同一Entrez时取唯一基因计数。若同一Entrez同时满足显著上调和下调，则从方向分析排除，combined保留。每个比较的原始特征到ID关系和折叠后的基因集均在 `02_gene_sets`。
- ID转换及GO关系来自本地 org.Hs.eg.db 3.22.0；GO.db 3.22.0的GO来源日期与其一致（2025-07-22）。名称、版本和元数据已缓存。KEGG来自本次下载的官方人类REST数据并单独记录下载时间和SHA256；不用 org.Hs.eg.db 中陈旧的KEGG通路表。

## 检验和筛选

使用标准单侧超几何过度代表分析（ORA），由R基础函数 `phyper` 实现，不依赖未安装的clusterProfiler。定义 N=注释背景基因数，n=注释前景基因数，K=背景中该条目的基因数，k=前景与该条目交集；P=P(X≥k)。

每组至少5个有注释前景基因；检验背景中10≤K≤500的条目。所有合格条目都参与BH校正，包括零交集条目（P=1），避免只校正命中条目。

- padj_family：单一细胞×方向×数据库内的BH FDR。
- padj_analysis_database：在同一种分析、同一种数据库内，合并所有细胞和三个方向的全部检验后BH校正。GO_BP、GO_CC、GO_MF和KEGG分开，主分析与敏感性分析分开。
- 主要结果和图使用 padj_analysis_database<0.05，且k≥3。后者为展示/解释过滤，不改变原始检验P值。
- GeneRatio=k/n，BgRatio=K/N，FoldEnrichment=(k/n)/(K/N)。表中保留交集Entrez和基因名称。
- ORA不等同于GSEA，也不能将上调基因富集直接解释为通路被激活。当前未做基因长度/表达检出概率偏倚校正；结果是探索性功能解释。

## 图和可复现性

每个比较一份4页PDF（BP、CC、MF、KEGG），每页分up/down/combined，展示严格筛选后全局FDR最小的最多8个条目；同时提供各页PNG。无显著结果或基因不足时显示明确说明。

{len(validation)}项校验全部通过，包括背景/前景计数、P值范围、分析组完整性、完整检验数量以及12个实际检验的独立精确组合数复核。PDF渲染和视觉检查记录见 `05_reports`。原DESeq2拟合与统计结果未修改。

## 参考来源

- [Bioconductor org.Hs.eg.db](https://bioconductor.org/packages/org.Hs.eg.db/)
- [Bioconductor GO.db](https://bioconductor.org/packages/GO.db/)
- [KEGG官方API](https://www.kegg.jp/kegg/rest/keggapi.html)
'''
 (ROOT/'05_reports/REPORT.md').write_text(report,encoding='utf-8')
 guide='''# GO / KEGG 富集分析入口

先看 [结果报告](05_reports/REPORT.md) 和 `03_results/enrichment_summary.csv`。

| 目录 | 内容 |
| --- | --- |
| 00_scripts | 下载、富集、绘图、校验报告及一键运行脚本 |
| 01_databases | GO/KEGG缓存、Entrez转换、歧义/未映射基因、来源和版本 |
| 02_gene_sets | 各比较的实际前景/背景候选基因集、方向冲突和ID转换覆盖 |
| 03_results | 汇总、全部检验、显著条目、主分析/敏感性共有条目及分细胞结果 |
| 04_figures | 各比较4页PDF及BP/CC/MF/KEGG的PNG、整体概览图 |
| 05_reports | 中文报告、软件版本、校验和运行日志 |

主分析目录为 `primary_baseline_SLE_vs_HD`，未治疗敏感性分析为 `sensitivity_naive_SLE_vs_HD`；LDG_CD14_plus / minus 单独保存。

在原data目录完整重跑（本地数据库已缓存，无需再次下载）：

```powershell
.\\enrichment_analysis\\00_scripts\\run_pipeline.ps1
```

只重绘：

```powershell
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' enrichment_analysis/00_scripts/plot_enrichment.R
```

`*_all.csv.gz` 是压缩CSV，包含全部被检验条目（含P=1）；每个比较的 `significant_terms.csv` 是可直接查看的显著结果。没有显著条目时CSV仅保留表头，状态及背景人数见 `status_and_background.csv`。

`padj_family` 是单组FDR，`padj_analysis_database` 是跨细胞/方向的严格FDR。图和主要显著表使用后者；详情与阈值见报告。图的条目与交集基因保存在 `03_results/plot_top8_terms.csv`。
'''
 (ROOT/'README.md').write_text(guide,encoding='utf-8')
 manifest=[]
 for path in sorted((ROOT/'00_scripts').glob('*')):
  if path.is_file():manifest.append({'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 pd.DataFrame(manifest).to_csv(ROOT/'05_reports/code_manifest_sha256.csv',index=False)
 print(f'PASS: {len(validation)} checks; {total} tests; {len(sig)} significant results; {status}')

if __name__=='__main__':main()

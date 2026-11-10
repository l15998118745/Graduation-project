"""Prepare a concise report and verify GO/KEGG enrichment outputs."""
from pathlib import Path
import hashlib,json
import pandas as pd

ROOT=Path(__file__).resolve().parent.parent
E=ROOT/'07_enrichment'
def tbl(d):
 return '\n'.join(['| '+' | '.join(d.columns)+' |','| '+' | '.join(['---']*len(d.columns))+' |']+['| '+' | '.join(map(str,r))+' |' for r in d.itertuples(index=False,name=None)])

def main():
 allres=pd.read_csv(E/'all_enrichment_tests.csv.gz')
 bg=pd.read_csv(E/'background_and_mapping_summary.csv')
 status=pd.read_csv(E/'case_status.csv')
 counts=pd.read_csv(E/'hypothesis_counts.csv')
 provenance=json.loads((E/'enrichment_provenance.json').read_text(encoding='utf-8'))
 checks=[('all_60_comparisons',len(status)==60 and (status.status=='OK').all()),
         ('global_FDR_is_valid',allres.padj_global.between(0,1).all()),
         ('reported_hypergeom_has_hits',allres.Count.gt(0).all()),
         ('named_entrez_hits',allres.geneID.notna().all()),
         ('background_bigger_than_hits',allres.QueryMapped.ge(allres.Count).all()),
         ('pathway_tests_recorded',counts.tested_hypotheses.gt(0).all())]
 validation=pd.DataFrame(checks,columns=['check','passed']);validation.to_csv(E/'validation_checks.csv',index=False)
 if not validation.passed.all():raise RuntimeError('Enrichment validation failed')
 sig=allres[allres.padj_global<.05].copy()
 sig.to_csv(E/'significant_terms_long.csv',index=False,encoding='utf-8-sig')
 summary=allres.groupby(['analysis','database','ontology','cell','direction']).agg(
  tested_categories=('TestedTerms','max'),terms_with_hits=('ID','size'),significant_global=('padj_global',lambda x:int((x<.05).sum())),
  background=('UniverseAnnotated','max'),query_mapped=('QueryMapped','max')).reset_index()
 summary.to_csv(E/'term_counts_by_cell.csv',index=False)
 recurrent=sig.groupby(['analysis','database','ontology','ID','Description']).agg(
  cell_types=('cell','nunique'),directions=('direction',lambda x:'/'.join(sorted(set(x)))),
  min_global_FDR=('padj_global','min'),cells=('cell',lambda x:'; '.join(sorted(set(x))))).reset_index().sort_values(['analysis','database','min_global_FDR'])
 recurrent.to_csv(E/'recurrent_significant_terms.csv',index=False,encoding='utf-8-sig')
 # Freeze exact downloaded and locally installed knowledgebase versions.
 for path in [E/'00_cache/go-basic.obo',E/'00_cache/kegg_pathways_hsa.tsv',E/'00_cache/kegg_gene_pathway_hsa.tsv']:
  h=hashlib.sha256(path.read_bytes()).hexdigest()
  provenance.setdefault('sha256_files',{})[path.name]=h
 (E/'enrichment_provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2),encoding='utf-8')
 counts_by=sig.groupby(['analysis','database']).size().rename('global_FDR05_terms_with_hits').reset_index()
 top=sig.sort_values('padj_global').groupby(['analysis','database']).head(10)
 topdisp=top[['analysis','database','ontology','cell','direction','ID','Description','Count','padj_global']].copy()
 topdisp.padj_global=topdisp.padj_global.map(lambda x:f'{x:.3g}')
 bgdisp=bg.groupby('analysis').agg(
  median_tested_protein_coding=('tested_protein_coding_ensembl','median'),median_Entrez_background=('mapped_background_entrez','median'),
  median_GO_annotated_background=('GO_annotated_background','median'),median_KEGG_annotated_background=('KEGG_annotated_background','median')).reset_index()
 rpt=f'''# GO 与 KEGG 富集分析报告

## 结果概览

每项比较按上调/下调分别做单侧超几何 ORA；输入为基因水平分析通过全局 DE FDR < 0.05、|log2FC| >= 1 的 **protein-coding 基因**。背景为同一分析/细胞内有非缺失 DESeq2 p 值且 protein-coding 的基因，经人类 Ensembl->Entrez 转换后，在相应数据库中有注释的 ID。低背景或未注释基因不加入该数据库背景。每组 pathway/GO term 至少10个、至多500个背景基因。每个 cell/direction 的局部 BH 与每个数据库、每种疾病分析中跨全部细胞和上下调 term 检验的全局 BH 同时输出；下列显著结果按全局 FDR < 0.05 定义。超几何全测试量包含无重叠 term，但为控制文件体积，长表仅保留 Count>0 的 term。分类多重检验总数见 hypothesis_counts.csv。

按上下调拆分后共运行 120 个 gene-list 对照；有数据库映射的背景与 DEG 数按比较汇总于 background_and_mapping_summary.csv。

{tbl(bgdisp)}

## 全局显著结果数量

{tbl(counts_by)}

基因列表较小且跨细胞类型统一校正严格，若某类没有全局显著 term，其结果目录仍含该比较所有有重叠 term 的表格，图中会注明没有通过全局 FDR 的结果。不要把空图解释为没有任何生物过程参与。

## 全局 FDR 最小的代表 term

{tbl(topdisp)}

以上取每个分析/数据库排名靠前的代表项；更多结果见显著长表和跨细胞重复 term 表。相似 GO term 经常共享基因，不能视为彼此独立的机制。KEGG pathway 图中的富集仅反映成员基因超额出现，不直接代表通路活性增强/减弱。

## 可复核文件

- `GO_global_FDR05.csv` 和 `KEGG_global_FDR05.csv`：跨细胞类型校正后显著项。
- `GO_all_tests.csv.gz` 和 `KEGG_all_tests.csv.gz`：有重叠的完整 term-level 结果。
- `significant_terms_long.csv`：两库的显著结果合并长表。
- `term_counts_by_cell.csv`：每比较、每方向、每 ontology 的检验/命中/显著 term 数和背景。
- `recurrent_significant_terms.csv`：跨多个细胞/方向出现的显著 term。
- `background_and_mapping_summary.csv`、`gene_id_mapping_all_comparisons.csv.gz`：每个比较的背景和映射审计。
- `00_cache/`：GO 本体、KEGG API 原始缓存及版本/哈希。
- `validation_checks.csv`：本次校验。

## 方法与版本

{provenance['GO_ontology']}; {provenance['GO_gene_annotations']}; {provenance['KEGG']}。算法为标准单侧超几何过度富集检验，方向为上尾；全局 BH 对每个疾病分析中的 GO 全 ontology 或 KEGG 所有细胞、方向和 term 统一调整。local FDR 仅供本比较筛选参考。

注意：这是基于 OrgDb 当前安装版本 2025 年注释、GO basic OBO 和下载时 KEGG 通路表的富集结果。Ensembl ID 通过本地 OrgDb 映射为 Entrez；PAR_Y 使用基础 Ensembl ID 查询后保留原 ID。映射比例和逐基因映射表均已保存。富集不包含 lncRNA、假基因或其他非 protein-coding 条目。主/敏感性分析各自校正，没有跨两种研究问题合并调整。

**整体富集不能单独证明因果机制。** GO 和 KEGG 版本、基因映射、基因集筛选门槛或上下调合并方式变化，都可能改变结果；报告同时保留全检验表和数据库快照，便于复现。
'''
 (E/'ENRICHMENT_REPORT.md').write_text(rpt,encoding='utf-8')
 print(f'PASS: {len(validation)} checks; {len(status)} comparisons; {len(allres)} hit terms; {len(sig)} globally significant terms.')
 print(counts_by.to_string(index=False))

if __name__=='__main__':main()

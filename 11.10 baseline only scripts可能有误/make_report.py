"""Run after summarize_results.R in the baseline-only analysis."""
from pathlib import Path
import hashlib
import re
import pandas as pd

script_root = Path(__file__).resolve().parent
root = script_root.parent / '.work'
data_root = script_root.parent.parent
manifest = []
for path in [*data_root.glob('gencode*.csv'),data_root/'metadata anno_final_aligned',script_root/'run_deseq2.R',script_root/'summarize_results.R',script_root/'make_report.py']:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(1024*1024),b''):
            h.update(block)
    manifest.append({'file':path.name,'bytes':path.stat().st_size,'sha256':h.hexdigest()})
pd.DataFrame(manifest).to_csv(root/'input_code_manifest_sha256.csv',index=False)
s = pd.read_csv(root / 'analysis_summary_with_global_FDR.csv')
m = pd.read_csv(root / 'metadata_audited.csv')
q = pd.read_csv(root / 'library_qc_all.csv')
g = pd.read_csv(root / 'DEG_global_FDR05_absLFC1.csv')
v = pd.read_csv(root / 'validation_checks.csv')
full = pd.read_csv(root/'all_tests_global_FDR.csv.gz')
integrity = []
for row in s.itertuples(index=False):
    if row.status != 'OK':
        continue
    slug = re.sub(r'[^A-Za-z0-9]+','_',row.cell.replace('-','_minus').replace('+','_plus'))
    sample_file = root/row.analysis/slug/'samples.csv'
    md = pd.read_csv(sample_file)
    integrity.append(not md.patient.duplicated().any() and md.identity_ok.all() and md.qc_pass.all())
extra = pd.DataFrame([
    {'check':'all_30_primary_patient_and_QC_checks','passed':len(integrity)==30 and all(integrity)},
    {'check':'all_30_primary_result_keys_unique','passed':full.groupby(['analysis','cell']).ngroups==30 and not full.duplicated(['analysis','cell','Geneid']).any()},
])
v = pd.concat([v[~v.check.isin(extra.check)],extra],ignore_index=True)
v.to_csv(root/'validation_checks.csv',index=False)
if not v.passed.all():
    raise RuntimeError('Integrity check failed; inspect validation_checks.csv')

def table(df):
    return '\n'.join(['| ' + ' | '.join(map(str, df.columns)) + ' |',
                      '| ' + ' | '.join(['---'] * len(df.columns)) + ' |'] +
                     ['| ' + ' | '.join(map(str, row)) + ' |' for row in df.itertuples(index=False, name=None)])

primary = s[s.analysis == 'primary_baseline_SLE_vs_HD'].copy()
display = primary[['cell', 'status', 'HD', 'SLE', 'genes_retained', 'global_up', 'global_down', 'global_DEG']].copy()
for col in ['global_up','global_down','global_DEG']:
    display[col] = display[col].fillna(0).astype(int)
display = display.sort_values('global_DEG',ascending=False)
display.to_csv(root/'primary_baseline_summary.csv',index=False)
lines = [
    '# 原始计数预处理和 DESeq2 结果报告', '',
    '## 分析状态', '',
    f'计划30项 baseline SLE vs HD 主分析，结果表记录{s.shape[0]}项；成功{(s.status == "OK").sum()}项。详见 validation_checks.csv。', '',
    table(v), '',
    '## 预处理', '',
    f'- 原始矩阵60,649个基因、{len(q):,}个文库；其中{q.has_metadata.sum():,}个具有注释。',
    f'- 无注释文库{(~q.has_metadata).sum()}个，未自行推断分组。',
    f'- 有注释文库中低信息量筛除{(~m.qc_pass).sum()}条；身份映射不明确{(~m.identity_ok).sum()}条；两者可能重叠。',
    f'- 同时通过两项检查的记录{(m.qc_pass & m.identity_ok).sum():,}条，随后按分析人群、细胞和每患者一条记录继续筛选。',
    '- 所有过滤均使用原始整数计数；不使用已有小数表达矩阵进行DESeq2检验。', '',
    '## 比较和判定', '',
    '主分析：仅纳入 Timepoint="baseline" 的 SLE，与全部 HD 比较；其他 SLE 时点（包括 "-"）均不分析。每种细胞内每患者仅保留一条记录。模型控制性别及由ID日期提取的月份代理批次。', '',
    '下表采用主分析内跨基因、跨细胞类型的全局BH FDR<0.05，并额外要求 |未收缩log2FC|≥1。up/down均指SLE相对HD。', '',
    table(display), '',
    f'主分析中共{g[g.analysis == "primary_baseline_SLE_vs_HD"].Geneid.nunique():,}个不同基因在至少一个细胞/组分满足上述筛选条件（含EV和Plasma）。', '',
    '## 解释边界', '',
    '- EV、Plasma是非细胞组分；28种细胞和2种组分的结果应分别解释。',
    '- 身份冲突采用整组保守排除，未自动纠正患者编号；这会损失样本、改变分析人群。',
    '- 月份不是经过确认的实验批次；年龄、真实批次、测序质量等信息缺失，不能排除混杂。',
    '- HD 的 Timepoint 均记录为 "-"，表示健康对照不适用随访时点；SLE 中 Timepoint="-" 及所有非 baseline 随访样本均已排除。',
    '- 没有匹配GENCODE版本的GTF，未做基因符号/基因类型映射，也未作GO/KEGG富集；不能把全部行直接称为蛋白编码mRNA。',
    '- 用户已确认每列为分选细胞群的计数文库，不是单细胞；当前统计单位适用，无需单细胞pseudobulk汇总。',
    '- 未进行LFC收缩，效应大小特别是低计数基因应谨慎解读。未显著不代表不存在差异。', '',
    '## 其他问题的可行性', '',
    'paired_analysis_feasibility.csv提供按细胞和治疗分层的baseline–6M配对数。清洁数据多数细胞belimumab仅4对、Telitacicept仅2对；治疗效应暂不作为主要分析。', '',
    '后续可在SLE内部分析SLEDAI连续关联、肾脏等器官受累，先明确混杂和患者身份。不能把HD缺失的SLEDAI填0，也不能根据文库列数推断细胞丰度。', '',
    '方法、运行命令、所有排除规则及结果字段说明见README.md；完整模型和标准化计数按比较保存。',
]
(root/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(display.to_string(index=False))
print('Report:',root/'REPORT.md')

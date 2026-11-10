"""Group generated files without moving the three original input files."""
from pathlib import Path
import argparse
import csv
import hashlib
import os

ROOT = Path(__file__).resolve().parent
if ROOT.name == '00_scripts':
    ROOT = ROOT.parent
CATEGORIES = ['00_scripts','01_annotation','02_quality_control','03_results','04_figures','05_models','06_reports']
ANALYSES = ['primary_baseline_SLE_vs_HD']

def organize(source=None):
    source = Path(source).resolve() if source else ROOT
    for name in CATEGORIES:
        (ROOT/name).mkdir(exist_ok=True)
    moves = []
    def move(src, dst):
        if not src.exists() or src.resolve() == dst.resolve():
            return
        if not src.resolve().is_relative_to(ROOT) or not dst.resolve().is_relative_to(ROOT):
            raise ValueError('Path outside analysis directory')
        dst.parent.mkdir(parents=True, exist_ok=True)
        size=src.stat().st_size
        os.replace(src,dst)
        moves.append((str(src.relative_to(ROOT)),str(dst.relative_to(ROOT)),size))
    qc_names={'identity_exclusions.csv','library_qc_all.csv','metadata_audited.csv','raw_without_metadata.csv'}
    for src in list(source.iterdir()):
        if src.is_file():
            if src.name=='README.md' and source==ROOT:
                continue
            if src.name in {'organize_outputs.py','annotate_genes.py','plot_annotated.R'} or src.suffix in {'.R','.py','.ps1'}:
                dst=ROOT/'00_scripts'/src.name
            elif src.name in qc_names: dst=ROOT/'02_quality_control'/src.name
            elif src.name=='aligned_raw_counts.rds': dst=ROOT/'05_models'/src.name
            elif src.suffix in {'.png','.pdf'}: dst=ROOT/'04_figures'/'overview'/src.name
            elif src.suffix=='.log': dst=ROOT/'06_reports'/'logs'/src.name
            elif src.suffix in {'.md','.txt'} or 'manifest' in src.name or src.name=='validation_checks.csv': dst=ROOT/'06_reports'/src.name
            elif src.name.endswith(('.csv','.csv.gz')): dst=ROOT/'03_results'/src.name
            else: continue
            move(src,dst)
    for an in ANALYSES:
        old=source/an
        if not old.exists(): continue
        for cell in list(old.iterdir()):
            if not cell.is_dir():continue
            for src in list(cell.iterdir()):
                if src.suffix=='.rds':dst=ROOT/'05_models'/an/cell.name/src.name
                elif src.name=='diagnostics.pdf':dst=ROOT/'04_figures'/'diagnostics'/an/(cell.name+'.pdf')
                elif src.name in {'samples.csv','batch_by_condition.csv','gene_filter.csv','repeated_patient_exclusions.csv','PCA_scores.csv'}:dst=ROOT/'02_quality_control'/an/cell.name/src.name
                else:dst=ROOT/'03_results'/an/cell.name/src.name
                move(src,dst)
            cell.rmdir()
        old.rmdir()
    for an in ANALYSES:
        for src in (ROOT/'03_results'/an).glob('*/PCA_scores.csv'):
            move(src,ROOT/'02_quality_control'/an/src.parent.name/src.name)
    log=ROOT/'06_reports'/'file_migration.csv'
    with log.open('a',encoding='utf-8',newline='') as f:
        w=csv.writer(f)
        if f.tell()==0:w.writerow(['old_path','new_path','bytes'])
        w.writerows(moves)
    print(f'Organized {len(moves)} files')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source');args=parser.parse_args()
    organize(args.source)

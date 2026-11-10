"""Annotate every gene-indexed CSV using UCSC's frozen GENCODE v38 table."""
from pathlib import Path
import csv,gzip,re,json,hashlib,datetime,urllib.request,socket
import pandas as pd

ROOT=Path(__file__).resolve().parent
if ROOT.name=='00_scripts': ROOT=ROOT.parent
DB=ROOT/'01_annotation'
URL='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeAttrsV38.txt.gz'
GTF=DB/'wgEncodeGencodeAttrsV38.txt.gz'
COLS=['gene_symbol','gene_type','annotation_match','annotation_source']

def sha256(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def make_annotation():
 DB.mkdir(exist_ok=True)
 if not GTF.exists() or GTF.stat().st_size==0:
  socket.setdefaulttimeout(60)
  part=GTF.with_suffix('.gz.part')
  urllib.request.urlretrieve(URL,part);part.replace(GTF)
 db=pd.read_csv(GTF,sep='\t',header=None,usecols=[0,1,2],names=['Geneid','gene_symbol','gene_type'],dtype=str,keep_default_na=False).drop_duplicates()
 if db.Geneid.duplicated().any():raise ValueError('Duplicate GTF gene IDs')
 # Preserve PAR_Y while removing only the version number.
 db['stable_id']=db.Geneid.str.replace(r'\.\d+(?=_|$)','',regex=True)
 if db.stable_id.duplicated().any():raise ValueError('Ambiguous unversioned mapping')
 raw=next(ROOT.parent.glob('gencode*.csv'))
 ids=pd.read_csv(raw,usecols=[0]).iloc[:,0]
 exact=db.set_index('Geneid');stable=db.set_index('stable_id')
 mapped=[]
 for gid in ids:
  key=re.sub(r'\.\d+(?=_|$)','',gid)
  if gid in exact.index:r=exact.loc[gid];match='exact_version'
  elif gid.endswith('_PAR_Y') and gid.removesuffix('_PAR_Y') in exact.index:
   r=exact.loc[gid.removesuffix('_PAR_Y')];match='PAR_Y_base_gene_exact_version'
  elif key in stable.index:r=stable.loc[key];match='stable_id_only'
  else:r=None;match='unmapped'
  mapped.append({'Geneid':gid,**{c:(r[c] if r is not None else '') for c in COLS[:2]},'annotation_match':match,'annotation_source':'GENCODE v38 / UCSC hg38'})
 ann=pd.DataFrame(mapped)
 ann.to_csv(DB/'gene_annotation.csv',index=False,encoding='utf-8-sig')
 ann[ann.annotation_match=='unmapped'].to_csv(DB/'unmapped_genes.csv',index=False)
 ann.groupby('gene_type').size().rename('n_genes').to_csv(DB/'gene_type_counts.csv')
 provenance={'source_url':URL,'schema_url':URL.replace('.txt.gz','.sql'),'release':'GENCODE v38','assembly':'GRCh38.p13 / hg38','retrieved_or_checked_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'database_sha256':sha256(GTF),'raw_genes':len(ann),'database_genes':len(db),'mapping_counts':ann.annotation_match.value_counts().to_dict(),'notes':'Gene names and biotypes from frozen GENCODE v38 UCSC attributes. PAR_Y copies share the exact-version base gene name, retain their original ID, and are marked separately; no chromosome coordinates are inferred. Names are not asserted to be current HGNC approved symbols.'}
 (DB/'annotation_provenance.json').write_text(json.dumps(provenance,indent=2,ensure_ascii=False),encoding='utf-8')
 return ann,provenance

def main():
 ann,provenance=make_annotation();ann=ann.set_index('Geneid')
 changed=[];checks=[]
 for parent in [ROOT/'03_results',ROOT/'02_quality_control']:
  for path in sorted(parent.rglob('*')):
   if not path.is_file() or not path.name.endswith(('.csv','.csv.gz')):continue
   head=pd.read_csv(path,nrows=0)
   if 'Geneid' not in head.columns:continue
   df=pd.read_csv(path);old=df.drop(columns=COLS,errors='ignore')
   joined=ann.reindex(df.Geneid)[COLS].reset_index(drop=True)
   result=pd.concat([old.reset_index(drop=True),joined],axis=1)
   # In-place annotation must not change row order, cardinality, or ANY original value.
   pd.testing.assert_frame_equal(old.reset_index(drop=True),result[old.columns],check_exact=True)
   result.to_csv(path,index=False,encoding='utf-8-sig',compression='gzip' if path.name.endswith('.gz') else None)
   changed.append({'file':str(path.relative_to(ROOT)),'rows':len(result),'named':int(result.gene_symbol.notna().sum())})
 pd.DataFrame(changed).to_csv(DB/'annotated_file_inventory.csv',index=False)
 print(json.dumps(provenance,ensure_ascii=True,indent=2));print('Annotated files:',len(changed))

if __name__=='__main__':main()

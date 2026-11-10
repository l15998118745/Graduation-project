"""Cache the two small, official KEGG human pathway tables used by ORA."""
from pathlib import Path
from urllib.request import Request,urlopen
import hashlib,json,time,datetime

ROOT=Path(__file__).resolve().parent.parent
CACHE=ROOT/'07_enrichment'/'00_cache';CACHE.mkdir(parents=True,exist_ok=True)
ENDPOINTS={'kegg_pathways_hsa.tsv':'https://rest.kegg.jp/list/pathway/hsa','kegg_gene_pathway_hsa.tsv':'https://rest.kegg.jp/link/pathway/hsa'}
manifest=[]
for name,url in ENDPOINTS.items():
 p=CACHE/name
 if not p.exists() or p.stat().st_size==0:
  req=Request(url,headers={'User-Agent':'Academic-RNAseq-enrichment/1.0'})
  with urlopen(req,timeout=45) as r:
   data=r.read()
  if len(data)<1000:raise RuntimeError(f'Incomplete KEGG response: {name}')
  p.write_bytes(data);time.sleep(1.1)
 raw=p.read_bytes()
 if len(raw.splitlines())<100:raise RuntimeError(f'Unexpected KEGG table: {name}')
 manifest.append({'file':name,'url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()})
(CACHE/'kegg_provenance.json').write_text(json.dumps({'api':'KEGG REST','tables':manifest,'license_note':'KEGG API use is limited to academic users belonging to academic institutions; this project is an academic thesis analysis.'},indent=2),encoding='utf-8')
print(json.dumps(manifest,indent=2))

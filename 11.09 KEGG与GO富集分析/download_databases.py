"""Download and freeze GO/KEGG resources; subsequent runs reuse local files."""
from pathlib import Path
import urllib.request,socket,tarfile,shutil,json,hashlib,datetime
ROOT=Path(__file__).resolve().parent.parent
DB=ROOT/'01_databases';DB.mkdir(parents=True,exist_ok=True)
SOURCES={
 'kegg_hsa_pathway_names.tsv':'https://rest.kegg.jp/list/pathway/hsa',
 'kegg_hsa_gene_pathway.tsv':'https://rest.kegg.jp/link/pathway/hsa',
 'kegg_hsa_entrez.tsv':'https://rest.kegg.jp/conv/ncbi-geneid/hsa',
 'GO.db_3.22.0.tar.gz':'https://bioconductor.posit.co/packages/3.22/data/annotation/src/contrib/GO.db_3.22.0.tar.gz',
}
socket.setdefaulttimeout(30)
for name,url in SOURCES.items():
 path=DB/name
 if not path.exists() or not path.stat().st_size:
  part=path.with_suffix(path.suffix+'.part');urllib.request.urlretrieve(url,part);part.replace(path)
if not (DB/'GO.sqlite').exists():
 with tarfile.open(DB/'GO.db_3.22.0.tar.gz') as tf:
  name=next(n for n in tf.getnames() if n.endswith('/GO.sqlite'))
  with tf.extractfile(name) as src,(DB/'GO.sqlite').open('wb') as dst:shutil.copyfileobj(src,dst)
manifest=[]
for name,url in SOURCES.items():
 path=DB/name
 manifest.append({'file':name,'source_url':url,'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
target=DB/'download_manifest.json'
if not target.exists():
 target.write_text(json.dumps({'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':manifest},indent=2),encoding='utf-8')
print('All GO/KEGG database files are cached')

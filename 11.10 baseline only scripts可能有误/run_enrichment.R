#!/usr/bin/env Rscript
# GO and KEGG over-representation analysis (ORA); run from the data directory.
if(.Platform$OS.type=='windows')Sys.setlocale('LC_ALL','Chinese_China.utf8')
suppressPackageStartupMessages({library(data.table);library(AnnotationDbi);library(org.Hs.eg.db);library(ggplot2)})
root<-'deseq2_baseline_only_analysis';out<-file.path(root,'07_enrichment')
dir.create(file.path(out,'GO_results'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,'KEGG_results'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,'GO_figures'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,'KEGG_figures'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,'00_cache'),recursive=TRUE,showWarnings=FALSE)
if(!requireNamespace('KEGGREST',quietly=TRUE))stop('KEGGREST is required')
cache<-file.path(out,'00_cache');obo<-file.path(cache,'go-basic.obo')
stopifnot(file.exists(obo),file.info(obo)$size>1000000)
sha256<-function(path){con<-file(path,'rb');on.exit(close(con));digest::digest(con,algo='sha256')}
if(!requireNamespace('digest',quietly=TRUE))stop('R package digest is required to record database checksums')
if(!requireNamespace('jsonlite',quietly=TRUE))stop('R package jsonlite is required to record database provenance')
message('Reading Gene Ontology basic ontology')
lines<-readLines(obo,warn=FALSE,encoding='UTF-8');sections<-which(lines=='[Term]');go<-list()
for(i in seq_along(sections)){
 st<-sections[i]+1L;en<-if(i<length(sections))sections[i+1L]-1L else length(lines)
 b<-lines[st:en];id<-sub('^id: ','',b[grepl('^id: GO:',b)][1]);
 if(is.na(id)||!nzchar(id))next
 nm<-sub('^name: ','',b[grepl('^name: ',b)][1]);ns<-sub('^namespace: ','',b[grepl('^namespace: ',b)][1])
 obs<-any(grepl('^is_obsolete: true',b))
 if(!obs&&nzchar(nm)&&nzchar(ns))go[[id]]<-c(name=nm,namespace=ns)
}
go<-rbindlist(lapply(names(go),function(id)data.table(GOID=id,Description=unname(go[[id]]['name']),namespace=unname(go[[id]]['namespace']))))
paths<-fread(file.path(cache,'kegg_pathways_hsa.tsv'),header=FALSE,col.names=c('pathway','Description'))
links<-fread(file.path(cache,'kegg_gene_pathway_hsa.tsv'),header=FALSE,col.names=c('gene','pathway'))
paths[,ID:=sub('^path:','',pathway)];paths[,pathway:=NULL]
links[,ENTREZID:=sub('^hsa:','',gene)];links[,ID:=sub('^path:','',pathway)]
links<-unique(links[,.(ENTREZID,ID)])
paths<-unique(paths[,.(ID,Description)])
kegg_sets<-split(links$ENTREZID,links$ID)
kegg_sets<-lapply(kegg_sets,unique)
all_deg<-fread(file.path(root,'03_results','DEG_global_FDR05_absLFC1.csv'))
s<-fread(file.path(root,'03_results','analysis_summary_with_global_FDR.csv'))
analysis_names<-c(primary_baseline_SLE_vs_HD='Baseline-only SLE vs HD')
strip_version<-function(x)sub('\\.[0-9]+$','',sub('_PAR_Y$','',x))
map_ensembl<-function(ids){
 ids<-unique(strip_version(ids));
 if(!length(ids))return(data.table(ensembl=character(),ENTREZID=character()))
 z<-mapIds(org.Hs.eg.db,keys=ids,keytype='ENSEMBL',column='ENTREZID',multiVals='first');
 data.table(ensembl=names(z),ENTREZID=unname(as.character(z)))
}
map_rows<-list();universe_records<-list();ora<-list();hypotheses<-list();case_log<-list();kegg_ids<-sort(unique(links$ENTREZID))
message('Mapping tested protein-coding backgrounds and DEGs to Entrez IDs')
for(i in seq_len(nrow(s))){
 if(s$status[i]!='OK')next
 an<-s$analysis[i];cell<-s$cell[i];tag<-gsub('[^A-Za-z0-9]+','_',gsub('+','_plus',gsub('-','_minus',cell,fixed=TRUE),fixed=TRUE))
 f<-file.path(root,'03_results',an,tag,'DESeq2_all_genes.csv')
 r<-fread(f);tested<-r[gene_type=='protein_coding' & !is.na(pvalue),Geneid]
  cell0<-cell
  hit<-all_deg[analysis==an & cell==cell0 & gene_type=='protein_coding']
 if(!length(tested)){case_log[[length(case_log)+1L]]<-data.table(analysis=an,cell=cell,status='no_tested_protein_coding_genes');next}
 u_map<-map_ensembl(tested);d_map<-map_ensembl(hit$Geneid)
 # Add a readable gene label to each ID mapping while retaining the Ensembl ID.
 names_by_stable<-unique(data.table(ensembl=strip_version(r$Geneid),symbol=r$gene_symbol))
 u_map<-merge(u_map,names_by_stable,by='ensembl',all.x=TRUE,sort=FALSE)
 names_hit<-unique(data.table(ensembl=strip_version(hit$Geneid),symbol=hit$gene_symbol))
 d_map<-merge(d_map,names_hit,by='ensembl',all.x=TRUE,sort=FALSE)
 u_map<-u_map[!is.na(ENTREZID)&nzchar(ENTREZID)];d_map<-d_map[!is.na(ENTREZID)&nzchar(ENTREZID)]
 u_ids<-unique(u_map$ENTREZID);d_ids<-unique(d_map$ENTREZID)
 # retain only IDs with annotations in the appropriate resource.
 kegg_u<-intersect(u_ids,kegg_ids);kegg_d<-intersect(d_ids,kegg_ids)
 go_raw<-AnnotationDbi::select(org.Hs.eg.db,keys=u_ids,keytype='ENTREZID',columns=c('GOALL'))
 go_raw<-unique(as.data.table(go_raw)[!is.na(GOALL),.(ENTREZID=as.character(ENTREZID),GOID=GOALL)])
 go_raw<-merge(go_raw,go,by='GOID',all=FALSE)
 go_sets<-split(go_raw$ENTREZID,go_raw$GOID);go_sets<-lapply(go_sets,unique)
 go_bg<-unique(go_raw$ENTREZID)
 # Report the Ensembl-to-Entrez conversion for the tested universe and selected hits.
 map_rows[[length(map_rows)+1L]]<-rbindlist(list(
   u_map[,.(analysis=an,cell=cell0,set='background',ensembl,ENTREZID,symbol)],
   d_map[,.(analysis=an,cell=cell0,set='DEG',ensembl,ENTREZID,symbol)]),fill=TRUE)
 universe_records[[length(universe_records)+1L]]<-data.table(analysis=an,cell=cell,
  tested_protein_coding_ensembl=length(tested),mapped_background_entrez=length(u_ids),GO_annotated_background=length(go_bg),KEGG_annotated_background=length(kegg_u),
  up_protein_coding_ensembl=uniqueN(hit[log2FoldChange>0]$Geneid),down_protein_coding_ensembl=uniqueN(hit[log2FoldChange<0]$Geneid),mapped_DEG_entrez=length(d_ids))
 for(direction in c('up','down')){
  raw_hits<-hit[if(direction=='up')log2FoldChange>0 else log2FoldChange<0]
  hmap<-map_ensembl(raw_hits$Geneid)
  hmap<-merge(hmap,names_by_stable,by='ensembl',all.x=TRUE,sort=FALSE)
  genes<-intersect(unique(hmap$ENTREZID[!is.na(hmap$ENTREZID)]),go_bg)
  make_ora<-function(sets,univ,query,term_table,db,ontology=''){
   sets<-lapply(sets,function(x)intersect(unique(as.character(x)),univ))
   sizes<-lengths(sets);valid<-sizes>=10 & sizes<=500
   ids<-names(sets)[valid];if(!length(ids)){hypotheses[[paste(an,db,sep='|')]]<<-c(hypotheses[[paste(an,db,sep='|')]],0L);return(data.table())}
   q<-intersect(unique(query),univ);N<-length(univ);nq<-length(q);k<-vapply(sets[ids],function(x)sum(x %in% q),integer(1));M<-sizes[valid]
   pv<-phyper(k-1,M,N-M,nq,lower.tail=FALSE);hypotheses[[paste(an,db,sep='|')]]<<-c(hypotheses[[paste(an,db,sep='|')]],length(ids))
   # Include zero-overlap categories in BH; omit only their bulky display rows.
   keep<-k>0;if(!any(keep))return(data.table())
   tt<-term_table[match(ids[keep],term_table$ID)]
   count_ids<-lapply(ids[keep],function(id)intersect(sets[[id]],q))
   ens<-lapply(count_ids,function(e)paste(unique(hmap[ENTREZID %in% e]$ensembl),collapse='/'))
   syms<-lapply(count_ids,function(e){z<-unique(hmap[ENTREZID %in% e & !is.na(symbol)&nzchar(symbol)]$symbol);if(length(z))paste(z,collapse='/')else paste(e,collapse='/')})
   data.table(analysis=an,cell=cell,direction=direction,database=db,ontology=ontology,ID=ids[keep],Description=tt$Description,
      Count=k[keep],GeneRatio=paste0(k[keep],'/',nq),BgRatio=paste0(M[keep],'/',N),
      RichFactor=k[keep]/M[keep],FoldEnrichment=(k[keep]/nq)/(M[keep]/N),pvalue=pv[keep],
      p.adjust_set=p.adjust(pv,method='BH')[keep],geneID=unlist(syms),ensemblID=unlist(ens),
      UniverseAnnotated=N,QueryMapped=nq,TestedTerms=length(ids))
  }
  for(ont in c('BP','CC','MF')){
   ns<-switch(ont,BP='biological_process',CC='cellular_component',MF='molecular_function')
   go_ont<-go_raw[namespace==ns]
   gs<-split(go_ont$ENTREZID,go_ont$GOID);gs<-lapply(gs,unique)
   go_bg_ont<-unique(go_ont$ENTREZID)
   terms<-go[match(names(gs),GOID),.(ID=GOID,Description)]
   if(length(gs))ora[[length(ora)+1L]]<-make_ora(gs,go_bg_ont,genes,terms,'GO',ont)
  }
  khits<-intersect(unique(hmap$ENTREZID[!is.na(hmap$ENTREZID)]),kegg_u)
  kt<-paths[ID %in% names(kegg_sets),.(ID,Description)]
  ora[[length(ora)+1L]]<-make_ora(kegg_sets,kegg_u,khits,kt,'KEGG','KEGG')
 }
 case_log[[length(case_log)+1L]]<-data.table(analysis=an,cell=cell,status='OK',tested_protein_coding=length(tested),mapped_background=length(u_ids),GO_background=length(go_bg),KEGG_background=length(kegg_u),mapped_DEG=length(d_ids))
 message(an,' | ',cell,' | background ',length(tested),' tested protein-coding genes')
}
if(length(map_rows))fwrite(rbindlist(map_rows,fill=TRUE),file.path(out,'00_cache','gene_id_mapping_all_comparisons.csv.gz'))
if(!length(ora))stop('No enrichment tests were generated')
res<-rbindlist(ora,fill=TRUE)
testing_table<-rbindlist(lapply(names(hypotheses),function(k){v<-strsplit(k,'\\|')[[1]];data.table(analysis=v[1],database=v[2],tested_hypotheses=sum(hypotheses[[k]]))}))
totals<-testing_table[,.(total_hypotheses=sum(tested_hypotheses)),by=.(analysis,database)]
res[,padj_global:=p.adjust(pvalue,method='BH',n=totals[.BY$analysis==analysis & .BY$database==database,total_hypotheses][1]),by=.(analysis,database)]
setorder(res,analysis,database,padj_global,pvalue)
# One combined result plus typed, easy-to-browse tables.
fwrite(res,file.path(out,'all_enrichment_tests.csv.gz'))
for(db in c('GO','KEGG')){
 sub<-res[database==db]
 fwrite(sub,file.path(out,paste0(db,'_all_tests.csv.gz')))
 fwrite(sub[!is.na(padj_global)&padj_global<.05],file.path(out,paste0(db,'_global_FDR05.csv')),bom=TRUE)
 sig<-sub[!is.na(padj_global)&padj_global<.05]
 if(nrow(sig))fwrite(sig[,.(analysis,cell,direction,ontology,ID,Description,Count,GeneRatio,BgRatio,RichFactor,FoldEnrichment,pvalue,p.adjust_set,padj_global,geneID,ensemblID,UniverseAnnotated,QueryMapped,TestedTerms)],file.path(out,paste0(db,'_significant_named_terms.csv')),bom=TRUE)
}
# A consistent per-cell four-panel file: GO BP, CC, MF and KEGG; rows show both directions.
plot_data<-res[!is.na(padj_global)&padj_global<.05]
for(an in unique(s$analysis))for(cell in unique(s$cell)){
 tag<-gsub('[^A-Za-z0-9]+','_',gsub('+','_plus',gsub('-','_minus',cell,fixed=TRUE),fixed=TRUE))
 for(db in c('GO','KEGG')){
  od<-file.path(out,if(db=='GO')'GO_results' else 'KEGG_results',an,tag);dir.create(od,recursive=TRUE,showWarnings=FALSE)
  cell0<-cell;x<-res[analysis==an&cell==cell0&database==db]
  fwrite(x,file.path(od,paste0(db,'_all_tests.csv')),bom=TRUE)
  fwrite(x[!is.na(padj_global)&padj_global<.05],file.path(od,paste0(db,'_significant.csv')),bom=TRUE)
 }
 for(db in c('GO','KEGG')){
  cell0<-cell;x<-plot_data[analysis==an&cell==cell0&database==db]
  if(db=='GO')x<-x[ontology %in% c('BP','CC','MF')]
  fig<-file.path(out,if(db=='GO')'GO_figures' else 'KEGG_figures',an);dir.create(fig,recursive=TRUE,showWarnings=FALSE)
  if(nrow(x)){
   x[,label:=paste0(Description,' (',direction,')')];x<-x[order(padj_global)][1:min(.N,30)]
   p<-ggplot(x,aes(FoldEnrichment,reorder(label,FoldEnrichment)))+geom_point(aes(size=Count,color=-log10(pmax(padj_global,1e-300))))+
    facet_grid(ontology~direction,scales='free_y',space='free_y')+scale_color_viridis_c(option='D',name='-log10(global FDR)')+theme_bw(base_size=9)+
    theme(axis.text.y=element_text(size=7),strip.text=element_text(size=8),panel.grid.minor=element_blank())+
    labs(title=paste(db,analysis_names[[an]],cell,sep=' | '),x='Fold enrichment',y=NULL,size='Gene count')
  }else{
   p<-ggplot()+annotate('text',x=0,y=0,label='No pathway / term passes global FDR < 0.05\nAll tested categories are retained in the tables.',size=5)+
    xlim(-1,1)+ylim(-1,1)+theme_void()+labs(title=paste(db,analysis_names[[an]],cell,sep=' | '))
  }
  ggsave(file.path(fig,paste0(tag,'.pdf')),p,width=12,height=ifelse(db=='GO',8,5))
  ggsave(file.path(fig,paste0(tag,'.png')),p,width=12,height=ifelse(db=='GO',8,5),dpi=180)
 }
}
fwrite(rbindlist(universe_records,fill=TRUE),file.path(out,'background_and_mapping_summary.csv'))
fwrite(rbindlist(case_log,fill=TRUE),file.path(out,'case_status.csv'))
fwrite(testing_table,file.path(out,'hypothesis_counts.csv'))
prov<-list(GO_ontology='GO basic ontology from OBO Foundry',GO_ontology_file='00_cache/go-basic.obo',GO_ontology_sha256=sha256(obo),
 GO_gene_annotations=paste0('org.Hs.eg.db ',as.character(packageVersion('org.Hs.eg.db')),'; GOALL; evidence codes inclusive'),
 Entrez_mapping=paste0('org.Hs.eg.db ',as.character(packageVersion('org.Hs.eg.db')),'; ENSEMBL -> ENTREZID; multiVals=first; versions/PAR_Y suffix removed only for mapping'),
 KEGG='KEGG REST API human pathways; cached responses and retrieval times in 00_cache/kegg_provenance.json',
 input='global FDR < 0.05 and abs(log2FC) >= 1; protein_coding only; up/down tested separately',
 background='Per analysis and cell: all genes with non-NA DESeq2 pvalue and gene_type=protein_coding, mapped to ENTREZID and annotated in the respective database',
 test='One-sided hypergeometric over-representation; GOALL propagates to parent GO terms; per-set BH plus global BH across cell, direction, and ontology terms within each analysis and database',
 size_limits='Term/background annotation size 10-500; zero-overlap hypotheses included in global BH but omitted from bulky row tables',
 significance='Global BH adjusted p < 0.05',date_UTC=as.character(Sys.time()),packages=c(R=R.version.string,AnnotationDbi=as.character(packageVersion('AnnotationDbi')),org.Hs.eg.db=as.character(packageVersion('org.Hs.eg.db')),KEGGREST=as.character(packageVersion('KEGGREST'))))
jsonlite::write_json(prov,file.path(out,'enrichment_provenance.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
message('ENRICHMENT COMPLETE')

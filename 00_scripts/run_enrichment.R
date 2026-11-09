#!/usr/bin/env Rscript
# From data/: Rscript enrichment_analysis/00_scripts/run_enrichment.R
if(.Platform$OS.type=='windows')Sys.setlocale('LC_ALL','Chinese_China.utf8')
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(org.Hs.eg.db);library(AnnotationDbi);library(DBI);library(RSQLite)})
root<-'enrichment_analysis';db<-file.path(root,'01_databases')
for(d in c('02_gene_sets','03_results','04_figures','05_reports'))dir.create(file.path(root,d),showWarnings=FALSE,recursive=TRUE)
slug<-function(x)gsub('[^A-Za-z0-9]+','_',gsub('+','_plus',gsub('-','_minus',x,fixed=TRUE),fixed=TRUE))
src<-'deseq2_analysis/03_results/all_tests_global_FDR.csv.gz'
a<-as.data.table(read.csv(gzfile(src)))
stopifnot(!anyDuplicated(a[,.(analysis,cell,Geneid)]),all(c('gene_symbol','padj_global') %in% names(a)))
a[,stable_id:=sub('_PAR_Y$','',sub('\\.[0-9]+(?=_|$)','',Geneid,perl=TRUE))]
valid<-intersect(unique(a$stable_id),keys(org.Hs.eg.db,keytype='ENSEMBL'))
mp<-as.data.table(AnnotationDbi::select(org.Hs.eg.db,keys=valid,keytype='ENSEMBL',columns=c('ENTREZID','SYMBOL','GENENAME')))
mp<-unique(mp[!is.na(ENTREZID)]);mp[,n_entrez:=uniqueN(ENTREZID),by=ENSEMBL]
fwrite(mp,file.path(db,'ensembl_entrez_all_mappings.csv'),bom=TRUE)
amb<-unique(mp[n_entrez>1,.(ENSEMBL)]);fwrite(amb,file.path(db,'ambiguous_ensembl_excluded.csv'))
mp<-unique(mp[n_entrez==1,.(stable_id=ENSEMBL,ENTREZID=as.character(ENTREZID),ncbi_symbol=SYMBOL,ncbi_gene_name=GENENAME)])
stopifnot(!anyDuplicated(mp$stable_id))
a<-merge(a,mp,by='stable_id',all.x=TRUE,sort=FALSE)
fwrite(a[,.(analysis,cell,Geneid,stable_id,gene_symbol,ENTREZID,ncbi_symbol,pvalue,padj_global,log2FoldChange)],file.path(root,'02_gene_sets','mapped_gene_features.csv.gz'))
fwrite(unique(a[is.na(ENTREZID),.(Geneid,stable_id,gene_symbol)]),file.path(db,'unmapped_tested_features.csv'),bom=TRUE)
metadata_org<-AnnotationDbi::metadata(org.Hs.eg.db);fwrite(metadata_org,file.path(db,'org_Hs_eg_db_metadata.csv'))
writeLines(capture.output(sessionInfo()),file.path(root,'05_reports','sessionInfo.txt'))
# GOALL includes ancestor propagation, not merely directly assigned GO terms.
go<-as.data.table(DBI::dbGetQuery(AnnotationDbi::dbconn(org.Hs.eg.db),'SELECT DISTINCT genes.gene_id AS ENTREZID, go_all.go_id AS term_id, go_all.ontology AS ontology FROM genes JOIN go_all ON genes._id=go_all._id'))
go[,ENTREZID:=as.character(ENTREZID)];go<-go[ENTREZID %in% mp$ENTREZID]
con<-DBI::dbConnect(RSQLite::SQLite(),file.path(db,'GO.sqlite'),flags=RSQLite::SQLITE_RO)
terms<-as.data.table(DBI::dbGetQuery(con,'SELECT go_id AS term_id, term AS Description, ontology FROM go_term'))
fwrite(as.data.table(DBI::dbGetQuery(con,'SELECT * FROM metadata')),file.path(db,'GO_database_metadata.csv'));DBI::dbDisconnect(con)
go<-merge(go,terms,by=c('term_id','ontology'))
fwrite(go,file.path(db,'GOALL_human_term_gene_used.csv.gz'))
kn<-fread(file.path(db,'kegg_hsa_pathway_names.tsv'),header=FALSE,col.names=c('term_id','Description'))
kn[,term_id:=sub('^path:','',term_id)];kn[,Description:=sub(' - Homo sapiens \\(human\\)$','',Description)]
kl<-fread(file.path(db,'kegg_hsa_gene_pathway.tsv'),header=FALSE,col.names=c('kegg_gene','term_id'))
kl[,term_id:=sub('^path:','',term_id)]
kc<-fread(file.path(db,'kegg_hsa_entrez.tsv'),header=FALSE,col.names=c('kegg_gene','ENTREZID'))
kc[,ENTREZID:=sub('^ncbi-geneid:','',ENTREZID)]
kl<-merge(merge(kl,kc,by='kegg_gene'),kn,by='term_id')
kl<-unique(kl[,.(term_id,Description,ENTREZID)])
fwrite(kl,file.path(db,'KEGG_human_term_gene_used.csv.gz'))
links<-list(GO_BP=unique(go[ontology=='BP',.(term_id,Description,ENTREZID)]),GO_CC=unique(go[ontology=='CC',.(term_id,Description,ENTREZID)]),GO_MF=unique(go[ontology=='MF',.(term_id,Description,ENTREZID)]),KEGG=kl)
ans<-list();summaries<-list();mapping_summary<-list();idx<-0L
comparisons<-unique(a[,.(analysis,cell)])
for(i in seq_len(nrow(comparisons))){
 an<-comparisons$analysis[i];ce<-comparisons$cell[i];tag<-slug(ce);message(an,' | ',ce)
 x<-a[analysis==an & cell==ce & !is.na(pvalue)&!is.na(padj_global)]
 y<-x[!is.na(ENTREZID)]
 y[,up:=padj_global<.05 & log2FoldChange>=1];y[,down:=padj_global<.05 & log2FoldChange<= -1]
 collapsed<-y[,.(up=any(up),down=any(down),n_input_features=.N,gene_symbols=paste(sort(unique(gene_symbol)),collapse='/'),original_gene_ids=paste(sort(unique(Geneid)),collapse='/')),by=ENTREZID]
 collapsed[,conflicting_DE_directions:=up&down]
 mapping_summary[[i]]<-data.frame(analysis=an,cell=ce,eligible_features=nrow(x),mapped_features=nrow(y),unique_entrez=nrow(collapsed),conflicting_entrez=sum(collapsed$conflicting_DE_directions))
 setdir<-file.path(root,'02_gene_sets',an);dir.create(setdir,recursive=TRUE,showWarnings=FALSE)
 fwrite(collapsed,file.path(setdir,paste0(tag,'.csv')),bom=TRUE)
 symbols<-setNames(collapsed$gene_symbols,collapsed$ENTREZID)
 for(database in names(links)){
  lk<-links[[database]];universe<-intersect(collapsed$ENTREZID,lk$ENTREZID);N<-length(universe)
  bg<-lk[ENTREZID %in% universe]
  sizes<-bg[,.(K=uniqueN(ENTREZID)),by=.(term_id,Description)][K>=10 & K<=500]
  for(direction in c('up','down','combined')){
   selected<-if(direction=='combined')collapsed[up|down,ENTREZID] else if(direction=='up')collapsed[up&!conflicting_DE_directions,ENTREZID]else collapsed[down&!conflicting_DE_directions,ENTREZID]
   fg<-intersect(selected,universe);n<-length(fg)
   status<-if(n<5)'insufficient_foreground'else if(nrow(sizes)==0)'no_testable_terms'else 'OK'
   idx<-idx+1L
   summaries[[idx]]<-data.frame(analysis=an,cell=ce,database=database,direction=direction,universe_N=N,selected_entrez=length(selected),annotated_foreground_n=n,tested_terms=if(status=='OK')nrow(sizes)else 0L,status=status)
   if(status!='OK')next
   hits<-bg[ENTREZID %in% fg,.(k=uniqueN(ENTREZID),overlap_entrez=paste(sort(unique(ENTREZID)),collapse='/'),overlap_gene_symbols=paste(symbols[sort(unique(ENTREZID))],collapse=';')),by=term_id]
   r<-merge(sizes,hits,by='term_id',all.x=TRUE);r[is.na(k),`:=`(k=0L,overlap_entrez='',overlap_gene_symbols='')]
   # Include zero-overlap eligible terms (P=1) in every BH family.
   r[,`:=`(universe_N=N,foreground_n=n)]
   r[,pvalue:=phyper(k-1,K,N-K,n,lower.tail=FALSE)]
   r[,padj_family:=p.adjust(pvalue,'BH')]
   r[,`:=`(GeneRatio=k/n,BgRatio=K/N,FoldEnrichment=(k/n)/(K/N),analysis=an,cell=ce,database=database,direction=direction)]
   ans[[idx]]<-r
  }
 }
}
all<-rbindlist(ans,fill=TRUE)
all[,padj_analysis_database:=p.adjust(pvalue,'BH'),by=.(analysis,database)]
setcolorder(all,c('analysis','cell','database','direction','term_id','Description','k','foreground_n','K','universe_N','GeneRatio','BgRatio','FoldEnrichment','pvalue','padj_family','padj_analysis_database','overlap_entrez','overlap_gene_symbols'))
fwrite(all,file.path(root,'03_results','all_enrichment_tests.csv.gz'))
sig<-all[padj_analysis_database<.05 & k>=3]
fwrite(sig[order(padj_analysis_database)],file.path(root,'03_results','significant_terms_global_FDR05.csv'),bom=TRUE)
fwrite(rbindlist(mapping_summary),file.path(root,'02_gene_sets','ID_mapping_summary.csv'))
sm<-rbindlist(summaries)
counts<-all[,.(significant_family=sum(padj_family<.05 & k>=3),significant_global=sum(padj_analysis_database<.05 & k>=3)),by=.(analysis,cell,database,direction)]
sm<-merge(sm,counts,by=c('analysis','cell','database','direction'),all.x=TRUE)
sm[is.na(significant_family),`:=`(significant_family=0L,significant_global=0L)]
fwrite(sm,file.path(root,'03_results','enrichment_summary.csv'),bom=TRUE)
for(i in seq_len(nrow(comparisons))){
 an<-comparisons$analysis[i];ce<-comparisons$cell[i];folder<-file.path(root,'03_results',an,slug(ce));dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 for(d in names(links))fwrite(all[analysis==an&cell==ce&database==d],file.path(folder,paste0(d,'_all.csv.gz')))
 fwrite(sig[analysis==an&cell==ce][order(padj_analysis_database)],file.path(folder,'significant_terms.csv'),bom=TRUE)
 fwrite(sm[analysis==an&cell==ce],file.path(folder,'status_and_background.csv'),bom=TRUE)
}
stopifnot(nrow(sm)==720,!anyDuplicated(all[,.(analysis,cell,database,direction,term_id)]),all(all$k<=all$foreground_n),all(all$k<=all$K),all(all$foreground_n<=all$universe_N))
writeLines(capture.output(sessionInfo()),file.path(root,'05_reports','sessionInfo.txt'))
message('ENRICHMENT COMPLETE: ',nrow(all),' tests; ',nrow(sig),' globally significant term/group results')

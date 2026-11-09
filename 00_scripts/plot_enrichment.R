if(.Platform$OS.type=='windows')Sys.setlocale('LC_ALL','Chinese_China.utf8')
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
root<-'enrichment_analysis'
slug<-function(x)gsub('[^A-Za-z0-9]+','_',gsub('+','_plus',gsub('-','_minus',x,fixed=TRUE),fixed=TRUE))
s<-fread(file.path(root,'03_results','enrichment_summary.csv'))
sig<-fread(file.path(root,'03_results','significant_terms_global_FDR05.csv'))
cmp<-unique(s[,.(analysis,cell)])
chosen<-list();inventory<-list()
for(i in seq_len(nrow(cmp))){
 an<-cmp$analysis[i];ce<-cmp$cell[i];tag<-slug(ce)
 folder<-file.path(root,'04_figures',an);dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 pdf(file.path(folder,paste0(tag,'_GO_KEGG.pdf')),width=13,height=14)
 for(db in c('GO_BP','GO_CC','GO_MF','KEGG')){
  z<-copy(sig[analysis==an&cell==ce&database==db]);setorder(z,padj_analysis_database,-FoldEnrichment)
  z<-z[,head(.SD,8),by=direction]
  if(nrow(z))chosen[[paste(an,ce,db)]]<-z
  z[,display_key:=paste(direction,Description,term_id,sep='___')]
  missing<-setdiff(c('up','down','combined'),z$direction)
  empty<-data.table(direction=missing,display_key=paste0(missing,'___No significant terms'))
  if(length(missing)){
   for(j in seq_along(missing)){
    st<-s[analysis==an&cell==ce&database==db&direction==missing[j],status]
    if(length(st)&&st!='OK')empty$display_key[j]<-paste0(missing[j],'___Insufficient annotated genes / no testable terms')
   }
  }
  levels<-rev(c(z$display_key,empty$display_key))
  z[,display_key:=factor(display_key,levels=levels)];empty[,display_key:=factor(display_key,levels=levels)]
  z[,direction:=factor(direction,levels=c('up','down','combined'))];empty[,direction:=factor(direction,levels=c('up','down','combined'))]
  label<-function(x)vapply(x,function(y)paste(strwrap(gsub('___',' | ',sub('^[^_]+___','',y)),width=64),collapse='\n'),character(1))
  p<-ggplot(z,aes(GeneRatio,display_key))+geom_point(aes(size=k,color=-log10(pmax(padj_analysis_database,1e-300))),alpha=.85)+
   geom_text(data=empty,aes(x=.1,y=display_key,label='No significant terms / insufficient genes'),inherit.aes=FALSE,size=3.5,color='grey45')+
   facet_wrap(~direction,ncol=1,scales='free_y',drop=FALSE)+scale_y_discrete(labels=label)+
   scale_color_gradient(low='#316B9F',high='#C63C35')+scale_size_continuous(range=c(3,8))+
   theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),strip.background=element_rect(fill='#EDF1F5'),legend.position='bottom',plot.margin=margin(10,20,10,10))+
   labs(title=paste(ce,'-',db),subtitle=if(an=='primary_baseline_SLE_vs_HD')'Baseline / non-follow-up SLE vs HD'else'Treatment-naive SLE vs HD',
    x='Annotated foreground gene ratio (k / n)',y=NULL,size='Overlap genes',color='-log10(global FDR)',
    caption='Top 8 per direction among global BH FDR < 0.05 with >=3 overlap genes.\nGlobal correction: all cells and directions within each analysis and database. Up/down refer to SLE vs HD.')
  print(p)
  ggsave(file.path(folder,paste0(tag,'_',db,'.png')),p,width=13,height=14,dpi=140)
 }
 dev.off()
 inventory[[i]]<-data.frame(analysis=an,cell=ce,pdf=file.path('04_figures',an,paste0(tag,'_GO_KEGG.pdf')),pages=4)
}
fwrite(rbindlist(chosen,fill=TRUE),file.path(root,'03_results','plot_top8_terms.csv'),bom=TRUE)
fwrite(rbindlist(inventory),file.path(root,'05_reports','figure_inventory.csv'))
# Overview separates directions; counts are descriptive, not biological effect sizes.
s[,direction:=factor(direction,levels=c('up','down','combined'))]
p<-ggplot(s,aes(database,cell,fill=log10(significant_global+1)))+geom_tile(color='white')+
 facet_grid(analysis~direction)+scale_fill_gradient(low='white',high='#216B9B')+theme_bw(base_size=10)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())+
 labs(title='GO / KEGG enrichment overview',x=NULL,y=NULL,fill='log10(terms+1)',caption='Globally significant terms with >=3 overlap genes. Counts depend on annotation coverage and gene-list size.')
ggsave(file.path(root,'04_figures','enrichment_overview.png'),p,width=15,height=15,dpi=150)
message('ENRICHMENT FIGURES COMPLETE')

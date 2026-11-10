# Regenerate figures from saved results; no differential-expression model refitting.
if(.Platform$OS.type=='windows') Sys.setlocale('LC_ALL','Chinese_China.utf8')
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(DESeq2)})
root<-'deseq2_baseline_only_analysis'
ann<-fread(file.path(root,'01_annotation','gene_annotation.csv'))
s<-fread(file.path(root,'03_results','analysis_summary.csv'))
slug<-function(x)gsub('[^A-Za-z0-9]+','_',gsub('+','_plus',gsub('-','_minus',x,fixed=TRUE),fixed=TRUE))
pick10<-function(r){
 x<-copy(r[!is.na(padj_global)&padj_global<.05&abs(log2FoldChange)>=1 & is.finite(log2FoldChange)])
 setorderv(x,c('padj_global','log2FoldChange','Geneid'),c(1,-1,1))
 # Five per direction by smallest global FDR; replenish from the other direction.
 y<-rbind(head(x[log2FoldChange>0],5),head(x[log2FoldChange<0],5))
 if(nrow(y)<10){
  selected_ids<-y[['Geneid']]
  y<-rbind(y,head(x[!Geneid %in% selected_ids],10-nrow(y)))
 }
 setorder(y,padj_global,Geneid);y[,label_rank:=seq_len(.N)]
 y[,plot_label:=ifelse(!is.na(gene_symbol)&nzchar(gene_symbol),gene_symbol,Geneid)]
 if(nrow(y))y[grepl('_PAR_Y$',Geneid),plot_label:=paste0(plot_label,' [Y-PAR]')]
 y
}
# Deterministic label spacing, with connector lines: no optional plotting dependency.
space_labels<-function(labels,xrange,yrange){
 if(!nrow(labels))return(labels)
 labels[,side:=ifelse(log2FoldChange>0,'up','down')]
 labels[,label_x:=ifelse(log2FoldChange>0,1,-1)*xrange*1.12]
 labels[,label_y:=y]
 for(si in c('up','down')){
  ix<-which(labels$side==si);if(!length(ix))next
  ix<-ix[order(labels$y[ix])];v<-labels$y[ix]
  gap<-max(yrange*.085,.15)
  for(j in seq_along(v))if(j>1)v[j]<-max(v[j],v[j-1]+gap)
  if(max(v)>yrange*1.04)v<-v-(max(v)-yrange*1.04)
  labels$label_y[ix]<-v
 }
 labels
}
all_labels<-list();inventory<-list()
for(i in seq_len(nrow(s))){
 if(s$status[i]!='OK')next
 an<-s$analysis[i];cell<-s$cell[i];tag<-slug(cell)
 message(an,' | ',cell)
 rd<-file.path(root,'03_results',an,tag);models<-file.path(root,'05_models',an,tag)
 vd<-file.path(root,'04_figures','volcano',an);dd<-file.path(root,'04_figures','diagnostics',an)
 dir.create(vd,recursive=TRUE,showWarnings=FALSE);dir.create(dd,recursive=TRUE,showWarnings=FALSE)
 r<-fread(file.path(rd,'DESeq2_all_genes.csv'))
 r[,y:=-log10(pmax(padj_global,1e-300))]
 r[,volcano_group:=ifelse(is.na(padj_global),'not_testable',ifelse(padj_global<.05&abs(log2FoldChange)>=1,ifelse(log2FoldChange>0,'SLE higher','SLE lower'),'not_selected'))]
 labels<-pick10(r)
 saved<-copy(labels);saved[,analysis:=an];saved[,cell:=cell]
 fwrite(saved[,!c('y','volcano_group'),with=FALSE],file.path(rd,'volcano_top10_genes.csv'),bom=TRUE)
 all_labels[[paste(an,tag)]]<-saved
 xmax<-max(abs(r$log2FoldChange),na.rm=TRUE);ymax<-max(r$y[is.finite(r$y)],na.rm=TRUE);if(!is.finite(ymax)||ymax<=0)ymax<-1
 labels<-space_labels(labels,xmax,ymax)
 title_an<-if(an=='primary_baseline_SLE_vs_HD')'Baseline SLE vs HD' else 'Treatment-naive SLE vs HD'
 p<-ggplot(r[is.finite(log2FoldChange)&is.finite(y)],aes(log2FoldChange,y))+geom_point(aes(color=volcano_group),alpha=.55,size=.8)+
  scale_color_manual(values=c('SLE higher'='#CD4438','SLE lower'='#2673B8','not_selected'='#B8BDC5','not_testable'='#E1E1E1'),drop=FALSE)+
  geom_vline(xintercept=c(-1,1),linetype=2,color='grey60')+geom_hline(yintercept=-log10(.05),linetype=2,color='grey60')+
  theme_bw(base_size=12)+theme(legend.position='bottom',panel.grid.minor=element_blank())+
  labs(title=cell,subtitle=title_an,x='log2 fold change (SLE / HD)',y='-log10(global BH FDR)',color=NULL,
   caption='Labels: up to 5 higher + 5 lower genes, ranked by global FDR; |log2FC| >= 1.\nGENCODE v38 gene names; unshrunk effects. Missing FDR values are omitted.')+
  coord_cartesian(xlim=c(-xmax*1.65,xmax*1.65),ylim=c(0,ymax*1.13),clip='off')
 if(nrow(labels))p<-p+geom_segment(data=labels,aes(x=log2FoldChange,y=y,xend=label_x,yend=label_y),inherit.aes=FALSE,color='grey45',linewidth=.3)+
  geom_label(data=labels,aes(x=label_x,y=label_y,label=plot_label),inherit.aes=FALSE,size=3.2,fill='white',linewidth=.15)
 ggsave(file.path(vd,paste0(tag,'.pdf')),p,width=12,height=8)
 ggsave(file.path(vd,paste0(tag,'.png')),p,width=12,height=8,dpi=180)
 # Gene symbols remain secondary annotations; original versioned row names are intact.
 dds<-readRDS(file.path(models,'dds.rds'));ar<-as.data.frame(ann[match(rownames(dds),ann$Geneid)])
 for(nm in setdiff(names(ar),'Geneid'))rowData(dds)[[nm]]<-ar[[nm]]
 saveRDS(dds,file.path(models,'dds.rds'))
 nc<-readRDS(file.path(models,'normalized_counts.rds'));attr(nc,'gene_annotation')<-ar;saveRDS(nc,file.path(models,'normalized_counts.rds'))
 pp<-fread(file.path(root,'02_quality_control',an,tag,'PCA_scores.csv'));pp$batch<-factor(pp$batch)
 pdf(file.path(dd,paste0(tag,'.pdf')),width=12,height=8)
 plotDispEsts(dds,main=paste(cell,'- dispersion'))
 ma<-ggplot(r[baseMean>0],aes(log10(baseMean),log2FoldChange))+geom_point(aes(color=volcano_group),alpha=.5,size=.6)+
  scale_color_manual(values=c('SLE higher'='#CD4438','SLE lower'='#2673B8','not_selected'='#B8BDC5','not_testable'='#E1E1E1'))+theme_bw()+theme(legend.position='bottom')+
  labs(title=paste(cell,'- MA plot'),x='log10(mean normalized count)',y='log2 fold change',color=NULL)
 # All ten selected names are also supplied as a readable figure-side key.
 key<-if(nrow(saved))paste(sprintf('%s: log2FC %.2f',saved$plot_label,saved$log2FoldChange),collapse='\n')else'No genes pass the label criteria.'
 ma<-ma+labs(subtitle=title_an,caption=paste('Top labelled genes (GENCODE v38):\n',key))
 print(ma)
 print(ggplot(pp,aes(PC1,PC2,color=condition,shape=batch))+geom_point(size=2.5)+scale_shape_manual(values=rep(0:25,length.out=nlevels(pp$batch)))+theme_bw()+labs(title=paste(cell,'- PCA'),subtitle=title_an))
 # Include the same named volcano in the diagnostic PDF as well.
 print(p);dev.off()
 inventory[[length(inventory)+1]]<-data.frame(analysis=an,cell=cell,labels=nrow(saved),pdf=file.path('04_figures','volcano',an,paste0(tag,'.pdf')))
}
fwrite(rbindlist(all_labels,fill=TRUE),file.path(root,'03_results','volcano_top10_all_comparisons.csv'),bom=TRUE)
fwrite(rbindlist(inventory),file.path(root,'06_reports','figure_inventory.csv'))
rawfile<-file.path(root,'05_models','aligned_raw_counts.rds')
raw<-readRDS(rawfile);raw$gene_annotation<-as.data.frame(ann[match(rownames(raw$counts),ann$Geneid)]);saveRDS(raw,rawfile,compress=FALSE)
message('ANNOTATED FIGURES COMPLETE')

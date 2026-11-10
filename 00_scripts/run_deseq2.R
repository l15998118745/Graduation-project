#!/usr/bin/env Rscript
# Run from the data directory: Rscript deseq2_baseline_only_analysis/00_scripts/run_deseq2.R
options(stringsAsFactors=FALSE)
if(.Platform$OS.type=='windows') Sys.setlocale('LC_ALL','Chinese_China.utf8')
suppressPackageStartupMessages({library(DESeq2);library(data.table);library(ggplot2)})
set.seed(20261108)
out <- 'deseq2_baseline_only_analysis/.work'; dir.create(out,showWarnings=FALSE,recursive=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
rawfile <- list.files(pattern='^gencode.*csv$',full.names=TRUE)


stopifnot(length(rawfile)==1)
m <- as.data.frame(fread('metadata anno_final_aligned',na.strings=NULL))
stopifnot(!anyDuplicated(m$id))
m$gender <- tolower(m$Gender)
m$condition <- m$BinaryGroup
m$batch <- ifelse(grepl('^20[0-9]{6}',m$id),substr(m$id,1,6),'unknown')
m$sample_conflict <- ave(m$patient,m$sample,FUN=function(x) length(unique(x)))>1
m$patient_cross_group <- ave(m$condition,m$patient,FUN=function(x) length(unique(x)))>1
m$identity_ok <- !m$sample_conflict & !m$patient_cross_group & !m$patient %in% c('','#N/A','-')
fwrite(m[!m$identity_ok,],file.path(out,'identity_exclusions.csv'))
message('Reading raw integer counts')




z <- fread(rawfile,check.names=FALSE)
genes <- z[[1]]; stopifnot(!anyDuplicated(genes),!anyDuplicated(names(z)[-1]))
rawids <- names(z)[-1]
qc <- data.frame(id=rawids,total_counts=vapply(z[,-1],sum,numeric(1)),detected_genes=vapply(z[,-1],function(x)sum(x>0),integer(1)))
# Conservative gross-failure screens, not universal sequencing quality thresholds.
qc$qc_pass <- qc$total_counts>=10000 & qc$detected_genes>=1000
qc$has_metadata <- qc$id %in% m$id
fwrite(qc,file.path(out,'library_qc_all.csv'))
fwrite(data.frame(id=setdiff(rawids,m$id)),file.path(out,'raw_without_metadata.csv'))
stopifnot(all(m$id %in% rawids))
m <- merge(m,qc,by='id',sort=FALSE)
fwrite(m,file.path(out,'metadata_audited.csv'))
mat <- as.matrix(z[,m$id,with=FALSE]); rm(z);gc()
stopifnot(!anyNA(mat),all(mat>=0),all(mat==floor(mat)),max(mat)<.Machine$integer.max)
storage.mode(mat)<-'integer';rownames(mat)<-genes;colnames(mat)<-m$id
saveRDS(list(counts=mat,coldata=m),file.path(out,'aligned_raw_counts.rds'),compress=FALSE)
p<-ggplot(m,aes(total_counts,detected_genes,color=condition))+geom_point(alpha=.5)+scale_x_log10()+theme_bw()
ggsave(file.path(out,'library_QC.png'),p,width=8,height=5,dpi=160)
summary_rows<-list(); all_p<-list(); skipped<-list()





run_one<-function(md,cell,analysis){
  tag<-gsub('[^A-Za-z0-9]+','_',gsub('+','_plus',gsub('-','_minus',cell,fixed=TRUE),fixed=TRUE)); folder<-file.path(out,analysis,tag);dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  md<-md[md$cell==cell & md$identity_ok & md$qc_pass,]
  # One independent observation per patient and cell. Prefer earlier dated sample;
  # missing dates sort last. Never sum biological visits or distinct specimens.
  md$date_key<-ifelse(grepl('^20[0-9]{6}',md$id),substr(md$id,1,8),'99999999')
  md<-md[order(md$date_key,md$id),]; discarded<-md[duplicated(md$patient),]
  fwrite(discarded,file.path(folder,'repeated_patient_exclusions.csv'))
  md<-md[!duplicated(md$patient),]
  md$condition<-factor(md$condition,levels=c('HD','SLE'));md$gender<-factor(md$gender);md$batch<-factor(md$batch)
  n<-table(md$condition)
  if(any(n<3)){return(data.frame(analysis=analysis,cell=cell,status='fewer_than_3_per_group',HD=n[1],SLE=n[2]))}
  # Candidate technical batch is the YYYYMM prefix, NOT a validated run identifier.
  terms<-c(if(nlevels(md$gender)>1)'gender',if(nlevels(md$batch)>1)'batch','condition')
  form<-reformulate(terms); mm<-model.matrix(form,md)
  if(qr(mm)$rank<ncol(mm) || nrow(mm)-ncol(mm)<3){
    fwrite(md,file.path(folder,'samples.csv'))
    return(data.frame(analysis=analysis,cell=cell,status='confounded_or_insufficient_residual_df',HD=n[1],SLE=n[2]))
  }
  fwrite(as.data.frame(table(md$batch,md$condition)),file.path(folder,'batch_by_condition.csv'))
  cnt<-mat[,md$id,drop=FALSE]; min_n<-min(n)
  keep<-rowSums(cnt>=10)>=min_n
  fwrite(data.frame(Geneid=genes,kept=keep),file.path(folder,'gene_filter.csv'))
  cnt<-cnt[keep,,drop=FALSE];rownames(md)<-md$id
  dds<-DESeqDataSetFromMatrix(cnt,md,form)




  # poscounts accommodates sparse counts; input remains unnormalized integer counts.
  cache<-file.path(folder,'dds.rds')
  if(Sys.getenv('DESEQ2_RESUME')=='1' && file.exists(cache)){
    fitted<-readRDS(cache)
    stopifnot(identical(counts(fitted),cnt),identical(deparse(design(fitted)),deparse(form)))
    # Verify every model covariate as well as the complete count matrix.
    for(nm in all.vars(form)) stopifnot(identical(as.character(colData(fitted)[[nm]]),as.character(md[[nm]])))
    dds<-fitted
  } else dds<-DESeq(dds,sfType='poscounts',quiet=TRUE,parallel=FALSE)
  res<-results(dds,contrast=c('condition','SLE','HD'),alpha=.05)
  r<-data.frame(Geneid=rownames(res),as.data.frame(res))
  r$ensembl_unversioned<-sub('\\.[0-9]+(?=_|$)','',r$Geneid,perl=TRUE)
  r$padj_global<-NA_real_
  r$direction<-ifelse(!is.na(r$padj)&r$padj<.05 & abs(r$log2FoldChange)>=1,ifelse(r$log2FoldChange>0,'SLE_up','SLE_down'),'not_selected')
  fwrite(r[order(r$padj),],file.path(folder,'DESeq2_all_genes.csv'))
  fwrite(r[r$direction!='not_selected',],file.path(folder,'DEG_FDR05_absLFC1.csv'))
  md$size_factor<-sizeFactors(dds);fwrite(md,file.path(folder,'samples.csv'))
  saveRDS(dds,file.path(folder,'dds.rds'))
  saveRDS(counts(dds,normalized=TRUE),file.path(folder,'normalized_counts.rds'))
  pdf(file.path(folder,'diagnostics.pdf'),width=8,height=6)
  plotDispEsts(dds);plotMA(res,main=paste(analysis,cell),ylim=c(-6,6))




  # log2 normalized counts for visualization only, never used for hypothesis tests.
  logc<-log2(counts(dds,normalized=TRUE)+1)
  vv<-apply(logc,1,var); top<-head(order(vv,decreasing=TRUE),500)
  pc<-prcomp(t(logc[top,,drop=FALSE])); pv<-100*pc$sdev^2/sum(pc$sdev^2)
  pp<-data.frame(id=md$id,PC1=pc$x[,1],PC2=pc$x[,2],condition=md$condition,batch=md$batch)
  print(ggplot(pp,aes(PC1,PC2,color=condition,shape=batch))+geom_point(size=2)+scale_shape_manual(values=rep(0:25,length.out=nlevels(md$batch)))+theme_bw()+labs(title=cell,x=sprintf('PC1 %.1f%%',pv[1]),y=sprintf('PC2 %.1f%%',pv[2])))
  print(ggplot(r,aes(log2FoldChange,-log10(pmax(padj,1e-300)),color=direction))+geom_point(alpha=.4,size=.5,na.rm=TRUE)+theme_bw()+labs(title=paste(analysis,cell)))
  dev.off();fwrite(pp,file.path(folder,'PCA_scores.csv'))
  all_p[[paste(analysis,tag,sep='/')]]<<-data.frame(analysis=analysis,cell=cell,Geneid=r$Geneid,pvalue=r$pvalue,log2FoldChange=r$log2FoldChange,padj_within_cell=r$padj)
  data.frame(analysis=analysis,cell=cell,status='OK',HD=unname(n[1]),SLE=unname(n[2]),genes_retained=sum(keep),tested_pvalues=sum(!is.na(r$pvalue)),FDR05=sum(r$padj<.05,na.rm=TRUE),SLE_up=sum(r$direction=='SLE_up'),SLE_down=sum(r$direction=='SLE_down'),design=paste(deparse(form),collapse=''),note=ifelse(cell %in% c('EV','Plasma'),'noncell_fraction_exploratory','cell_type'))
}
for(analysis in 'primary_baseline_SLE_vs_HD'){
  # Primary analysis: all healthy controls versus baseline SLE only.
  # Timepoint is not applicable to HD (recorded as "-"); SLE visits with an
  # unknown or non-baseline timepoint are deliberately excluded.
  eligible<-m$condition=='HD' | (m$condition=='SLE' & m$Timepoint=='baseline')
  md<-m[eligible,]
  for(cell in sort(unique(md$cell))){
    message(analysis,' | ',cell)
    ans<-tryCatch(run_one(md,cell,analysis),error=function(e)data.frame(analysis=analysis,cell=cell,status=paste('ERROR',conditionMessage(e))))
    summary_rows[[length(summary_rows)+1]]<-ans
    fwrite(rbindlist(summary_rows,fill=TRUE),file.path(out,'analysis_summary.csv'))
  }
}
if(any(grepl('^ERROR',rbindlist(summary_rows,fill=TRUE)$status))) stop('One or more comparisons failed; inspect analysis_summary.csv')
if(length(all_p)){
  a<-rbindlist(all_p)
  # Additional family-wide BH across all genes/cell types WITHIN each analysis.
  a[,padj_global:=p.adjust(pvalue,method='BH'),by=analysis]
  fwrite(a,file.path(out,'all_tests_global_FDR.csv.gz'))
  fwrite(a[!is.na(padj_global)&padj_global<.05&abs(log2FoldChange)>=1],file.path(out,'DEG_global_FDR05_absLFC1.csv'))
}
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
message('COMPLETE')

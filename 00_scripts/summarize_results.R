# Run AFTER run_deseq2.R, from the data directory.
if(.Platform$OS.type=='windows') Sys.setlocale('LC_ALL','Chinese_China.utf8')
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
out<-'deseq2_baseline_only_analysis/.work'
m<-fread(file.path(out,'metadata_audited.csv'))
b<-m[identity_ok==TRUE & qc_pass==TRUE]
pairs<-b[condition=='SLE' & Treatment!='-',.(n_pairs=sum(vapply(split(Timepoint,patient),function(t) all(c('baseline','6M') %in% t),logical(1)))),by=.(cell,Treatment)]
fwrite(pairs,file.path(out,'paired_analysis_feasibility.csv'))
s<-fread(file.path(out,'analysis_summary.csv'))
a<-as.data.table(read.csv(gzfile(file.path(out,'all_tests_global_FDR.csv.gz'))))
g<-a[!is.na(padj_global)&padj_global<.05&abs(log2FoldChange)>=1]
global<-g[,.(global_DEG=.N,global_up=sum(log2FoldChange>0),global_down=sum(log2FoldChange<0)),by=.(analysis,cell)]
s<-merge(s,global,by=c('analysis','cell'),all.x=TRUE)
for(c in c('global_DEG','global_up','global_down'))set(s,which(is.na(s[[c]])&s$status=='OK'),c,0L)
fwrite(s,file.path(out,'analysis_summary_with_global_FDR.csv'))
ggsave(file.path(out,'DEG_counts_global_FDR.png'),ggplot(s[status=='OK'],aes(global_DEG,reorder(cell,global_DEG),fill=analysis))+geom_col(position='dodge')+theme_bw()+labs(x='DEGs: global BH FDR < 0.05 and |unshrunk log2FC| >= 1',y=NULL),width=11,height=10,dpi=160)
recurrent<-g[,.(n_cell_types=uniqueN(cell),n_up=sum(log2FoldChange>0),n_down=sum(log2FoldChange<0)),by=.(analysis,Geneid)][order(analysis,-n_cell_types)]
fwrite(recurrent,file.path(out,'recurrent_genes_across_cells.csv'))
checks<-data.frame(check=c('no_model_errors','30_primary_comparisons_accounted_for','all_fitted_have_min3_each_group','all_fitted_have_global_counts'),passed=c(!any(grepl('^ERROR',s$status)),nrow(s)==30 && uniqueN(s$analysis)==1 && unique(s$analysis)=='primary_baseline_SLE_vs_HD',all(s[status=='OK',HD>=3&SLE>=3]),!anyNA(s[status=='OK',global_DEG])))
fwrite(checks,file.path(out,'validation_checks.csv'))
print(s[,.(analysis,cell,status,HD,SLE,global_DEG)])
print(checks)
# Populate the global-FDR column in each comparison's own result files as well.
for(i in seq_len(nrow(s))){
 if(s$status[i]!='OK') next
 an<-s$analysis[i];ce<-s$cell[i]
 folder<-file.path(out,an,gsub('[^A-Za-z0-9]+','_',gsub('+','_plus',gsub('-','_minus',ce,fixed=TRUE),fixed=TRUE)))
 rr<-fread(file.path(folder,'DESeq2_all_genes.csv'))
 vals<-a[analysis==an & cell==ce]
 stopifnot(all(rr$Geneid %in% vals$Geneid))
 rr$padj_global<-vals$padj_global[match(rr$Geneid,vals$Geneid)]
 fwrite(rr,file.path(folder,'DESeq2_all_genes.csv'))
 fwrite(rr[direction!='not_selected'],file.path(folder,'DEG_FDR05_absLFC1.csv'))
}

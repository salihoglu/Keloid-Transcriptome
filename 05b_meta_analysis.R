source('scripts/common.R')
dir_create('results/meta'); dir_create('reports')
manifest <- get_manifest()
merged <- read_tsv_safe('results/deg/multi_dataset_gene_level_merged_all_available.tsv')
if (!nrow(merged)) merged <- read_tsv_safe('results/deg/multi_dataset_gene_level_merged.tsv')
if (!nrow(merged)) { write_tsv(data.frame(), 'results/meta/keloid_meta_analysis_all_genes.tsv'); writeLines('# Meta-analysis\n\nNo input DEG results.', 'reports/meta_analysis_report.md'); quit(save='no') }
if (!('accession' %in% names(merged))) merged$accession <- sub('_.*$','', merged$dataset)
merged$include_in_meta_analysis <- manifest_flag(merged$accession, 'include_in_meta_analysis', default=TRUE, manifest=manifest)
meta_input <- merged[merged$include_in_meta_analysis,]
write_tsv(unique(meta_input[,intersect(c('dataset','accession','intended_role','intended_role_manifest','include_in_meta_analysis'), names(meta_input)), drop=FALSE]), 'results/meta/meta_analysis_dataset_inclusion.tsv')
if (!nrow(meta_input)) { write_tsv(data.frame(), 'results/meta/keloid_meta_analysis_all_genes.tsv'); writeLines('# Meta-analysis\n\nNo manifest-eligible meta-analysis input DEG results.', 'reports/meta_analysis_report.md'); quit(save='no') }
meta <- do.call(rbind, lapply(split(meta_input, meta_input$gene_symbol), function(x) {
  p <- pmax(pmin(as.numeric(x$P.Value),1), .Machine$double.xmin)
  lfc <- as.numeric(x$logFC)
  z <- qnorm(p/2, lower.tail=FALSE) * sign(lfc)
  k <- sum(is.finite(z))
  zc <- sum(z, na.rm=TRUE)/sqrt(max(1,k))

  q <- if (k > 1) sum((lfc - mean(lfc, na.rm=TRUE))^2, na.rm=TRUE) else NA_real_
  i2 <- if (k > 1 && is.finite(q) && q > 0) max(0, (q - (k-1))/q) else NA_real_
  data.frame(gene_symbol=x$gene_symbol[1], n_datasets=k, mean_logFC=mean(lfc, na.rm=TRUE), median_logFC=median(lfc, na.rm=TRUE), stouffer_z=zc, meta_p=2*pnorm(abs(zc), lower.tail=FALSE), n_up=sum(lfc>0, na.rm=TRUE), n_down=sum(lfc<0, na.rm=TRUE), direction=ifelse(mean(lfc,na.rm=TRUE)>0,'up','down'), direction_consistency=pmax(sum(lfc>0,na.rm=TRUE),sum(lfc<0,na.rm=TRUE))/max(1,k), descriptive_I2=i2)
}))
meta$meta_fdr <- p.adjust(meta$meta_p, 'BH')
meta <- meta[order(meta$meta_fdr, -meta$direction_consistency, -meta$n_datasets),]
write_tsv(meta, 'results/meta/keloid_meta_analysis_all_genes.tsv')
sig <- meta[meta$meta_fdr<0.05 & abs(meta$mean_logFC)>=0.5 & meta$n_datasets>=2 & meta$direction_consistency>=0.67,]
write_tsv(sig, 'results/meta/keloid_meta_analysis_consensus_signature.tsv')
writeLines(c('# Meta-analysis report','',paste('Manifest-eligible DEG datasets:', 

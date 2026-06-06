source('scripts/common.R')
dir_create('results/deg'); dir_create('reports')
asg <- read_tsv_safe('data/metadata/sample_group_assignment.tsv')
manifest <- get_manifest()
files <- list.files('data/processed', pattern='_expression_preprocessed.tsv.gz$', full.names=TRUE)
files <- files[!grepl('_SUPP', files)]
run_deg <- function(f, asg, manifest) {
  dataset <- sub('_expression_preprocessed.tsv.gz$','',basename(f)); acc <- strsplit(dataset,'_')[[1]][1]
  role <- manifest_role(acc, manifest)[1]
  include_primary <- manifest_flag(acc, 'include_in_primary_deg', default=TRUE, manifest=manifest)[1]
  include_meta <- manifest_flag(acc, 'include_in_meta_analysis', default=TRUE, manifest=manifest)[1]
  include_context <- manifest_flag(acc, 'include_in_contextual_validation', default=FALSE, manifest=manifest)[1]
  analyze_all_literature <- Sys.getenv('KELOID_ANALYZE_ALL_LITERATURE_GEO', unset='1') %in% c('1','true','TRUE','yes','YES')
  active_bulk_route <- include_primary || include_meta || (analyze_all_literature && include_context)
  if (!active_bulk_route) {
    return(list(summary=data.frame(dataset=dataset, accession=acc, intended_role=role, status='skipped', skip_reason='manifest_excludes_all_bulk_deg_routes', n_keloid=NA, n_control=NA, n_significant=NA), deg=data.frame()))
  }
  exprdf <- read_tsv_safe(f); if (!nrow(exprdf)) stop('empty expression')
  genes <- exprdf$gene_symbol; expr <- as.matrix(exprdf[, setdiff(names(exprdf),'gene_symbol'), drop=FALSE]); storage.mode(expr) <- 'numeric'; rownames(expr) <- genes


  if (acc == 'GSE90051') {
    y <- expr
    n_pairs <- ncol(y)
    if (n_pairs < 3) {
      return(list(summary=data.frame(dataset=dataset, accession=acc, intended_role=role, status='skipped', skip_reason='insufficient paired two-channel arrays', n_keloid=0, n_control=0, n_significant=NA), deg=data.frame()))
    }
    if (requireNamespace('limma', quietly=TRUE)) {
      design <- matrix(1, nrow=n_pairs, ncol=1); colnames(design) <- 'Intercept'
      fit <- limma::lmFit(y, design); fit <- limma::eBayes(fit)
      tab <- limma::topTable(fit, coef='Intercept', number=Inf, sort.by='P')
      deg <- data.frame(gene_symbol=rownames(tab), gene=rownames(tab), logFC=tab$logFC, AveExpr=tab$AveExpr, t=tab$t, P.Value=tab$P.Value, adj.P.Val=tab$adj.P.Val, B=tab$B, stringsAsFactors=FALSE)
    } else {
      p <- apply(y,1,function(v) tryCatch(t.test(v, mu=0)$p.value, error=function(e) NA_real_))
      lfc <- rowMeans(y, na.rm=TRUE)
      deg <- data.frame(gene_symbol=rownames(y), gene=rownames(y), logFC=lfc, P.Value=p, adj.P.Val=p.adjust(p,'BH'))
    }
    deg$dataset <- dataset; deg$accession <- acc; deg$intended_role <- role; deg$direction <- ifelse(deg$logFC>0,'up','down')
    sig <- deg[!is.na(deg$adj.P.Val) & deg$adj.P.Val < 0.05 & abs(deg$logFC) >= 1,]
    write_tsv(data.frame(sample_col=colnames(y), sample_id=colnames(y), assigned_group='paired_keloid_control_logratio', exclusion_reason='two_channel_ratio_model', sample_type='bulk_tissue'), file.path('results/deg', paste0(dataset, '_sample_mapping_used.tsv')))
    write_tsv(deg, file.path('results/deg', paste0(dataset, '_all_genes_DE_results.tsv')))
    write_tsv(sig, file.path('results/deg', paste0(dataset, '_significant_DEGs_FDR005_logFC1.tsv')))
    writeLines(sig$gene_symbol[sig$logFC>0], file.path('results/deg', paste0(dataset, '_upregulated_genes.txt')))
    writeLines(sig$gene_symbol[sig$logFC<0], file.path('results/deg', paste0(dataset, '_downregulated_genes.txt')))
    write_tsv(data.frame(dataset=dataset, accession=acc, intended_role=role, formula='~ 1; one-sample test of paired log-ratio vs 0', n_keloid=n_pairs, n_control=n_pairs, low_power=(n_pairs<4)), file.path('results/deg', paste0(dataset, '_model_design_report.tsv')))
    return(list(summary=data.frame(dataset=dataset, accession=acc, intended_role=role, status='complete', skip_reason='', n_keloid=n_pairs, n_control=n_pairs, n_significant=nrow(sig), sample_use_policy='paired_logratio_all_pairs', effective_balanced_n_per_group=n_pairs, cohort_balance_ratio=1), deg=deg))
  }

  map <- resolve_sample_groups(colnames(expr), acc, asg)
  map <- map[map$assigned_group %in% c('keloid','control') & !(tolower(map$exclusion_reason) %in% c('treated_or_perturbed')),]
  if (!('sample_type' %in% names(map))) map$sample_type <- NA_character_
  map$sample_type[is.na(map$sample_type) | map$sample_type==''] <- mapply(infer_sample_type, map$sample_col[is.na(map$sample_type) | map$sample_type==''], USE.NAMES=FALSE)
  if (!nrow(map) || sum(map$assigned_group=='keloid') < 2 || sum(map$assigned_group=='control') < 2) {
    write_tsv(map, file.path('results/deg', paste0(dataset, '_sample_mapping_failed_or_insufficient.tsv')))
    return(list(summary=data.frame(dataset=dataset, accession=acc, intended_role=role, status='skipped', skip_reason='insufficient curated keloid/control groups after robust sample matching', n_keloid=sum(map$assigned_group=='keloid', na.rm=TRUE), n_control=sum(map$assigned_group=='control', na.rm=TRUE), n_significant=NA), deg=data.frame()))
  }
  map_unbalanced <- map
  map <- balance_case_control_map(map, dataset=dataset, out_dir='results/deg')

  write_tsv(map_unbalanced, file.path('results/deg', paste0(dataset, '_sample_mapping_before_balancing.tsv')))
  write_tsv(map, file.path('results/deg', paste0(dataset, '_sample_mapping_used.tsv')))
  y <- expr[, map$sample_col, drop=FALSE]
  group <- factor(map$assigned_group, levels=c('control','keloid'))
  model_terms <- '~ group'
  if (requireNamespace('limma', quietly=TRUE)) {
    st <- factor(map$sample_type)
    use_st <- length(unique(st)) > 1 && all(table(group, st) > 0)
    if (use_st) { design <- model.matrix(~ st + group); coef_name <- 'groupkeloid'; model_terms <- '~ sample_type + group' } else { design <- model.matrix(~ group); coef_name <- 'groupkeloid' }
    fit <- limma::lmFit(y, design); fit <- limma::eBayes(fit)
    tab <- limma::topTable(fit, coef=coef_name, number=Inf, sort.by='P')
    deg <- data.frame(gene_symbol=rownames(tab), gene=rownames(tab), logFC=tab$logFC, AveExpr=tab$AveExpr, t=tab$t, P.Value=tab$P.Value, adj.P.Val=tab$adj.P.Val, B=tab$B, stringsAsFactors=FALSE)
  } else {
    p <- apply(y,1,function(v) tryCatch(t.test(v[group=='keloid'], v[group=='control'])$p.value, error=function(e) NA_real_))
    lfc <- rowMeans(y[,group=='keloid',drop=FALSE], na.rm=TRUE)-rowMeans(y[,group=='control',drop=FALSE], na.rm=TRUE)
    deg <- data.frame(gene_symbol=rownames(y), gene=rownames(y), logFC=lfc, P.Value=p, adj.P.Val=p.adjust(p,'BH'))
  }
  deg$dataset <- dataset; deg$accession <- acc; deg$intended_role <- role; deg$direction <- ifelse(deg$logFC>0,'up','down')
  sig <- deg[!is.na(deg$adj.P.Val) & deg$adj.P.Val < 0.05 & abs(deg$logFC) >= 1,]
  write_tsv(deg, file.path('results/deg', paste0(dataset, '_all_genes_DE_results.tsv')))
  write_tsv(sig, file.path('results/deg', paste0(dataset, '_significant_DEGs_FDR005_logFC1.tsv')))
  writeLines(sig$gene_symbol[sig$logFC>0], file.path('results/deg', paste0(dataset, '_upregulated_genes.txt')))
  writeLines(sig$gene_symbol[sig$logFC<0], file.path('results/deg', paste0(dataset, '_downregulated_genes.txt')))
  write_tsv(data.frame(dataset=dataset, accession=acc, intended_role=role, formula=model_terms, n_keloid=sum(group=='keloid'), n_control=sum(group=='control'), balanced_equal_n=(sum(group=='keloid')==sum(group=='control')), sample_use_policy=unique(map$sample_use_policy)[1], effective_balanced_n_per_group=min(sum(group=='keloid'), sum(group=='control')), cohort_balance_ratio=round(min(sum(group=='keloid'), sum(group=='control'))/max(sum(group=='keloid'), sum(group=='control')),4), low_power=(sum(group=='keloid')<4 | sum(group=='control')<4)), file.path('results/deg', paste0(dataset, '_model_design_report.tsv')))
  list(summary=data.frame(dataset=dataset, accession=acc, intended_role=role, status='complete', skip_reason='', n_keloid=sum(group=='keloid'), n_control=sum(group=='control'), n_significant=nrow(sig), sample_use_policy=unique(map$sample_use_policy)[1], effective_balanced_n_per_group=min(sum(group=='keloid'), sum(group=='control')), cohort_balance_ratio=round(min(sum(group=='keloid'), sum(group=='control'))/max(sum(group=='keloid'), sum(group=='control')),4)), deg=deg)
}
summ <- list()
for (f in files) {
  res <- tryCatch(run_deg(f, asg, manifest), error=function(e){ dataset<-sub('_expression_preprocessed.tsv.gz$','',basename(f)); acc<-strsplit(dataset,'_')[[1]][1]; log_msg(dataset, 'DEG error:', conditionMessage(e)); list(summary=data.frame(dataset=dataset, accession=acc, intended_role=manifest_role(acc, manifest)[1], status='failed', skip_reason=conditionMessage(e), n_keloid=NA, n_control=NA, n_significant=NA), deg=data.frame()) })
  summ[[length(summ)+1]] <- res$summary
  log_msg(res$summary$dataset, ':', ifelse(res$summary$status=='complete', paste('DEG complete; significant=', res$summary$n_significant), paste('skipped/failed', res$summary$skip_reason)))
}
write_tsv(union_rbind(summ), 'results/deg/DEG_summary.tsv')

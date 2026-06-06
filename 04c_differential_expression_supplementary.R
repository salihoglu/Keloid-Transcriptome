source('scripts/common.R')
dir_create('results/deg')
asg <- read_tsv_safe('data/metadata/sample_group_assignment.tsv')
manifest <- get_manifest()
files <- list.files('data/processed', pattern='_SUPP_.*_expression_preprocessed.tsv.gz$', full.names=TRUE)
map_by_colname <- function(cols, acc) {

  out <- data.frame(sample_col=cols, sample_id=cols, assigned_group=NA_character_, exclusion_reason='', sample_type=NA_character_, stringsAsFactors=FALSE)
  for (i in seq_along(cols)) {
    cur0 <- curate_processed_column_group(acc, cols[i])
    if (!is.na(cur0[['group']]) && nzchar(cur0[['group']])) {
      out$assigned_group[i] <- cur0[['group']]; out$exclusion_reason[i] <- cur0[['reason']]; out$sample_type[i] <- cur0[['sample_type']]; next
    }
    cur <- curate_known_accession_sample(acc, title=cols[i], source=cols[i], characteristics='', sample_id=cols[i])
    if (!is.na(cur[['group']]) && nzchar(cur[['group']])) {
      out$assigned_group[i] <- cur[['group']]; out$exclusion_reason[i] <- cur[['reason']]; out$sample_type[i] <- cur[['sample_type']]
    }
  }
  out
}

normalize_supp_mapping <- function(map, cols, acc, asg) {
  if (!is.data.frame(map)) map <- data.frame()
  for (cc in c('sample_col','sample_id','assigned_group','exclusion_reason','sample_type')) {
    if (!(cc %in% names(map))) map[[cc]] <- NA_character_
  }
  map$sample_col <- as.character(map$sample_col)
  map$sample_id <- as.character(map$sample_id)
  map$assigned_group <- as.character(map$assigned_group)
  map$exclusion_reason <- as.character(map$exclusion_reason)
  map$sample_type <- as.character(map$sample_type)

  if (!length(map$sample_col) || any(!cols %in% map$sample_col)) {
    restored <- data.frame(sample_col=cols, sample_id=cols, assigned_group=NA_character_, exclusion_reason='', sample_type=NA_character_, stringsAsFactors=FALSE)
    if (nrow(map)) {
      hit <- match(restored$sample_col, map$sample_col)
      ok <- !is.na(hit)
      restored$sample_id[ok] <- map$sample_id[hit[ok]]
      restored$assigned_group[ok] <- map$assigned_group[hit[ok]]
      restored$exclusion_reason[ok] <- map$exclusion_reason[hit[ok]]
      restored$sample_type[ok] <- map$sample_type[hit[ok]]
    }
    map <- restored
  }
  map
}

run_supp_deg <- function(f) {
  dataset <- sub('_expression_preprocessed.tsv.gz$','',basename(f)); acc <- strsplit(dataset,'_')[[1]][1]
  role <- manifest_role(acc, manifest)[1]
  include_primary <- manifest_flag(acc, 'include_in_primary_deg', default=TRUE, manifest=manifest)[1]
  include_meta <- manifest_flag(acc, 'include_in_meta_analysis', default=TRUE, manifest=manifest)[1]
  include_context <- manifest_flag(acc, 'include_in_contextual_validation', default=FALSE, manifest=manifest)[1]
  analyze_all_literature <- Sys.getenv('KELOID_ANALYZE_ALL_LITERATURE_GEO', unset='1') %in% c('1','true','TRUE','yes','YES')
  active_bulk_route <- include_primary || include_meta || (analyze_all_literature && include_context)
  if (!active_bulk_route) {
    return(data.frame(dataset=dataset, accession=acc, intended_role=role, status='skipped', skip_reason='manifest_excludes_all_bulk_deg_routes', n_keloid=NA, n_control=NA, n_significant=NA, formula=''))
  }
  exprdf <- read_tsv_safe(f); if (!nrow(exprdf)) stop('empty expression')
  genes <- exprdf$gene_symbol; expr <- as.matrix(exprdf[, setdiff(names(exprdf),'gene_symbol'), drop=FALSE]); storage.mode(expr) <- 'numeric'; rownames(expr) <- genes

  map <- normalize_supp_mapping(map_by_colname(colnames(expr), acc), colnames(expr), acc, asg)
  if (sum(map$assigned_group=='keloid', na.rm=TRUE) < 2 || sum(map$assigned_group=='control', na.rm=TRUE) < 2) {
    m2 <- resolve_sample_groups(colnames(expr), acc, asg)
    if (sum(m2$assigned_group=='keloid', na.rm=TRUE) >= 2 && sum(m2$assigned_group=='control', na.rm=TRUE) >= 2) map <- normalize_supp_mapping(m2, colnames(expr), acc, asg)
  }
  if (sum(map$assigned_group=='keloid', na.rm=TRUE) < 2 || sum(map$assigned_group=='control', na.rm=TRUE) < 2) {
    mf <- force_manifest_sample_groups(colnames(expr), acc, asg, manifest)
    if (sum(mf$assigned_group=='keloid', na.rm=TRUE) >= 2 && sum(mf$assigned_group=='control', na.rm=TRUE) >= 2) map <- normalize_supp_mapping(mf, colnames(expr), acc, asg)
  }

  bad_reason <- grepl('treated|perturbed|hypertrophic|nonlesional|perilesional|contextual_not_healthy', tolower(map$exclusion_reason))
  map <- map[map$assigned_group %in% c('keloid','control') & !bad_reason, , drop=FALSE]
  if (!nrow(map) || sum(map$assigned_group=='keloid') < 2 || sum(map$assigned_group=='control') < 2) {
    write_tsv(map, file.path('results/deg', paste0(dataset, '_supplementary_sample_mapping_failed_or_insufficient.tsv')))
    return(data.frame(dataset=dataset, accession=acc, intended_role=role, status='skipped', skip_reason='insufficient curated keloid/control groups after supplementary sample matching', n_keloid=sum(map$assigned_group=='keloid', na.rm=TRUE), n_control=sum(map$assigned_group=='control', na.rm=TRUE), n_significant=NA, formula=''))
  }
  map_unbalanced <- map
  map <- balance_case_control_map(map, dataset=dataset, out_dir='results/deg')

  write_tsv(map_unbalanced, file.path('results/deg', paste0(dataset, '_supplementary_sample_mapping_before_balancing.tsv')))
  write_tsv(map, file.path('results/deg', paste0(dataset, '_supplementary_sample_mapping.tsv')))
  y <- expr[, map$sample_col, drop=FALSE]
  group <- factor(map$assigned_group, levels=c('control','keloid'))
  formula_used <- '~ group'
  if (requireNamespace('limma', quietly=TRUE)) {
    design <- model.matrix(~ group); fit <- limma::lmFit(y, design); fit <- limma::eBayes(fit)
    tab <- limma::topTable(fit, coef='groupkeloid', number=Inf, sort.by='P')
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
  write_tsv(data.frame(dataset=dataset, accession=acc, intended_role=role, formula=formula_used, n_keloid=sum(group=='keloid'), n_control=sum(group=='control'), balanced_equal_n=(sum(group=='keloid')==sum(group=='control')), sample_use_policy=unique(map$sample_use_policy)[1], effective_balanced_n_per_group=min(sum(group=='keloid'), sum(group=='control')), cohort_balance_ratio=round(min(sum(group=='keloid'), sum(group=='control'))/max(sum(group=='keloid'), sum(group=='control')),4), low_power=(sum(group=='keloid')<4 | sum(group=='control')<4)), file.path('results/deg', paste0(dataset, '_model_design_report.tsv')))
  data.frame(dataset=dataset, accession=acc, intended_role=role, status='complete', skip_reason='', n_keloid=sum(group=='keloid'), n_control=sum(group=='control'), n_significant=nrow(sig), formula=formula_used, sample_use_policy=unique(map$sample_use_policy)[1], effective_balanced_n_per_group=min(sum(group=='keloid'), sum(group=='control')), cohort_balance_ratio=round(min(sum(group=='keloid'), sum(group=='control'))/max(sum(group=='keloid'), sum(group=='control')),4))
}
summ <- list()
for (f in files) {
  res <- tryCatch(run_supp_deg(f), error=function(e){ dataset<-sub('_expression_preprocessed.tsv.gz$','',basename(f)); acc<-strsplit(dataset,'_')[[1]][1]; log_msg(dataset, 'supplementary DEG error:', conditionMessage(e)); data.frame(dataset=dataset, accession=acc, intended_role=manifest_role(acc, manifest)[1], status='failed', skip_reason=conditionMessage(e), n_keloid=NA, n_control=NA, n_significant=NA, formula='') })
  summ[[length(summ)+1]] <- res
}
write_tsv(union_rbind(summ), 'results/deg/DEG_summary_supplementary.tsv')
old <- read_tsv_safe('results/deg/DEG_summary.tsv')
all <- union_rbind(c(list(old), summ))
write_tsv(all, 'results/deg/DEG_summary_expanded_all.tsv')
log_msg('Supplementary DEG completed for', length(files), 'supplementary datasets; completed=', sum(union_rbind(summ)$status=='complete', na.rm=TRUE))

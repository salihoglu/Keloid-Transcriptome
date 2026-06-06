source('scripts/common.R')
dir_create('data/processed'); dir_create('data/metadata'); dir_create('reports')
args <- commandArgs(trailingOnly=TRUE)
exclude <- character()
if ('--exclude-accessions' %in% args) exclude <- strsplit(args[which(args=='--exclude-accessions')+1], ',')[[1]]
exclude <- unique(exclude)
respect_manifest <- '--respect-inclusion-decisions' %in% args
manifest_cfg <- get_manifest()
analyze_all_literature <- Sys.getenv('KELOID_ANALYZE_ALL_LITERATURE_GEO', unset='1') %in% c('1','true','TRUE','yes','YES')
allowed_flags <- truthy(manifest_cfg$include_in_primary_deg) | truthy(manifest_cfg$include_in_meta_analysis)
if (analyze_all_literature) allowed_flags <- allowed_flags | truthy(manifest_cfg$include_in_contextual_validation)
allowed_accessions <- unique(as.character(manifest_cfg$accession[allowed_flags]))
files <- list.files('data/raw/geo', pattern='_ExpressionSet\\.rds$', full.names=TRUE)
log_msg('Starting preprocessing script: v45 runtime-safe empty-ExpressionSet/supplementary-route preprocessing')
log_msg('Found expression inputs to preprocess:', paste(basename(files), collapse=', '))
aud <- list()
for (f in files) {
  base <- sub('_ExpressionSet\\.rds$','',basename(f)); acc <- strsplit(base,'_')[[1]][1]
  if (acc %in% exclude) {
    log_msg('Skipping', basename(f), ': excluded accession')
    aud[[length(aud)+1]] <- data.frame(dataset=base, accession=acc, status='skipped', reason='excluded_accession', stringsAsFactors=FALSE)
    next
  }
  if (respect_manifest && length(allowed_accessions) && !(acc %in% allowed_accessions)) {
    log_msg('Skipping', basename(f), ': accession not in active v34 bulk/contextual analysis route')
    aud[[length(aud)+1]] <- data.frame(dataset=base, accession=acc, status='skipped', reason='not_in_active_v34_bulk_contextual_route', stringsAsFactors=FALSE)
    next
  }
  skip_contextual_direct <- Sys.getenv('KELOID_SKIP_DIRECT_CONTEXTUAL_ESET_PREPROCESS', unset='true') %in% c('1','true','TRUE','yes','YES')
  if (skip_contextual_direct && nrow(manifest_cfg) && acc %in% manifest_cfg$accession) {
    row <- manifest_cfg[match(acc, manifest_cfg$accession), , drop=FALSE]
    is_primary <- truthy(row$include_in_primary_deg) | truthy(row$include_in_meta_analysis)
    is_contextual_only <- truthy(row$include_in_contextual_validation) & !is_primary
    if (isTRUE(is_contextual_only)) {
      log_msg('Skipping direct ExpressionSet preprocessing for', basename(f), ': contextual/sensitivity accession routed to supplementary/import audit')
      aud[[length(aud)+1]] <- data.frame(dataset=base, accession=acc, status='skipped', reason='contextual_direct_eset_skipped_v45_supplementary_or_audit_route', stringsAsFactors=FALSE)
      next
    }
  }
  status <- 'failed'; reason <- ''
  tryCatch({
    if (!requireNamespace('Biobase', quietly=TRUE)) stop('Biobase not installed')
    log_msg('Preprocessing', base, 'from', basename(f))
    eset <- readRDS(f)
    expr <- Biobase::exprs(eset)
    feat <- as.data.frame(Biobase::fData(eset), check.names=FALSE)
    # v37 runtime fix: several GEO Series Matrix ExpressionSet objects for processed
    # RNA-seq accessions contain 0 expression rows even though their supplementary files
    # hold the usable count/FPKM/TPM matrix. Treat this as an audited skip, not an error.
    if (is.null(expr) || length(expr) == 0 || nrow(expr) == 0 || ncol(expr) == 0) {
      status <- 'skipped'
      reason <- paste0('empty_or_zero_row_ExpressionSet; nrow=', ifelse(is.null(expr), NA_integer_, nrow(expr)),
                       '; ncol=', ifelse(is.null(expr), NA_integer_, ncol(expr)),
                       '; supplementary_expression_import_required')
      write_tsv(data.frame(dataset=base, accession=acc,
                           symbol_source='none_empty_expression_set', symbol_score=0,
                           feature_id_fallback=FALSE, n_features_written=0,
                           n_samples=ifelse(is.null(expr), 0, ncol(expr)), output='',
                           action='skipped_to_supplementary_route', stringsAsFactors=FALSE),
                file.path('data/metadata', paste0(base, '_preprocessing_symbol_source_audit.tsv')))
      log_msg('SKIPPING preprocessing', basename(f), ':', reason)
      next
    }
    expr <- as.matrix(expr)
    expr <- expr[rowSums(is.finite(suppressWarnings(matrix(as.numeric(expr), nrow=nrow(expr))))) >= 1, , drop=FALSE]
    if (nrow(expr) == 0) {
      status <- 'skipped'
      reason <- 'ExpressionSet has no finite expression rows after numeric coercion; supplementary_expression_import_required'
      write_tsv(data.frame(dataset=base, accession=acc,
                           symbol_source='none_empty_numeric_expression_set', symbol_score=0,
                           feature_id_fallback=FALSE, n_features_written=0,
                           n_samples=ncol(expr), output='',
                           action='skipped_to_supplementary_route', stringsAsFactors=FALSE),
                file.path('data/metadata', paste0(base, '_preprocessing_symbol_source_audit.tsv')))
      log_msg('SKIPPING preprocessing', basename(f), ':', reason)
      next
    }
    feat <- as.data.frame(Biobase::fData(eset), check.names=FALSE)
    expr <- maybe_log2_transform(expr)
    platform_id <- strsplit(base, '_')[[1]][2]
    chosen <- choose_gene_symbol(expr, feat, platform_id=platform_id)
    feature_fallback <- FALSE
    if (!length(chosen$symbol)) {

      rid <- rownames(expr)
      if (is.null(rid) || length(rid) != nrow(expr) || mean(!is.na(rid) & nzchar(rid)) < 0.50) {
        rid <- paste0('FEATURE_', seq_len(nrow(expr)))
      }
      chosen <- list(symbol=make.unique(clean_symbol(rid)), source='feature_id_fallback_no_gene_symbol_review_required', score=0)
      feature_fallback <- TRUE
      log_msg('WARNING preprocessing', basename(f), ': no usable gene-symbol column; using feature-id fallback; dataset will be audited downstream')
    }
    mat <- collapse_by_gene(expr, chosen$symbol)
    if (nrow(mat) == 0) {
 
      rid <- chosen$symbol
      if (length(rid) != nrow(expr) || mean(!is.na(rid) & nzchar(rid)) < 0.50) rid <- rownames(expr)
      if (is.null(rid) || length(rid) != nrow(expr)) rid <- paste0('FEATURE_', seq_len(nrow(expr)))
      rid <- make.unique(clean_symbol(rid))
      keep <- rowSums(is.finite(expr)) >= max(2, ceiling(0.25 * ncol(expr)))
      mat <- expr[keep, , drop=FALSE]
      rownames(mat) <- rid[keep]
      feature_fallback <- TRUE
      chosen$source <- paste0(chosen$source, ';feature_level_no_hgnc_collapse')
      log_msg('WARNING preprocessing', basename(f), ': no mapped gene symbols after collapsing; writing feature-level matrix for downstream audit')
    }
    out <- data.frame(gene_symbol=rownames(mat), mat, check.names=FALSE)
    outpath <- file.path('data/processed', paste0(base,'_expression_preprocessed.tsv.gz'))
    write_tsv(out, outpath)
    write_tsv(data.frame(dataset=base, accession=acc, platform_id=platform_id,
                         symbol_source=chosen$source, symbol_score=round(chosen$score, 4),
                         feature_id_fallback=feature_fallback,
                         n_features_written=nrow(mat), n_samples=ncol(mat),
                         output=basename(outpath), stringsAsFactors=FALSE),
              file.path('data/metadata', paste0(base, '_preprocessing_symbol_source_audit.tsv')))
    status <- 'processed'; reason <- paste0('features=', nrow(mat), ';samples=', ncol(mat), ';symbol_source=', chosen$source, ';feature_fallback=', feature_fallback)
  }, error=function(e){ reason <<- conditionMessage(e); log_msg('ERROR preprocessing', basename(f), ':', reason) })
  aud[[length(aud)+1]] <- data.frame(dataset=base, accession=acc, status=status, reason=reason, stringsAsFactors=FALSE)
}
manifest <- union_rbind(aud)
if (!nrow(manifest)) manifest <- data.frame(dataset=character(), accession=character(), status=character(), reason=character())
write_tsv(manifest, 'data/metadata/preprocessing_manifest.tsv')
writeLines(c('# Preprocessing log','',capture.output(print(manifest))), 'reports/preprocessing_log.txt')
processed <- if (nrow(manifest)) sum(manifest$status=='processed', na.rm=TRUE) else 0
skipped <- if (nrow(manifest)) sum(manifest$status=='skipped', na.rm=TRUE) else 0
failed <- if (nrow(manifest)) sum(manifest$status=='failed', na.rm=TRUE) else 0
log_msg('Preprocessing completed with fail-soft policy. Processed=', processed, '; skipped=', skipped, '; failed=', failed)

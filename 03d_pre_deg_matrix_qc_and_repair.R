#!/usr/bin/env Rscript

source('scripts/common.R')
dir_create('reports'); dir_create('data/processed/quarantined'); dir_create('results/qc/pre_deg_matrix_qc')

min_genes <- as.integer(Sys.getenv('KELOID_MATRIX_QC_MIN_GENES', '5000'))
min_samples <- as.integer(Sys.getenv('KELOID_MATRIX_QC_MIN_SAMPLES', '4'))
min_symbol_fraction <- as.numeric(Sys.getenv('KELOID_MATRIX_QC_MIN_SYMBOL_FRACTION', '0.70'))
max_ercc_fraction <- as.numeric(Sys.getenv('KELOID_MATRIX_QC_MAX_ERCC_FRACTION', '0.05'))
max_lowvar_fraction <- as.numeric(Sys.getenv('KELOID_MATRIX_QC_MAX_LOWVAR_FRACTION', '0.90'))
action <- tolower(Sys.getenv('KELOID_MATRIX_QC_ACTION', 'report_only'))
report_only <- action %in% c('report_only','report-only','none')

infer_scale <- function(mat) {
  vals <- as.numeric(mat[is.finite(mat)])
  if (!length(vals)) return(list(scale='empty', transform='none', q01=NA, q50=NA, q99=NA, max=NA, min=NA))
  qs <- suppressWarnings(quantile(vals, c(.01,.5,.99), na.rm=TRUE))
  mx <- max(vals, na.rm=TRUE); mn <- min(vals, na.rm=TRUE)
  int_like <- mean(abs(vals - round(vals)) < 1e-8, na.rm=TRUE) > 0.95
  transform <- 'none'
  scale <- 'log_like'
  if (is.finite(mx) && is.finite(qs[3]) && mn >= 0 && (qs[3] > 30 || mx > 100)) {
    transform <- 'log2_plus_1'
    scale <- if (int_like) 'count_like_raw' else 'continuous_unlogged'
  }
  list(scale=scale, transform=transform, q01=unname(qs[1]), q50=unname(qs[2]), q99=unname(qs[3]), max=mx, min=mn)
}

repair_matrix <- function(f) {
  x <- read_tsv_safe(f)
  dataset <- sub('_expression_preprocessed.tsv.gz$', '', basename(f))
  acc <- strsplit(dataset, '_')[[1]][1]
  if (!nrow(x) || !('gene_symbol' %in% names(x))) {
    return(list(summary=data.frame(dataset=dataset, accession=acc, file=basename(f), status='fail', action='none', reasons='empty_or_missing_gene_symbol'), repaired=NULL))
  }
  sample_cols <- setdiff(names(x), 'gene_symbol')
  mat0 <- as.matrix(x[, sample_cols, drop=FALSE]); suppressWarnings(storage.mode(mat0) <- 'numeric')
  raw_symbols <- as.character(x$gene_symbol)
  symbols <- strict_hgnc_symbol(map_ensembl_to_symbol(raw_symbols))
  ercc_fraction <- mean(grepl('^ERCC', raw_symbols, ignore.case=TRUE), na.rm=TRUE)
  interpretable_fraction_raw <- mean(!is.na(symbols) & nzchar(symbols), na.rm=TRUE)
  keep <- !is.na(symbols) & nzchar(symbols) & rowSums(is.finite(mat0)) >= max(2, ceiling(0.25*ncol(mat0)))
  mat <- mat0[keep, , drop=FALSE]
  sym <- symbols[keep]
  if (nrow(mat) > 0) {
    rownames(mat) <- sym
    # Collapse duplicates by maximum variance, preserving a representative probe.
    mat <- collapse_by_gene(mat, rownames(mat), method='max_variance')
  }
  sc <- infer_scale(mat)
  if (nrow(mat) > 0 && sc$transform == 'log2_plus_1') mat <- log2(pmax(mat, 0) + 1)
  lowvar_fraction <- if (nrow(mat) > 0) mean(apply(mat, 1, stats::var, na.rm=TRUE) <= 1e-12, na.rm=TRUE) else 1
  reasons <- character(); status <- 'pass'
  if (length(sample_cols) < min_samples) { reasons <- c(reasons, paste0('too_few_samples_', length(sample_cols))); status <- 'fail' }
  if (nrow(mat) < min_genes) { reasons <- c(reasons, paste0('too_few_interpretable_genes_after_repair_', nrow(mat), '_lt_', min_genes)); status <- 'fail' }
  if (is.finite(interpretable_fraction_raw) && interpretable_fraction_raw < min_symbol_fraction) { reasons <- c(reasons, sprintf('low_official_symbol_fraction_raw_%.3f', interpretable_fraction_raw)); status <- 'fail' }
  if (is.finite(ercc_fraction) && ercc_fraction > max_ercc_fraction) { reasons <- c(reasons, sprintf('high_ERCC_fraction_%.3f', ercc_fraction)); status <- 'fail' }
  if (is.finite(lowvar_fraction) && lowvar_fraction > max_lowvar_fraction) { reasons <- c(reasons, sprintf('too_many_low_variance_genes_%.3f', lowvar_fraction)); status <- 'fail' }
  if (!length(reasons)) reasons <- 'passed_preDEG_matrix_qc'
  repaired <- NULL
  if (nrow(mat) > 0) {
    repaired <- data.frame(gene_symbol=rownames(mat), mat, check.names=FALSE)
  }
  summ <- data.frame(dataset=dataset, accession=acc, file=basename(f), status=status,
                     n_rows_input=nrow(x), n_samples=length(sample_cols),
                     official_symbol_fraction_raw=round(interpretable_fraction_raw,4),
                     ercc_fraction=round(ercc_fraction,4), n_genes_after_repair=nrow(mat),
                     lowvar_fraction=round(lowvar_fraction,4), scale=sc$scale,
                     transform_applied=sc$transform, q01=sc$q01, q50=sc$q50, q99=sc$q99,
                     min=sc$min, max=sc$max, reasons=paste(reasons, collapse=';'),
                     action='none', stringsAsFactors=FALSE)
  list(summary=summ, repaired=repaired)
}

files <- sort(list.files('data/processed', pattern='_expression_preprocessed.tsv.gz$', full.names=TRUE, recursive=FALSE))
if (!length(files)) {
  write_tsv(data.frame(), 'results/qc/pre_deg_matrix_qc/pre_deg_matrix_qc_summary.tsv')
  writeLines('# Pre-DEG matrix QC and repair\n\nNo processed expression matrices found.', 'reports/PRE_DEG_MATRIX_QC_AND_REPAIR_REPORT.md')
  quit(save='no')
}

summaries <- list()
for (f in files) {
  res <- tryCatch(repair_matrix(f), error=function(e) {
    dataset <- sub('_expression_preprocessed.tsv.gz$', '', basename(f)); acc <- strsplit(dataset, '_')[[1]][1]
    list(summary=data.frame(dataset=dataset, accession=acc, file=basename(f), status='fail', action='none', reasons=conditionMessage(e), stringsAsFactors=FALSE), repaired=NULL)
  })
  s <- res$summary
  if (s$status == 'pass' && !is.null(res$repaired)) {
    write_tsv(res$repaired, f)
    s$action <- ifelse(s$transform_applied == 'log2_plus_1', 'repaired_symbols_collapsed_and_log2_transformed', 'repaired_symbols_collapsed')
  } else if (s$status == 'fail' && file.exists(f)) {
    if (report_only) {
      s$action <- 'reported_only'
    } else {
      dst <- file.path('data/processed/quarantined', sub('_expression_preprocessed.tsv.gz$', '_expression_preprocessed.QUARANTINED.tsv.gz', basename(f)))
      ok <- file.rename(f, dst)
      s$action <- ifelse(ok, paste0('quarantined_to_', basename(dst)), 'quarantine_failed')
    }
  }
  summaries[[length(summaries)+1]] <- s
}
qc <- union_rbind(summaries)
write_tsv(qc, 'results/qc/pre_deg_matrix_qc/pre_deg_matrix_qc_summary.tsv')
write_tsv(qc[qc$status=='fail',,drop=FALSE], 'results/qc/pre_deg_matrix_qc/pre_deg_matrix_qc_failed.tsv')
write_tsv(qc[qc$status=='pass',,drop=FALSE], 'results/qc/pre_deg_matrix_qc/pre_deg_matrix_qc_passed.tsv')
report <- c('# Pre-DEG matrix QC and repair','',
            paste('Processed matrices audited:', nrow(qc)),
            paste('Passed/repaired:', sum(qc$status=='pass', na.rm=TRUE)),
            paste('Failed/quarantined:', sum(qc$status=='fail', na.rm=TRUE)),
            paste('Official HGNC whitelist available:', official_hgnc_available()),'',
            'This gate is intentionally before PCA and DEG. It prevents annotation descriptors, ERCC-only tables, unlogged expression scales and non-informative matrices from generating misleading volcano plots or consensus signals.', '',
            'Outputs:',
            '- results/qc/pre_deg_matrix_qc/pre_deg_matrix_qc_summary.tsv',
            '- data/processed/quarantined/ for unsafe matrices')
writeLines(report, 'reports/PRE_DEG_MATRIX_QC_AND_REPAIR_REPORT.md')
log_msg('Pre-DEG matrix QC/repair complete. passed=', sum(qc$status=='pass', na.rm=TRUE), ' failed=', sum(qc$status=='fail', na.rm=TRUE), ' report_only=', report_only)

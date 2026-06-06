#!/usr/bin/env Rscript

source('scripts/common.R')
suppressPackageStartupMessages({
  if (!requireNamespace('ggplot2', quietly=TRUE)) stop('ggplot2 is required')
})
dir_create('results/qc/pre_deg_pca')
dir_create('reports')

asg <- read_tsv_safe('data/metadata/sample_group_assignment.tsv')
manifest <- get_manifest()
files <- list.files('data/processed', pattern='_expression_preprocessed.tsv.gz$', full.names=TRUE)
files <- sort(files)

num_safe <- function(x) suppressWarnings(as.numeric(x))
scale_label <- function(mat) {
  qs <- suppressWarnings(quantile(as.numeric(mat), c(0, .01, .25, .5, .75, .99, 1), na.rm=TRUE))
  if (!all(is.finite(qs))) return('unknown')
  if (qs[7] > 100 || qs[6] > 50) return('nonlog_or_count_like')
  if (qs[1] < -5 && qs[7] > 20) return('mixed_or_suspicious_scale')
  'log_like'
}
robust_var <- function(mat) {
  v <- apply(mat, 1, stats::var, na.rm=TRUE)
  v[is.finite(v)]
}
plot_safe_name <- function(x) gsub('[^A-Za-z0-9._-]+', '_', x)

write_pca <- function(expr, map, dataset, acc, out_prefix) {
  y <- expr

  keep <- rowSums(is.finite(y)) >= max(3, floor(ncol(y) * 0.5))
  y <- y[keep,,drop=FALSE]
  rv <- robust_var(y)
  if (!length(rv)) stop('no finite gene-level variance for PCA')
  y <- y[names(sort(rv, decreasing=TRUE))[seq_len(min(length(rv), 5000))],,drop=FALSE]
  y <- y[apply(y, 1, function(z) stats::var(z, na.rm=TRUE) > 0),,drop=FALSE]
  if (nrow(y) < 2 || ncol(y) < 3) stop('insufficient matrix dimensions for PCA')

  if (anyNA(y)) {
    med <- apply(y, 1, stats::median, na.rm=TRUE)
    idx <- which(is.na(y), arr.ind=TRUE)
    y[idx] <- med[idx[,1]]
  }
  pcs <- stats::prcomp(t(y), center=TRUE, scale.=FALSE)
  ve <- round(100 * (pcs$sdev^2 / sum(pcs$sdev^2))[1:2], 1)
  pdat <- data.frame(sample_col=rownames(pcs$x), PC1=pcs$x[,1], PC2=pcs$x[,2], stringsAsFactors=FALSE)
  pdat <- merge(pdat, map, by='sample_col', all.x=TRUE, sort=FALSE)
  if (!'assigned_group' %in% names(pdat)) pdat$assigned_group <- 'unknown'
  if (!'sample_type' %in% names(pdat)) pdat$sample_type <- 'unspecified'
  pdat$assigned_group[is.na(pdat$assigned_group) | pdat$assigned_group==''] <- 'unassigned'
  pdat$sample_type[is.na(pdat$sample_type) | pdat$sample_type==''] <- 'unspecified'
  p <- ggplot2::ggplot(pdat, ggplot2::aes(PC1, PC2, shape=sample_type, color=assigned_group, label=sample_col)) +
    ggplot2::geom_point(size=3, alpha=.88) +
    ggplot2::labs(title=paste0(acc, ' pre-DEG PCA'), subtitle=dataset,
                  x=paste0('PC1 (', ve[1], '%)'), y=paste0('PC2 (', ve[2], '%)'),
                  color='Group', shape='Sample type') +
    ggplot2::theme_classic(base_size=11) +
    ggplot2::theme(plot.title=ggplot2::element_text(face='bold'), legend.position='right')
  ggplot2::ggsave(paste0(out_prefix, '_PCA.pdf'), p, width=6.8, height=5.2, device=grDevices::cairo_pdf, bg='white')
  ggplot2::ggsave(paste0(out_prefix, '_PCA.png'), p, width=6.8, height=5.2, dpi=300, bg='white')
  write_tsv(pdat, paste0(out_prefix, '_PCA_coordinates.tsv'))
  invisible(pdat)
}

write_corr_heatmap <- function(expr, dataset, acc, out_prefix) {
  y <- expr
  keep <- apply(y, 1, function(z) stats::var(z, na.rm=TRUE) > 0)
  y <- y[keep,,drop=FALSE]
  if (ncol(y) < 3 || nrow(y) < 2) return(FALSE)
  if (anyNA(y)) {
    med <- apply(y, 1, stats::median, na.rm=TRUE)
    idx <- which(is.na(y), arr.ind=TRUE); y[idx] <- med[idx[,1]]
  }
  cm <- suppressWarnings(stats::cor(y, use='pairwise.complete.obs', method='pearson'))
  long <- as.data.frame(as.table(cm), stringsAsFactors=FALSE)
  names(long) <- c('sample_1','sample_2','correlation')
  p <- ggplot2::ggplot(long, ggplot2::aes(sample_1, sample_2, fill=correlation)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(limits=c(-1,1), midpoint=0) +
    ggplot2::labs(title=paste0(acc, ' sample correlation'), subtitle=dataset, x=NULL, y=NULL, fill='r') +
    ggplot2::theme_minimal(base_size=8) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1, vjust=.5), panel.grid=ggplot2::element_blank())
  ggplot2::ggsave(paste0(out_prefix, '_sample_correlation.pdf'), p, width=6.8, height=6.2, device=grDevices::cairo_pdf, bg='white')
  ggplot2::ggsave(paste0(out_prefix, '_sample_correlation.png'), p, width=6.8, height=6.2, dpi=300, bg='white')
  write_tsv(long, paste0(out_prefix, '_sample_correlation_long.tsv'))
  TRUE
}

write_distribution <- function(expr, dataset, acc, out_prefix) {
  if (ncol(expr) < 1 || nrow(expr) < 10) return(FALSE)
  vals <- as.data.frame(expr[seq_len(min(nrow(expr), 10000)),,drop=FALSE], check.names=FALSE)
  long <- utils::stack(vals)
  names(long) <- c('expression','sample_col')
  long$expression <- num_safe(long$expression)
  long <- long[is.finite(long$expression),]
  if (!nrow(long)) return(FALSE)
  p <- ggplot2::ggplot(long, ggplot2::aes(sample_col, expression)) +
    ggplot2::geom_boxplot(outlier.size=.25, linewidth=.25) +
    ggplot2::labs(title=paste0(acc, ' expression distribution'), subtitle=dataset, x=NULL, y='Expression') +
    ggplot2::theme_classic(base_size=8) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1, vjust=.5), plot.title=ggplot2::element_text(face='bold'))
  ggplot2::ggsave(paste0(out_prefix, '_boxplot.pdf'), p, width=7.2, height=4.8, device=grDevices::cairo_pdf, bg='white')
  ggplot2::ggsave(paste0(out_prefix, '_boxplot.png'), p, width=7.2, height=4.8, dpi=300, bg='white')
  TRUE
}

qc <- list()
for (f in files) {
  dataset <- sub('_expression_preprocessed.tsv.gz$','', basename(f))
  acc <- strsplit(dataset, '_')[[1]][1]
  prefix <- file.path('results/qc/pre_deg_pca', plot_safe_name(dataset))
  msg <- ''
  status <- 'complete'
  tryCatch({
    x <- read_tsv_safe(f)
    if (!nrow(x) || !('gene_symbol' %in% names(x))) stop('empty or missing gene_symbol')
    genes <- as.character(x$gene_symbol)
    expr <- as.matrix(x[, setdiff(names(x), 'gene_symbol'), drop=FALSE]); storage.mode(expr) <- 'numeric'; rownames(expr) <- genes
    map <- resolve_sample_groups(colnames(expr), acc, asg)
    if (!'sample_col' %in% names(map)) map$sample_col <- colnames(expr)
    map <- map[match(colnames(expr), map$sample_col),,drop=FALSE]
    if (!nrow(map)) map <- data.frame(sample_col=colnames(expr), assigned_group='unassigned', sample_type='unspecified')
    if (!'sample_type' %in% names(map)) map$sample_type <- 'unspecified'
    map$sample_type[is.na(map$sample_type) | map$sample_type==''] <- mapply(infer_sample_type, map$sample_col[is.na(map$sample_type) | map$sample_type==''], USE.NAMES=FALSE)
    write_pca(expr, map, dataset, acc, prefix)
    write_corr_heatmap(expr, dataset, acc, prefix)
    write_distribution(expr, dataset, acc, prefix)
    v <- robust_var(expr)
    q <- suppressWarnings(quantile(as.numeric(expr), c(.01,.5,.99), na.rm=TRUE))
    qc[[length(qc)+1]] <- data.frame(dataset=dataset, accession=acc, status=status,
      n_genes=nrow(expr), n_samples=ncol(expr), n_keloid=sum(map$assigned_group=='keloid', na.rm=TRUE),
      n_control=sum(map$assigned_group=='control', na.rm=TRUE), n_unassigned=sum(!(map$assigned_group %in% c('keloid','control')), na.rm=TRUE),
      scale_status=scale_label(expr), q01=q[1], median=q[2], q99=q[3], zero_variance_genes=sum(v==0, na.rm=TRUE),
      pca_pdf=paste0(prefix, '_PCA.pdf'), notes=msg, stringsAsFactors=FALSE)
  }, error=function(e) {
    qc[[length(qc)+1]] <- data.frame(dataset=dataset, accession=acc, status='failed', n_genes=NA, n_samples=NA,
      n_keloid=NA, n_control=NA, n_unassigned=NA, scale_status='unknown', q01=NA, median=NA, q99=NA,
      zero_variance_genes=NA, pca_pdf='', notes=conditionMessage(e), stringsAsFactors=FALSE)
    log_msg('Pre-DEG QC failed for ', dataset, ': ', conditionMessage(e))
  })
}
qc_df <- union_rbind(qc)
write_tsv(qc_df, 'results/qc/pre_deg_pca/pre_deg_qc_summary.tsv')
writeLines(c('# Pre-DEG PCA and sample QC report','',
             'This report is generated before differential-expression modelling. It is intended to detect sample-label mismatch, outliers, scale problems and low-information matrices before volcano plots or consensus analysis are interpreted.', '',
             paste('Processed expression matrices evaluated:', nrow(qc_df)),
             paste('QC failures:', sum(qc_df$status!='complete', na.rm=TRUE)), '',
             'Primary outputs:',
             '- results/qc/pre_deg_pca/*_PCA.pdf / .png',
             '- results/qc/pre_deg_pca/*_sample_correlation.pdf / .png',
             '- results/qc/pre_deg_pca/*_boxplot.pdf / .png',
             '- results/qc/pre_deg_pca/pre_deg_qc_summary.tsv'),
           'reports/PRE_DEG_QC_PCA_REPORT.md')
log_msg('Pre-DEG QC/PCA completed for ', nrow(qc_df), ' processed matrices')

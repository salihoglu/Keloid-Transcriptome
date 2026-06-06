source('scripts/common.R')

dir_create('results/ppi'); dir_create('results/figures'); dir_create('reports'); dir_create('data/raw/stringdb')

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
first_existing <- function(paths) { for (p in paths) if (file.exists(p) && file.info(p)$size > 0) return(p); '' }

load_consensus <- function() {
  p <- first_existing(c('results/meta/all_eligible_quality_weighted_consensus_high_confidence.tsv','results/meta/balanced_consensus_signature.tsv','results/meta/layered_consensus_signature.tsv','results/deg/multi_dataset_gene_level_merged.tsv'))
  if (!nzchar(p)) return(data.frame())
  x <- read_tsv_safe(p)
  if (!nrow(x) || !('gene_symbol' %in% names(x))) return(data.frame())
  x$gene_symbol <- strict_hgnc_symbol(x$gene_symbol)
  x <- x[!is.na(x$gene_symbol) & nzchar(x$gene_symbol), , drop=FALSE]
  score_col <- intersect(c('v42_consensus_score','consensus_score','weighted_signed_evidence','mean_neglog10p'), names(x))[1]
  if (is.na(score_col)) { x$evidence_score <- 1 } else { x$evidence_score <- abs(as_num(x[[score_col]])); x$evidence_score[!is.finite(x$evidence_score)] <- 0 }
  dir <- rep('', nrow(x))
  for (dc in intersect(c('weighted_direction','consensus_direction','direction'), names(x))) dir[dir==''] <- tolower(as.character(x[[dc]][dir=='']))
  if ('weighted_mean_logFC' %in% names(x)) dir[dir==''] <- ifelse(as_num(x$weighted_mean_logFC[dir=='']) >= 0, 'up', 'down')
  if ('mean_logFC' %in% names(x)) dir[dir==''] <- ifelse(as_num(x$mean_logFC[dir=='']) >= 0, 'up', 'down')
  x$direction <- dir
  x <- x[order(-x$evidence_score), ]
  x[!duplicated(x$gene_symbol), c('gene_symbol','evidence_score','direction'), drop=FALSE]
}

curated_edges <- function() {
  data.frame(protein1=c('TGFB2','TGFB2','TGFBR1','TGFBR1','COL1A1','COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL6A1','FN1','FN1','FN1','POSTN','POSTN','SPARC','SPARC','BGN','SERPINH1','LOXL2','LOXL2','TGFBI','TGFBI','COL1A1','COL3A1','MMP2','ITGB1','ITGA5'),
             protein2=c('TGFBR1','SMAD3','SMAD2','SMAD3','COL1A2','COL3A1','COL3A1','COL5A1','COL5A2','COL6A1','FN1','ITGA5','ITGB1','COL1A1','COL1A1','FN1','COL1A1','COL3A1','COL1A1','COL1A1','COL1A1','COL3A1','ITGB1','FN1','SERPINH1','SERPINH1','COL1A1','PTK2','PTK2'),
             combined_score=c(880,820,850,850,900,890,880,860,830,810,760,850,840,810,820,780,800,790,760,760,780,770,760,740,760,750,700,700,700),
             edge_source='curated_fibrosis_ECM_prior', stringsAsFactors=FALSE)
}

cons <- load_consensus()
if (!nrow(cons)) {
  write_tsv(data.frame(protein1=character(), protein2=character(), combined_score=numeric(), edge_source=character()), 'results/ppi/string_interactions.tsv')
  write_tsv(data.frame(status='skipped_no_input'), 'results/ppi/ppi_audit.tsv')
  log_msg('PPI/network skipped because no consensus/DEG input was available.'); quit(save='no')
}
anchors <- c('TGFB2','TGFBR1','SMAD2','SMAD3','COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL6A1','POSTN','SPARC','LOXL2','SERPINH1','FN1','TGFBI','BGN','ITGA5','ITGB1','MMP2','PTK2')
input_genes <- unique(c(head(cons$gene_symbol, 200), anchors))

edges <- data.frame(); mode <- 'curated_fallback'
if (requireNamespace('STRINGdb', quietly=TRUE)) {
  sdb <- tryCatch(STRINGdb::STRINGdb$new(version='12.0', species=9606, score_threshold=400, input_directory='data/raw/stringdb'), error=function(e) NULL)
  if (!is.null(sdb)) {
    mapped <- tryCatch(sdb$map(data.frame(gene_symbol=input_genes), 'gene_symbol', removeUnmappedRows=TRUE), error=function(e) data.frame())
    if (nrow(mapped) >= 3) {
      inter <- tryCatch(sdb$get_interactions(mapped$STRING_id), error=function(e) data.frame())
      if (nrow(inter)) {
        revmap <- mapped$gene_symbol; names(revmap) <- mapped$STRING_id
        edges <- data.frame(protein1=unname(revmap[inter$from]), protein2=unname(revmap[inter$to]), combined_score=as_num(inter$combined_score), edge_source='STRINGdb', stringsAsFactors=FALSE)
        edges <- edges[!is.na(edges$protein1) & !is.na(edges$protein2) & edges$protein1 != edges$protein2, , drop=FALSE]
        mode <- 'STRINGdb'
      }
    }
  }
}
if (!nrow(edges)) {
  ce <- curated_edges()
  keep <- ce$protein1 %in% input_genes | ce$protein2 %in% input_genes
  edges <- ce[keep, , drop=FALSE]
  if (!nrow(edges)) edges <- ce
}
edges$key <- apply(edges[, c('protein1','protein2')], 1, function(z) paste(sort(z), collapse='__'))
edges <- edges[!duplicated(edges$key), c('protein1','protein2','combined_score','edge_source'), drop=FALSE]
write_tsv(edges, 'results/ppi/string_interactions.tsv')
write_tsv(edges, 'results/ppi/ppi_network_edges.tsv')

deg <- table(c(edges$protein1, edges$protein2))
hubs <- data.frame(gene_symbol=names(deg), degree=as.integer(deg), stringsAsFactors=FALSE)
hubs$evidence_score <- cons$evidence_score[match(hubs$gene_symbol, cons$gene_symbol)]
hubs$evidence_score[!is.finite(hubs$evidence_score)] <- 0
hubs$network_score <- hubs$degree * log1p(hubs$evidence_score)
hubs$direction <- cons$direction[match(hubs$gene_symbol, cons$gene_symbol)]
hubs <- hubs[order(-hubs$network_score, -hubs$degree), ]
write_tsv(hubs, 'results/ppi/hub_genes.tsv')
write_tsv(data.frame(status='completed', network_mode=mode, input_genes=length(input_genes), interactions=nrow(edges), hub_genes=nrow(hubs), note=ifelse(mode=='STRINGdb','STRINGdb interactions retrieved locally','STRINGdb unavailable/empty; curated fibrosis/ECM PPI-prior fallback used and explicitly annotated')), 'results/ppi/ppi_audit.tsv')

if (requireNamespace('ggplot2', quietly=TRUE) && nrow(hubs)) {
  suppressPackageStartupMessages(library(ggplot2))
  h <- head(hubs, 20); h$gene_symbol <- factor(h$gene_symbol, levels=rev(h$gene_symbol))
  p <- ggplot(h, aes(x=network_score, y=gene_symbol, size=degree)) + geom_point() + theme_bw(base_size=10) + labs(title='Network-prioritized fibrosis/ECM hub genes', x='Network score: degree × log1p(signature evidence)', y=NULL)
  ggsave('results/figures/Figure_ppi_network_hubs.png', p, width=7, height=5, dpi=300)
  ggsave('results/figures/Figure_ppi_network_hubs.pdf', p, width=7, height=5)
}


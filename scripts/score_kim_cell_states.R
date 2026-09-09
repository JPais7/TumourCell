args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop('usage: Rterm --file=score_kim_cell_states.R --args <counts.csv> <metadata.csv> <gtf> <out_dir>')
counts_file <- args[[1]]; metadata_file <- args[[2]]; gtf_file <- args[[3]]; out_dir <- args[[4]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
counts <- read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
genes <- counts[[1]]; x <- as.matrix(counts[-1]); rownames(x) <- genes
lib <- colSums(x); logcpm <- log2(sweep(x, 2, pmax(lib, 1), '/') * 1e6 + 1)

ann <- new.env(hash = TRUE, parent = emptyenv())
con <- file(gtf_file, open = 'r'); on.exit(close(con), add = TRUE)
repeat {
  z <- readLines(con, n = 100000L)
  if (!length(z)) break
  z <- z[substr(z, 1, 1) != '#']
  gid <- sub('.*gene_id "([^"]+)".*', '\\1', z)
  gname <- sub('.*gene_name "([^"]+)".*', '\\1', z)
  ok <- grepl('gene_id "', z) & grepl('gene_name "', z) & gid != z & gname != z
  for (i in which(ok)) if (!exists(gid[i], ann, inherits = FALSE)) assign(gid[i], gname[i], ann)
}
symbols <- vapply(genes, function(g) if (exists(g, ann, inherits = FALSE)) get(g, ann) else NA_character_, character(1))
write.csv(data.frame(gene_id = genes, gene_symbol = symbols), file.path(out_dir, 'kim_gene_id_symbol_map.csv'), row.names = FALSE, quote = FALSE)

sets <- list(
  proliferation = c('MKI67','TOP2A','PCNA','TYMS','MCM2','MCM5','MCM6','CCNB1','CCNB2','CDK1'),
  EMT = c('VIM','ZEB1','ZEB2','SNAI1','SNAI2','TWIST1','FN1','ITGA5','MMP2','MMP9'),
  hypoxia = c('HIF1A','VEGFA','CA9','LDHA','SLC2A1','ENO1','NDRG1','ADM','DDIT4'),
  interferon = c('ISG15','IFI6','IFIT1','IFIT2','IFIT3','MX1','OAS1','OAS2','STAT1','IRF7'),
  MYC = c('MYC','MCM2','MCM5','MCM6','LDHA','ODC1','NPM1',' nucleolin'),
  mTOR = c('MTOR','RPTOR','EIF4EBP1','RPS6KB1','AKT1','PIK3CA'),
  angiogenesis = c('VEGFA','KDR','FLT1','ANGPT1','ANGPT2','PECAM1','VWF'),
  resistance = c('ABCB1','ABCC1','ABCG2','ALDH1A1','SOX2','NANOG','PROM1','EPCAM')
)
score <- matrix(NA_real_, nrow = length(sets), ncol = ncol(logcpm), dimnames = list(names(sets), colnames(logcpm)))
for (s in names(sets)) {
  ix <- which(!is.na(symbols) & symbols %in% sets[[s]])
  if (length(ix)) score[s, ] <- colMeans(logcpm[ix, , drop = FALSE])
}
scores <- data.frame(cell_id = colnames(logcpm), t(score), check.names = FALSE)
scores <- merge(scores, meta, by = 'cell_id', all.x = TRUE, sort = FALSE)
write.csv(scores, file.path(out_dir, 'kim_cell_state_scores.csv'), row.names = FALSE, quote = FALSE)
writeLines(c('Scores are mean log2(CPM+1) over predefined marker genes found in the GRCh38 GTF.',
             'Scores are exploratory, not validated malignant-state labels.',
             'MYC score is limited to genes present in the supplied marker list and annotation.'),
           file.path(out_dir, 'README_scores.md'))

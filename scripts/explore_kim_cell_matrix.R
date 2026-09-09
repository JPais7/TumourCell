args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop('usage: Rterm --file=explore_kim_cell_matrix.R --args <counts.csv> <metadata.csv> <out_dir>')
counts_file <- args[[1]]; metadata_file <- args[[2]]; out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
counts <- read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
genes <- counts[[1]]; x <- as.matrix(counts[-1]); rownames(x) <- genes
colnames(x) <- sub('\\.bam$', '', colnames(x))
meta <- meta[match(colnames(x), meta$cell_id), , drop = FALSE]
if (anyNA(meta$cell_id)) stop('metadata does not cover all matrix columns')

lib <- colSums(x); cpm <- sweep(x, 2, pmax(lib, 1), '/') * 1e6
logcpm <- log2(cpm + 1)
v <- apply(logcpm, 1, var)
keep <- order(v, decreasing = TRUE)[seq_len(min(1000, nrow(logcpm)))]
pca <- prcomp(t(logcpm[keep, , drop = FALSE]), center = TRUE, scale. = TRUE)
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
emb <- data.frame(cell_id = colnames(x), PC1 = pca$x[,1], PC2 = pca$x[,2],
                  PC3 = pca$x[,3], PC4 = pca$x[,4])
emb <- cbind(emb, meta[, c('run_accession','clinical_id','candidate_timepoint',
                           'timepoint_confidence','sample_alias')])
set.seed(441)
k <- min(4, nrow(emb) - 1)
km <- kmeans(pca$x[, seq_len(min(10, ncol(pca$x))), drop = FALSE], centers = k)
emb$cluster <- paste0('C', km$cluster)
write.csv(emb, file.path(out_dir, 'kim_pilot_pca_clusters.csv'), row.names = FALSE, quote = FALSE)
png(file.path(out_dir, 'kim_pilot_pca_clusters.png'), width = 1200, height = 900, res = 150)
cols <- as.integer(factor(emb$cluster))
plot(emb$PC1, emb$PC2, col = cols, pch = 19, xlab = 'PC1', ylab = 'PC2',
     main = 'Kim pilot PCA: cell libraries')
text(emb$PC1, emb$PC2, labels = emb$clinical_id, pos = 3, cex = 0.55)
legend('topright', legend = sort(unique(emb$cluster)), col = seq_along(sort(unique(emb$cluster))), pch = 19)
dev.off()
write.csv(data.frame(component = paste0('PC', seq_along(var_explained)),
                     variance_fraction = var_explained),
          file.path(out_dir, 'kim_pilot_pca_variance.csv'), row.names = FALSE)

group_summary <- aggregate(cbind(PC1, PC2, PC3, PC4) ~ clinical_id + candidate_timepoint,
                           data = emb, FUN = mean)
write.csv(group_summary, file.path(out_dir, 'kim_pilot_group_summary.csv'), row.names = FALSE)

top_by_cluster <- lapply(sort(unique(emb$cluster)), function(cl) {
  a <- which(emb$cluster == cl); b <- which(emb$cluster != cl)
  score <- rowMeans(logcpm[, a, drop = FALSE]) - rowMeans(logcpm[, b, drop = FALSE])
  ix <- order(score, decreasing = TRUE)[seq_len(min(50, length(score)))]
  data.frame(cluster = cl, rank = seq_along(ix), gene_id = rownames(logcpm)[ix], score = score[ix])
})
write.csv(do.call(rbind, top_by_cluster), file.path(out_dir, 'kim_pilot_top_genes_by_cluster.csv'), row.names = FALSE)

cond <- paste(meta$clinical_id, meta$candidate_timepoint, sep = '__')
condition_tables <- lapply(sort(unique(cond)), function(g) {
  a <- which(cond == g); b <- which(cond != g)
  if (length(a) < 2 || length(b) < 2) return(NULL)
  score <- rowMeans(logcpm[, a, drop = FALSE]) - rowMeans(logcpm[, b, drop = FALSE])
  ix <- order(abs(score), decreasing = TRUE)[seq_len(min(100, length(score)))]
  data.frame(group = g, rank = seq_along(ix), gene_id = rownames(logcpm)[ix], score = score[ix])
})
write.csv(do.call(rbind, condition_tables), file.path(out_dir, 'kim_pilot_top_genes_by_condition.csv'), row.names = FALSE)
write.csv(data.frame(cell_id = colnames(x), library_size = lib, detected_genes = colSums(x > 0), meta),
          file.path(out_dir, 'kim_pilot_cell_qc.csv'), row.names = FALSE, quote = FALSE)
writeLines(c('PCA and k-means exploration generated with R base.',
             'The matrix is gene-by-cell-library; no UMAP package was assumed.',
             'Timepoints remain candidate labels and clone identity remains unresolved.',
             'Gene IDs are retained as supplied by the GRCh38 GTF; no symbol conversion was assumed.'),
           file.path(out_dir, 'README.md'))

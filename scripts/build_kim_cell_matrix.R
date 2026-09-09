args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop('usage: Rterm --file=build_kim_cell_matrix.R --args <counts.csv> <manifest.csv> <out_dir>')
counts_file <- args[[1]]
manifest_file <- args[[2]]
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
run <- sub('\\.bam$', '', names(counts)[-1])
idx <- match(run, manifest$run_accession)
if (anyNA(idx)) stop('count matrix contains runs absent from manifest')

cell_id <- make.unique(paste(manifest$clinical_id[idx], manifest$sample_alias[idx], sep = '__'))
names(counts)[-1] <- cell_id
metadata <- manifest[idx, c('run_accession','clinical_id','candidate_timepoint',
                            'timepoint_confidence','classification_basis',
                            'sample_accession','sample_alias','library_strategy')]
metadata$cell_id <- cell_id
metadata <- metadata[, c('cell_id', setdiff(names(metadata), 'cell_id'))]

write.csv(counts, file.path(out_dir, 'kim_pilot_gene_by_cell_counts.csv'),
          row.names = FALSE, quote = FALSE)
write.csv(metadata, file.path(out_dir, 'kim_pilot_cell_metadata.csv'),
          row.names = FALSE, quote = FALSE)
writeLines(c(
  'Interpretation: each RNA-seq run is treated as one cell/nucleus library because the ENA sample_alias identifies a patient-specific cell.',
  'Timepoints remain candidate labels: 0cell -> pre_candidate, 2cell -> mid_candidate, OPcell -> unknown.',
  'Clone identity is not inferred by this file; it must be linked from the DNA-cell tables/clone assignments.'
), file.path(out_dir, 'README.md'))

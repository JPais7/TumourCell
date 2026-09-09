args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop('usage: Rscript quantify_kim_pilot.R <fastq_dir> <index_base> <gtf> <out_dir>')

fastq_dir <- normalizePath(args[[1]], mustWork = TRUE)
index_base <- args[[2]]
gtf <- normalizePath(args[[3]], mustWork = TRUE)
out_dir <- args[[4]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, 'bam'), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(Rsubread))
paired1 <- sort(list.files(fastq_dir, pattern = '_1\\.fastq\\.gz$', full.names = TRUE))
paired2 <- sub('_1\\.fastq\\.gz$', '_2.fastq.gz', paired1)
if (any(!file.exists(paired2))) stop('missing mate FASTQ')
single <- sort(list.files(fastq_dir, pattern = '^[^_]+\\.fastq\\.gz$',
                          full.names = TRUE))
fq1 <- c(paired1, single)
fq2 <- c(paired2, rep(NA_character_, length(single)))
if (!length(fq1)) stop('no FASTQ files found')

sample <- sub('_1\\.fastq\\.gz$', '', basename(fq1))
sample <- sub('\\.fastq\\.gz$', '', sample)
bam <- file.path(out_dir, 'bam', paste0(sample, '.bam'))

for (i in seq_along(fq1)) {
  if (!file.exists(bam[[i]])) {
    if (is.na(fq2[[i]])) {
      align(index = index_base, readfile1 = fq1[[i]], output_file = bam[[i]],
            type = 'rna', nthreads = 4, input_format = 'gzFASTQ',
            output_format = 'BAM', unique = TRUE)
    } else {
      align(index = index_base, readfile1 = fq1[[i]], readfile2 = fq2[[i]],
            output_file = bam[[i]], type = 'rna', nthreads = 4,
            input_format = 'gzFASTQ', output_format = 'BAM', unique = TRUE)
    }
  }
}

fc <- featureCounts(files = bam, annot.ext = gtf, isGTFAnnotationFile = TRUE,
                    GTF.featureType = 'exon', GTF.attrType = 'gene_id',
                    isPairedEnd = c(rep(TRUE, length(paired1)), rep(FALSE, length(single))),
                    countReadPairs = FALSE,
                    nthreads = 4, strandSpecific = 0)
counts <- as.data.frame(fc$counts, check.names = FALSE)
counts <- cbind(gene_id = rownames(counts), counts)
write.csv(counts, file.path(out_dir, 'kim_pilot_gene_by_library_counts.csv'),
          row.names = FALSE, quote = FALSE)
write.csv(fc$stat, file.path(out_dir, 'kim_pilot_featurecounts_summary.csv'),
          quote = FALSE)
write.csv(data.frame(library = sample, r1 = fq1, r2 = fq2, bam = bam),
          file.path(out_dir, 'kim_pilot_alignment_manifest.csv'), row.names = FALSE)

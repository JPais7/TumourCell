args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop('usage: Rterm --file=paired_kim_pilot_analysis.R --args <scores.csv> <featurecounts.csv> <out_dir>')
scores_file <- args[[1]]; fc_file <- args[[2]]; out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
scores <- read.csv(scores_file, stringsAsFactors = FALSE)
fc <- read.csv(fc_file, check.names = FALSE, stringsAsFactors = FALSE)
status <- fc[fc$Status == 'Assigned', , drop = FALSE]
run <- sub('\\.bam$', '', names(fc)[-(1:2)])
assigned <- as.numeric(status[1, -(1:2)])
total <- as.numeric(fc[fc$Status == 'Assigned', -(1:2)]) + as.numeric(fc[fc$Status == 'Unassigned_Unmapped', -(1:2)])
qc <- data.frame(run_accession = run, assigned = assigned, total = total,
                 assigned_fraction = assigned / pmax(total, 1))
qc$qc_flag <- ifelse(qc$assigned_fraction >= 0.10, 'pass_exploratory', 'low_assignment_flag')
write.csv(qc, file.path(out_dir, 'kim_pilot_library_qc.csv'), row.names = FALSE, quote = FALSE)

sc <- c('proliferation','EMT','hypoxia','interferon','MYC','mTOR','angiogenesis','resistance')
z <- merge(scores, qc, by = 'run_accession', all.x = TRUE)
write.csv(z, file.path(out_dir, 'kim_pilot_scores_with_qc.csv'), row.names = FALSE, quote = FALSE)

pre <- z[z$candidate_timepoint == 'pre_candidate', ]
mid <- z[z$candidate_timepoint == 'mid_candidate', ]
keys <- intersect(pre$clinical_id, mid$clinical_id)
diffs <- do.call(rbind, lapply(keys, function(p) {
  a <- pre[pre$clinical_id == p, ]; b <- mid[mid$clinical_id == p, ]
  data.frame(clinical_id = p, n_pre = nrow(a), n_mid = nrow(b),
             patient_clonal_outcome = NA_character_,
             setNames(as.list(colMeans(b[sc], na.rm = TRUE) - colMeans(a[sc], na.rm = TRUE)), paste0('delta_', sc)))
}))
link <- file.path(dirname(scores_file), '..', 'kim_pilot_cell_matrix', 'kim_pilot_cell_clone_linkage.csv')
if (file.exists(link)) {
  lk <- read.csv(link, stringsAsFactors = FALSE)
  out <- unique(lk[, c('clinical_id','patient_clonal_outcome')])
  diffs$patient_clonal_outcome <- out$patient_clonal_outcome[match(diffs$clinical_id, out$clinical_id)]
}
write.csv(diffs, file.path(out_dir, 'kim_pilot_paired_mid_minus_pre.csv'), row.names = FALSE, quote = FALSE)

overall <- data.frame(score = sc,
                      mean_delta = vapply(sc, function(s) mean(diffs[[paste0('delta_', s)]], na.rm = TRUE), numeric(1)),
                      n_patients = nrow(diffs))
write.csv(overall, file.path(out_dir, 'kim_pilot_paired_effect_summary.csv'), row.names = FALSE, quote = FALSE)
writeLines(c('QC flag threshold: assigned fraction >= 0.10 is exploratory pass; this is not a validated biological cutoff.',
             'Timepoint comparison is paired at patient level and uses mean cell-library scores per patient/timepoint.',
             'Only three patients have both candidate timepoints; inferential power is insufficient for confirmatory claims.'),
           file.path(out_dir, 'README.md'))

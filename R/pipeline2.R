#' pipeline2
#' Pipeline 2 is best used when the main goal is to identify potential biomarkers, rather than characterize the original
#' composition.
#'
#' @family pipeline2
#'
#' @import DescTools
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta dataframe with columns is_control, sample_type
#' @param blocklist Vector of known previously identified contaminant features
#' @param control_name Character name of controls within metadata
#' @param remove_if Threshold for number of steps feature must be identified as potential contaminant to be removed from final cleaned count matrix. Default set to 1.
#' @param step2_threshold Threshold value for prevalence method of decontam
#' @param seed Random seed
#' @return List object with contaminant matrix, decontaminated count matrix, character string of all removed contaminants,
#' and filtering loss (FL) statistic
#' @export

pipeline2 = function(counts, meta, blocklist, control_name, remove_if = 1,
                     step2_threshold = 0.5, seed = 42) {

  set.seed(seed)

  # Step 0: Check if there is multiple batches
  if (length(unique(meta$batch))<2) {
    stop('Only one batch detected. Use pipeline1.')
  }

  # Step 0: W2W check
  well2well(counts, meta, control_name = control_name, seed = seed)

  # Step 1: Remove features that showed different abundance in different batches
  ## ancombc comparison across batches

  s1_res = step1(counts, meta)

  # Step 2: Remove features that are differentially abundant in negative controls
  ## decontam

  s2_res = step2(counts, meta, step2_threshold)

  # Step 4: Remove known 'blocklist' of contaminants

  s4_res = step4(counts, blocklist)

  # Create dataframe indicating TRUE if contaminant and FALSE if not tagged
  res = data.frame('feature' = colnames(counts))

  res = data.frame('feature' = colnames(counts),
                   'step1' = ifelse(colnames(counts) %in% s1_res, TRUE, FALSE),
                   'step2' = ifelse(colnames(counts) %in% s2_res, TRUE, FALSE),
                   'step4' = ifelse(colnames(counts) %in% s4_res, TRUE, FALSE))

  # return column with summed cases where feature was true

  # transpose to same as counts matrix
  res2 = res
  res2$remove = rowSums(res2[,-c(1)])
  res2 = t(res2)

  # remove features above specified threshold from original counts frame
  counts_rem = rbind(counts, res2['remove',])
  rownames(counts_rem)[nrow(counts_rem)] = 'remove'
  final_counts = counts[,counts_rem['remove',]<remove_if]

  removed = setdiff(colnames(counts), colnames(final_counts))

  # determine filtering loss value

  FL = FL(counts, final_counts)

  # Create deliverable

  deliv = list('contaminant_id' = res,
               'decontaminated_count' = final_counts,
               'removed' = removed,
               'filtering_loss' = FL,
               'pipeline' = 'pipeline2')

  return(deliv)
}

#' step1
#'
#' Run step 1 of pipeline 2 to identify features that are differentially abundant
#' between batches
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta Matrix with columns is_control, sample_type, and batch
#'
#' @family pipeline2
#'
#' @return  Vector of features tagged as contaminants
#' @export

step1 = function(counts, meta) {

  phyloseq = wrap_phyloseq(counts, meta)

  # run differential analysis
  suppressWarnings({
    s1 = ANCOMBC::ancombc(data = phyloseq, assay_name = "counts",  # use package ANCOMBC
                 group = "batch", p_adj_method = "BH", lib_cut = 0,
                 formula = "batch",
                 struc_zero = TRUE, neg_lb = FALSE,
                 tol = 1e-5, max_iter = 100, conserve = FALSE,
                 alpha = 0.05, global = TRUE)
  })
  # create results matrix
  s1_res = do.call(cbind, s1$res)

  # identify column for diff results
  col = ncol(s1_res)

  # return indices for which differentially abundant across batches
  ind = which(s1_res[,col]==TRUE)

  # return list of tagged contaminant features
  return(s1_res[ind,1])
}

#' step2
#'
#' Run step 2 of pipeline 2 to identify features that are expressed higher in negative
#' controls and lower in samples
#'
#' @family pipeline2
#'
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta Matrix with columns is_control, sample_type, and batch
#' @param threshold Threshold value for prevalence method of decontam
#' @return Vector of features tagged as contaminants
#' @export

step2 = function(counts, meta, threshold) {

  # subset to only batches that contain negative controls

  # create phyloseq object
  phyloseq =  wrap_phyloseq(counts, meta)

  # run decontam prevalence method
  s2_res = decontam::isContaminant(phyloseq, method="prevalence", neg="is_control", threshold=threshold)

  # return indices for which features identified as contaminant by decontam prevalence method
  ind = which(s2_res$contaminant)

  # return list of tagged contaminant features
  return(rownames(s2_res[ind,]))

}

#' step4
#'
#' Run step 4 of pipeline 2 to identify features previously identified as contaminants
#' based on blocklist. Requires blocklist items of the lowest taxonomy.
#'
#' @family pipeline2
#'
#' @import DescTools
#'
#' @param blocklist List of known previously identified contaminant features
#' @param counts Count matrix with samples as rows and features as counts
#' @param meta Matrix with columns is_control, sample_type, and batch
#' @return Vector of features tagged as contaminants
#' @export

step4 = function(counts, blocklist) {
  allTaxa = colnames(counts)

  if (class(blocklist)!='character') {
    warning('blocklist input must be formatted as a character vector')
    break
  }

  for (i in length(blocklist)) {
    blocklist = paste0('%', blocklist)
    blocklist = paste0(blocklist, '%')
  }

  return(allTaxa[(tolower(allTaxa) %like any% tolower(blocklist))])
}

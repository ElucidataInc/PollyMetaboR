#' make_group_summary_from_peak_detailed
#'
#' It converts the peak detailed to group summary format of elmaven
#'
#' @param maven_output_df Provide peak detailed format of elmaven
#' @param quant_type The type of qunatification to be used which is from peak detailed format columns
#' @return A dataframe of group summary format
#' @import dplyr
#' @export
make_group_summary_from_peak_detailed <- function(maven_output_df = NULL, quant_type = "peakAreaCorrected"){
  
  message("Make Group Summary From Peak Detailed Started...")
  
  if (identical(maven_output_df, NULL)){
    warning("maven_output_df is NULL dataframe")
    return (NULL)
  }
  
  if (nrow(maven_output_df) < 1){
    warning("Not a valid dataframe, requires atleast one row for processing")
    return (NULL)
  }
  
  required_cols <- c("groupId", "compound", "compoundId", "formula", "sample", "peakMz", "rt")
  
  diff_cols <- setdiff(required_cols, colnames(maven_output_df))
  if (length(diff_cols) > 0){
    warning(paste("The following columns are missing from maven_output_df :", paste(diff_cols, collapse = ", "), sep = " "))
    return (NULL)
  }   
  
  if (!(quant_type %in% colnames(maven_output_df))){
    warning(c("The ", quant_type, " is not present in maven_output_df column"))
    return (NULL)
  }    
  
  default_cols_vec <- c("label", "metaGroupId", "groupId", "goodPeakCount", "medMz", "medRt", 
                        "maxQuality", "adductName", "isotopeLabel", "compound", "compoundId", 
                        "formula", "expectedRtDiff", "ppmDiff", "parent")
  
  empty_cols_vec <- c("label", "goodPeakCount", "maxQuality", "adductName", "isotopeLabel", 
                      "expectedRtDiff", "ppmDiff", "parent")
  
  
  
  summarize_fun_list <- list(compound = ~unique(compound)[[1]], compoundId = ~unique(compoundId)[[1]], 
                             formula = ~unique(formula)[[1]], medMz = ~round(mean(peakMz), 6), 
                             medRt = ~round(mean(rt), 6))
  
  peakml_cols <- c('peakML_label', 'peakML_probability')  
  if (all(peakml_cols %in% colnames(maven_output_df))){
    default_cols_vec <- c(default_cols_vec, peakml_cols)
    summarize_fun_list <- c(summarize_fun_list, list(peakML_label = ~unique(peakML_label)[[1]], 
                                                     peakML_probability = ~unique(peakML_probability)[[1]]))
  }    
  
  samples_vec <-  unique(maven_output_df[, "sample"])  
  summarize_peak_detailed <- dplyr::summarize_at(dplyr::group_by_at(maven_output_df, "groupId"), "sample", 
                                                 summarize_fun_list)    
  
  interm_empty_df <- data.frame(matrix(ncol = length(empty_cols_vec), nrow = 1), stringsAsFactors = FALSE, check.names = FALSE)
  interm_empty_df[1, ] <- ""
  colnames(interm_empty_df) <- empty_cols_vec    
  
  summarize_peak_detailed <- data.frame(summarize_peak_detailed, metaGroupId = summarize_peak_detailed$groupId, 
                                        interm_empty_df, stringsAsFactors = FALSE, check.names = FALSE)
  peaks_to_wide_df <- reshape::cast(maven_output_df, groupId ~ sample, value = quant_type, fill = 0)
  
  group_summary_df <- merge(summarize_peak_detailed, peaks_to_wide_df, by = "groupId")
  group_summary_df <- group_summary_df[, c(default_cols_vec, samples_vec)]    
  
  message("Make Group Summary From Peak Detailed Completed...")
  
  return(group_summary_df)
}
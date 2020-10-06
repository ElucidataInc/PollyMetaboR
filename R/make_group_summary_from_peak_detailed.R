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
  require(dplyr)
  
  if (identical(maven_output_df, NULL)){
    warning("maven_output_df is NULL dataframe")
    return (NULL)
  }
  
  if (nrow(maven_output_df) < 1){
    warning("Not a valid dataframe, requires atleast one row for processing")
    return (NULL)
  }
  
  if (!(quant_type %in% colnames(maven_output_df))){
    warning(c("The ", quant_type, " is not present in maven_output_df column"))
    return (NULL)
  }
  
  group_summary_df <- data.frame(check.names = FALSE, stringsAsFactors = FALSE)
  for (groupId_val in unique(maven_output_df$groupId)){
    feature_df <- dplyr::filter(maven_output_df, groupId == groupId_val)
    compound <- as.character(unique(feature_df$compound))
    compoundId <- as.character(unique(feature_df$compoundId))
    formula <- as.character(unique(feature_df$formula))
    samples_keys <- as.character(feature_df$sample)
    samples_values <- as.numeric(feature_df[[quant_type]])
    medMz <- mean(feature_df$peakMz)
    medRt <- round(mean(feature_df$rt),3)
    default_cols_vec <- c("label" = "", "metaGroupId" = groupId_val, "groupId"= groupId_val, 
                          "goodPeakCount" = "", "medMz" = medMz, "medRt" = medRt, 
                          "maxQuality" = "", "adductName" = "", "isotopeLabel" = "", 
                          "compound" = compound, "compoundId" = compoundId, 
                          "formula" = formula, "expectedRtDiff" = "", "ppmDiff" = "",
                          "parent" = "")
    samples_vec <- samples_values
    names(samples_vec) <- make.unique(samples_keys)
    group_summary_cols <- c(default_cols_vec, samples_vec)
    interm_df <- as.data.frame(matrix(,ncol=length(group_summary_cols), nrow=0))
    names(interm_df)<-names(group_summary_cols)
    interm_df[1,names(group_summary_cols)] <- group_summary_cols
    group_summary_df <- dplyr::bind_rows(group_summary_df, interm_df)
  }
  
  numeric_cols <- c("metaGroupId", "groupId", "medMz", "medRt", as.character(unique(maven_output_df$sample)))
  group_summary_df[numeric_cols] <- sapply(group_summary_df[numeric_cols],as.numeric)
  
  message("Make Group Summary From Peak Detailed Completed...")
  
  return(group_summary_df)
}
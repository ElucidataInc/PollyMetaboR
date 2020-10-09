#' make_group_summary_from_metab_ident_format
#'
#' It convert metabolite identification step format to group summary format of elmaven
#'
#' @param metab_identified_df The output of metabolite identification step
#' @return The group summary format of elmaven
#' @examples
#' make_group_summary_from_metab_ident_format(metab_identified_df)
#' @import dplyr
#' @export
make_group_summary_from_metab_ident_format <- function(metab_identified_df = NULL){
  message("Make Group Summary From Metab Ident Format Started...")
  require(dplyr)
  
  if (identical(metab_identified_df, NULL)){
    warning("The metab_identified_df is NULL")
    return (NULL)
  }   
  
  if ((!class(metab_identified_df)=='data.frame')){
    warning("The metab_identified_df is not a dataframe")
    return (NULL)
  }
  
  default_cols  <- c('mz', 'mzmin', 'mzmax', 'rt', 'rtmin', 'rtmax', 'groupId', 
                     'isotopes', 'adduct', 'pcgroup', 'adduct_type', 'basemass',
                     'isotope_id', 'isotope_type', 'feature_group', 'compound', 
                     'id', 'formula', 'mass')
  utp_output_cols <- colnames(metab_identified_df)
  
  if (!all(default_cols %in% utp_output_cols)){
    warning(c("The metab_identified_df should have the following columns :", paste0(default_cols, collapse = ", ")))
    return (NULL)  
  }
  
  samples_vec <- utp_output_cols[!utp_output_cols %in% default_cols]
  subset_cols <- c("groupId", "mz", "rt", "compound",  "id", "formula", samples_vec)
  maven_output_format <- metab_identified_df[, subset_cols, drop = FALSE]
  maven_output_format  <- maven_output_format %>% dplyr::rename_at(vars(c("mz","rt", "id")), ~ c("medMz", "medRt", "compoundId"))
  maven_output_format  <- maven_output_format %>% dplyr::mutate(label = "", goodPeakCount= "",maxQuality="", adductName="", isotopeLabel="",expectedRtDiff="",ppmDiff="",parent="")
  maven_output_format$metaGroupId <- maven_output_format$groupId
  maven_default_cols <- c("label", "metaGroupId", "groupId", "goodPeakCount", "medMz",
                          "medRt", "maxQuality", "adductName", "isotopeLabel", "compound",
                          "compoundId", "formula", "expectedRtDiff", "ppmDiff", "parent")
  maven_cols_order <- c(maven_default_cols, samples_vec)
  maven_output_format <- maven_output_format[maven_cols_order]
  cmpd_nan_bool <- is.na(maven_output_format$compound)
  identified_df <- maven_output_format[!cmpd_nan_bool,]
  unidentified_df <- maven_output_format[cmpd_nan_bool,]
  unidentified_df$compound <- paste(unidentified_df$medMz, unidentified_df$medRt, sep = "@")
  unidentified_df$compoundId <- paste(unidentified_df$medMz, unidentified_df$medRt, sep = "@")
  unidentified_df$formula <- ""
  maven_output_format <- rbind(identified_df, unidentified_df)
  
  # Replace NA's in samples with zero
  maven_output_format[samples_vec][is.na(maven_output_format[samples_vec])] <- 0
  
  message("Make Group Summary From Metab Ident Format Completed...")
  
  return (maven_output_format)
}
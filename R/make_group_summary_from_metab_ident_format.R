#' make_group_summary_from_metab_ident_format
#'
#' It converts metabolites identification step or restructured camera output format to group summary format of elmaven
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
  
  metab_default_cols  <- c('mz', 'mzmin', 'mzmax', 'rt', 'rtmin', 'rtmax', 'groupId',
                           'identified', 'isotopes', 'adduct', 'pcgroup', 'adduct_type',
                           'basemass', 'isotope_id', 'isotope_type', 'feature_group',
                           'compound', 'id', 'formula', 'mass', 'rt_db', 'identification_score')
  
  maven_default_cols <- c("label", "metaGroupId", "groupId", "goodPeakCount", "medMz",
                          "medRt", "maxQuality", "adductName", "isotopeLabel", "compound",
                          "compoundId", "formula", "expectedRtDiff", "ppmDiff", "parent")
  
  subset_cols <- c("groupId", "mz", "rt", "compound",  "id", "formula", "feature_group",
                   "adduct_type", "isotope_type", "basemass")    
  
  utp_output_cols <- colnames(metab_identified_df)
  if (!all(metab_default_cols %in% utp_output_cols)){
    if (!all(metab_default_cols[1:16] %in% utp_output_cols)){
      warning(c("The metab_identified_df should have at least the following columns : ",
                paste0(metab_default_cols[1:16], collapse = ", "),
                "\nThe following columns are optional which are from metabolite identification step : ",
                paste0(metab_default_cols[17:20], collapse = ", ")))
      return (NULL)
    }
    
    if (!("compound" %in% utp_output_cols)){
      metab_identified_df$compound <- NA
    }
    
    if (!("id" %in% utp_output_cols)){
      metab_identified_df$id <- NA
    }
    
    if (!("formula" %in% utp_output_cols)){
      metab_identified_df$formula <- NA
    }      
  }
  
  if ("identification_score" %in% utp_output_cols){
    maven_default_cols <- c(maven_default_cols, "identification_score")
    subset_cols <- c(subset_cols, "identification_score")
  }
  
  samples_vec <- utp_output_cols[!utp_output_cols %in% metab_default_cols]
  filter_cols <- c(subset_cols, samples_vec)
  maven_output_format <- metab_identified_df[, filter_cols, drop = FALSE]
  maven_output_format  <- dplyr::rename_at(maven_output_format, dplyr::vars(c("mz","rt", "id")), ~ c("medMz", "medRt", "compoundId"))
  maven_output_format  <- dplyr::mutate(maven_output_format, label = "", goodPeakCount= "",
                                        maxQuality="", adductName="", isotopeLabel="",
                                        expectedRtDiff="", ppmDiff="", parent="")
  maven_output_format$metaGroupId <- maven_output_format$feature_group
  maven_output_format$adductName <- maven_output_format$adduct_type
  maven_output_format$isotopeLabel <- maven_output_format$isotope_type
  maven_output_format$parent <- maven_output_format$basemass    
  
  maven_cols_order <- c(maven_default_cols, samples_vec)
  maven_output_format <- maven_output_format[, maven_cols_order, drop = FALSE]
  cmpd_nan_bool <- is.na(maven_output_format$compound)
  identified_df <- maven_output_format[!cmpd_nan_bool, , drop = FALSE]
  unidentified_df <- maven_output_format[cmpd_nan_bool, , drop = FALSE]
  unidentified_df$compound <- paste(unidentified_df$medMz, unidentified_df$medRt, sep = "@")
  unidentified_df$compoundId <- paste(unidentified_df$medMz, unidentified_df$medRt, sep = "@")
  unidentified_df$formula <- ""
  maven_output_format <- rbind(identified_df, unidentified_df)
  
  # Replace NA's in samples with zero
  maven_output_format[samples_vec][is.na(maven_output_format[samples_vec])] <- 0
  maven_output_format[maven_default_cols][is.na(maven_output_format[maven_default_cols])] <- ""
  
  message("Make Group Summary From Metab Ident Format Completed...")
  
  return (maven_output_format)
}
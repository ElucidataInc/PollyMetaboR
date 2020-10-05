#' merge_identified_with_restructured_camera
#'
#' It merges the identified metabolites data with restructured camera output data.
#'
#' @param metab_identified_df The identified metabolites dataframe
#' @param restructured_camera_df The restructured camera output
#' @return The dataframe with merged identified and restructured camera data
#' @examples
#' merge_identified_with_restructured_camera(metab_identified_df, restructured_camera_df)
#' @import dplyr
#' @export
merge_identified_with_restructured_camera <- function(metab_identified_df = NULL, restructured_camera_df = NULL) {
  message("Merge Identified With Restructured Camera Started...")
  require(dplyr)
  
  if (identical(metab_identified_df, NULL)){
    warning("The metab_identified_df is NULL")
    return (NULL)
  }   
  
  if ((!class(metab_identified_df)=='data.frame')){
    warning("The metab_identified_df is not a dataframe")
    return (NULL)
  }
  
  if (identical(restructured_camera_df, NULL)){
    warning("The restructured_camera_df is NULL")
    return (NULL)
  }
  
  if ((!class(restructured_camera_df)=='data.frame')){
    warning("The restructured_camera_df is not a dataframe")
    return (NULL)
  }
  
  common_cols <- c('mz', 'rt', 'basemass', 'groupId')
  
  if (!all(common_cols %in% colnames(metab_identified_df))){
    warning(c("The following columns should be present in metab_identified_df: ", paste0(common_cols, collapse = ", ")))
    return (NULL)
  }
  
  if (!all(common_cols %in% colnames(restructured_camera_df))){
    warning(c("The following columns should be present in restructured_camera_df: ", paste0(common_cols, collapse = ", ")))
    return (NULL)
  }
  
  camera_results_cols <- colnames(restructured_camera_df)
  drop_cols <- camera_results_cols[!(camera_results_cols %in% common_cols)]
  metab_identified_df <- metab_identified_df[,!(colnames(metab_identified_df) %in% drop_cols)]
  merged_ident_annot_df <- dplyr::full_join(restructured_camera_df, metab_identified_df, by = common_cols)
  
  message("Merge Identified With Restructured Camera Completed...")
  
  return(merged_ident_annot_df)
}
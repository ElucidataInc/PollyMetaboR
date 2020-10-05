#' get_feature_group_representative
#'
#' It returns a dataframe with only representative features from each feature groups.
#'
#' @param restructured_camera_df The restructured camera output
#' @param polarity The polarity of the data
#' @return The dataframe with only representative features
#' @examples
#' get_feature_group_representative(restructured_camera_df, polarity = "positive")
#' @import dplyr
#' @export
get_feature_group_representative <- function(restructured_camera_df = NULL, polarity = "positive"){
  message("Get Feature Group Representative Started...")
  require(dplyr)
  
  if (identical(restructured_camera_df, NULL)){
    warning("The restructured_camera_df is NULL")
    return (NULL)
  }
  
  if ((!class(restructured_camera_df)=='data.frame')){
    warning("The restructured_camera_df is not a dataframe")
    return (NULL)
  }
  
  final_grouped_camera_result_cols <- colnames(restructured_camera_df)
  identifier_cols <- c('mz', 'mzmin', 'mzmax', 'rt', 'rtmin', 'rtmax', 'groupId', 
                       'isotopes', 'adduct', 'pcgroup', 'adduct_type', 'basemass',
                       'isotope_id', 'isotope_type', 'feature_group')
  
  all_samples_vec <- final_grouped_camera_result_cols[!(final_grouped_camera_result_cols %in% identifier_cols)]
  
  feature_group_representative_df <- data.frame()
  
  for (ftr_grp in sort(unique(restructured_camera_df$feature_group))){
    feature_group_df <- dplyr::filter(restructured_camera_df, feature_group == ftr_grp)
    isotopes_vec <- feature_group_df$isotope_type[!is.na(feature_group_df$isotope_type)]
    drop_isotopes_vec <- isotopes_vec[!(isotopes_vec %in% 'M')]
    filtered_feature_group_df <- dplyr::filter(feature_group_df, !isotope_type %in% drop_isotopes_vec)
    
    if (polarity == "positive"){
      M_H_bool_check <- any((filtered_feature_group_df$adduct_type == '[M+H]+') == TRUE)
      if (M_H_bool_check == TRUE){
        filtered_feature_group_df <- dplyr::filter(filtered_feature_group_df, adduct_type == '[M+H]+')
      }
    } else {
      M_H_bool_check <- any((filtered_feature_group_df$adduct_type == '[M-H]-') == TRUE)
      if (M_H_bool_check == TRUE){
        filtered_feature_group_df <- dplyr::filter(filtered_feature_group_df, adduct_type == '[M-H]-')
      }
    }
    
    filtered_feature_group_samples_df <- filtered_feature_group_df[all_samples_vec]
    row_sum_vec <- (rowSums(filtered_feature_group_samples_df, na.rm=TRUE))
    max_value_index_vec <- which(row_sum_vec == max(row_sum_vec))
    interm_df <-  filtered_feature_group_df[max_value_index_vec[1],]
    feature_group_representative_df <- rbind(feature_group_representative_df, interm_df)
    
  }
  
  order_feature_group_representative_df <- feature_group_representative_df
  
  message("Get Feature Group Representative Completed...")
  
  return(order_feature_group_representative_df)
}
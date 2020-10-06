#' separate_identified_from_unidentified
#'
#' It separates identified features from unidentified features in elmaven
#'
#' @param maven_output_df The peak detailed format of maven output
#' @return A list of two dataframes having identified and unidentified features information
#' @examples
#' separate_identified_from_unidentified(maven_output_df)
#' @export
separate_identified_from_unidentified <- function(maven_output_df = NULL){
  message("Separate Identified From Unidentified Started...")
  
  if (identical(maven_output_df, NULL)){
    warning("The maven_output_df is NULL")
    return (NULL)  
  }
  
  if ((!class(maven_output_df)=='data.frame')){
    warning("The maven_output_df is not a dataframe")
    return (NULL)
  }    
  
  if (!("compound" %in% colnames(maven_output_df))){
    warning("The compound column is not present in maven_output_df")
    return (NULL)  
  }
  
  mz_rt_bool <- grepl('[0-9]+@[0-9]+', maven_output_df$compound)
  untar_data <- maven_output_df[mz_rt_bool, , drop = FALSE]
  
  if (nrow(untar_data) >= 1){
    rownames(untar_data) <- 1:nrow(untar_data)
  } 
  else {
    untar_data <- NULL
  }  
  
  tar_data <- maven_output_df[!mz_rt_bool, , drop = FALSE]
  if (nrow(tar_data) >= 1){
    rownames(tar_data) <- 1:length(rownames(tar_data))
  }
  else {
    tar_data <- NULL
  }  
  
  message("Separate Identified From Unidentified Completed...")
  
  return (list(unidentified = untar_data, identified = tar_data))
}
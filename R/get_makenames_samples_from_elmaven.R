#' get_makenames_samples_from_elmaven
#'
#' It returns a named vector of samplenames where make.names are assigned to samples.
#'
#' @param maven_output_df The peak detailed format of elmaven
#' @return The named vector of samplenames using make.names
#' @examples
#' get_makenames_samples_from_elmaven(maven_output_df)
#' @export
get_makenames_samples_from_elmaven <- function(maven_output_df = NULL){ 
  message("Get Makenames Samples From Elmaven Started...")
  
  if (identical(maven_output_df, NULL)){
    warning("The maven_output_df is NULL")
  }
  
  if ((!class(maven_output_df)=='data.frame')){
    warning("The maven_output_df is not a dataframe")
    return (NULL)
  }    
  
  if (!('sample' %in% colnames(maven_output_df))){
    warning("The sample column is not present in maven_output_df dataframe")
    return (NULL)
  }
  
  maven_samples <- as.vector(unique(maven_output_df$sample))
  names(maven_samples) <- make.names(maven_samples, unique = TRUE)
  
  message("Get Makenames Samples From Elmaven Completed...")
  
  return (maven_samples)
}
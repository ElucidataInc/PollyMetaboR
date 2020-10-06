#' make_valid_elmaven_samplenames_by_makenames
#'
#' It returns a dataframe with only representative features from each feature groups.
#'
#' @param maven_output_df The peak detailed format of elmaven
#' @param metadata_df The sample to cohort mapping file where samples should be present in first column
#' @return The list of two dataframes with replaced samples names by make.names
#' @examples
#' make_valid_elmaven_samplenames_by_makenames(maven_output_df, metadata_df)
#' @import dplyr
#' @export
make_valid_elmaven_samplenames_by_makenames <- function(maven_output_df = NULL, metadata_df = NULL){
  message("Make Valid Elmaven Samplenames by Makenames Started...")
  require(dplyr)
  
  if (identical(maven_output_df, NULL) && identical(metadata_df, NULL)){
    warning("Both maven_output_df and metadata_df are NULL")
  }
  
  ###Make Maven Samples Valid Names###
  if(!identical(maven_output_df, NULL)){
    
    if ((!class(maven_output_df)=='data.frame')){
      warning("The maven_output_df is not a dataframe")
      return (NULL)
    }
    
    if (!("sample" %in% colnames(maven_output_df))){
      warning("The sample column is not present in maven_output_df dataframe")
      return (NULL)   
    }
    
    maven_samples <- as.vector(unique(maven_output_df$sample))
    maven_samples_valid_names <- make.names(maven_samples, unique = TRUE)
    for (maven_sample_index in 1:length(maven_samples)){
      maven_sample_name <- as.character(maven_samples[maven_sample_index])
      maven_sample_valid_name <- as.character(maven_samples_valid_names[maven_sample_index])
      maven_output_df <- maven_output_df %>% dplyr::mutate(sample = replace(as.character(sample), sample == maven_sample_name, maven_sample_valid_name))
    }  
  }
  
  ###Make Metadata Samples Valid Names###
  if(!identical(metadata_df, NULL)){
    
    if ((!class(metadata_df)=='data.frame')){
      warning("The metadata_df is not a dataframe")
      return (NULL)
    }
    
    sample_col <- colnames(metadata_df)[1]
    metadata_samples <- as.vector(unique(metadata_df[[sample_col]]))
    metadata_samples_valid_names <- make.names(metadata_samples, unique = TRUE)
    for (metadata_sample_index in 1:length(metadata_samples)){
      metadata_sample_name <- as.character(metadata_samples[metadata_sample_index])
      metadata_sample_valid_name <- as.character(metadata_samples_valid_names[metadata_sample_index])
      metadata_df <- metadata_df %>% dplyr::mutate(!! sym(sample_col) := replace(as.character(!! sym(sample_col)), !! sym(sample_col) == metadata_sample_name, metadata_sample_valid_name))
    }
  }
  
  message("Make Valid Elmaven Samplenames by Makenames Completed...")
  
  return(list(maven_data = maven_output_df, metadata = metadata_df))
}
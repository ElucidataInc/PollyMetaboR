#' replace_names_by_names_vector
#'
#' It returns a dataframe after replacing names using names vector.
#'
#' @param data_df The data where names to be replaced
#' @param names_vector The named vector used for replacing names 
#' @param replace_in the column where names are present (it can be a column or "by_columns")
#' @return The named vector of samplenames using make.names
#' @examples
#' replace_names_by_names_vector(data_df, names_vector, replace_in = "sample")
#' @import dplyr
#' @export
replace_names_by_names_vector <- function(data_df = NULL, names_vector = NULL, replace_in = NULL){
  message("Replace Names By Names Vector Started...")
  require(dplyr)
  
  if (identical(data_df, NULL)){
    warning("The data_df is NULL")
  }
  
  if ((!class(data_df)=='data.frame')){
    warning("The data_df is not a dataframe")
    return (NULL)
  }
  
  if (identical(names_vector, NULL)){
    warning("The names_vector is NULL")
    return (NULL)
  }
  
  if (identical(names(names_vector), NULL)){
    warning("The names are not defined for names_vector")
    return (NULL)
  }
  
  if (replace_in == 'by_columns'){
    data_cols <- colnames(data_df)
    common_names <- intersect(names(names_vector), data_cols)
    
    if (length(common_names) == 0){
      message("No common names to replace")
      return (data_df)
    }
    
    common_names_vector <- names_vector[common_names]
    data_df <- data_df %>% dplyr::rename_at(vars(names(names_vector)), ~ unname(names_vector))
  }
  else {  
    if (!(replace_in %in% colnames(data_df))){
      warning(c(replace_in, " not in found data_df columns"))
      return (NULL)
    }
    
    data_names <- unique(data_df[[replace_in]])
    common_names <- intersect(names(names_vector), data_names)
    
    if (length(common_names) == 0){
      message("No common names to replace")
      return (NULL)
    }  
    
    for (names_index in 1:length(names_vector)){
      replace_by_name <- names(names_vector)[names_index]
      replace_with_name <- unname(names_vector)[names_index]
      data_df <- data_df %>% dplyr::mutate(!! sym(replace_in) := replace(as.character(!! sym(replace_in)), !! sym(replace_in) == replace_by_name, replace_with_name))
    } 
  }
  
  message("Replace Names By Names Vector Completed...")
  
  return(data_df)  
}
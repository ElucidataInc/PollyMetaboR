#' smartbind_dataframes
#'
#' It combines two dataframe having same columns
#'
#' @param df_1 The dataframe 1
#' @param df_2 The dataframe 2
#' @return The combined dataframe
#' @examples
#' smartbind_dataframes(df_1, df_2)
#' @export
smartbind_dataframes <- function(df_1 = NULL, df_2 = NULL){
      
  message("Smartbind Dataframe Started...")

  if (identical(df_1, NULL) && identical(df_2, NULL)){
    warning("Both dataframes are NULL\n")
    combined_df <- NULL
  } 
  else if ((class(df_1) == 'data.frame') && (class(df_2) == 'data.frame')){
    combined_df <- gtools::smartbind(df_1, df_2)
    rownames(combined_df) <- 1:length(rownames(combined_df))
  }
  else if ((class(df_1) == 'data.frame') && (class(df_2) != 'data.frame')){
    combined_df <- df_1
  }
  else if ((class(df_1) != 'data.frame') && (class(df_2) == 'data.frame')){
    combined_df <- df_2
  }
  else{
    warning("Please provide vaild dataframes\n")
    combined_df <- NULL
  }

  message("Smartbind Dataframe Completed...")
  
  return (combined_df)    
}
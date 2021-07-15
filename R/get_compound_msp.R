#' get_compound_msp
#'
#' It converts the string text of a single compound msp into the list of msp data for a compound
#'
#' @param element_msp A string of a single compound msp
#' @return A list of a single compound msp
#' @examples
#' get_compound_msp(element_msp)
#' @export
get_compound_msp <- function(element_msp = NULL){
  message("Get Compound MSP Started...")  
  require(stringi)  
  require(stringr)
  
  if (identical(element_msp, NULL)){
    warning("The element_msp is NULL")
    return (NULL)
  }
  
  if (!identical(class(element_msp), "character")){
    warning("The element_msp is not a character")
    return (NULL)
  }
  
  compound_msp <- list()
  element_msp <- stringi::stri_split_lines(element_msp)[[1]]
  names_bool <- sapply(element_msp, function(x) grepl(": ", x, fixed = TRUE))
  names_elements <- element_msp[names_bool]
  if (length(names_elements) > 0){
    names_split_list <- lapply(names_elements, function(x){
      name_list <- list()
      name_split <- stringr::str_trim(stringr::str_split(x, ": ", 2)[[1]])
      num_val <- suppressWarnings(as.numeric(name_split[2]))
      if (!is.na(num_val)){
        name_list[[name_split[1]]] <- num_val
      }
      else {
        name_list[[name_split[1]]] <- name_split[2]
      }
      return (name_list) 
    })
    if (length(names_split_list) > 0){
      for (name_element in names_split_list){
        compound_msp <- c(compound_msp, name_element)
      }
    }
  }
  
  peaks_elements <- element_msp[!names_bool]
  if (length(peaks_elements) > 0){                        
    peaks_split_list <- lapply(peaks_elements, function(x) {
      peak_list <- list()
      peak_split <- stringr::str_trim(base::strsplit(x, "(\\s+|\\t+)(?=(?:[^\"]*\"[^\"]*\")*(?![^\"]*\"))", perl = TRUE)[[1]])
      return (peak_split)
    })
    peaks_split_list <- peaks_split_list[lengths(peaks_split_list) > 0 & peaks_split_list != ""]
    if (length(peaks_split_list) > 0){                      
      if (length(peaks_split_list) > 1){
        peaks_df <- as.data.frame(stringi::stri_list2matrix(peaks_split_list, byrow=TRUE), stringsAsFactors = FALSE)
      }
      else {
        peaks_df <- as.data.frame(t(data.frame(peaks_split_list)), stringsAsFactors = FALSE)
      }
      elements_len <- ncol(peaks_df)
      if (elements_len > 1){
        extra_cols_len <- elements_len - 2
        if (extra_cols_len > 1){ extra_cols <- paste("other_info", c(1:extra_cols_len), sep = "_")}
        else if (extra_cols_len == 1){ extra_cols <- "other_info"}
        else {extra_cols <- NULL}
        spectrum_cols <- c("mz", "intensity")
        colnames(peaks_df)<- c(spectrum_cols, extra_cols)
        peaks_df[, spectrum_cols] <- suppressWarnings(sapply(peaks_df[, spectrum_cols], as.numeric))
        row.names(peaks_df) <- NULL   
        compound_msp[["pspectrum"]] <- peaks_df
      }
    }  
  }
  
  message("Get Compound MSP Completed...") 
  
  return (compound_msp)            
}
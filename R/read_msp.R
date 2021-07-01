#' read_msp
#'
#' It reads the msp format.
#'
#' @param msp_path The path to msp file
#' @param compound A vector of features/compounds to search in the msp file which is optional
#' @param exact_match Search exact match (TRUE) or do partial matching (FALSE) of compound names in the msp file
#' @return A list of msp data
#' @examples
#' read_msp(msp_path)
#' @export
read_msp <- function(msp_path = NULL, compound = NULL, exact_match = TRUE){
  message("Read MSP Started...")
  require(readr)  
  require(stringi)  
  require(stringr)  
  
  if (identical(msp_path, NULL)){
    warning("The msp_path is NULL")
    return (NULL)
  }
  
  if (!file.exists(msp_path)){
    warning(paste("The", sQuote(msp_path), "is not a valid path", sep = " "))
    return (NULL)
  }
  
  if (!identical(tools::file_ext(msp_path), "msp")){
    warning(paste("The", sQuote(msp_path), "does not end with '.msp'", sep = " "))
    return (NULL)      
  }
  
  if (!identical(compound, NULL)){
    if (identical(exact_match, NULL)){ exact_match <- TRUE}  
  }
  
  id_list <- list()
  msp_list <- list()
  
  data_as_str <- readr::read_file(msp_path)
  compound_data_list <- strsplit(gsub("Name: ","~Name: ", data_as_str, ignore.case = T), split = "~")[[1]]
  if (!identical(compound, NULL)){
    compound_data_list <- compound_data_list[sapply(compound_data_list, function(x) any(sapply(compound, function(y) any(grepl(y, x, fixed = TRUE)))))]
  }
  compound_data_list <- compound_data_list[lengths(compound_data_list) > 0 & compound_data_list != ""] 
  
  get_compound_msp <- function(element_msp = NULL){
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
    
    return (compound_msp)             
  }                                                                                           
  
  name_index <- 0                                                                                       
  for (compound_data in compound_data_list){
    compound_msp <- get_compound_msp(compound_data)
    if (length(compound_msp) > 0){
      if ("Name" %in% names(compound_msp)){
        if (identical(compound, NULL)){
          name_index <- name_index + 1
          msp_list[[name_index]] <- compound_msp
          id_list[[name_index]] <- compound_msp$Name 
        }
        else {
          
          if (identical(exact_match, FALSE)){ 
            name_index <- name_index + 1
            msp_list[[name_index]] <- compound_msp
            id_list[[name_index]] <- compound_msp$Name  
          }
          else {
            if (any(compound_msp$Name %in% compound)){
              name_index <- name_index + 1
              msp_list[[name_index]] <- compound_msp
              id_list[[name_index]] <- compound_msp$Name  
            }
          }
        }
      }
    }
  }
  
  message("Read MSP Completed...")                                                                                      
  
  return (list(compound = id_list, msp = msp_list))                                                                                           
}
#' get_elements_from_msp
#'
#' It returns the list of elements found in the msp file using the compound key
#'
#' @param msp_path The path to msp file
#' @param compound_key A key name from msp file where features/compounds are present
#' @param load_on_memory Load whole text into a list if TRUE else read lines one by one to save memory
#' @return A vector of features/compounds
#' @examples
#' get_elements_from_msp(msp_path, compound_key)
#' @export
get_elements_from_msp <- function(msp_path = NULL, compound_key = NULL, load_on_memory = TRUE){
  message("Get Elements From MSP Started...")
  require(readr)
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
  
  if (identical(compound_key, NULL)){
    warning("The compound_key is NULL")
    return (NULL)
  }
  
  load_on_memory <- as.logical(load_on_memory)
  if (is.na(load_on_memory)){
    warning("The load_on_memory should be TRUE/FALSE ")
    return (NULL)
  } 
  
  if (identical(load_on_memory, TRUE)){
    data_as_readlines <- readr::read_lines(msp_path)
    compound_key_list <- data_as_readlines[sapply(data_as_readlines, function(x) grepl(pattern = paste(compound_key, ":", sep = ""), x = x, ignore.case = TRUE, fixed = FALSE))]                                              
    compound_values  = as.character(sapply(compound_key_list, function(x) stringr::str_trim(stringr::str_split(x, ":", 2)[[1]][2])))
  } else {
    compound_values <- c()
    file_con  <- file(msp_path, open = "r")
    while (length(one_line <- readLines(file_con, n = 1, warn = FALSE)) > 0){
      if(grepl(pattern = paste(compound_key, ":", sep = ""), x = one_line, ignore.case = TRUE, fixed = FALSE)){
        element_val <- stringr::str_trim(stringr::str_split(one_line, ":", 2)[[1]])[2]
        compound_values <- c(compound_values, element_val)                                 
      }
    }
    close(file_con)  
  }    
  
  message("Get Elements From MSP Completed...")
  
  return (compound_values) 
}
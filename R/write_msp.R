#' write_msp
#'
#' It writes the compounds msp data into a text file with .msp extention
#'
#' @param msp_data A list of compounds msp data which is the output of read_msp
#' @param save_path File path to save the data as .msp format
#' @return The file saved as .msp format
#' @examples
#' write_msp(msp_data, save_path)
#' @export
write_msp <- function(msp_data = NULL, save_path = NULL){
  message("Write MSP Started...")
  require(readr)
  
  if (identical(msp_data, NULL)){
    warning("The msp_data is NULL")
    return (NULL)
  }
  
  if (!all(c("compound", "msp") %in% names(msp_data))){
    warning("The msp_data should contain compound and msp names")
    return (NULL)
  }
  
  if (length(msp_data$msp) < 1){
    warning("The msp_data should have at least one element in msp")
    return (NULL) 
  }
  
  if (identical(save_path, NULL)){
    warning("The save_path is NULL, using msp_data.msp as default name")
    save_path <- "msp_data.msp"
  }
  
  dir_path <- dirname(save_path)  
  if (!dir.exists(dir_path)){
    warning(paste0("The ", sQuote(dir_path), " does not exist"))
    return (NULL)   
  }
  
  if (!identical(tools::file_ext(save_path), "msp")){
    save_path <- paste0(save_path, ".msp")    
  }
  
  correct_path <- TRUE  
  tryCatch(readr::write_file("", save_path), error = function(e) { correct_path <<- FALSE})
  if (!correct_path){ return (NULL)} 
  
  msp_string <- NULL
  for (compound_ind in 1:length(msp_data$msp)){
    compound_msp_str <- suppressMessages(PollyMetaboR::compound_msp_to_string(msp_data$msp[[compound_ind]]))
    if (compound_ind == 1){
      msp_string <- compound_msp_str
    }
    else {
      msp_string <- paste(msp_string, compound_msp_str, sep = "\r\n\r\n")
    }
  }
  
  readr::write_file(msp_string, save_path)
  message(paste0("The msp file is saved at: ", sQuote(save_path)))
  
  message("Write MSP Completed...")
  
}
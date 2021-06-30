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
  require(data.table)
  require(dplyr)
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
  
  file_con  <- file(msp_path, open = "r")
  msp_list <- list()
  name_list <- list()
  id_list <- list()
  name_index <- 0
  name_bool <- FALSE
  num_peaks <- 0
  count_peaks <- 0
  while (length(one_line <- readLines(file_con, n = 1, warn = FALSE)) > 0){
    try({
      if (identical(compound, NULL)){
        if (grepl("Name: ", one_line, fixed = FALSE, ignore.case = TRUE)){
          name_list <- list()
          spectrum_list <- NULL
          name_index <- name_index + 1
          name_bool <- TRUE
          num_peaks <- 0
          count_peaks <- 0
        }
      }
      else {
        compound <- stringr::str_trim(compound)
        if (any(sapply(compound, function(x) grepl(x, one_line, fixed = TRUE))) && grepl("Name: ", one_line, fixed = FALSE, ignore.case = TRUE)){
          if (identical(exact_match, FALSE)){
            name_list <- list()
            spectrum_list <- NULL
            name_index <- name_index + 1
            name_bool <- TRUE
            num_peaks <- 0
            count_peaks <- 0
          }
          else {
            split_line <- stringr::str_trim(stringr::str_split(one_line, ": ", 2)[[1]])
            if (any(split_line[2] %in% compound)){
              name_list <- list()
              spectrum_list <- NULL
              name_index <- name_index + 1
              name_bool <- TRUE
              num_peaks <- 0
              count_peaks <- 0
            }
          }
        }
      }  
      
      if (name_bool){
        if (grepl(": ", one_line, fixed = TRUE)){
          split_line <- stringr::str_trim(stringr::str_split(one_line, ": ", 2)[[1]])
          num_val <- suppressWarnings(as.numeric(split_line[2]))
          if (!is.na(num_val)){
            name_list[[split_line[1]]] <- num_val
            if (grepl(paste(c("Num Peaks", "NumPeaks"), collapse = "|"), split_line[1],  fixed = FALSE, ignore.case = TRUE)){
              num_peaks <- num_val
            }
          }
          else{ name_list[[split_line[1]]] <- split_line[2] }
        }
        else{
          if (num_peaks > 0){
            split_line <- stringr::str_trim(base::strsplit(one_line, "(\\s+|\\t+)(?=(?:[^\"]*\"[^\"]*\")*(?![^\"]*\"))", perl = TRUE)[[1]]) 
            elements_len <- length(split_line)
            if (elements_len > 1){
              extra_cols_len <- elements_len - 2
              if (extra_cols_len > 1){ extra_cols <- paste("other_info", c(1:extra_cols_len), sep = "_")}
              else if (extra_cols_len == 1){ extra_cols <- "other_info"}
              else {extra_cols <- NULL}
              spectrum_cols <- c("mz", "intensity")
              names(split_line) <- c(spectrum_cols, extra_cols)
              interm_df <- as.data.frame(dplyr::bind_rows(split_line), stringsAsFactors = FALSE)
              interm_df[, spectrum_cols] <- sapply(interm_df[, spectrum_cols], as.numeric)           
              spectrum_list <- as.data.frame(data.table::rbindlist(list(spectrum_list, interm_df), fill = TRUE), stringsAsFactors = FALSE)
              name_list[["pspectrum"]] <- spectrum_list
            }
            count_peaks <- count_peaks + 1 
            if (count_peaks >= num_peaks){ name_bool <- FALSE}
          }
        }
        msp_list[[name_index]] <- name_list
        id_list[[name_index]] <- name_list[[1]]
      }
      
    }, silent = TRUE)                     
  }
  close(file_con)
  
  message("Read MSP Completed...")
  
  return (list(compound = id_list, msp = msp_list))
}
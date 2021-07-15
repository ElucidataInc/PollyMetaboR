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
  
  name_index <- 0                                                                                       
  for (compound_data in compound_data_list){
    compound_msp <- suppressMessages(PollyMetaboR::get_compound_msp(compound_data))
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
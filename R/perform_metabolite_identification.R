# library(data.table)
# library(dplyr)
# library(future)
# 
# calculate_delta <- function(mz_source, mz_tolerence_unit = "ppm", mz_tolerence = 20){
#     delta <- NULL
#     if (mz_tolerence_unit == 'ppm'){
#         delta <- mz_source * mz_tolerence * 10^(-6)
#     }
#     else if (mz_tolerence_unit == 'Da'){
#         delta <- mz_tolerence
#     }
#     else {
#         warning("Please select valid mass tolerence unit (ppm or Da)")
#         return()
#     }
# 
#     return(delta)
# }
# 
# mz_search_comp_data <- function(mz_source, comp_data, delta){
#     if (is.null(delta)){
#         warning("Please input valid delta")
#         return()
#     }
#     mzmin <- mz_source - delta
#     mzmax <- mz_source + delta
#     comp_data_filter <- filter(comp_data, mzMass >= mzmin & mzMass <= mzmax)
# 
#     return(comp_data_filter)
# }
# 

#' mz_match_check
#'
#' 
#'
#' @param mz_source mz source value
#' @param mz_target mz target value
#' @param mz_tolerence_unit The mz tolerance unit (ppm or Da)
#' @param mz_tolerence Value of mz tolerence
#' @return Return TRUE if mz source and mz target are within mz tolerence 
#' @examples
#' mz_match_check(mz_source, mz_target, mz_tolerence_unit = "ppm", mz_tolerence = 20)
#' @export
mz_match_check <- function(mz_source, mz_target, mz_tolerence_unit = "ppm", mz_tolerence = 20){
  message("MZ Match Check Started...")
  
  mz_delta_bool_check <- FALSE
  delta <- calculate_delta(mz_source, mz_tolerence_unit = mz_tolerence_unit, mz_tolerence = mz_tolerence)
  if (!is.null(delta)){
    mz_diff <- abs(mz_source - mz_target)
    if (delta >= mz_diff){
      mz_delta_bool_check <- TRUE
    }
  }
  
  message("MZ Match Check Completed...")
  
  return(mz_delta_bool_check)
}

#' perform_metabolite_identification
#'
#' It restructures the camera output to group features
#'
#' @param mz_data A dataframe or numeric vector of mz values
#' @param comp_data The compound database used for identification
#' @param mz_colname A mz column name present in mz_data dataframe
#' @param mz_tolerence_unit The mz tolerance unit (ppm or Da)
#' @param mz_tolerence Value of mz tolerence
#' @param numcores Number of cores used for processing
#' @return A dataframe with identified metabolites
#' @examples
#' perform_metabolite_identification(mz_data, comp_data, mz_colname = 'basemass',
#'                                   mz_tolerence_unit = "ppm", mz_tolerence = 20, numcores = 2)
#' @import data.table dplyr
#' @export
perform_metabolite_identification <- function(mz_data = NULL,  comp_data = NULL, 
                                              mz_colname = 'basemass', mz_tolerence_unit = "ppm",
                                              mz_tolerence = 20, numcores = 2){
  message("Perform Metabolite Identification Started...")
  require(data.table)
  require(dplyr)
  
  if(identical(mz_data, NULL)){
    warning("No mz_data was given")
    return (NULL)
  }
    
  if ((!class(mz_data)=='data.frame') & (!class(mz_data)=='numeric')){
    warning("Please input valid mz_data (a dataframe containing 'mz data' or a vector of 'mz values')")
    return (NULL)
  }

  if(identical(comp_data, NULL)){
    warning("No comp_data was given")
    return (NULL)
  }
  
  if ((!class(comp_data)=='data.frame') | (("mzMass" %in% colnames(comp_data))==FALSE)){
    warning("Please input valid compound database (a dataframe containing 'mzMass' column)")
    return (NULL)
  }
    
  if (class(mz_data)=='numeric'){
    mz_data <- data.frame(mz_index = 1: length(mz_data), mz_colname = mz_data)
    mz_data <- setnames(mz_data, "mz_colname", mz_colname)
  }
    
  if (!(mz_colname %in% colnames(mz_data))){
    warning(c("The ", mz_colname, " is not present in mz_data columns which should have basemass mz values"))
    return (NULL)
  }
    
  if (identical(mz_tolerence_unit, NULL)){
    warning("The mz_tolerence_unit is NULL")
    return (NULL)  
  }    
  
  if (!(mz_tolerence_unit %in% c("ppm", "Da"))){
    warning("The mz_tolerence_unit is not valid, please choose from ppm and Da")
    return (NULL)  
  }
    
  mz_identify_metab <- function(mz_row){
    mz_source <- as.numeric(mz_row[[mz_colname]])
    delta <- calculate_delta(mz_source, mz_tolerence_unit = mz_tolerence_unit, mz_tolerence = mz_tolerence)
    comp_data_filter <- mz_search_comp_data(mz_source, comp_data, delta)
    if (nrow(comp_data_filter) > 0){
      interm_df <- cbind(mz_row, comp_data_filter, row.names = NULL)   
    }
    else {
      comp_data_cols <- colnames(comp_data_filter)
      interm_df <- mz_row
      interm_df[,comp_data_cols] <- NA   
      }
    return(interm_df)
  }
    
  run_split_df <- function(split_mz_data_element){
    metabolite_identified_split_df <- data.frame()
    for (row_index in 1:nrow(split_mz_data_element)){
      mz_row <- split_mz_data_element[row_index,]
      interm_df <- mz_identify_metab(mz_row)
      metabolite_identified_split_df <- dplyr::bind_rows(metabolite_identified_split_df, interm_df)
      rownames(metabolite_identified_split_df) <- NULL
      }
    return(metabolite_identified_split_df)
  }
    
  ###Split mz data into chunks of 1000 features each###
  split_mz_data_list <- split(mz_data, (as.numeric(rownames(mz_data))-1) %/% 1000)
    
  ###Run each future run on only specified core (numcores)###
  cores_split_mz_data_list <- split(split_mz_data_list, ((1:length(split_mz_data_list))-1) %/% numcores)
  
  metabolite_identified_df <- data.frame()
  for (i in names(cores_split_mz_data_list)){
    identify_metab_future <- lapply(cores_split_mz_data_list[[i]],
                                    function(split_mz_data_element) future::future({run_split_df(split_mz_data_element)}))
    metabolite_identified_list <- lapply(identify_metab_future, value) # grab the results
    interm_core_split_future_df <- data.table::rbindlist(metabolite_identified_list, fill=TRUE)
    metabolite_identified_df <- rbind(metabolite_identified_df, interm_core_split_future_df)  
  }
  
  metabolite_identified_df <- as.data.frame(metabolite_identified_df, stringsAsFactors = FALSE)
    
  message("Perform Metabolite Identification Completed...")
    
  return(metabolite_identified_df)
}
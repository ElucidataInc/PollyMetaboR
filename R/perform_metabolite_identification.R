#' calculate_mz_delta
#'
#' It calculates delta of particular mz value
#'
#' @param mz_source mz source numeric value
#' @param mz_tolerence_unit The mz tolerance unit (ppm or Da)
#' @param mz_tolerence Value of mz tolerence
#' @return A numeric value of mz delta 
#' @examples
#' calculate_mz_delta(mz_source, mz_target, mz_tolerence_unit = "ppm", mz_tolerence = 20)
#' @export
calculate_mz_delta <- function(mz_source = NULL, mz_tolerence_unit = "ppm", mz_tolerence = 20){
  message("Calculate MZ Delta Started...")
  
  if (identical(mz_source, NULL)){
    warning("The mz_source is NULL")
    return (NULL)
  }
  
  delta <- NULL
  if (mz_tolerence_unit == 'ppm'){
    delta <- mz_source * mz_tolerence * 10^(-6) 
  }
  else if (mz_tolerence_unit == 'Da'){
    delta <- mz_tolerence
  }
  else {
    warning("Please select valid mass tolerence unit (ppm or Da)")
    return (NULL)
  }
  
  message("Calculate MZ Delta Completed...")
  
  return(delta)
}

#' mz_search_with_comp_data
#'
#' It searches the mz within compound database
#'
#' @param mz_source mz source numeric value
#' @param comp_data reference compound database used for mz seraching with compulsory mzMass column
#' @param mz_tolerence_unit The mz tolerance unit (ppm or Da)
#' @param mz_tolerence Value of mz tolerence
#' @param rt_source rt source numeric value in minutes
#' @param rt_tolerence Value of rt tolerence in minutes
#' @return A dataframe with mapped metabolites 
#' @examples
#' mz_search_with_comp_data(mz_source, comp_data, mz_target, mz_tolerence_unit = "ppm", mz_tolerence = 20)
#' @export
mz_search_with_comp_data <- function(mz_source, comp_data, mz_tolerence_unit = "ppm", 
                                     mz_tolerence = 20, rt_source = NULL, rt_tolerence = NULL){
  message("MZ Search With Comp Data Started...")

  if (identical(mz_source, NULL)){
    warning("The mz_source is NULL")
    return (NULL)
  }  
  
  if ((!class(comp_data)=='data.frame') | (("mass" %in% colnames(comp_data))==FALSE)){
    warning("Please input valid compound database (a dataframe containing 'mass' column)")
    return (NULL)
  }
  
  if (identical(mz_tolerence_unit, NULL)){
    warning("The mz_tolerence_unit is NULL")
    return (NULL)  
  }
  
  if (identical(mz_tolerence, NULL)){
    warning("The mz_tolerence is NULL")
    return (NULL)  
  }
  
  delta <- PollyMetaboR::calculate_mz_delta(mz_source, mz_tolerence_unit = mz_tolerence_unit, mz_tolerence = mz_tolerence)
  
  if (identical(delta, NULL)){
    warning("Please input valid parameters")
    return (NULL)
  }
  
  mzmin <- mz_source - delta
  mzmax <- mz_source + delta
  comp_data_filter <- dplyr::filter(comp_data, mass >= mzmin & mass <= mzmax)
  
  if (nrow(comp_data_filter) >=1){
    comp_data_filter$identification_type <- "mz"
  }   
  
  if ("rt" %in% colnames(comp_data_filter)){
    comp_data_filter <- dplyr::rename(comp_data_filter, rt_db = rt)
    if (!identical(rt_source, NULL) && !identical(rt_tolerence, NULL)){
      message("Matching comp_data with both mz and rt")
      
      comp_data_filter$rt_db <- as.numeric(comp_data_filter$rt_db)
      
      rt_bool <- !is.na(comp_data_filter$rt_db)
      comp_data_with_rt <- comp_data_filter[rt_bool, , drop = FALSE]
      rtmin <- rt_source - rt_tolerence
      rtmax <- rt_source + rt_tolerence
      comp_data_filter_rt <- dplyr::filter(comp_data_with_rt, rt_db >= rtmin & rt_db <= rtmax)
      
      if (nrow(comp_data_filter_rt) >= 1){
        comp_data_filter <- comp_data_filter_rt
        comp_data_filter$identification_type <- "mz and rt"  
      } else {
        comp_data_filter <- comp_data_filter[!rt_bool, , drop = FALSE]
      } 
    } else {
      message("Matching comp_data with only mz")
    }
  }
  
  message("MZ Search With Comp Data Completed...")
  
  return(comp_data_filter)
}

#' mz_match_check
#'
#' It checks if mz_target is within the delta of mz_source
#'
#' @param mz_source mz source numeric value
#' @param mz_target mz target numeric value
#' @param mz_tolerence_unit The mz tolerance unit (ppm or Da)
#' @param mz_tolerence Value of mz tolerence
#' @return Return TRUE if mz source and mz target are within mz tolerence 
#' @examples
#' mz_match_check(mz_source, mz_target, mz_tolerence_unit = "ppm", mz_tolerence = 20)
#' @export
mz_match_check <- function(mz_source = NULL, mz_target = NULL, mz_tolerence_unit = "ppm", mz_tolerence = 20){
  message("MZ Match Check Started...")
  
  if (identical(mz_source, NULL)){
    warning("The mz_source is NULL")
    return (NULL)
  }
  
  if (identical(mz_target, NULL)){
    warning("The mz_target is NULL")
    return (NULL)      
  }
  
  mz_delta_bool_check <- FALSE
  delta <- PollyMetaboR::calculate_mz_delta(mz_source, mz_tolerence_unit = mz_tolerence_unit, mz_tolerence = mz_tolerence)
  
  if (identical(delta, NULL)){
    warning("Please input valid parameters")
    return (NULL)
  }
  
  mz_diff <- abs(mz_source - mz_target)
  if (delta >= mz_diff){
    mz_delta_bool_check <- TRUE
  }
  
  message("MZ Match Check Completed...")
  
  return(mz_delta_bool_check)
}

#' calc_mass_from_formula_comp_data
#'
#' It calculates monoisotopic mass from molecular formula within compound database
#'
#' @param comp_data The compound database having atleast formula column
#' @return The compound database with calculated mass for the compounds
#' @examples
#' calc_mass_from_formula_comp_data(comp_data)
#' @export
calc_mass_from_formula_comp_data <- function(comp_data = NULL){
  message("Calc Mass From Formula Comb Data Started...")
  
  if (identical(comp_data, NULL)){
    warning("The comp_data is NULL")  
  }
  
  if ((!class(comp_data)=='data.frame') | (!any(c("mass", "formula") %in% colnames(comp_data)))){
    warning("Please input valid compound database (a dataframe containing 'mass' or 'formula' column)")
    return (NULL)
  }    
  
  if (!("mass" %in% colnames(comp_data))){
    comp_data$mass <- NA
  } 
  
  if ("formula" %in% colnames(comp_data)){
    get_mass_from_formula <- function(x){
      exact_mass <- as.numeric(x[["mass"]])
      if (is.na(exact_mass)){
        comp_formula <- stringr::str_trim(as.character(x[["formula"]]))
        if (!is.na(comp_formula)){ 
          if (!(comp_formula == "")){
            try(x[["mass"]] <- as.numeric(Rdisop::getMolecule(formula = comp_formula, maxisotopes = 1)$exactmass), silent = TRUE)
          }
        }
      }
      return (x)
    }    
    
    comp_data <- as.data.frame(t(apply(comp_data, 1, function(x) get_mass_from_formula(x))), check.names = FALSE, stringsAsFactors = FALSE) 
    comp_data$mass <- as.numeric(comp_data$mass)
  }
  
  message("Calc Mass From Formula Comb Data Completed...")
  
  return (comp_data)
}

#' perform_metabolite_identification
#'
#' It restructures the camera output to group features
#'
#' @param mz_data A dataframe or numeric vector of mz values
#' @param comp_data The compound database used for identification
#' @param mz_colname The mz column name present in mz_data dataframe
#' @param mz_tolerence_unit The mz tolerance unit (ppm or Da)
#' @param mz_tolerence Value of mz tolerence
#' @param rt_tolerence Value of rt tolerence in minutes
#' @param rt_colname The rt column name present in mz_data dataframe
#' @param wrap_comp Wrap compounds in same row detected for same feature (mz and rt)
#' @param numcores Number of cores used for processing
#' @return A dataframe with identified metabolites
#' @examples
#' perform_metabolite_identification(mz_data, comp_data, mz_colname = 'basemass',
#'                                   mz_tolerence_unit = "ppm", mz_tolerence = 20)
#' @export
perform_metabolite_identification <- function(mz_data = NULL,  comp_data = NULL, mz_colname = 'basemass',
                                              mz_tolerence_unit = "ppm", mz_tolerence = 20, rt_tolerence = NULL,
                                              rt_colname = 'rt', wrap_comp = FALSE, numcores = 2){
  message("Perform Metabolite Identification Started...")
  
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
  
  if ((!class(comp_data)=='data.frame') | (!any(c("mass", "formula") %in% colnames(comp_data)))){
    warning("Please input valid compound database (a dataframe containing 'mass' or 'formula' column)")
    return (NULL)
  }
  
  if (class(mz_data)=='numeric'){
    mz_data <- data.frame(mz_index = 1:length(mz_data), mz_colname = mz_data)
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
  
  if (!("mass" %in% colnames(comp_data))){
    comp_data$mass <- NA
  }
  
  if (!identical(rt_tolerence, NULL)){
    if (!(rt_colname %in% colnames(mz_data))){
      warning(c("The ", rt_colname, " is not present in mz_data columns which should have rt values"))
      return (NULL)
    }      
  }  
  
  comp_data <- PollyMetaboR::calc_mass_from_formula_comp_data(comp_data)
  
  mz_identify_metab <- function(mz_row){
    mz_source <- as.numeric(mz_row[[mz_colname]])
    rt_source <- as.numeric(mz_row[[rt_colname]])
    comp_data_filter <- suppressMessages(PollyMetaboR::mz_search_with_comp_data(mz_source, comp_data, mz_tolerence_unit = mz_tolerence_unit, mz_tolerence = mz_tolerence, rt_source = rt_source, rt_tolerence = rt_tolerence))
    if (nrow(comp_data_filter) > 0){
      if (identical(wrap_comp, TRUE)){
        comp_data_filter_unique <- data.frame(stringsAsFactors = FALSE, check.names = FALSE)
        for (coln in colnames(comp_data_filter)){
          if (length(unique(comp_data_filter[, coln])) > 1){  
            comp_data_filter_unique[1, coln] <- paste(as.character(comp_data_filter[, coln]), collapse = ";")
          }
          else {
            comp_data_filter_unique[1, coln] <- unique(comp_data_filter[, coln])[1]
          }  
        }
        interm_df <- cbind(mz_row, comp_data_filter_unique, row.names = NULL)    
      }
      else {
        interm_df <- cbind(mz_row, comp_data_filter, row.names = NULL)   
      }  
    }
    else {
      comp_data_cols <- colnames(comp_data_filter)
      interm_df <- mz_row
      interm_df[, comp_data_cols] <- NA   
    }
    return(interm_df)
  }

  run_split_df <- function(split_mz_data_element){
    metabolite_identified_split_df <- data.frame()
    for (row_index in 1:nrow(split_mz_data_element)){
      mz_row <- split_mz_data_element[row_index,]
      interm_df <- mz_identify_metab(mz_row)
      metabolite_identified_split_df <- data.table::rbindlist(list(metabolite_identified_split_df, interm_df), fill = TRUE)
      rownames(metabolite_identified_split_df) <- NULL
    }
    return(metabolite_identified_split_df)
  }
  
  ###Split mz data into chunks of 1000 features each###
  split_mz_data_list <- split(mz_data, (as.numeric(1:nrow(mz_data))-1) %/% 1000)
  
  ###Run each future run on only specified core (numcores)###
  cores_split_mz_data_list <- split(split_mz_data_list, ((1:length(split_mz_data_list))-1) %/% numcores)
  
  metabolite_identified_df <- data.frame()
  for (i in names(cores_split_mz_data_list)){
    identify_metab_future <- lapply(cores_split_mz_data_list[[i]],
                                    function(split_mz_data_element) future::future({run_split_df(split_mz_data_element)}))
    metabolite_identified_list <- lapply(identify_metab_future, future::value) # grab the results
    interm_core_split_future_df <- data.table::rbindlist(metabolite_identified_list, fill = TRUE)
    metabolite_identified_df <- rbind(metabolite_identified_df, interm_core_split_future_df)  
  }
  
  metabolite_identified_df <- as.data.frame(metabolite_identified_df, stringsAsFactors = FALSE)
  
  message("Perform Metabolite Identification Completed...")
  
  return(metabolite_identified_df)
}
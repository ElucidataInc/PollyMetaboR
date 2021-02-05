#' get_peak_detailed_from_emdb
#'
#' It extracts the peak detailed format of elmaven from emdb.
#'
#' @param emdb_path The path to emDB file
#' @param table_name The peak table name where desired metabolite is present
#' @return ggplot object or plotly object
#' @examples
#' get_peak_detailed_from_emdb(emdb_path, table_name)
#' @export
get_peak_detailed_from_emdb <- function(emdb_path, table_name) {
  message("Get Peak Detailed From EMDB Started...")

  if (identical(emdb_path, NULL)){
    warning("The emdb_path is NULL")
    return (NULL)
  }
  
  if (identical(table_name, NULL)){
    warning("The table_name is NULL")
    return (NULL)
  }    
  
  if (!file.exists(emdb_path)){
    warning(paste("The", sQuote(emdb_path), "is not a valid path", sep = " "))
    return (NULL)
  }
  
  con_db <- DBI::dbConnect(RSQLite::SQLite(), emdb_path)
  
  peaktables_db <- DBI::dbSendQuery(con_db, "SELECT DISTINCT table_name FROM peakgroups")
  peaktables_vec <- DBI::dbFetch(peaktables_db)[, 1]
  DBI::dbClearResult(peaktables_db)
  
  if (!(table_name %in% peaktables_vec)){
    warning(paste("The", table_name, "is not a valid peak table name", sep = " "))
    return (NULL)       
  }
  
  peakgroups_required_cols <- c('group_id', 'table_group_id')
  
  peaks_required_cols <- c('group_id', 'sample_id', 'peak_mz', 'mzmin', 'mzmax', 'rt', 'rtmin', 'rtmax', 'quality',
                           'peak_intensity', 'peak_area', 'peak_spline_area', 'peak_area_top', 'peak_area_corrected',
                           'peak_area_top_corrected', 'no_noise_obs', 'signal_baseline_ratio', 'from_blank_sample')
  
  merged_peaks_required_cols <- c('table_group_id', 'compound', 'compoundId', 'formula', 'sample', 
                                  'adductName', 'peak_mz', 'mzmin', 'mzmax', 'rt', 'rtmin', 'rtmax', 
                                  'quality', 'peak_intensity', 'peak_area', 'peak_spline_area', 'peak_area_top', 
                                  'peak_area_corrected', 'peak_area_top_corrected', 'no_noise_obs', 
                                  'signal_baseline_ratio', 'from_blank_sample')    
  
  peak_detailed_required_cols <- c('groupId', 'compound', 'compoundId', 'formula', 'sample', 'adductName', 'peakMz', 'mzmin',
                                   'mzmax', 'rt', 'rtmin', 'rtmax', 'quality', 'peakIntensity', 'peakArea', 'peakSplineArea',
                                   'peakAreaTop', 'peakAreaCorrected', 'peakAreaTopCorrected', 'noNoiseObs', 'signalBaseLineRatio',
                                   'fromBlankSample')  
  
  peakgroups_cols <- DBI::dbListFields(con_db, "peakgroups")
  peakml_cols <- c('peakML_label', 'peakML_probability')  
  if (all(peakml_cols %in% peakgroups_cols)){
    peakgroups_required_cols <- c(peakgroups_required_cols, peakml_cols)
    merged_peaks_required_cols <- c(merged_peaks_required_cols, peakml_cols)
    peak_detailed_required_cols <- c(peak_detailed_required_cols, peakml_cols)  
  }  
  
  table_name <- dQuote(table_name, q = getOption("UTF-8"))    
  peakgroups_condition <- paste0('SELECT ', paste(peakgroups_required_cols, collapse = ", "), ' FROM peakgroups WHERE table_name = ', table_name)
  peakgroups_db <- DBI::dbSendQuery(con_db, peakgroups_condition)
  peakgroups <- DBI::dbFetch(peakgroups_db)
  DBI::dbClearResult(peakgroups_db)
  
  group_id <- peakgroups[, "group_id"] 
  peaks_condition <- paste0('SELECT ', paste(peaks_required_cols, collapse = ", "), ' FROM peaks WHERE group_id IN ', paste0("(",paste0(group_id, collapse = ", "), ")"))
  peaks_db <- DBI::dbSendQuery(con_db, peaks_condition)
  peaks <- DBI::dbFetch(peaks_db)
  DBI::dbClearResult(peaks_db)
  
  samples <- DBI::dbReadTable(con_db, 'samples')
  DBI::dbDisconnect(con_db)
  
  peak_detailed_df <- dplyr::summarize_at(dplyr::group_by_at(peaks, "group_id"), "sample_id", list(mean_mz = ~round(mean(peak_mz), 6), mean_rt = ~round(mean(rt), 6)))    
  peak_detailed_df <- merge(peak_detailed_df, peaks, by = "group_id")
  peak_detailed_df$compound <- paste(peak_detailed_df$mean_mz,'@', peak_detailed_df$mean_rt, sep = "")
  peak_detailed_df$compoundId <- peak_detailed_df$compound 
  peak_detailed_df$formula <- ''
  peak_detailed_df$adductName <- ''
  
  peak_detailed_df <- merge(peak_detailed_df, samples[, c('sample_id', 'name')], by.x = 'sample_id', by.y = 'sample_id')        
  peak_detailed_df$sample <- peak_detailed_df$name
  peak_detailed_df$name <- NULL
  peak_detailed_df$sample_id <-NULL
  
  peak_detailed_df <- merge(peak_detailed_df, peakgroups, by = 'group_id')
  peak_detailed_df <- dplyr::rename_at(peak_detailed_df, dplyr::vars(merged_peaks_required_cols), ~ peak_detailed_required_cols)
  peak_detailed_df <- peak_detailed_df[, peak_detailed_required_cols]
  
  peak_detailed_df <- peak_detailed_df[with(peak_detailed_df, order(groupId, sample)), ]  
  
  message("Get Peak Detailed From EMDB Completed...")
  
  return (peak_detailed_df)
}
#' plot_eic_from_emdb
#'
#' It plots the extracted ion chromatogram (EIC) for specific groupid from emdb.
#'
#' @param emdb_path The path to emDB file
#' @param table_group_id The group id of desired metabolite
#' @param table_name The peak table name where desired metabolite is present
#' @param x_label Label x-axis
#' @param y_label Label x-axis
#' @param title_label Title of the plot
#' @param legend_label Label for the lagends
#' @param x_label_size Text size of x label
#' @param y_label_size Text size of y label
#' @param title_label_size Text size of text label
#' @param legend_label_size Text size of legend label
#' @param interactive Make plot interactive using plotly
#' @return ggplot object or plotly object
#' @examples
#' plot_eic_from_emdb <- function(emdb_path, table_group_id, table_name)
#' @export
plot_eic_from_emdb <- function(emdb_path = NULL, table_group_id = NULL, table_name = NULL, 
                               x_label = "RT", y_label = "Intensity", title_label = "", 
                               legend_label = "Sample", x_label_size = 15, y_label_size = 15, 
                               title_label_size = 15, legend_label_size = 15, interactive = T){
  message("Plot EIC From EMDB Started...")
  require(DBI)
  require(RSQLite)
  require(stringr)
  
  if (identical(emdb_path, NULL)){
    warning("The emdb_path is NULL")
    return (NULL)
  }
  
  if (identical(table_group_id, NULL)){
    warning("The table_group_id is NULL")
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
  
  con_db <- DBI::dbConnect(RSQLite::SQLite(), dbname = emdb_path)
  peaktables_db <- DBI::dbSendQuery(con_db, "SELECT DISTINCT table_name FROM peakgroups")
  peaktables_vec <- DBI::dbFetch(peaktables_db)[, 1]
  DBI::dbClearResult(peaktables_db)
  
  if (!(table_name %in% peaktables_vec)){
    warning(paste("The", sQuote(table_name), "is not a valid peak table name", sep = " "))
    return (NULL)       
  }
  
  table_name <- dQuote(table_name, q = getOption("UTF-8"))
  groupid_condition <- paste0('SELECT * FROM peakgroups WHERE table_group_id = ', table_group_id, ' AND table_name = ', table_name)
  peak_groups <- data.frame()
  try(expr = {
    peak_groups_db <- DBI::dbSendQuery(con_db, groupid_condition)
    peak_groups <- DBI::dbFetch(peak_groups_db)
    DBI::dbClearResult(peak_groups_db)
  }, silent = T)
  
  if (nrow(peak_groups) < 1){
    warning(paste("The", sQuote(table_group_id), "is not a valid table group id", sep = " "))
    return (NULL) 
  }
  
  compound_name <- peak_groups[, "compound_name"]
  group_id <- peak_groups[, "group_id"]
  peaks_condition <- paste0('SELECT * FROM peaks WHERE group_id = ', group_id)
  peaks <- data.frame()
  try(expr = {
    peaks_db <- DBI::dbSendQuery(con_db, peaks_condition)
    peaks <- DBI::dbFetch(peaks_db)
    DBI::dbClearResult(peaks_db)
  }, silent = TRUE)
  
  if (nrow(peaks) < 1){
    warning(paste("The", sQuote(table_group_id), "is not a valid table group id", sep = " "))
    return (NULL) 
  }    
  
  if (all(is.na(peaks$eic_rt)) | all(is.na(peaks$eic_intensity))){
    warning("The emDB file does not have eic raw data")
    return (NULL)   
  }
  
  samples_df <- DBI::dbReadTable(con_db, "samples")[, c("sample_id", "name")]
  colnames(samples_df) <- c("sample_id", "sample_name")
  peaks <- merge(samples_df, peaks, by = "sample_id")
  rtmin <- min(peaks$rtmin)
  rtmax <- max(peaks$rtmax)
  
  DBI::dbDisconnect(con_db)
  
  mean_mz <- round(mean(peaks$peak_mz), 6) 
  mean_rt <- round(mean(peaks$rt), 6)
  mz_at_rt <- paste(mean_mz, mean_rt, sep = "@")
  
  intensity_data <- data.frame()
  for (row_id in 1:nrow(peaks)){
    sample_name <- peaks$sample_name[row_id]
    int_data <- as.numeric(stringr::str_split(peaks$eic_intensity[row_id], pattern = ",")[[1]])
    rt_data <- as.numeric(stringr::str_split(peaks$eic_rt[row_id], pattern = ",")[[1]])
    interm_data <- data.frame(sample = sample_name, intensity = int_data, rt = rt_data, stringsAsFactors = FALSE)
    intensity_data <- rbind(intensity_data, interm_data)
  }
  
  if (identical(title_label, "") | identical(title_label, NULL)){
    if (identical(compound_name, "") | is.na(compound_name)){
      title_label <- mz_at_rt
    }
    else {
      title_label <- paste(compound_name, mz_at_rt, sep = "\n")
    }         
  }
  
  p <- PollyMetaboR::plot_eic(intensity_data, rt_min = rtmin, rt_max = rtmax, x_label = x_label, y_label = y_label, title_label = title_label, 
                              legend_label = legend_label, x_label_size = x_label_size, y_label_size = y_label_size, 
                              title_label_size = title_label_size, legend_label_size = legend_label_size, interactive = interactive)
  
  message("Plot EIC From EMDB Completed...")
  
  return (p)
}
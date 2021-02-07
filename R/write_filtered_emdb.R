#' write_filtered_emdb
#'
#' It filters the emDB by groupIds and saves in new emDB file.
#'
#' @param emdb_path The path to emDB file
#' @param table_name The peak table name where desired metabolite is present
#' @param filtered_groups_df A dataframe with filtered groupids
#' @param save_path The path to saved filtered emDB file
#' @return A file is saved at save_path
#' @examples
#' get_peak_detailed_from_emdb(emdb_path, filtered_groups_df, table_name)
#' @export
write_filtered_emdb <- function(emdb_path = NULL, table_name = NULL, 
                                filtered_groups_df = NULL, save_path = NULL){
  message("Write Filtered EMDB Started...")
  
  if (identical(emdb_path, NULL)){
    warning("The emdb_path is NULL")
    return (NULL)
  }
  
  if (identical(table_name, NULL)){
    warning("The table_name is NULL")
    return (NULL)
  }    
  
  if (identical(filtered_groups_df, NULL)){
    warning("The filtered_groups_df is NULL")
    return (NULL)
  }    
  
  if (!file.exists(emdb_path)){
    warning(paste("The", sQuote(emdb_path), "is not a valid path", sep = " "))
    return (NULL)
  }
  
  if (!identical(tools::file_ext(emdb_path), "emDB")){
    warning(paste("The", sQuote(emdb_path), "does not end with '.emDB'", sep = " "))
    return (NULL)      
  }
  
  if (!(identical(as.character(class(filtered_groups_df)), "data.frame"))){
    warning("The filtered_groups_df is not a dataframe")
    return (NULL)      
  } 
  
  if(!("groupId" %in% colnames(filtered_groups_df))){
    warning("The 'groupId' column is not present in filtered_groups_df")  
    return (NULL)
  }
  
  if(nrow(filtered_groups_df) < 1){
    warning("The filtered_groups_df is empty dataframe")   
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
  
  if (identical(save_path, NULL) | identical(save_path, "")){
    new_emdb_path <- paste0(tools::file_path_sans_ext(emdb_path), "_filtered.emDB")
  }
  else{
    if (dir.exists(save_path)){
      new_emdb_path <- paste0(tools::file_path_sans_ext(emdb_path), "_filtered.emDB")    
    }
    else {
      if (!identical(tools::file_ext(save_path), "emDB")){
        filname <- tools::file_path_sans_ext(save_path)
        new_emdb_path <- paste0(filname, ".emDB")
      }
      else{
        new_emdb_path <- save_path 
      }
    }  
  }  
  
  if (file.exists(new_emdb_path)) {
    file.remove(new_emdb_path)
  }
  
  new_emdb <- DBI::dbConnect(RSQLite::SQLite(), dbname = new_emdb_path)
  message(paste0("The filtered emDB file is saved at : ", sQuote(new_emdb_path, q = getOption("UTF-8"))))
  
  alignment_rts <- DBI::dbReadTable(con_db, 'alignment_rts')
  DBI::dbWriteTable(new_emdb, "alignment_rts", alignment_rts, overwrite = TRUE)
  
  samples <- DBI::dbReadTable(con_db, 'samples')
  DBI::dbWriteTable(new_emdb, "samples", samples, overwrite = TRUE)      
  
  user_settings_global_db <- DBI::dbSendQuery(con_db, "SELECT * FROM user_settings WHERE domain = 'global'")
  user_settings_global <- DBI::dbFetch(user_settings_global_db)
  DBI::dbClearResult(user_settings_global_db)
  DBI::dbWriteTable(new_emdb, "user_settings", user_settings_global, overwrite = TRUE)
  
  table_name <- dQuote(table_name, q = getOption("UTF-8"))
  groupids_db <- DBI::dbSendQuery(con_db, paste0('SELECT DISTINCT table_group_id FROM peakgroups WHERE table_name = ', table_name))
  groupids_vec <- DBI::dbFetch(groupids_db)[, 1]
  DBI::dbClearResult(groupids_db)    
  
  filtered_group_ids <- unique(filtered_groups_df[, 'groupId'])  
  common_groupids <- intersect(groupids_vec, filtered_group_ids)  
  if (length(common_groupids) < 1){
    warning("No common groupId found in emdb and filtered groups data")
    return (NULL) 
  }
  
  for (groupid in common_groupids){
    peakgroups_condition <- paste0('SELECT * FROM peakgroups WHERE table_group_id = ', groupid, ' AND table_name = ', table_name)
    peakgroups_db <- DBI::dbSendQuery(con_db, peakgroups_condition)
    peakgroups <- DBI::dbFetch(peakgroups_db)
    DBI::dbClearResult(peakgroups_db)
    DBI::dbWriteTable(new_emdb, "peakgroups", peakgroups, append = TRUE, overwrite = FALSE)    
    
    compound_id <- dQuote(peakgroups[, "compound_id"], q = getOption("UTF-8"))
    compound_name <- dQuote(peakgroups[, "compound_name"], q = getOption("UTF-8"))
    compound_db <- dQuote(peakgroups[, "compound_db"], q = getOption("UTF-8"))
    compounds_condition <- paste0('SELECT * FROM compounds WHERE compound_id = ', compound_id, ' AND name = ', compound_name, ' AND db_name = ', compound_db)
    compounds_db <- DBI::dbSendQuery(con_db, compounds_condition)
    compounds <- DBI::dbFetch(compounds_db)
    DBI::dbClearResult(compounds_db)      
    DBI::dbWriteTable(new_emdb, "compounds", compounds, append = TRUE, overwrite = FALSE)
    
    group_id <- peakgroups[, "group_id"]
    peaks_condition <- paste0('SELECT * FROM peaks WHERE group_id = ', group_id)
    peaks_db <- DBI::dbSendQuery(con_db, peaks_condition)
    peaks <- DBI::dbFetch(peaks_db)
    DBI::dbClearResult(peaks_db)
    DBI::dbWriteTable(new_emdb, "peaks", peaks, append = TRUE, overwrite = FALSE)      
    
    user_settings_condition <- paste0('SELECT * FROM user_settings WHERE domain = ', group_id)
    user_settings_db <- DBI::dbSendQuery(con_db, user_settings_condition)
    user_settings <- DBI::dbFetch(user_settings_db)
    DBI::dbClearResult(user_settings_db)
    DBI::dbWriteTable(new_emdb, "user_settings", user_settings, append = TRUE, overwrite = FALSE)  
  }
  
  message("Write Filtered EMDB Completed...")
  
}
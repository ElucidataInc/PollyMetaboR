#' update_label_in_emdb
#'
#' It updates the label of feature in emDB file.
#'
#' @param emdb_path The path to emDB file
#' @param table_group_id The group id of desired metabolite
#' @param table_name The peak table name where desired metabolite is present
#' @param label_group_id A label of that feature, select from (good, bad, maybe)
#' @return Update label in emDB file
#' @examples
#' update_label_in_emdb(emdb_path, table_group_id, table_name, label_group_id)
#' @export
update_label_in_emdb <- function(emdb_path = NULL, table_group_id = NULL,
                                 table_name = NULL, label_group_id = NULL){
  message("Update Label In EMDB Started...")
  
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
  
  if (identical(label_group_id, NULL)){
    warning("The label_group_id is NULL")
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
  con_db <- DBI::dbConnect(RSQLite::SQLite(), dbname = emdb_path)
  peaktables_db <- DBI::dbSendQuery(con_db, "SELECT DISTINCT table_name FROM peakgroups")
  peaktables_vec <- DBI::dbFetch(peaktables_db)[, 1]
  DBI::dbClearResult(peaktables_db)  
  if (!(table_name %in% peaktables_vec)){
    warning(paste0("The ", sQuote(table_name), " is not a valid peak table name, please select from ", "(",
                   paste0(sQuote(peaktables_vec), collapse = ", "), ")"))
    return (NULL)       
  }
  
  table_name <- dQuote(table_name, q = getOption("UTF-8"))
  
  groupids_db <- DBI::dbSendQuery(con_db, paste0('SELECT DISTINCT table_group_id FROM peakgroups WHERE table_name = ', table_name))
  groupids_vec <- DBI::dbFetch(groupids_db)[, 1]
  DBI::dbClearResult(groupids_db)    
  if (!(table_group_id %in% groupids_vec)){
    warning(paste("The", sQuote(table_group_id), "is not a valid table group id", sep = " "))
    return (NULL) 
  }
  
  if (identical(label_group_id, "good")){
    label <- "g"
    peakML_label <- "Signal"  
    peakML_label_id <- 1 
  }
  else if (identical(label_group_id, "bad")){
    label <- "b"
    peakML_label <- "Noise"  
    peakML_label_id <- 0      
  }
  else if (identical(label_group_id, "maybe")){
    label <- ""
    peakML_label <- "May be good signal"  
    peakML_label_id <- 5      
  }
  else{
    warning(paste0("Select valid label from ","(", paste0(sQuote(c("good", "bad", "maybe")), collapse = ", "), ")"))
    return (NULL)  
  }         
  
  peakml_cols <- c("peakML_label", "peakML_label_id")  
  peakgroups_cols <- DBI::dbListFields(con_db, "peakgroups")
  label <- dQuote(label, q = getOption("UTF-8"))    
  peakML_label <- dQuote(peakML_label, q = getOption("UTF-8"))    
  if (all(peakml_cols %in% peakgroups_cols)){ 
    label_condition <- paste0("UPDATE peakgroups SET label = ", label, ", peakML_label = ", peakML_label,
                              ", peakML_label_id =", peakML_label_id," WHERE table_group_id = ", table_group_id, 
                              " AND table_name = ", table_name)
  }
  else{
    label_condition <- paste0("UPDATE peakgroups SET label = ", label," WHERE table_group_id = ", table_group_id, 
                              " AND table_name = ", table_name) 
  }   
  
  execute <- DBI::dbExecute(con_db, label_condition)
  
  message("Update Label In EMDB Completed...")
}
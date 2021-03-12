#' restructure_camera_annotations
#'
#' It restructures the camera output to group features
#'
#' @param annotated_peaks_df The camera output
#' @param polarity The polarity of the data
#' @return A dataframe with restructure camera output
#' @examples
#' restructure_camera_annotations(annotated_peaks_df, polarity = 'positive')
#' @export
restructure_camera_annotations <- function(annotated_peaks_df = NULL, polarity = 'positive') {
  message("Restructure CAMERA Annotations Started...")
  
  if (identical(annotated_peaks_df, NULL)){
    warning("The annotated_peaks_df is NULL")
    return (NULL)  
  }
  
  if (!identical(as.character(class(annotated_peaks_df)), "data.frame")){
    warning("The annotated_peaks_df parameter is not a dataframe")
    return (NULL)  
  }
  
  if (identical(polarity, NULL)){
    warning("Polarity is NULL")
    return (NULL)  
  }
  
  if (!(polarity %in% c("positive", "negative"))){
    warning("Please select polarity from positive and negative")
    return (NULL)  
  }    
  
  annotated_peaks_df <- transform(annotated_peaks_df, pcgroup = as.numeric(pcgroup))
  H_mass <- 1.007825
  # Replace emplty values with "NA"
  annotated_peaks_df[annotated_peaks_df==""] <- NA
  # Sorting of table by "pcgroup" number
  camera_result_df_sorted = annotated_peaks_df[order(annotated_peaks_df$pcgroup),]
  pcgroup_vec <- unique(camera_result_df_sorted$pcgroup) 
  all_isotope_adduct_grouped_result_df <- data.frame()
  all_unannotated_camera_df <- data.frame()
  count = 1
  
  for (pcgroup_id in pcgroup_vec){
    pcgroup_camera_result_df_sorted <- dplyr::filter(camera_result_df_sorted, pcgroup == pcgroup_id)
    
    # Filter cells by NAs in columns "isotopes" and "adduct"
    isna_isotopes <- is.na(pcgroup_camera_result_df_sorted$isotopes)
    isna_adducts <- is.na(pcgroup_camera_result_df_sorted$adduct)                   
    filter_isotopes_adducts_df <- pcgroup_camera_result_df_sorted[!(isna_isotopes) | !(isna_adducts),]
    
    if (nrow(filter_isotopes_adducts_df) > 0){
      ## Add basemass column
      get_basemass_with_adduct <- function(adduct_str = ''){
        basemass_with_adduct_list = list()
        adduct_type_vec <- vector()
        basemass_vec <- numeric()
        if (is.na(adduct_str) | (is.nan(adduct_str))) {
          basemass_with_adduct_list[['adduct_type']] <- adduct_type_vec
          basemass_with_adduct_list[['basemass']] <- basemass_vec
          return (basemass_with_adduct_list)}
        separate_adducts_str <- strsplit(gsub("(\\[.*?\\])(.*?)\\ ", "\\1\\2\\_", as.character(adduct_str)), " ")[[1]]
        for (i_index in 1:length(separate_adducts_str)){
          basemass_adduct_pair <- strsplit(separate_adducts_str[i_index], "_")[[1]]
          adduct_type <- as.character(basemass_adduct_pair[1])
          basemass <- as.numeric(basemass_adduct_pair[2])
          if (!is.na(basemass)){
            adduct_type_vec <- c(adduct_type_vec, adduct_type)
            basemass_vec <- c(basemass_vec, basemass)
          }
        }
        basemass_with_adduct_list[['adduct_type']] <- adduct_type_vec
        basemass_with_adduct_list[['basemass']] <- basemass_vec
        
        return (basemass_with_adduct_list)
      }
      
      # Make dataframe for basemass column
      basemass_camera_df = data.frame()
      for (i_row in 1:nrow(filter_isotopes_adducts_df)){
        row_df <- filter_isotopes_adducts_df[i_row,]
        adduct <- as.character(row_df$adduct)
        basemass_adduct_vec = get_basemass_with_adduct(adduct)
        if (length(basemass_adduct_vec[[1]])==0){
          row_df$adduct_type <- NA
          row_df$basemass <- NA
          basemass_camera_df <- rbind(basemass_camera_df, row_df)
        }
        else{
          for (i_index in 1: length(basemass_adduct_vec[[1]])){
            row_df$adduct_type <- basemass_adduct_vec[[1]][i_index]
            row_df$basemass <- basemass_adduct_vec[[2]][i_index]
            basemass_camera_df <- rbind(basemass_camera_df, row_df)
          }
        }
      }
      
      # Add isotope_id and isotope_type columns
      basemass_camera_df_cp <- basemass_camera_df
      split_isotopes_str <- function(str_to_split = '', iso_par = ''){
        str_split_vec = vector()
        if (is.na(str_to_split) | (is.nan(str_to_split))) return (str_split_vec)
        str_split_vec <- stringr::str_extract_all(as.character(str_to_split), "\\[(.*?)\\]")[[1]]
        str_split_vec <- substring(str_split_vec, 2, nchar(str_split_vec)-1)
        if (iso_par == 'isotope_id'){ return (str_split_vec[1]) }
        else if (iso_par == 'isotope_type'){return (str_split_vec[2])}
        else {return (str_split_vec)}
      }
      
      # Add new columns to the basemass dataframe
      isotope_id_list <- as.numeric(apply(basemass_camera_df['isotopes'],1, function(isotopes) split_isotopes_str(isotopes, iso_par='isotope_id')))
      isotope_type_list <- as.character(apply(basemass_camera_df['isotopes'],1, function(isotopes) split_isotopes_str(isotopes, iso_par='isotope_type')))
      
      if ((length(isotope_id_list) > 0) & (length(isotope_type_list) > 0)){
        basemass_camera_df_cp$isotope_id <- isotope_id_list
        basemass_camera_df_cp$isotope_type <- isotope_type_list
      } else {
        basemass_camera_df_cp$isotope_id <- NA
        basemass_camera_df_cp$isotope_type <- NA
      }
      
      basemass_camera_df_cp[basemass_camera_df_cp=='logical(0)'] <- NA                          
      basemass_camera_df_cp$isotope_id <- as.numeric(unlist(basemass_camera_df_cp$isotope_id))                           
      basemass_camera_df_cp$isotope_type <- as.character(unlist(basemass_camera_df_cp$isotope_type))                               
      
      #Group features having adducts with their isotopes
      all_basemass_vec <- unique(basemass_camera_df_cp[!(is.na(basemass_camera_df_cp$basemass)),]$basemass)
      cluster_adduct_isotopes_df <- data.frame()
      for (basemass_val in all_basemass_vec){
        mass_df <- dplyr::filter(basemass_camera_df_cp, basemass==basemass_val)
        mass_df$feature_group <- count
        isotope_id_vec <- as.numeric(mass_df$isotope_id)
        isotope_id_vec <- isotope_id_vec[!is.na(isotope_id_vec)]
        cluster_adduct_isotopes_df <- rbind(cluster_adduct_isotopes_df, mass_df)
        for (isotope_id_val in isotope_id_vec){
          isotope_id_df <- unique(dplyr::filter(basemass_camera_df_cp, isotope_id==isotope_id_val, basemass %in% c(basemass_val, NA)))
          adduct_type_val <- dplyr::filter(isotope_id_df, isotope_type == "M")$adduct_type[[1]]
          isotope_id_df$basemass <- basemass_val
          isotope_id_df$adduct_type <- adduct_type_val
          isotope_id_df$feature_group <- count
          cluster_adduct_isotopes_df <- unique(rbind(cluster_adduct_isotopes_df, isotope_id_df))
        }
        count <- count+1
      }
      
      ## Group features having only isotopes without adducts
      adduct_isotope_id_vec <- as.numeric(cluster_adduct_isotopes_df$isotope_id)
      adduct_isotope_id_vec <- unique(adduct_isotope_id_vec[!is.na(adduct_isotope_id_vec)])
      
      all_isotope_id_vec = as.numeric(basemass_camera_df_cp$isotope_id)
      all_isotope_id_vec <- unique(all_isotope_id_vec[!is.na(all_isotope_id_vec)])
      
      only_rest_isotopes_id_vec <- setdiff(all_isotope_id_vec, adduct_isotope_id_vec)
      
      rest_isotope_df <- data.frame()
      
      for (isotope_id_val in only_rest_isotopes_id_vec){
        interm_isotope_df <- dplyr::filter(basemass_camera_df_cp, isotope_id == isotope_id_val)
        parent_isotope_df <- dplyr::filter(interm_isotope_df, isotope_type == 'M')
        parent_isotope_mz <- parent_isotope_df$mz
        if (polarity=='positive'){
          interm_isotope_df$adduct_type <- '[M+H]+'
          interm_isotope_df$basemass <- parent_isotope_mz - H_mass
        } else {
          interm_isotope_df$adduct_type <- '[M-H]-'
          interm_isotope_df$basemass <- parent_isotope_mz + H_mass  
        }
        
        interm_isotope_df$feature_group <- count
        rest_isotope_df <- rbind(rest_isotope_df, interm_isotope_df)
        count <- count + 1
        
      }
      # Annotated dataframe
      isotope_adduct_grouped_result_df <- dplyr::bind_rows(cluster_adduct_isotopes_df, rest_isotope_df)
      
      
      # Combine pcgroup unannotated dataframe to global dataframe
      all_isotope_adduct_grouped_result_df <- dplyr::bind_rows(all_isotope_adduct_grouped_result_df, isotope_adduct_grouped_result_df)  
      all_isotope_adduct_grouped_result_df <- all_isotope_adduct_grouped_result_df[with(all_isotope_adduct_grouped_result_df, 
                                                                                        order(pcgroup, feature_group, adduct_type, 
                                                                                              isotope_id, isotope_type)), ]
    }
    
    #Provide feature_group id to features without adducts and isotopes
    # Making un-annotated dataframe            
    pcgroup_unannotated_camera_df <- pcgroup_camera_result_df_sorted[isna_isotopes & isna_adducts, ]
    #print (unannotated_camera_df)
    if (nrow(pcgroup_unannotated_camera_df) > 0){
      unannotated_camera_df <- pcgroup_unannotated_camera_df
      new_cols_add <- c('adduct_type', 'basemass' ,'isotope_id', 'isotope_type','feature_group')
      unannotated_camera_df[, new_cols_add] <- NA
      unannotated_initial_count <- count
      count <- count + length(unannotated_camera_df$groupId)
      unannotated_camera_df$feature_group <- seq(from = unannotated_initial_count, to = count - 1, by = 1)
      
      if (polarity=='positive'){
        unannotated_camera_df$adduct_type <- '[M+H]+'
        unannotated_camera_df$basemass <- unannotated_camera_df$mz - H_mass
      } else {
        unannotated_camera_df$adduct_type <- '[M-H]-'
        unannotated_camera_df$basemass <- unannotated_camera_df$mz + H_mass  
      }
      
      # Combine pcgroup unannotated dataframe to global dataframe
      all_unannotated_camera_df <- dplyr::bind_rows(all_unannotated_camera_df, unannotated_camera_df)
      all_unannotated_camera_df <- all_unannotated_camera_df[with(all_unannotated_camera_df, 
                                                                  order(pcgroup, feature_group, adduct_type, 
                                                                        isotope_id, isotope_type)), ]
    }
  }
  # Combining all-annotated and all-un-annotated dataframe
  final_grouped_camera_result_df <- dplyr::bind_rows(all_isotope_adduct_grouped_result_df, all_unannotated_camera_df)
  final_grouped_camera_result_df <- final_grouped_camera_result_df[with(final_grouped_camera_result_df, 
                                                                        order(pcgroup, feature_group, adduct_type, 
                                                                              isotope_id, isotope_type)), ]
  
  if (nrow(all_isotope_adduct_grouped_result_df) > 0) {row.names(all_isotope_adduct_grouped_result_df) <- NULL }
  if (nrow(all_unannotated_camera_df) > 0) {row.names(all_unannotated_camera_df) <- NULL }
  if (nrow(final_grouped_camera_result_df) > 0) {row.names(final_grouped_camera_result_df) <- NULL }
  
  output_dfs <- list("annotated" = all_isotope_adduct_grouped_result_df, "unannotated" = all_unannotated_camera_df, "combined" = final_grouped_camera_result_df)
  
  message("Restructure CAMERA Annotations Completed...")
  
  return(output_dfs)
}
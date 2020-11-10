#' create_xcms_object_from_elmaven
#'
#' It create xcms object from peak detailed format of El-MAVEN
#'
#' @param maven_output_df The peak detailed format of El-MAVEN output
#' @param polarity The polarity of the data
#' @return xcms object
#' @examples
#' create_xcms_object_from_elmaven(maven_output_df)
#' @import xcms dplyr
#' @export
create_xcms_object_from_elmaven <- function(maven_output_df = NULL, polarity = NULL){
  message("Create XCMS Object From ElMaven Started...")
  require(xcms)
  require(dplyr)
  
  if (identical(maven_output_df, NULL)){
    warning("The maven_output_df is NULL")
    return (NULL)
  }
  
  if (!identical(as.character(class(maven_output_df)), "data.frame")){
    warning("The maven_output_df parameter is not a dataframe")
    return (NULL)  
  }
  
  required_cols <- c('peakMz', 'mzmin', 'mzmax', 'rt', 'rtmin', 'rtmax', 'groupId', 'sample', 'peakIntensity', 'peakAreaCorrected')  
  if (!all(required_cols %in% colnames(maven_output_df))){
    warning(c("The maven_output_df dataframe should have the following columns: ", paste0(required_cols, collapse = ", ")))
    return (NULL)
  }
  
  if (identical(polarity, NULL)){
    polarity <- character()
  }
  
  maven_output_df$peakMz <- as.numeric(maven_output_df$peakMz)
  maven_output_df <- dplyr::filter(maven_output_df, !is.na(peakMz) & peakMz != 0)  
  samples_vec <- unique(maven_output_df$sample)
  
  # Make peaks df
  maven_output_df$Index <- 1:length(rownames(maven_output_df))
  
  elmaven_peaks_df <- data.frame(maven_output_df, check.names = FALSE, stringsAsFactors = FALSE)
  
  for (sample_index in 1:length(samples_vec)){
    sample_val <- as.character(samples_vec[sample_index])
    elmaven_peaks_df <- dplyr::mutate(elmaven_peaks_df, sample = replace(as.character(sample), sample == sample_val, sample_index))
  }
  elmaven_peaks_df$sample <- as.numeric(elmaven_peaks_df$sample)
  elmaven_peaks_df <-  dplyr::rename_at(elmaven_peaks_df, dplyr::vars(c('peakMz', 'peakIntensity', 'peakAreaCorrected')), ~ c('mz', 'into', 'maxo'))
  
  drops <- c('compound', 'compound', 'compoundId', 'formula')
  elmaven_peaks_df <- elmaven_peaks_df[ , !(names(elmaven_peaks_df) %in% drops)]
  elmaven_peaks_df <- elmaven_peaks_df[c('mz', 'mzmin', 'mzmax', 'rt', 'rtmin', 'rtmax', 
                                         'into', 'maxo', 'Index', 'groupId', 'sample')]
  
  # Make metadata
  metadata_df <- data.frame("sample" = samples_vec)
  metadata_df$class <- "sample"
  row.names(metadata_df) <- metadata_df$sample
  metadata_df$sample <- NULL
  
  # Make groupId list
  groupIdx_list <- list()
  
  # Make groups df
  elmaven_groups_df <- data.frame()
  ident_groups_vec <- vector()  
  index_count <- 0
  for (groupId_val in unique(maven_output_df$groupId)){
    index_count <- index_count + 1
    feature_df <- dplyr::filter(maven_output_df, groupId == groupId_val, sample %in% samples_vec)
    compound_val <- as.character(unique(feature_df$compound))
    samples_index_v <- feature_df$Index
    groupIdx_list[[index_count]] <- samples_index_v
    
    mz_rt_bool <- grepl('[0-9]+@[0-9]+', compound_val)
    if (mz_rt_bool) {
      group_mz_rt <- strsplit(compound_val, "@")
      group_mzmed <- as.numeric(group_mz_rt[[1]][1])
      group_rtmed <- as.numeric(group_mz_rt[[1]][2])
      identified_val <- 0  
    } else {
      ident_groups_vec <- c(ident_groups_vec, groupId_val)  
      group_mzmed <- round(mean(feature_df$peakMz), 6)
      group_rtmed <- round(mean(feature_df$rt), 3)
      identified_val <- 1  
    }
    
    mzmin <- min(feature_df$mzmin)
    mzmax <- max(feature_df$mzmax)
    rtmin <- min(feature_df$rtmin)
    rtmax <- max(feature_df$rtmax)
    interm_df <- data.frame("groupId" = groupId_val, identified = identified_val,
                            "mzmed" = group_mzmed, "mzmin" = mzmin, "mzmax" = mzmax, 
                            "rtmed" = group_rtmed, "rtmin" = rtmin, "rtmax" = rtmax)
    elmaven_groups_df <- rbind(elmaven_groups_df, interm_df)
  }
  elmaven_groups_df <- elmaven_groups_df[c('mzmed', 'mzmin', 'mzmax', 'rtmed', 'rtmin' ,'rtmax', 'groupId', 'identified')]
  
  # Make xcms object
  xcms_obj <- new("xcmsSet",
                  peaks = as.matrix(elmaven_peaks_df),
                  groups = as.matrix(elmaven_groups_df),
                  groupidx = groupIdx_list,
                  phenoData = metadata_df,
                  filepaths = as.vector(row.names(metadata_df)), 
                  polarity = as.character(polarity))
  
  message("Create XCMS Object From ElMaven Completed...")
  
  return(xcms_obj)
}
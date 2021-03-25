#' create_xcms_object_from_elmaven
#'
#' It create xcms object from peak detailed format of El-MAVEN
#'
#' @param maven_output_df The peak detailed format of El-MAVEN output
#' @param polarity The polarity of the data
#' @return xcms object
#' @examples
#' create_xcms_object_from_elmaven(maven_output_df)
#' @import dplyr xcms
#' @export
create_xcms_object_from_elmaven <- function(maven_output_df = NULL, polarity = NULL){
  message("Create XCMS Object From ElMaven Started...")
  require(dplyr)  
  require(xcms)
  
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
  
  # Make peaks
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
  
  # Make groups  
  get_mz_and_rt <- function(compound, peakMz, rt){
    compound <- unique(compound)  
    mz_rt_bool <- grepl('[0-9]+@[0-9]+', compound)
    if (mz_rt_bool) {
      group_mz_rt <- strsplit(compound, "@")
      group_mzmed <- round(as.numeric(group_mz_rt[[1]][1]), 6)
      group_rtmed <- round(as.numeric(group_mz_rt[[1]][2]), 6)
      identified_val <- 0  
    } else {
      group_mzmed <- round(mean(peakMz), 6)
      group_rtmed <- round(mean(rt), 6)
      identified_val <- 1  
    } 
    return(list(group_mzmed, group_rtmed, identified_val))   
  }
  
  elmaven_groups_df  <- maven_output_df %>% dplyr::group_by(groupId) %>% dplyr::summarise(mzmed = get_mz_and_rt(compound, peakMz, rt)[[1]],
                                                                                          rtmed = get_mz_and_rt(compound, peakMz, rt)[[2]],
                                                                                          identified = get_mz_and_rt(compound, peakMz, rt)[[3]],
                                                                                          mzmin = round(min(mzmin), 6), mzmax = round(max(mzmax), 6),
                                                                                          rtmin = round(min(rtmin), 6), rtmax = round(max(rtmax), 6))
  
  elmaven_groups_df <- elmaven_groups_df[, c('mzmed', 'mzmin', 'mzmax', 'rtmed', 'rtmin' ,'rtmax', 'groupId', 'identified')]
  
  # Make groupidx    
  groupidx <- elmaven_groups_df[, "groupId", drop = FALSE]
  groupidx$idx <- row.names(groupidx)  
  groupidx <- merge(groupidx, maven_output_df[, c("groupId", "Index"), drop = FALSE], by = "groupId")
  groupidx  <- groupidx %>% dplyr::group_by(idx) %>% dplyr::summarise(sample_index = list(Index))
  groupidx$idx <- as.numeric(groupidx$idx)  
  groupidx <- groupidx[with(groupidx, order(as.numeric(idx))), ]  
  groupidx <- as.list(stats::setNames(groupidx$sample_index, groupidx$idx))
  names(groupidx)<- groupidx$idx
  
  # Create xcms object
  xcms_obj <- new("xcmsSet",
                  peaks = as.matrix(elmaven_peaks_df),
                  groups = as.matrix(elmaven_groups_df),
                  groupidx = groupidx,
                  phenoData = metadata_df,
                  filepaths = as.vector(row.names(metadata_df)), 
                  polarity = as.character(polarity))
  
  message("Create XCMS Object From ElMaven Completed...")
  
  return(xcms_obj)
}
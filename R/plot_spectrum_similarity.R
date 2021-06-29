#' plot_spectrum_similarity
#'
#' It plots the query spectra w.r.t the reference spectra
#'
#' @param spec_query The query spectrum dataframe with mz and intensity columns
#' @param spec_ref The reference spectrum dataframe with mz and intensity columns
#' @param mz_tolerence Value of mz tolerence
#' @param mz_tolerence_unit The mz tolerence unit (ppm or Da)
#' @param annotate_mz Annotate mz values on spectrum, select from TRUE/FALSE
#' @param annotate_mz_cutoff The intensity cutoff in percentage w.r.t base peak (highest intensity peak) used to annotate mz peaks
#' @param intensity_cutoff The intensity cutoff in percentage w.r.t base peak (highest intensity peak) used to filter data
#' @param mz_precision The numeric precision of mz which is equal to the number of digits to show on mz annotation
#' @param mz_range The mz range used to calculate score
#' @param x_label Label x-axis
#' @param y_label Label x-axis
#' @param title_label Title of the plot
#' @param top_label Specify top label
#' @param bottom_label Specify bottom label
#' @param interactive Make plot interactive using plotly
#' @return ggplot object or plotly object
#' @examples
#' plot_spectrum_similarity(spec_query, spec_ref, mz_tolerence = 10, mz_tolerence_unit = "ppm")
#' @import plotly ggplot2
#' @export
plot_spectrum_similarity <- function (spec_query = NULL, spec_ref = NULL,  mz_tolerence = 20,
                                      mz_tolerence_unit = "ppm", annotate_mz = FALSE, annotate_mz_cutoff = NULL,
                                      intensity_cutoff = NULL, mz_precision = NULL, mz_range = NULL,
                                      x_label = NULL, y_label = NULL, title_label = NULL,
                                      top_label = NULL, bottom_label = NULL, interactive = FALSE){
  message("Plot Spectrum Similarity Started...")
  require(plotly)
  require(ggplot2)
  require(ggrepel)
  
  if(identical(spec_query, NULL)){
    warning("No spec_query was given")
    return (NULL)
  }
  
  if ((!class(spec_query)=='data.frame') | (!any(c("mz", "intensity") %in% colnames(spec_query)))){
    warning("Please input valid spec_query (a dataframe containing 'mz' and 'intensity' columns)")
    return (NULL)
  }
  
  if (nrow(spec_query) < 1){
    warning("The spec_query is an empty dataframe") 
    return (NULL)    
  } 
  
  if(identical(spec_ref, NULL)){
    warning("No spec_ref was given")
    return (NULL)
  }
  
  if ((!class(spec_ref)=='data.frame') | (!any(c("mz", "intensity") %in% colnames(spec_ref)))){
    warning("Please input valid spec_ref (a dataframe containing 'mz' and 'intensity' columns)")
    return (NULL)
  }
  
  if (nrow(spec_ref) < 1){
    warning("The spec_ref is an empty dataframe") 
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
  
  if (identical(mz_tolerence, NULL)){
    warning("The mz_tolerence is NULL")
    return (NULL)  
  }
  
  if(!identical(mz_tolerence, NULL)){
    mz_tolerence <- as.numeric(mz_tolerence)
    if (is.na(mz_tolerence)){
      warning("The mz_tolerence is not a numeric value") 
      return (NULL)
    }  
  }
  
  if(!identical(annotate_mz_cutoff, NULL)){
    annotate_mz_cutoff <- as.numeric(annotate_mz_cutoff)
    if (is.na(annotate_mz_cutoff)){
      warning("The annotate_mz_cutoff is not a numeric value") 
      return (NULL)
    }  
  }    
  
  if(!identical(intensity_cutoff, NULL)){
    intensity_cutoff <- as.numeric(intensity_cutoff)
    if (is.na(intensity_cutoff)){
      warning("The intensity_cutoff is not a numeric value") 
      return (NULL)
    }  
  }    
  
  if(!identical(mz_precision, NULL)){
    mz_precision <- as.numeric(mz_precision)
    if (is.na(mz_precision)){
      warning("The mz_precision is not a numeric value") 
      return (NULL)
    }  
  }
  
  if(!identical(mz_range, NULL)){
    if (!is.numeric(mz_range)){
      warning("The mz_range is not a numeric vector") 
      return (NULL)
    }
    
    if (length(mz_range) != 2){
      warning("The mz_range is not a numeric vector of two elements") 
      return (NULL)
    }      
  }  
  
  if (identical(x_label, NULL)){
    x_label <- "mz"
  }                                              
  if (identical(y_label, NULL)){
    y_label <- "intensity(%)"
  }                                                                                           
  if (identical(title_label, NULL)){
    title_label <- ""
  }                                              
  if (identical(top_label, NULL)){
    top_label <- "Query Spectra"
  }                                             
  if (identical(bottom_label, NULL)){
    bottom_label <- "Reference Spectra"
  }    
  
  query_tmp <- data.frame(mz = spec_query[, 1], intensity = spec_query[, 2])
  if (!identical(mz_range, NULL)){ query_tmp <- subset(query_tmp, mz >= mz_range[1] & mz <= mz_range[2])}  
  if (nrow(query_tmp) < 1){
    warning("No mz values found in specified mz range in query_tmp data")
    return (NULL)   
  }
  query_tmp$normalized <- round(((query_tmp$intensity/max(query_tmp$intensity)) * 100), 2)
  if(!identical(intensity_cutoff, NULL)){ query_tmp <- subset(query_tmp, normalized >= intensity_cutoff)}  
  query_score <- data.frame(mz = query_tmp$mz, intensity = query_tmp$normalized)
  query_score$matched_color <- "black"
  
  ref_tmp <- data.frame(mz = spec_ref[, 1], intensity = spec_ref[, 2])
  if (!identical(mz_range, NULL)){ ref_tmp <- subset(ref_tmp, mz >= mz_range[1] & mz <= mz_range[2])}  
  if (nrow(ref_tmp) < 1){
    warning("No mz values found in specified mz range in ref_tmp data")
    return (NULL)   
  }
  ref_tmp$normalized <- round(((ref_tmp$intensity/max(ref_tmp$intensity)) * 100), 2)
  if(!identical(intensity_cutoff, NULL)){ ref_tmp <- subset(ref_tmp, normalized >= intensity_cutoff)}  
  ref_score <- data.frame(mz = ref_tmp$mz, intensity = -ref_tmp$normalized)
  ref_score$matched_color <- "red"   
  ref_score$matched_type <- "mismatched"
  
  query_score$delta <- suppressMessages(sapply(query_score[, 1], function(x) PollyMetaboR::calculate_mz_delta(x, mz_tolerence_unit = mz_tolerence_unit, mz_tolerence = mz_tolerence)))
  query_score$mzmin <- query_score$mz - query_score$delta
  query_score$mzmax <- query_score$mz + query_score$delta
  
  for (i in 1:nrow(ref_score)) {
    query_matched_df <- subset(query_score, mzmin <= ref_score[i, "mz"] & mzmax >= ref_score[i, "mz"])  
    if (nrow(query_matched_df) >= 1) {
      ref_score[i, "matched_color"] <- "blue"
      ref_score[i, "matched_type"] <- "matched"
    }
  }
  
  if (nchar(top_label) > 0){
    query_score$data_type <- paste("data : ", top_label, sep = "")  
  } else{
    query_score$data_type <- ""  
  }                                            
  
  if (nchar(bottom_label) > 0){
    ref_score$data_type <- paste("data : ", bottom_label, "<br>", "matched_type : ", ref_score$matched_type, sep = "")  
  } else{
    ref_score$data_type <- paste("matched_type : ", ref_score$matched_type, sep = "")  
  }
  
  if(!identical(mz_precision, NULL)){ query_score[, "mz_annotate"] <- round(query_score$mz, mz_precision)}
  else { query_score[, "mz_annotate"] <- query_score$mz}
  
  if(!identical(annotate_mz_cutoff, NULL)){
    query_score[query_score[, "intensity"] < annotate_mz_cutoff, "mz_annotate"] <- ""
  }
  
  p <- ggplot() + 
    geom_bar(data = query_score, aes(x = mz, y = intensity, fill = matched_color, 
                                     text = data_type, width = 1), 
             stat = "identity", position = "dodge", orientation = "x") +
    geom_bar(data = ref_score, aes(x = mz, y = intensity, fill = matched_color,
                                   text = data_type, width = 1),
             stat = "identity", position = "dodge",  orientation = "x") +
    ggtitle(title_label) + 
    labs(x = x_label, y = y_label) +
    geom_hline(yintercept = 0, size = 0.4) +                                    
    ggsci::scale_color_aaas() + scale_fill_identity() +
    theme(axis.line = element_line(size = 1, colour = "black"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_blank(), panel.background = element_blank(), 
          plot.title = element_text(colour = "black", size = 18, 
                                    face = "plain", hjust = 0.5), 
          axis.title = element_text(colour = "black",margin = 1,
                                    size = 14, face = "plain"), 
          axis.text.x = element_text(colour = "black", size = 10, angle = 90, hjust = 1, 
                                     margin = unit(c(0.2,0.2, 0.1, 0.1), "cm"), face = "plain"),
          axis.text.y = element_text(colour = "black", size = 10, 
                                     margin = unit(c(0.2, 0.2, 0.1, 0.1), "cm"), face = "plain"),
          axis.ticks.length = unit(0.25, "cm"), legend.position = "none") +
    ggplot2::annotate(geom = "text",  x=Inf, y = Inf, label = top_label, vjust=1, hjust=1) + 
    ggplot2::annotate(geom = "text",  x= Inf, y = -Inf, label = bottom_label, vjust=-1, hjust=1)
  
  if (!identical(interactive, TRUE)){ 
    if (identical(annotate_mz, TRUE)){ p <- p + ggrepel::geom_text_repel(data = query_score, aes(x = mz, y = intensity, text = data_type, label = mz_annotate))}     
  }
  else {
    if (identical(annotate_mz, TRUE)){ p <- p + geom_text(data = query_score, aes(x = mz, y = intensity, text = data_type, label = mz_annotate))}                                            
    p <- ggplotly(p, tooltip = c("x", "y", "text")) %>% layout(hovermode = TRUE) %>% 
      plotly::config(displaylogo = FALSE,
                     modeBarButtons = list(list("zoom2d"),
                                           list("select2d"),
                                           list("lasso2d"),
                                           list("autoScale2d"),
                                           list("resetScale2d"),
                                           list("pan2d"),
                                           list("zoomIn2d"), 
                                           list("zoomOut2d"),
                                           list("hoverClosestCartesian"),
                                           list('toImage')), 
                     mathjax = 'cdn')  
  }                                            
  
  message("Plot Spectrum Similarity Completed...")
  
  return (p)
}
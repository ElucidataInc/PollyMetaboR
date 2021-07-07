#' plot_mass_spectrum
#'
#' It plots the mass spectrum
#'
#' @param ms_spectrum The query spectrum dataframe with mz and intensity columns
#' @param annotate_mz Annotate mz values on spectrum, select from TRUE/FALSE
#' @param annotate_mz_cutoff The intensity cutoff in percentage w.r.t base peak (highest intensity peak) used to annotate mz peaks
#' @param intensity_cutoff The intensity cutoff in percentage w.r.t base peak (highest intensity peak) used to filter data
#' @param mz_precision The numeric precision of mz which is equal to the number of digits to show on mz annotation
#' @param show_rel_intensity Show relative intensity in percentage (TRUE) on y axis which is calculated w.r.t base peak
#' @param mz_range The mz range used to calculate score
#' @param start_from_zero Force x and y axes to start from zero (TRUE/FALSE)
#' @param x_label Label x-axis
#' @param y_label Label x-axis
#' @param title_label Title of the plot
#' @param interactive Make plot interactive using plotly
#' @return ggplot object or plotly object
#' @examples
#' plot_mass_spectrum(ms_spectrum)
#' @import plotly ggplot2
#' @export
plot_mass_spectrum <- function (ms_spectrum = NULL, annotate_mz = FALSE, annotate_mz_cutoff = NULL, 
                                intensity_cutoff = NULL, mz_precision = NULL, show_rel_intensity = FALSE,
                                mz_range = NULL, start_from_zero = TRUE, x_label = NULL, y_label = NULL, 
                                title_label = NULL, interactive = FALSE){
  message("Plot Mass Spectrum Started...")
  require(plotly)
  require(ggplot2)
  require(ggrepel)  
  
  if(identical(ms_spectrum, NULL)){
    warning("No ms_spectrum was given")
    return (NULL)
  }
  
  if ((!class(ms_spectrum)=='data.frame') | (!any(c("mz", "intensity") %in% colnames(ms_spectrum)))){
    warning("Please input valid ms_spectrum (a dataframe containing 'mz' and 'intensity' columns)")
    return (NULL)
  }
  
  if (nrow(ms_spectrum) < 1){
    warning("The ms_spectrum is an empty dataframe") 
    return (NULL)    
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
  
  if (identical(show_rel_intensity, TRUE)){ 
    intensity_type <- "rel_intensity"
    if (identical(y_label, NULL)){ y_label <- "intensity(%)"}     
  }
  else {
    intensity_type <- "intensity"
    if (identical(y_label, NULL)){ y_label <- "intensity"}    
  }
  
  if (identical(x_label, NULL)){
    x_label <- "mz"
  }                                              
  
  if (identical(title_label, NULL)){
    title_label <- ""
  }                                              
  
  if (!identical(mz_range, NULL)){ ms_spectrum <- subset(ms_spectrum, mz >= mz_range[1] & mz <= mz_range[2])}
  
  if (nrow(ms_spectrum) < 1){
    warning(paste0("No data present in the mz range ", "(", paste(mz_range, collapse = ", "), ")"))
    return (NULL)    
  }  
  
  if(!identical(mz_precision, NULL)){ ms_spectrum[, "mz_annotate"] <- round(ms_spectrum$mz, mz_precision)}
  else { ms_spectrum[, "mz_annotate"] <- ms_spectrum$mz}
  
  ms_spectrum[, "rel_intensity"] <- round(((ms_spectrum$intensity/max(ms_spectrum$intensity)) * 100), 2)    
  if(!identical(annotate_mz_cutoff, NULL)){
    ms_spectrum[ms_spectrum[, "rel_intensity"] < annotate_mz_cutoff, "mz_annotate"] <- ""
  }
  if(!identical(intensity_cutoff, NULL)){ ms_spectrum <- subset(ms_spectrum, rel_intensity >= intensity_cutoff)}
  
  ms_spectrum$text_label<-  paste0(paste0("mz: ", ms_spectrum$mz), "<br>", paste0("intensity: ", ms_spectrum$intensity),
                                   "<br>", paste0("relatve intensity(%): ", ms_spectrum$rel_intensity))
  
  p <- ggplot(data = ms_spectrum, aes_string(x = "mz", y = intensity_type, text = "text_label", label = "mz_annotate",  width = 0.1)) + 
    geom_bar(color = "red", fill = "red",
             stat = "identity", position = "dodge", orientation = "x") +
    labs(x = x_label, y = y_label, title = title_label) +
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
          axis.ticks.length = unit(0.25, "cm"), legend.position = "none")
  
  
  if (identical(start_from_zero, TRUE)){ p <- p + expand_limits(x = 0, y = 0)}  
  
  if (identical(interactive, FALSE)){
    if (identical(annotate_mz, TRUE)){ p <- p + ggrepel::geom_text_repel()}     
  }
  else {
    if (identical(annotate_mz, TRUE)){ p <- p + geom_text(position = ggplot2::position_stack(vjust = 1.05))}
    p <- ggplotly(p, tooltip = c("text")) %>% layout(hovermode = TRUE)  %>% 
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
  
  message("Plot Mass Spectrum Completed...")
  
  return (p)
}
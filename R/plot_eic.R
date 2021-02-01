#' plot_eic
#'
#' It plots the extracted ion chromatogram (EIC) for provided intensity data.
#'
#' @param intensity_data The intensity data with sample, intensity and rt columns
#' @param rt_min Make dashed line at rt_min value
#' @param rt_max Make dashed line at rt_max value
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
#' plot_eic(intensity_data, rt_min, rt_min)
#' @export
plot_eic <- function(intensity_data = NULL, rt_min = NULL, rt_max = NULL, x_label = "RT",
                     y_label = "Intensity", title_label = "", legend_label = "Sample",
                     x_label_size = 15, y_label_size = 15, title_label_size = 20, 
                     legend_label_size = 15, interactive = TRUE){
  message("Plot EIC Started...")
  require(dplyr)
  require(plotly)
  require(ggplot2)
  
  if (identical(intensity_data, NULL)){
    warning("The intensity_data is NULL")
    return (NULL)  
  }
  
  required_cols <- c("sample", "intensity", "rt")
  if (!all(required_cols %in% colnames(intensity_data))){
    warning(paste0("The intensity_data should contain the following columns : ", paste0(required_cols, collapse = ", ")))
    return (NULL)  
  }
  
  if (!identical(rt_min, NULL)){
    rt_min <- as.numeric(rt_min)
    if (is.na(rt_min)){
      warning("The rt_min is not a numeric value, ignoring rt_min line.")
      rt_min <- NULL
    }
  }
  
  if (!identical(rt_max, NULL)){
    rt_max <- as.numeric(rt_max)
    if (is.na(rt_max)){
      warning("The rt_max is not a numeric value, ignoring rt_max line.")
      rt_max <- NULL
    }
  }
  
  if (identical(x_label, NULL)){ x_label <- "" }
  
  if (identical(y_label, NULL)){ y_label <- "" }
  
  if (identical(title_label, NULL)){ title_label <- "" }
  
  if (identical(legend_label, NULL)){ legend_label <- "" } 
  
  intensity_data <- intensity_data[with(intensity_data, order(sample, rt, intensity)), ]
  sample_with_intmax_df <- intensity_data %>% dplyr::group_by(sample) %>% dplyr::summarise(int_max = max(intensity))
  ordered_samples <- sample_with_intmax_df[with(sample_with_intmax_df, order(-int_max)), ]$sample
  intensity_data$sample <- factor(intensity_data$sample, levels = ordered_samples)
  
  if (interactive == TRUE){
    p <- plot_ly(data = intensity_data) 
    
    if (!identical(rt_min, NULL)){
      p <- p %>% add_segments(x = rt_min, xend = rt_min, y = 0, yend = max(intensity_data$intensity), 
                              line = list(dash = "dash", color ="grey", width = 1), opacity = 0.7, showlegend = F)
    }
    
    if (!identical(rt_max, NULL)){
      p <- p %>% add_segments(x = rt_max, xend = rt_max, y = 0, yend = max(intensity_data$intensity), 
                              line = list(dash = "dash", color ="grey", width = 1),  opacity = 0.7, showlegend = F) 
    }
    
    p <- p %>% add_trace(x = ~rt, y = ~intensity, color = ~sample, text = ~sample, 
                         fill= 'tozeroy', type = 'scatter', mode = 'lines') %>%        
      layout(title = list(text = title_label, font = list(size = title_label_size), xref = "paper", yref = "paper", x = 0.5), 
             xaxis = list(title = x_label, titlefont = list(size = y_label_size), 
                          autotick = TRUE, ticks = "outside", tick0 = 0, ticklen = 5, showgrid = F, showline = TRUE), 
             yaxis = list(title = y_label, titlefont = list(size = x_label_size), rangemode = "tozero",
                          autotick = TRUE, ticks = "outside", tick0 = 0, ticklen = 5, showgrid = F, showline = TRUE),
             legend = list(font = list(size = 11)), showlegend = T, 
             margin = list(t = 100)) %>% 
      add_annotations(text = legend_label, xref = "paper", yref = "paper",
                      x=1.04, xanchor="left",
                      y=0.99, yanchor="bottom",
                      font = list(size = legend_label_size),
                      legendtitle = TRUE, showarrow = FALSE)  %>% 
      plotly::config(displaylogo = FALSE,
                     modeBarButtons = list(list("zoomIn2d"), list("zoomOut2d"), list('toImage')), 
                     mathjax = 'cdn')
  }
  else {
    p <- ggplot(intensity_data, aes(rt, intensity, color = sample, fill = sample)) +
      theme_linedraw() +
      geom_area(alpha = 0.5) +
      geom_vline(xintercept = rt_min, size = 0.3, linetype = 'dashed') +
      geom_vline(xintercept = rt_max, size = 0.3, linetype = 'dashed') +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(title = title_label, x = x_label, y = y_label, color = legend_label, fill = legend_label) + # x and y axis labels
      #ggsci::scale_color_aaas() + # filling the point colors
      theme(legend.position = "right", legend.direction = "vertical", # legend positioned at the bottom, horizantal direction,
            axis.line = element_line(size = 1, colour = "black"), # axis line of size 1 inch in black color
            panel.grid.major = element_blank(), # major grids included
            panel.grid.minor = element_blank(), # no minor grids
            panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
            plot.title = element_text(colour = "black", size = title_label_size, face = "bold", hjust = 0.5),
            axis.title.x = element_text(colour ="black", size = x_label_size, face = "bold"), # axis title
            axis.title.y = element_text(colour = "black", size = y_label_size, face = "bold"), # axis title   
            axis.text.x = element_text(colour = "black", size = 10, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), # x-axis text in fontsize 10
            axis.text.y = element_text(colour = "black", size = 10, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), # y-axis text in fontsize 10
            legend.text = element_text(size = 10, face = "bold"),
            legend.title = element_text(colour = "black", size = title_label_size, face = "bold"),
            axis.ticks.length = unit(0.25, "cm")) # ticks facing inward with 0.25cm length
  }
  
  message("Plot EIC Completed...")
  
  return (p)
}
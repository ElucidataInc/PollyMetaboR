#' plot_eic
#'
#' It plots the extracted ion chromatogram (EIC) for provided intensity data.
#'
#' @param intensity_data The intensity data with sample, intensity and rt columns
#' @param rt_min Make dashed line at rt_min value
#' @param rt_max Make dashed line at rt_max value
#' @param set_rt_range Set rt range to plot sliced eic
#' @param x_label Label x-axis
#' @param y_label Label x-axis
#' @param title_label Title of the plot
#' @param legend_label Label for the lagends
#' @param x_label_size Text size of x label
#' @param y_label_size Text size of y label
#' @param title_label_size Text size of text label
#' @param legend_label_size Text size of legend label
#' @param x_text_size Text size of x ticks
#' @param y_text_size Text size of y ticks
#' @param legend_text_size Text size of legends
#' @param interactive Make plot interactive using plotly
#' @return ggplot object or plotly object
#' @examples
#' plot_eic(intensity_data, rt_min, rt_min)
#' @import dplyr plotly ggplot2
#' @export
plot_eic <- function(intensity_data = NULL, rt_min = NULL, rt_max = NULL, set_rt_range = NULL,
                     x_label = "RT", y_label = "Intensity", title_label = "", legend_label = "Sample",
                     x_label_size = 16, y_label_size = 16, title_label_size = 16, 
                     legend_label_size = 16, x_text_size = 14, y_text_size = 14, 
                     legend_text_size = 14, interactive = TRUE){
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
  
  if (!identical(set_rt_range, NULL)){  
    set_rt_range <- as.numeric(set_rt_range)
    if (length(set_rt_range) != 2){
      warning("The set_rt_range is not a numeric vector of 2 values of rt range, considering all rt values")   
    }
    else {  
      if (any(is.na(set_rt_range))){
        warning("The set_rt_range is not a numeric vector, considering all rt values")  
      }
      else{  
        intensity_data <- dplyr::filter(intensity_data, rt >= set_rt_range[1], rt <= set_rt_range[2])
        if (nrow(intensity_data) < 1){
          warning("The filtered data has 0 rows, please use valid set_rt_range values")
          return (NULL)  
        }  
      }
    }  
  }
  
  intensity_data <- intensity_data[with(intensity_data, order(sample, rt, intensity)), ]
  sample_with_intmax_df <- intensity_data %>% dplyr::group_by(sample) %>% dplyr::summarise(int_max = max(intensity))
  ordered_samples <- sample_with_intmax_df[with(sample_with_intmax_df, order(-int_max)), ]$sample
  intensity_data$sample <- factor(intensity_data$sample, levels = ordered_samples)
  sample_by_maxint <- data.frame(intensity_data)
  
  if (!identical(rt_min, NULL)){
    sample_by_maxint <- dplyr::filter(sample_by_maxint, rt >= rt_min)
  }
  
  if (!identical(rt_max, NULL)){
    sample_by_maxint <- dplyr::filter(sample_by_maxint, rt <= rt_max)
  }    
  
  if(nrow(sample_by_maxint) < 1){
    warning("The rt_min or rt_max or both are out of rt range")
    return (NULL)
  }
  
  sample_by_maxint <- sample_by_maxint %>% dplyr::group_by(sample) %>% dplyr::slice(which.max(intensity)) 
  
  if (interactive == TRUE){
    p <- plot_ly() 
    
    if (!identical(rt_min, NULL)){
      p <- p %>% add_segments(x = rt_min, xend = rt_min, y = 0, yend = max(intensity_data$intensity), name = "rt min",
                              line = list(dash = "dash", color ="grey", width = 1), opacity = 0.8, showlegend = F)
    }
    
    if (!identical(rt_max, NULL)){
      p <- p %>% add_segments(x = rt_max, xend = rt_max, y = 0, yend = max(intensity_data$intensity), name = "rt max",
                              line = list(dash = "dash", color ="grey", width = 1),  opacity = 0.8, showlegend = F) 
    }
    
    for (sample_n in levels(intensity_data$sample)){
      eic_data <- dplyr::filter(intensity_data, sample %in% sample_n)
      maxint_data <- dplyr::filter(sample_by_maxint, sample %in% sample_n)
      p <- p %>% add_trace(x = eic_data$rt, y = eic_data$intensity, color = sample_n, text = sample_n, 
                           hoverinfo = 'x+y+text', fill= 'tozeroy', type = 'scatter', mode = 'lines', 
                           legendgroup = sample_n)
      p <- p %>% add_trace(x = maxint_data$rt, y = maxint_data$intensity, color = sample_n, text = sample_n, 
                           hoverinfo = 'x+y+text', type = 'scatter', mode = 'markers', marker = list(size = 18, opacity = 1), 
                           legendgroup = sample_n, showlegend = F)
    }  
    
    p <- p %>% 
      layout(title = list(text = title_label, font = list(size = title_label_size), xref = "paper", yref = "paper", x = 0.5), 
             xaxis = list(title = x_label, titlefont = list(size = y_label_size), tickfont = list(size = x_text_size), 
                          autotick = TRUE, ticks = "outside", tick0 = 0, ticklen = 5, showgrid = F, showline = TRUE), 
             yaxis = list(title = y_label, titlefont = list(size = x_label_size), tickfont = list(size = y_text_size), rangemode = "tozero",
                          autotick = TRUE, ticks = "outside", tick0 = 0, ticklen = 5, showgrid = F, showline = TRUE),
             legend = list(font = list(size = legend_text_size), tracegroupgap = 1), showlegend = T, 
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
      geom_line() +
      geom_area(position = 'identity', alpha = 0.5) +
      geom_vline(xintercept = rt_min, size = 0.3, linetype = 'dashed') +
      geom_vline(xintercept = rt_max, size = 0.3, linetype = 'dashed') +
      geom_point(data = sample_by_maxint, aes(rt, intensity, color = sample, fill = sample), size = 7) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(sample_by_maxint$intensity) + 0.03* max(sample_by_maxint$intensity))) +
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
            axis.text.x = element_text(colour = "black", size = x_text_size, margin=unit(c(0.2,0.2,0.1,0.1), "cm"), face = "bold"), # x-axis text in fontsize 10
            axis.text.y = element_text(colour = "black", size = y_text_size, margin=unit(c(0.2,0.2,0.1,0.1), "cm"), face = "bold"), # y-axis text in fontsize 10
            legend.text = element_text(size = legend_text_size, face = "bold"),
            legend.title = element_text(colour = "black", size = title_label_size, face = "bold"),
            axis.ticks.length = unit(0.25, "cm")) # ticks facing inward with 0.25cm length
  }
  
  message("Plot EIC Completed...")
  
  return (p)
}
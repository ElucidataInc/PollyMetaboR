#' plot_hist_elements_frequency
#'
#' It plots a histogram showing the ferquency of each element
#'
#' @param elements_vec The dataframe 1
#' @param frequency_type Show bar plot frequency by occurrence (by_occurrence) or by elements (by_elements)
#' @param plot_title The title of the plot
#' @param xaxis_title The title of the x-axis
#' @param yaxis_title The title of the y-axis
#' @return The plotly object of bar plot
#' @examples
#' plot_hist_elements_frequency(elements_vec, frequency_type = "by_occurrence")
#' @import plotly
#' @export
plot_hist_elements_frequency <- function(elements_vec = NULL, frequency_type = "by_occurrence", plot_title = NULL, xaxis_title = NULL, yaxis_title = NULL){
  message("Plot Histogram Elements Frequency Started...")
  require(plotly)
  
  if (identical(elements_vec, NULL)){
    warning("The elements_vec is NULL")  
  }  
  
  if (length(elements_vec) < 1){
    warning("The elements_vec is not a vector")
    return (NULL)
  }
  
  if (identical(frequency_type, NULL)){
    warning("The frequency_type is NULL")
    return (NULL)  
  }
  
  if (!(frequency_type %in% c("by_occurrence", "by_elements"))){
    warning("The frequency_type should be from by_occurrence and by_elements")
    return (NULL)  
  } 
  
  frequency_df <- data.frame(table(Element = elements_vec))  
  if (frequency_type %in% "by_occurrence"){
    frequency_df$Freq <- factor(frequency_df$Freq, levels = as.character(sort(unique(frequency_df$Freq))))
    p <- plot_ly(x = frequency_df$Freq, type = "histogram", 
                 marker = list(line = list(color = 'black', width = 0.5)))  %>% 
      layout(title = plot_title, yaxis = list(title = yaxis_title), xaxis = list(title = xaxis_title))
  }
  
  if (frequency_type %in% "by_elements"){
    frequency_df$Element <- factor(frequency_df$Element, levels = as.character(sort(unique(frequency_df$Element))))    
    p <- plot_ly(x = frequency_df$Freq, y = frequency_df$Element, type = "bar", orientation = 'h',
                 marker = list(line = list(color = 'black', width = 1.5)))  %>% 
      layout(title = plot_title, yaxis = list(title = yaxis_title), xaxis = list(title = xaxis_title))
  }
  
  message("Plot Histogram Elements Frequency Completed...")
  
  return(p)
}
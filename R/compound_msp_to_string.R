#' compound_msp_to_string
#'
#' It converts the list of elements of a single compound msp to a string text
#'
#' @param compound_msp A list of elements of a single compound msp
#' @return A string text of a single compound msp
#' @examples
#' compound_msp_to_string(compound_msp)
#' @export
compound_msp_to_string <- function(compound_msp = NULL){
  message("Compound MSP to String Started...") 
  
  if (identical(compound_msp, NULL)){
    warning("The compound_msp is NULL")
    return (NULL)
  }
  
  msp_str <- NULL
  for (id in names(compound_msp)){
    if (identical(id, "Name")){
      msp_str <- paste(id, compound_msp[[id]], sep = ": ")
    }
    else {
      if (!identical(id, "pspectrum")){
        msp_str <- paste(msp_str, paste(id, compound_msp[[id]], sep = ": "), sep = "\r\n")
      }
      else {
        pspectrum_str <- paste(apply(compound_msp$pspectrum, 1, function(x) paste(x, collapse = " ")), collapse = "\r\n")
        msp_str <- paste(msp_str, pspectrum_str, sep = "\r\n")
      }                                
    }
  }
  
  message("Compound MSP to String Completed...")
  
  return (msp_str)                                   
}
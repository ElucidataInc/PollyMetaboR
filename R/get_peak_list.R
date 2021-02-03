#' get_peak_list
#'
#' Extract all information from an xsAnnotate object. Returns a peaklist with annotated peaks. 
#' It is cloned from getPeaklist of CAMERA.
#'
#' @param object xsAnnotate object
#' @param intval Choose intensity values. Allowed values are into, maxo, intb, intf, maxf, area, depending on the feature detection algorithm used.
#' @return The dataframe of grouped peaktable
#' @examples
#' getPeaklist(object, intval="maxo")
#' @export
get_peak_list <- function(object, intval="maxo") {
  message("Get Peak List Started...")
  
  if (!identical(as.character(class(object)), "xsAnnotate")){
    warning(" The object is not the xsAnnotate object")
    return (NULL)  
  }
  
  if (! intval %in% colnames(peaks(object@xcmsSet))) {
    stop("unknown intensity value!")
  }
  
  #generate peaktable
  #Check if xcmsSet contains only one sample
  if(object@sample[1] == 1 & length(sampnames(object@xcmsSet)) == 1){
    #intval is here ignored since all intensity values are already contained
    peaktable <- object@groupInfo;
  }else {
    #Case of xcmsSet with multiple samples
    #Use groupInfo information and replace intensity values
    peaktable <- object@groupInfo;
    
    #get intensity values from xcmsSet
    grpval <- groupval(object@xcmsSet, value=intval);
    
    #get column range for replacement
    grpval.ncol <- ncol(grpval)
    start <- ncol(peaktable) - grpval.ncol +1;
    ende  <- start + grpval.ncol - 1; 
    
    peaktable[, start:ende] <- grpval;
  }
  
  #allocate variables for CAMERA output
  adduct   <- vector("character", nrow(object@groupInfo));
  isotopes <- vector("character", nrow(object@groupInfo));
  pcgroup  <- vector("character", nrow(object@groupInfo));
  
  #default polarity set to positive
  polarity <- "+";
  
  if(length(object@polarity) > 0){
    if(object@polarity == "negative"){
      polarity <- "-";
    }
  }
  
  #First isotope informationen and adduct informationen
  for(i in seq(along = isotopes)){
    #check if adduct annotation is present for peak i
    if(length(object@derivativeIons) > 0 && !(is.null(object@derivativeIons[[i]]))) {
      #Check if we have more than one annotation for peak i
      if(length(object@derivativeIons[[i]]) > 1) {
        #combine ion species name and rounded mass hypophysis
        names <- paste(object@derivativeIons[[i]][[1]]$name, round(object@derivativeIons[[i]][[1]]$mass, 6));
        for(ii in 2:length(object@derivativeIons[[i]])) {
          names <- paste(names, object@derivativeIons[[i]][[ii]]$name, round(object@derivativeIons[[i]][[ii]]$mass, 6));
        }
        #save name in vector adduct
        adduct[i] <- names;
      } else {
        #Only one annotation
        adduct[i] <- paste(object@derivativeIons[[i]][[1]]$name, round(object@derivativeIons[[i]][[1]]$mass, 6));
      }
    } else {
      #no annotation empty name
      adduct[i] <- ""; 
    }
    
    #Check if we have isotope informationen about peak i
    if(length(object@isotopes) > 0&& !is.null(object@isotopes[[i]])) {
      num.iso <- object@isotopes[[i]]$iso;
      #Which isotope peak is peak i?
      if(num.iso == 0){
        str.iso <- "[M]";
      } else { 
        str.iso <- paste("[M+", num.iso, "]", sep="")
      }
      #Multiple charged?
      if(object@isotopes[[i]]$charge > 1){
        isotopes[i] <- paste("[", object@isotopes[[i]]$y, "]", str.iso, object@isotopes[[i]]$charge, polarity, sep="");
      }else{
        isotopes[i] <- paste("[", object@isotopes[[i]]$y, "]", str.iso, polarity, sep="");
      }
    } else { 
      #No isotope informationen available
      isotopes[i] <- ""; 
    }
  }
  
  #Have we more than one pseudospectrum?
  if(length(object@pspectra) < 1){
    pcgroup <- 0;
  } else {
    for(i in seq(along = object@pspectra)){
      index <- object@pspectra[[i]];
      pcgroup[index] <- i;
    }
  }
  
  rownames(peaktable)<-NULL;#Bugfix for: In data.row.names(row.names, rowsi, i) :  some row.names duplicated:
  
  message("Get Peak List Completed...")
  
  return(invisible(data.frame(peaktable, isotopes, adduct, pcgroup, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)));
}
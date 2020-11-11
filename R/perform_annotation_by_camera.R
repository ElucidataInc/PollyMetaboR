#' perform_annotation_by_camera
#'
#' It performs annotation by CAMERA using different functions of it.
#'
#' @param xcms_object The xcms object having peaks information
#' @param polarity The polarity of the data
#' @param ppm General ppm error
#' @param mzabs General absolute error in m/z
#' @param adducts_rules A dataframe with different adducts rules
#' @param isotope_matrix four column m/z-diff and ratio Matrix, for matching isotopic peaks
#' @param cor_exp_th Threshold for intensity correlations across samples
#' @param pval p-value threshold for testing correlation of significance (0-1)
#' @param perfwhm percentage of FWHM (Full Width at Half Maximum) width used in groupFWHM function for grouping features
#' @param sigma multiplier of the standard deviation used in groupFWHM function for grouping features
#' @param maxcharge maximum ion charge
#' @param maxiso maximum number of expected isotopes (0-8)
#' @param minfrac The percentage number of samples, which must satisfy the C12/C13 rule for isotope annotation
#' @param multiplier If no ruleset is provided, calculate ruleset with max. number n of [nM+x] cluster ions
#' @param max_peaks How many peaks will be calculated in every thread using the parallel mode
#' @return A dataframe of annotated peaks
#' @examples
#' perform_annotation_by_camera(xcms_object, polarity = "positive", ppm = 5, mzabs = 0.01)
#' @import CAMERA
#' @export
perform_annotation_by_camera <- function(xcms_object = NULL, polarity = NULL, ppm = 10, mzabs = 0,
                                         adducts_rules = NULL, isotope_matrix = NULL, cor_exp_th = 0.75, 
                                         pval = 0.05, perfwhm = 0.6, sigma = 6, maxcharge = 3, 
                                         maxiso = 4, minfrac = 0.5, multiplier = 3, max_peaks = 100){
  message("Perform Annotation BY CAMERA Started...")
  require(CAMERA)
  
  if (identical(xcms_object, NULL)){
    warning("The xcms_object is NULL")
    return (NULL)  
  }
  
  if (!identical(as.character(class(xcms_object)), "xcmsSet")){
    warning(" The xcms_object is not the xcmsSet object")
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
  
  if (!identical(adducts_rules, NULL)){
    if (!identical(as.character(class(adducts_rules)), "data.frame")){
      warning("The adducts_rules parameter is not a dataframe")
      return (NULL)  
    }
    
    adducts_rules_required_cols <- c('name', 'nmol', 'charge', 'massdiff', 'oidscore', 'quasi', 'ips')
    if (!all(adducts_rules_required_cols %in% colnames(adducts_rules))){
      warning(c("The adducts_rules dataframe should have the following columns: ", paste0(adducts_rules_required_cols, collapse = ", ")))
      return (NULL)  
    }
  }
  
  if (identical(isotope_matrix, NULL)){  
    isotope_matrix <- PollyMetaboR::isotope_matrix
    C13_mass = 13.0033548378
    C12_mass = 12.00000000   
    C_isotopic_mz <- C13_mass - C12_mass
    C_isotopic_mz_delta <- PollyMetaboR::calculate_mz_delta(mz_source = C_isotopic_mz, mz_tolerence_unit = "ppm", mz_tolerence = ppm)
    isotope_matrix[, "mzmin"] <- round(C_isotopic_mz - (mzabs + C_isotopic_mz_delta), 6)
    isotope_matrix[, "mzmax"] <- round(C_isotopic_mz + (mzabs + C_isotopic_mz_delta), 6)    
  } else {
    if (!identical(as.character(class(isotope_matrix)), "matrix")){
      warning("The isotope_matrix parameter is not a matrix")
      return (NULL) 
    }
    
    if (!all(c("mzmin", "mzmax", "intmin", "intmax") %in% colnames(isotope_matrix))){
      warning("The isotope_matrix should have the following columns : mzmin, mzmax, intmin, intmax")
      return (NULL)        
    } 
  }
  
  an <- CAMERA::xsAnnotate(xcms_object, polarity = polarity)
  anF <- CAMERA::groupFWHM(an, sigma = sigma, perfwhm = perfwhm, intval = "maxo")
  anIC <- CAMERA::groupCorr(anF, calcCiS = FALSE, calcCaS = TRUE, cor_exp_th = cor_exp_th, pval = pval, intval = "maxo")
  anI <- CAMERA::findIsotopes(anIC, maxcharge = maxcharge, maxiso = maxiso, ppm = ppm, mzabs = mzabs, intval = "maxo", minfrac = minfrac, isotopeMatrix = isotope_matrix)
  anFA <- CAMERA::findAdducts(anI, ppm = ppm, mzabs = mzabs, multiplier = multiplier, polarity = polarity, rules = adducts_rules, max_peaks = max_peaks, intval = "maxo")
  annotated_peaks_df <- PollyMetaboR::get_peak_list(anFA, intval = "maxo")
  
  message("Perform Annotation BY CAMERA Completed...")
  
  return(annotated_peaks_df)
}
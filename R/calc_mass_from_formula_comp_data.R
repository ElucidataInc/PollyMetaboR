#' calc_mass_from_formula_comp_data
#'
#' It calculates monoisotopic mass from molecular formula within compound database
#'
#' @param comp_data The compound database having atleast formula column
#' @return The compound database with calculated mass for the compounds
#' @examples
#' calc_mass_from_formula_comp_data(comp_data)
#' @import stringr Rdisop
#' @export
calc_mass_from_formula_comp_data <- function(comp_data = NULL){
  message("Calc Mass From Formula Comb Data Started...")
  
  if (identical(comp_data, NULL)){
    warning("The comp_data is NULL")  
  }
  
  if ((!class(comp_data)=='data.frame') | (!any(c("mass", "formula") %in% colnames(comp_data)))){
    warning("Please input valid compound database (a dataframe containing 'mass' or 'formula' column)")
    return (NULL)
  }    
  
  if (!("mass" %in% colnames(comp_data))){
    comp_data$mass <- NA
  } 
  
  if ("formula" %in% colnames(comp_data)){
    get_mass_from_formula <- function(x){
      exact_mass <- as.numeric(x[["mass"]])
      if (is.na(exact_mass)){
        comp_formula <- stringr::str_trim(as.character(x[["formula"]]))  
        if (!(comp_formula == "")){
          try(x[["mass"]] <- as.numeric(Rdisop::getMolecule(formula = comp_formula, maxisotopes = 1)$exactmass), silent = TRUE)
        }
      }
      return (x)
    }    
    
    comp_data <- as.data.frame(t(apply(comp_data, 1, function(x) get_mass_from_formula(x))), check.names = FALSE, stringsAsFactors = FALSE) 
    comp_data$mass <- as.numeric(comp_data$mass)
  }
  
  message("Calc Mass From Formula Comb Data Completed...")
  
  return (comp_data)
}
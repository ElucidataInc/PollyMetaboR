#' calc_spectrum_similarity
#'
#' It calculates the spectrum similarity score using weighted dot product algorithm.
#'
#' @param spec_query The query spectrum dataframe with mz and intensity columns
#' @param spec_ref The reference spectrum dataframe with mz and intensity columns
#' @param mz_tolerence_unit The mz tolerence unit (ppm or Da)
#' @param mz_tolerence_unit Value of mz tolerence
#' @param mz_range The mz range used to calculate score
#' @return The similarity score
#' @examples
#' calc_spectrum_similarity <- function (spec_query, spec_ref, mz_tolerence_unit = "ppm", 
#'                                       mz_tolerence = 10, mz_range = c(0, 1600))
#' @export
calc_spectrum_similarity <- function (spec_query = NULL, spec_ref = NULL, mz_tolerence_unit = "ppm", 
                                      mz_tolerence = 10, mz_range = c(0, 1600)){
  message("Calc Spectrum Similarity Started...")
  
  if(identical(spec_query, NULL)){
    warning("No spec_query was given")
    return (NULL)
  }

  if ((!class(spec_query)=='data.frame') | (!any(c("mz", "intensity") %in% colnames(spec_query)))){
    warning("Please input valid spec_query (a dataframe containing 'mz' and 'intensity' columns)")
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
  
  if (identical(mz_range, NULL)){
    warning("The mz_range is NULL")
    return (NULL)  
  }    
  
  query_tmp <- data.frame(mz = spec_query[, 1], intensity = spec_query[, 2])
  query_tmp$normalized <- (query_tmp$intensity/max(query_tmp$intensity)) * 100
  query_tmp <- subset(query_tmp, query_tmp$mz >= mz_range[1] & query_tmp$mz <= mz_range[2])
  query_score <- data.frame(mz = query_tmp$mz, intensity = query_tmp$normalized)    
  ref_tmp <- data.frame(mz = spec_ref[, 1], intensity = spec_ref[, 2])
  ref_tmp$normalized <- (ref_tmp$intensity/max(ref_tmp$intensity)) * 100
  ref_tmp <- subset(ref_tmp, ref_tmp$mz >= mz_range[1] & ref_tmp$mz <= mz_range[2])
  ref_score <- data.frame(mz = ref_tmp$mz, intensity = ref_tmp$normalized)
  query_score$delta <- suppressMessages(sapply(query_score[, 1], function(x) PollyMetaboR::calculate_mz_delta(x, mz_tolerence_unit = mz_tolerence_unit, mz_tolerence = mz_tolerence)))
  ref_score$delta <- suppressMessages(sapply(ref_score[, 1], function(x) PollyMetaboR::calculate_mz_delta(x, mz_tolerence_unit = mz_tolerence_unit, mz_tolerence = mz_tolerence)))
  
  for (i in 1:nrow(ref_score)) {
    query_score[, 1][ref_score[, 1][i] >= query_score[, 1] - query_score[, "delta"] & ref_score[, 1][i] <= query_score[, 1] + query_score[, "delta"]] <- ref_score[, 1][i]
    ref_score[, 1][ref_score[, 1][i] >= ref_score[, 1] - ref_score[, "delta"] & ref_score[, 1][i] <= ref_score[, 1] + ref_score[, "delta"]] <- ref_score[, 1][i]  
  }
  
  query_score$delta <- NULL
  ref_score$delta <- NULL
  query_score <- query_score[order(query_score[, 1], -abs(query_score[, 2])), ]                                                    
  query_score <- query_score[!duplicated(query_score[, 1]), ]
  ref_score <- ref_score[order(ref_score[, 1], -abs(ref_score[, 2])), ]                                                    
  ref_score <- ref_score[!duplicated(ref_score[, 1]), ]
  
  alignment <- merge(query_score, ref_score, by = 1, all = TRUE)
  names(alignment) <- c("mz", "intensity_query", "intensity_ref")
  matched_score <- length(unique(alignment[rowSums(is.na(alignment[ , c(2,3)])) == 0, ][, 1])) / nrow(ref_score)
  alignment[, c(2, 3)][is.na(alignment[, c(2, 3)])] <- 0                                  
  alignment <- alignment[order(alignment[, 1], -abs(alignment[, 2])), ]                                                    
  alignment <- alignment[!duplicated(alignment[, 1]), ]
  
  calc_dot_product <- function(u, v){
    sum_u <- sum(u)  
    for (i in 1:length(u)) {
      if (v[i] == 0) {u[i] <- u[i] / 2}
      query_weight = 1.0 / (1.0 + ((u[i]) / (sum_u - 50)))
      u[i] <- query_weight * u[i]  
    }
    dot_score <- as.vector((u %*% v)/(sqrt(sum(u^2)) * sqrt(sum(v^2))))
    
    return (dot_score)  
  }
  
  dot_product_score <- calc_dot_product(alignment[, 2], alignment[, 3])
  reverse_dot_product_score <- calc_dot_product(alignment[, 3], alignment[, 2])                                           
  similarity_score <-  round((dot_product_score + reverse_dot_product_score + matched_score) / 3, 4)
  
  message("Calc Spectrum Similarity Completed...")
  
  return (similarity_score)
}
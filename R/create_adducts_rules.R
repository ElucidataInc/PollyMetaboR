#' create_adducts_rules
#'
#' It creates adducts rules.
#'
#' @param ionlist_df File of known charged ions, an example is found in CAMERA/lists/ions.csv
#' @param neutralloss_df File of known neutral losses, an example is found in CAMERA/lists/neutralloss.csv
#' @param neutraladdition_df File of known adducts, an example is found in CAMERA/lists/lists/neutraladdition.csv
#' @param polarity The polarity of the data
#' @param maxcharge 
#' @param mol Number of molecules (xM) included in the molecule charge
#' @param nion .
#' @param nnloss .
#' @param nnadd .
#' @param nh .
#' @param lib.loc Path to local R library
#' @return A dataframe of adducts rules
#' @examples
#' create_adducts_rules(ionlist_df, neutralloss_df, neutraladdition_df, polarity = 'positive')
#' @import CAMERA
#' @export
create_adducts_rules <- function(ionlist_df = NULL, neutralloss_df = NULL, neutraladdition_df = NULL,
                                 polarity = 'positive', maxcharge = 3, mol = 3, nion = 2,
                                 nnloss = 1, nnadd = 1, nh = 2, lib.loc = .libPaths()){
    message("Create Adducts Rules Started...")
    require(CAMERA)
    
    rule_object <- new("ruleSet")
    rule_object@ionlist <- ionlist_df
    rule_object@neutralloss <- neutralloss_df
    rule_object@neutraladdition <- neutraladdition_df
    rule_object <- CAMERA::setParams(rule_object, maxcharge = maxcharge, mol = mol, nion = nion,
                  nnloss = nnloss, nnadd = nnadd, nh = nh, polarity = polarity, lib.loc = lib.loc)
    rule_object <- CAMERA::generateRules(rule_object)
    rules <- rule_object@rules

    message("Create Adducts Rules Completed...")
    
    return(rules)
}